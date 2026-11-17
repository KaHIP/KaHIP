/******************************************************************************
 * jet_refinement.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <vector>

#include "jet_refinement.h"
#include "../quotient_graph_refinement/complete_boundary.h"

jet_refinement::jet_refinement() {
}

jet_refinement::~jet_refinement() {
}

namespace {

const int    JET_PATIENCE       = 12;    // iterations without a significant improvement
const int    JET_MAX_ITERATIONS = 100;
const double JET_PHI            = 0.999; // significant: below phi times the best cut

struct move { NodeID v; PartitionID to; };

// Work is proportional to the moved nodes and their neighborhoods: only nodes next to a move are
// evaluated again, the cut is updated by the change of the moved edges, and the best partition is
// kept as a journal of the nodes moved since it was seen.
class jet {
public:
        jet(graph_access & graph, PartitionID blocks) : G(graph), k(blocks), n(graph.number_of_nodes()),
                part(n), bw(blocks, 0), conn(blocks, 0),
                dest(n, 0), gain(n, 0), is_cand(n, 0), in_cand(n, 0), dirty(n, 0), locked(n, 0),
                bnd(n, 0), bpos(n, 0), wdeg(n, 0), moved_now(n, 0), old_part(n, 0), journaled(n, 0), best_part(n, 0), ever(n, 0) {
                forall_nodes(G, v) {
                        part[v] = G.getPartitionIndex(v);
                        bw[part[v]] += G.getNodeWeight(v);
                } endfor
                cut = 0;
                forall_nodes(G, v) {
                        forall_out_edges(G, e, v) {
                                wdeg[v] += G.getEdgeWeight(e);
                        } endfor
                        update_bnd(v);
                        mark_dirty(v);
                        forall_out_edges(G, e, v) {
                                if(part[G.getEdgeTarget(e)] != part[v]) cut += G.getEdgeWeight(e);
                        } endfor
                } endfor
                cut /= 2;
        }

        graph_access & G;
        const PartitionID k;
        const NodeID n;
        std::vector<PartitionID> part;
        std::vector<long long>   bw;
        EdgeWeight cut;

        long long max_weight() const { return *std::max_element(bw.begin(), bw.end()); }

        // one synchronous label propagation step with the afterburner
        void lp(double c) {
                for(NodeID v : dirty_list) {
                        dirty[v] = 0;
                        evaluate(v, c);
                        if(is_cand[v] && !in_cand[v]) { in_cand[v] = 1; cand.push_back(v); }
                }
                dirty_list.clear();
                size_t j = 0;
                for(NodeID v : cand) {
                        if(is_cand[v]) cand[j++] = v; else in_cand[v] = 0;
                }
                cand.resize(j);

                // a neighbor counts as moved if it is a candidate of higher rank
                std::vector<move> keep;
                for(NodeID v : cand) {
                        EdgeWeight g = 0;
                        forall_out_edges(G, e, v) {
                                NodeID u = G.getEdgeTarget(e);
                                PartitionID pu = part[u];
                                if(is_cand[u] && (gain[u] > gain[v] || (gain[u] == gain[v] && u < v))) pu = dest[u];
                                if(pu == dest[v])      g += G.getEdgeWeight(e);
                                else if(pu == part[v]) g -= G.getEdgeWeight(e);
                        } endfor
                        if(g >= 0) keep.push_back({v, dest[v]});
                }

                apply(keep);
                unlock_all();
                for(const move & m : keep) { locked[m.v] = 1; locked_list.push_back(m.v); }
        }

        // weak rebalancing: move the cheapest nodes (loss per weight) out of every overloaded block
        void rebalance(long long ub) {
                unlock_all();
                struct option { double score; NodeID v; PartitionID to; };
                std::vector< std::vector<option> > options(k);
                PartitionID lightest = 0;
                for(PartitionID b = 1; b < k; b++) if(bw[b] < bw[lightest]) lightest = b;

                for(NodeID v = 0; v < n; v++) {
                        PartitionID own = part[v];
                        if(bw[own] <= ub) continue;
                        NodeWeight w = G.getNodeWeight(v);
                        if(!bnd[v]) { // interior: all edges are lost, the lightest block takes it
                                if(lightest != own) options[own].push_back({(double)wdeg[v] / (double)std::max<NodeWeight>(w, 1), v, lightest});
                                continue;
                        }
                        gather(v);
                        PartitionID best = lightest;
                        EdgeWeight  best_conn = 0;
                        for(PartitionID b : touched) {
                                if(b == own || bw[b] + (long long)w > ub) continue;
                                if(conn[b] > best_conn || (conn[b] == best_conn && bw[b] < bw[best])) {
                                        best = b; best_conn = conn[b];
                                }
                        }
                        if(best != own) {
                                double loss = (double)(conn[own] - best_conn);
                                options[own].push_back({loss / (double)std::max<NodeWeight>(w, 1), v, best});
                        }
                        clear();
                }

                std::vector<long long> nbw = bw;
                std::vector<move> moves;
                for(PartitionID b = 0; b < k; b++) {
                        if(options[b].empty()) continue;
                        std::stable_sort(options[b].begin(), options[b].end(),
                                         [](const option & x, const option & y) { return x.score < y.score; });
                        for(const option & o : options[b]) {
                                if(nbw[b] <= ub) break;
                                NodeWeight w = G.getNodeWeight(o.v);
                                if(nbw[o.to] + (long long)w > ub) continue;
                                nbw[b] -= w; nbw[o.to] += w;
                                moves.push_back({o.v, o.to});
                        }
                }
                apply(moves);
        }

        // the partition is now the best one seen
        void mark_best() {
                for(NodeID v : journal) journaled[v] = 0;
                journal.clear();
        }
        // go back to the best partition seen, returns the nodes that differ from the start
        const std::vector<NodeID> & restore_best() {
                for(NodeID v : journal) { part[v] = best_part[v]; journaled[v] = 0; }
                journal.clear();
                return ever_list;
        }

private:
        std::vector<EdgeWeight>  conn;   // connectivity of the gathered node to each block
        std::vector<PartitionID> touched;
        std::vector<PartitionID> dest;
        std::vector<EdgeWeight>  gain;
        std::vector<char> is_cand, in_cand, dirty, locked, bnd;
        std::vector<NodeID> bpos, blist, cand, dirty_list, locked_list;
        std::vector<EdgeWeight> wdeg;   // weighted degree
        std::vector<char> moved_now;
        std::vector<PartitionID> old_part;
        std::vector<char> journaled;
        std::vector<PartitionID> best_part;
        std::vector<NodeID> journal;
        std::vector<char> ever;
        std::vector<NodeID> ever_list;

        void gather(NodeID v) {
                forall_out_edges(G, e, v) {
                        PartitionID b = part[G.getEdgeTarget(e)];
                        if(conn[b] == 0) touched.push_back(b);
                        conn[b] += G.getEdgeWeight(e);
                } endfor
        }
        void clear() {
                for(PartitionID b : touched) conn[b] = 0;
                touched.clear();
        }
        void mark_dirty(NodeID v) {
                if(!dirty[v]) { dirty[v] = 1; dirty_list.push_back(v); }
        }
        void unlock_all() {
                for(NodeID v : locked_list) { locked[v] = 0; mark_dirty(v); }
                locked_list.clear();
        }
        void update_bnd(NodeID v) {
                bool b = false;
                forall_out_edges(G, e, v) {
                        if(part[G.getEdgeTarget(e)] != part[v]) { b = true; break; }
                } endfor
                if(b && !bnd[v]) {
                        bnd[v] = 1; bpos[v] = blist.size(); blist.push_back(v);
                } else if(!b && bnd[v]) {
                        bnd[v] = 0;
                        NodeID last = blist.back();
                        blist[bpos[v]] = last; bpos[last] = bpos[v];
                        blist.pop_back();
                }
        }
        void evaluate(NodeID v, double c) {
                is_cand[v] = 0;
                if(!bnd[v] || locked[v]) return;
                PartitionID own = part[v];
                gather(v);
                PartitionID best = own;
                EdgeWeight  best_conn = -1;
                for(PartitionID b : touched) {
                        if(b == own) continue;
                        if(conn[b] > best_conn || (conn[b] == best_conn && bw[b] < bw[best])) {
                                best = b; best_conn = conn[b];
                        }
                }
                if(best != own) {
                        EdgeWeight own_conn = conn[own];
                        EdgeWeight g        = best_conn - own_conn;
                        if(g >= 0 || -g < (EdgeWeight)(c * own_conn)) {
                                dest[v] = best; gain[v] = g; is_cand[v] = 1;
                        }
                }
                clear();
        }
        // apply all moves at once: weights, cut change, journal, and the nodes to evaluate again
        void apply(const std::vector<move> & moves) {
                for(const move & m : moves) {
                        NodeID v = m.v;
                        moved_now[v] = 1; old_part[v] = part[v];
                        if(!journaled[v]) { journaled[v] = 1; best_part[v] = part[v]; journal.push_back(v); }
                        if(!ever[v])      { ever[v] = 1; ever_list.push_back(v); }
                        NodeWeight w = G.getNodeWeight(v);
                        bw[part[v]] -= w; bw[m.to] += w;
                        part[v] = m.to;
                }
                for(const move & m : moves) {
                        NodeID v = m.v;
                        forall_out_edges(G, e, v) {
                                NodeID u = G.getEdgeTarget(e);
                                if(moved_now[u] && u < v) continue; // count edges between two moved nodes once
                                PartitionID ou = moved_now[u] ? old_part[u] : part[u];
                                EdgeWeight w = G.getEdgeWeight(e);
                                cut += w * ((part[u] != part[v]) - (ou != old_part[v]));
                        } endfor
                }
                for(const move & m : moves) {
                        NodeID v = m.v;
                        mark_dirty(v); update_bnd(v);
                        forall_out_edges(G, e, v) {
                                NodeID u = G.getEdgeTarget(e);
                                mark_dirty(u); update_bnd(u);
                        } endfor
                }
                for(const move & m : moves) moved_now[m.v] = 0;
        }
};

}

EdgeWeight jet_refinement::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {
        const NodeID n = G.number_of_nodes();
        const PartitionID k = G.get_partition_count();
        if(n == 0 || k < 2) return 0;

        bool unit_weights = true;
        forall_nodes(G, v) {
                if(G.getNodeWeight(v) != 1) { unit_weights = false; break; }
        } endfor

        jet s(G, k);
        const long long ub      = (long long)config.upper_bound_partition;
        const long long allowed = std::max(ub, s.max_weight()); // never leave a worse balance than given
        const double    c       = unit_weights ? 0.25 : 0.75;

        const EdgeWeight start_cut = s.cut;
        EdgeWeight best_cut = start_cut;

        int patience = 0;
        for(int it = 0; it < JET_MAX_ITERATIONS && patience < JET_PATIENCE; it++) {
                if(s.max_weight() > ub) s.rebalance(ub);
                else                    s.lp(c);

                if(s.max_weight() <= allowed && s.cut < best_cut) {
                        if(s.cut < JET_PHI * best_cut) patience = 0; else patience++;
                        best_cut = s.cut;
                        s.mark_best();
                } else {
                        patience++;
                }
        }

        if(best_cut >= start_cut) return 0;

        // apply the difference to the start partition one node at a time
        const std::vector<NodeID> & changed = s.restore_best();
        for(NodeID v : changed) {
                PartitionID from = G.getPartitionIndex(v);
                PartitionID to   = s.part[v];
                if(from == to) continue;
                G.setPartitionIndex(v, to);
                boundary_pair pair;
                pair.k   = config.k;
                pair.lhs = from;
                pair.rhs = to;
                boundary.postMovedBoundaryNodeUpdates(v, &pair, true, true);

                NodeWeight w = G.getNodeWeight(v);
                boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
                boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
                boundary.setBlockWeight( from, boundary.getBlockWeight(from)-w);
                boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+w);
        }

        return start_cut - best_cut;
}
