#include "VieClus_adapter.h"

#include <VieClus_interface.h>

#include "tools/timer.h"

VieClus_adapter::VieClus_adapter() {
    VieClus::init(nullptr, nullptr);
}

VieClus_adapter::~VieClus_adapter() {
    VieClus::finalize();
}

VieClus_adapter& VieClus_adapter::instance() {
    static VieClus_adapter adapter;
    return adapter;
}

std::vector<int> VieClus_adapter::compute_modularity_clustering(graph_access &G, const PartitionConfig &partition_config) {
    VieClus::Graph graph{
            static_cast<int>(G.number_of_nodes()),
            G.UNSAFE_metis_style_xadj_array(),
            G.UNSAFE_metis_style_adjncy_array(),
            G.UNSAFE_metis_style_vwgt_array(),
            G.UNSAFE_metis_style_adjwgt_array()
    };
    int no_clusters = -1;
    std::vector<int> partition_map(G.number_of_nodes());

    timer vieclus_timer;
    double modularity = VieClus::run(graph,
            partition_config.vieclus_time_limit,
            partition_config.seed,
            &no_clusters,
            partition_map.data());
    std::cout << "modularity: " << modularity << " with " << no_clusters << " clusters in: " << vieclus_timer.elapsed() << " s" << std::endl;

    delete[] graph.xadj;
    delete[] graph.adjncy;
    delete[] graph.vwgt;
    delete[] graph.adjwgt;
    return partition_map;
}