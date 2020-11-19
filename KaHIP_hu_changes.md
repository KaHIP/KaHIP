# Changes and addition for KaHip_hu

*The motivation* is to extend ParHip so that it also perform mapping when doing refinement. Instead of the usual refinement where a vertex is moved between blocks if this move will improve the cut, now we also want to take the Processing Element's (PE) communication into account and move vertices
if this improves the communication costs. A description of the system is given as input and passed to the parallel refiment algorithm that uses the information when considering vertex moves. 
For the sequantial version, KaHip, several techniques to do that were introduced in the 
paper [High-Quality Hierarchical Process Mapping](https://arxiv.org/pdf/2001.07134.pdf).
On of those is label propagation and we will follow this as it fits better to the 
parallel partitioning scenario.

**Note**: In the present version, this is applicable only when the system has a **tree architecture**.

Below are the necessary additions to the code. These were performed mostly during September-October, 2020.


## Technical changes

The system description must be given to ParHip through the command line. As we only target near-homegeneous (PEs are identical in their specs, communication costs can vary)
tree system this can be done using two extra arguments: **`hierarchy_parameter_string`** 
and **`distance_parameter_string`**. 

Both parameters is a sequence of numbers separated by `:`. 
If hierarchy_parameter_string=a1:a2:a3:...:ai:...:ak this means that the system has k levels, the nodes in level `i-1` have `ai` "children" and in total we have a1\*a2\*...\*ak leaves. 
In this system description
only the leaves are PEs and intermediate tree nodes account for sockets, racks, compute nodes etc.
Similarly, distance_parameter_string=d1:d2:...:dk gives the communication costs, ie, the "distance"  between leaves/PEs: leaves that have the same "father" in the tree have cost d1, the same "grandfather" have cost d2 etc. In other words, the communication cost between two leaves is relevant to the level
of their least common ancestor in the tree. Naturally, (but not necessarily) the distances increases,
as PEs that belong to different compute nodes, for example, have higher communication costs 
than PEs in the same node.

Example: A system has 4 compute nodes, where each node has 2 sockets, each sockets 4 CPUs and each CPU 6 cores. The communication costs is 1 for PEs in the same core, 5 for PEs in the same CPU but different cores, 20 for PEs in the same socket but different CPU and 100 for PEs in different compute node.
Then, ParHIP will be called like 

`mpirun -n 4 parhip file.graph --hierarchy_parameter_string=6:4:2:4 --distance_parameter_string=1:5:20:100 --k=192 --preconfiguration= ...`

**Note**: `k` must be equal a1\*a2\*...\*ak.


## Algorithmic changes

The core algorithmic adaptations were done in the refinement step located in `parallel/parallel_src/lib/parallel_label_compress/parallel_label_compress.h`. Remember that ParHip uses the multilevel method:
it perform several coarsening steps that result is much smaller graph than the original, 
then it finds
a partition on this coarsest graph and finally uncoarsens the graph
back to the original in steps and in each steps moves vertices between blocks 
if this move will improve the cut. This sequence is called a `v-cycle`;
By default, performs 2 `v-cycles`.
That is, after the first refinement, the coarsening, partitioning and refinement steps are repeated
one more tim, now in the partitioning graph obtained from the first `v-cycle`.

By providing the system's description (with communication costs), during refinement we can minimize, not for edge cut, but for the communication costs (CoCo), a metric that better represents the application's costs. 
Now, we will move nodes using a label propagation technique.
During the refinement steps, in each PE independently, we go over its local verices and, for each vertex `v`,
we visit all of its neighbors and we 
gather all the blocks where each neighbor belong to, notated `R(v)`.
At this point, `v` belongs to block `P(v)`. Then, we check, for all blocks in `R(v)`, if reassigning
`v` to one of these blocks improves the communication costs. We find the assignment that offers the maximum improvement, we update `P(v)` if needed and continue to the next vertex.
After a certain number of local nodes are visited, PEs communicate in order to update the ghost nodes partition.
When all local nodes are visited, we repeat the same procedure. How many label propagation iterations
we will perform is controlled by the command line parameter `label_iterations_refinement`.
The whole procedure is repeated on every uncoarsening step until we get the original graph.


## Miscellaneous additions and changes

The information about the system description are held in a processor_tree, found in `parallel/parallel_src/lib/data_structure/processor_tree.h`. The processor tree is constructed 
using the command line parameters mentioned above. 
The main feature is that, given two leaves, it calculates their distance in the tree efficiently 
by using bit label operations.

In `parallel/parallel_src/lib/tools` we added two new classes. The `distributed_quality_metrics` 
calculate various metric related to mapping like CoCo, max and total dilation and congestion.
In `system_info.h` we offer a function to output the memory usage of the application at run-time 
mostly for debugging purposes.


## Examples

### Installation

Use or adpat the `compile_withcmake.sh` scripts and follow the standart KaHip instructions.

### Run ParHip

The following command runs ParHip with 8 PEs, for a random hyperbolic graph and partitions
it into 96 blocks for hierarchical system with 3 levels using the `fastsocial` configuration
and performing only one `v-cycle`.

>mpirun -n 8 ./deploy/parhip /meshes/rhg_n23_d16.graph --distance_parameter_string 1:40:1600 --preconfiguration fastsocial --k 96 --hierarchy_parameter_string 12:4:2 --num_vcycles 1

Similarly, to partition a mesh into 192 blocks using 12 PEs you can call (no `--num_vcycles` 
parameter was given so ParHip will perform 2 `v-cycles` by default).

>mpirun -n 12 ./deploy/parhip /meshes/hugetrace-00020.graph --distance_parameter_string 1:40:400 --preconfiguration fastmesh --k 192 --hierarchy_parameter_string 12:4:4


## Future plans

When tested, the procedure obtained promising results as it manages to lower CoCo 
without significant impact in the running time. Unfortunately, it turned out that ParHip cannot handle some graph families very well. Precisely, on generated graph with the R-mat and
Barabasi-Albert or some real life graphs like `com-friendster`, ParHip crushes with out of memory
errors or running time is prohibitively high. On the other hand, it can handle random
hyerbolic graphs efficiently.

The reasons for that are not completely understood. 
ParHip operates in the following phases: coarsen, partition, refine, coarsen, 
improve partition, refine.
One possibility is that we get out of memory errors on the second coarsening phase.
Because in this step we have a partition, contrary to the first coarsening phase,
edges with endpoints that belong to different blocks are not contracted. This greatly affects
the number of contractions, leaving a too big coarsest graph. Additionally, the intermediate
coarsened graphs that are kept in memory are also too big.

**Possible solution** proposed by Christian Schultz. In the second coarsening round, before
the second partition step, i.e., for the coarsest graph, check the graph's size (for example,
it could be a constant factor of the first coarsest graph). If the graph is too big,
store it with the partition. Let that be `Gc` and its partition `Pc`.
Then neglect the partition and continue coarsening for some rounds.
This allows us to contract edges that would not be contracted otherwise.
Perform partitioning in this coarsest graph and uncoarsen+refine until we get graph `Gc`.
Now, we have the `Gc` with `Pc` from before and the new graph with another partition.
Compare the two partition and keep the best of the two. Continue with uncoarsening+refinement
as usual until the original graph.