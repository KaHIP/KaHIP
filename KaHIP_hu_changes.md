# Changes and addition for KaHIP_hu

The motivation is to extend ParHIP so that it also perform mapping when doing refinement. Instead of the usual refinement where a vertex is moved if this move will improve the cut, now we also want to take the Processing Element's (PE) communication into account. A description of the system is given as input and passed to the parallel refiment algorithm that used the information when considering vertex moves. 

**Note**: In the present version, this is applicable only when the system has a **tree architecture**.

Below are the necessary additions to the code. These were performed mostly during September-October, 2020.

## Necassary technical changes

The system description whould somehow be given to ParHIP through the command line. As we only target near-homegeneous tree system this can be done using two extra arguments: **`hierarchy_parameter_string`** 
and **`distance_parameter_string`**. Both parameters is a sequence of numbers separated by `:`. 
If hierarchy_parameter_string=a1:a2:a3:...:ai:...:ak this means that the system has k levels, the nodes in level `i-1` have `ai` "children" and in total we have a1\*a2\*...\*ak leaves. 
In this system description
only the leaves are PEs and intermediate tree nodes account for sockets, racks, compute nodes etc.
Similarly, distance_parameter_string=d1:d2:...:dk gives the communication costs, ie, the "distance"  between leaves/PEs: leaves that have the same "father" in the tree have cost d1, the same "grandfather" have cost d2 etc. In other words, the communication cost between two leaves is relevant to the level
of their least common ancestor in the tree. Naturally, (but not necessarily) the distances increase
as PEs that belong to different compute node for example have higher communication costs than PE
in the same node.

Example: A system has 4 compute nodes, where each node has 2 sockets, each sockets 4 CPUs and each CPU 6 cores. The communication costs is 1 for PEs in the same core, 5 for PEs in the same CPU but different cores, 20 for PEs in the same socket but different CPU and 100 for PEs in different compute node.
Then, ParHIP will be called like 

`mpirun -n 4 parhip file.graph --hierarchy_parameter_string=6:4:2:4 --distance_parameter_string=1:5:20:100 --k=192 --preconfiguration= ...`

**Note**: `k` must be equal a1\*a2\*...\*ak.

## Algorithmic changes

The core algorithmic adaptations were done in the refinement step located in `parallel/parallel_src/lib/parallel_label_compress/parallel_label_compress.h`. Remember that ParHIP uses the multilevel method:
it perform several coarsening steps that result is much smaller graph than the original, then it finds
a partition on this coarsest graph using an integer programming solver and finally uncoarsens the graph
back to the original in steps and in each steps moves vertices if this move will improve the cut.

By providing the system's description (with communication costs) we can minimize, not for edge cut, but for the communication costs (CoCo), a metric that better represents the application's costs. 
Now, we will move nodes using a label propagation technique.
During the refinement steps, in each PE independently, we go over its local verices and, for each vertex `v`,
we visit all of its neighbors and we 
gather all the blocks where each neighbor belong to, notated `R(v)`.
At this point, `v` belongs to block `P(v)`. Then, we check, for all blocks in `R(v)`, if reassigning
`v` to one of these blocks would improve the communication costs. We find the assignement that offers the maximum improvement and we update `P(v)` and continue to the next.
Every a certain number of local nodes are visited, PEs communicate in order to update the ghost nodes partition.
When all local nodes are visited, we repeat the same procedure. How many label propagation iterations
we will perform is controlled by the command line parameter `label_iterations_refinement`.
The whole procedure is repeated on every uncoarsening step until we get the original graph.

## Miscellanous additions and changes

The information about the system description are held in a processor_tree, found in `parallel/parallel_src/lib/data_structure/processor_tree.h`. The processor tree is constructed 
using the command line parameters mentioned above. 
The main feature is that, given two leaves, it calculates their distance in the tree efficiently 
by using bit label operations.

In `parallel/parallel_src/lib/tools` we added two new classes. The `distributed_quality_metrics` 
calculate various metric related to mapping like CoCo, max and total dilation and congestion.
In `system_info.h` we offer a function to output the memory usage of the application at runtime 
mostly for debugging purposes.