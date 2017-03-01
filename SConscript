#/******************************************************************************
# * SConscript
# *
# * Source of KaHIP -- Karlsruhe High Quality Partitioning.
# *
# ******************************************************************************
# * Copyright (C) 2015 Christian Schulz <christian.schulz@kit.edu>
# *
# * This program is free software: you can redistribute it and/or modify it
# * under the terms of the GNU General Public License as published by the Free
# * Software Foundation, either version 2 of the License, or (at your option)
# * any later version.
# *
# * This program is distributed in the hope that it will be useful, but WITHOUT
# * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# * more details.
# *
# * You should have received a copy of the GNU General Public License along with
# * this program.  If not, see <http://www.gnu.org/licenses/>.
# *****************************************************************************/


# The main SConscript file for the code.
#
# We simply import the main environment and then define the targets.  This
# submodule contains a sequential matching and contraction code and tests for
# the code.
import platform
import sys

# Get the current platform.
SYSTEM = platform.uname()[0]

Import('env')

# Build a library from the code in lib/.
libkaffpa_files = [   'lib/data_structure/graph_hierarchy.cpp',
                      'lib/algorithms/strongly_connected_components.cpp',
                      'lib/algorithms/topological_sort.cpp',
                      'lib/algorithms/push_relabel.cpp',
                      'lib/io/graph_io.cpp',
                      'lib/tools/quality_metrics.cpp',
                      'lib/tools/random_functions.cpp',
                      'lib/tools/graph_extractor.cpp',
                      'lib/tools/misc.cpp',
                      'lib/tools/partition_snapshooter.cpp',
                      'lib/partition/graph_partitioner.cpp',
                      'lib/partition/w_cycles/wcycle_partitioner.cpp',
                      'lib/partition/coarsening/coarsening.cpp',
                      'lib/partition/coarsening/contraction.cpp',
                      'lib/partition/coarsening/edge_rating/edge_ratings.cpp',
                      'lib/partition/coarsening/matching/matching.cpp',
                      'lib/partition/coarsening/matching/random_matching.cpp',
                      'lib/partition/coarsening/matching/gpa/path.cpp',
                      'lib/partition/coarsening/matching/gpa/gpa_matching.cpp',
                      'lib/partition/coarsening/matching/gpa/path_set.cpp',
                      'lib/partition/coarsening/clustering/node_ordering.cpp',
                      'lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp',
                      'lib/partition/initial_partitioning/initial_partitioning.cpp',
                      'lib/partition/initial_partitioning/initial_partitioner.cpp',
                      'lib/partition/initial_partitioning/initial_partition_bipartition.cpp',
                      'lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp',
                      'lib/partition/initial_partitioning/bipartition.cpp',
                      'lib/partition/initial_partitioning/initial_node_separator.cpp',
                      'lib/partition/uncoarsening/uncoarsening.cpp',
                      'lib/partition/uncoarsening/separator/area_bfs.cpp',
                      'lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp',
                      'lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp',
                      'lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp',
                      'lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp',
                      'lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp',
                      'lib/partition/uncoarsening/refinement/mixed_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/refinement.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp',
                      'lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp',
                      'lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp',
                      'lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp',
                      'lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp', 
                      'lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp', 
                      'lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp', 
                      'lib/partition/uncoarsening/refinement/node_separators/greedy_ns_local_search.cpp', 
                      'lib/partition/uncoarsening/refinement/node_separators/fm_ns_local_search.cpp', 
                      'lib/partition/uncoarsening/refinement/node_separators/localized_fm_ns_local_search.cpp', 
                      'lib/algorithms/cycle_search.cpp',
                      'lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp',
                      'lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp'
                      ]

libkaffpa_parallel_async  = ['lib/parallel_mh/parallel_mh_async.cpp',
                             'lib/parallel_mh/population.cpp',
                             'lib/parallel_mh/galinier_combine/gal_combine.cpp',
                             'lib/parallel_mh/galinier_combine/construct_partition.cpp',
                             'lib/parallel_mh/exchange/exchanger.cpp',
                             'lib/tools/graph_communication.cpp',
                             'lib/tools/mpi_tools.cpp' ]

libmapping                = ['lib/mapping/local_search_mapping.cpp',
                             'lib/mapping/full_search_space.cpp',
                             'lib/mapping/full_search_space_pruned.cpp',
                             'lib/mapping/communication_graph_search_space.cpp',
                             'lib/mapping/fast_construct_mapping.cpp',
                             'lib/mapping/construct_distance_matrix.cpp',
                             'lib/mapping/mapping_algorithms.cpp',
                             'lib/mapping/construct_mapping.cpp' ]

if env['program'] == 'kaffpa':
        env.Append(CXXFLAGS = '-DMODE_KAFFPA')
        env.Append(CCFLAGS  = '-DMODE_KAFFPA')
        env.Program('kaffpa', ['app/kaffpa.cpp']+libkaffpa_files+libmapping, LIBS=['libargtable2','gomp'])

if env['program'] == 'evaluator':
        env.Append(CXXFLAGS = '-DMODE_EVALUATOR')
        env.Append(CCFLAGS  = '-DMODE_EVALUATOR')
        env.Program('evaluator', ['app/evaluator.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'node_separator':
        env.Append(CXXFLAGS = ' -DMODE_NODESEP')
        env.Append(CCFLAGS  = ' -DMODE_NODESEP')
        env.Program('node_separator', ['app/node_separator_ml.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'label_propagation':
        env.Append(CXXFLAGS = '-DMODE_LABELPROPAGATION')
        env.Append(CCFLAGS  = '-DMODE_LABELPROPAGATION')
        env.Program('label_propagation', ['app/label_propagation.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'partition_to_vertex_separator':
        env.Append(CXXFLAGS = '-DMODE_PARTITIONTOVERTEXSEPARATOR')
        env.Append(CCFLAGS  = '-DMODE_PARTITIONTOVERTEXSEPARATOR')
        env.Program('partition_to_vertex_separator', ['app/partition_to_vertex_separator.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'interfacetest':
        env['CXX'] = 'mpicxx'
        env.Append(CXXFLAGS = '-DMODE_KAFFPA')
        env.Append(CCFLAGS  = '-DMODE_KAFFPA')
        env.Program('interface_test', ['app/interface_test.cpp','interface/kaHIP_interface.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'improve_vertex_separator':
        env.Append(CXXFLAGS = '-DMODE_IMPROVEVERTEXSEPARATOR')
        env.Append(CCFLAGS  = '-DMODE_IMPROVEVERTEXSEPARATOR')
        env.Program('improve_vertex_separator', ['app/improve_vertex_separator.cpp']+libkaffpa_files, LIBS=['libargtable2','gomp'])

if env['program'] == 'kaffpaE':
        env.Append(CXXFLAGS = '-DMODE_KAFFPAE')
        env.Append(CCFLAGS  = '-DMODE_KAFFPAE')

        if SYSTEM == 'Darwin':
                env['CXX'] = 'openmpicxx'
        else:
                env['CXX'] = 'mpicxx'
        env.Program('kaffpaE', ['app/kaffpaE.cpp']+libkaffpa_files+libkaffpa_parallel_async, LIBS=['libargtable2','gomp'])

if env['program'] == 'graphchecker':
        env.Append(CXXFLAGS = '-DMODE_GRAPHCHECKER')
        env.Append(CCFLAGS  = '-DMODE_GRAPHCHECKER')
        env.Program('graphchecker', ['app/graphchecker.cpp'], LIBS=['libargtable2','gomp'])

if env['program'] == 'library':
        env.Append(CXXFLAGS = '-fPIC')
        env.Append(CCFLAGS  = '-fPIC')
        SConscript('interface/SConscript',exports='env')
