cmake_minimum_required(VERSION 3.10)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# --march=nativeflag
option(NONATIVEOPTIMIZATIONS "Disable --march=native optimizations" OFF)

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  if( NOT NONATIVEOPTIMIZATIONS )
          add_definitions(-march=native)
  endif()
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-unused-result COMPILER_SUPPORTS_UNUSED)
if(COMPILER_SUPPORTS_UNUSED)
  add_definitions(-Wno-unused-result)
endif()
CHECK_CXX_COMPILER_FLAG(-Wno-sign-compare COMPILER_SUPPORTS_NOSIGNCOMP)
if(COMPILER_SUPPORTS_NOSIGNCOMP)
  add_definitions(-Wno-sign-compare)
endif()


# Check dependencies
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()
find_library(LIB_METIS metis)
if(LIB_METIS)
  message(STATUS "Metis support detected")
  add_definitions("-DUSEMETIS")
  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include_directories(/usr/local/include)
  endif()
endif()

# 64 Bit option
option(64BITMODE "64 bit mode" OFF)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()


# Optionally disable all MPI-dependent targets
option(NOMPI "disable all targets that depend on MPI (kaffpaE, ParHIP)" OFF)

# ParHIP
option(PARHIP "build ParHIP" ON)
option(DETERMINISTIC_PARHIP "enforce deterministic computations in ParHIP" OFF)

# Look for MPI (needed for ParHIP and kaffpaE)
# Report which MPI we actually found,
# may need use MPI_HOME hint to get the proper one (eg, on openSUSE)
if(NOT NOMPI)
  # Always look for MPI since kaffpaE also requires it
  find_package(MPI REQUIRED)

  if(${MPI_C_FOUND})
    message("MPI detected (can use MPI_HOME hint to direct the detection...)")
    message(STATUS "MPI include: ${MPI_C_INCLUDE_DIRS}")
    message(STATUS "MPI library: ${MPI_C_LIBRARIES}")
  endif()

  if(PARHIP)
    message(STATUS "ParHIP build requested")
  else()
    message(STATUS "ParHIP build disabled")
  endif()
else()
  message(STATUS "Build without MPI dependency: kaffpaE and ParHIP disabled")
endif()


# tcmalloc
option(USE_TCMALLOC "if available, link against tcmalloc" OFF)


# ILP improver
option(USE_ILP "build local ILP improver - introduces dependency on Gurobi" OFF)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.2.2)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)

set(LIBKAFFPA_SOURCE_FILES
  lib/data_structure/graph_hierarchy.cpp
  lib/algorithms/strongly_connected_components.cpp
  lib/algorithms/topological_sort.cpp
  lib/algorithms/push_relabel.cpp
  lib/io/graph_io.cpp
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp
  lib/tools/graph_extractor.cpp
  lib/tools/misc.cpp
  lib/tools/partition_snapshooter.cpp
  lib/partition/graph_partitioner.cpp
  lib/partition/w_cycles/wcycle_partitioner.cpp
  lib/partition/coarsening/coarsening.cpp
  lib/partition/coarsening/contraction.cpp
  lib/partition/coarsening/edge_rating/edge_ratings.cpp
  lib/partition/coarsening/matching/matching.cpp
  lib/partition/coarsening/matching/random_matching.cpp
  lib/partition/coarsening/matching/gpa/path.cpp
  lib/partition/coarsening/matching/gpa/gpa_matching.cpp
  lib/partition/coarsening/matching/gpa/path_set.cpp
  lib/partition/coarsening/clustering/node_ordering.cpp
  lib/partition/coarsening/clustering/size_constraint_label_propagation.cpp
  lib/partition/initial_partitioning/initial_partitioning.cpp
  lib/partition/initial_partitioning/initial_partitioner.cpp
  lib/partition/initial_partitioning/initial_partition_bipartition.cpp
  lib/partition/initial_partitioning/initial_refinement/initial_refinement.cpp
  lib/partition/initial_partitioning/bipartition.cpp
  lib/partition/initial_partitioning/initial_node_separator.cpp
  lib/partition/uncoarsening/uncoarsening.cpp
  lib/partition/uncoarsening/separator/area_bfs.cpp
  lib/partition/uncoarsening/separator/vertex_separator_algorithm.cpp
  lib/partition/uncoarsening/separator/vertex_separator_flow_solver.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/problem_factory.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph.cpp
  lib/partition/uncoarsening/refinement/mixed_refinement.cpp
  lib/partition/uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.cpp
  lib/partition/uncoarsening/refinement/refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/two_way_flow_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/boundary_bfs.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/cut_flow_problem_solver.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/most_balanced_minimum_cuts/most_balanced_minimum_cuts.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/quotient_graph_scheduling.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/active_block_quotient_graph_scheduler.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/augmented_Qgraph_fabric.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.cpp
  lib/partition/uncoarsening/refinement/kway_graph_refinement/multitry_kway_fm.cpp
  lib/partition/uncoarsening/refinement/node_separators/greedy_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/fm_ns_local_search.cpp
  lib/partition/uncoarsening/refinement/node_separators/localized_fm_ns_local_search.cpp
  lib/algorithms/cycle_search.cpp
  lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.cpp
  lib/partition/uncoarsening/refinement/tabu_search/tabu_search.cpp
  extern/argtable3-3.2.2/argtable3.c)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})

if(NOT NOMPI)
  set(LIBKAFFPA_PARALLEL_SOURCE_FILES
    lib/parallel_mh/parallel_mh_async.cpp
    lib/parallel_mh/population.cpp
    lib/parallel_mh/galinier_combine/gal_combine.cpp
    lib/parallel_mh/galinier_combine/construct_partition.cpp
    lib/parallel_mh/exchange/exchanger.cpp
    lib/tools/graph_communication.cpp
    lib/tools/mpi_tools.cpp)
  add_library(libkaffpa_parallel OBJECT ${LIBKAFFPA_PARALLEL_SOURCE_FILES})
  target_include_directories(libkaffpa_parallel PUBLIC ${MPI_CXX_INCLUDE_PATH})
endif()

set(LIBMAPPING_SOURCE_FILES
  lib/mapping/local_search_mapping.cpp
  lib/mapping/full_search_space.cpp
  lib/mapping/full_search_space_pruned.cpp
  lib/mapping/communication_graph_search_space.cpp
  lib/mapping/fast_construct_mapping.cpp
  lib/mapping/construct_distance_matrix.cpp
  lib/mapping/mapping_algorithms.cpp
  lib/mapping/construct_mapping.cpp)
add_library(libmapping OBJECT ${LIBMAPPING_SOURCE_FILES})

set(LIBSPAC_SOURCE_FILES lib/spac/spac.cpp)
add_library(libspac OBJECT ${LIBSPAC_SOURCE_FILES})

set(NODE_ORDERING_SOURCE_FILES
  lib/node_ordering/min_degree_ordering.cpp
  lib/node_ordering/nested_dissection.cpp
  lib/node_ordering/ordering_tools.cpp
  lib/node_ordering/reductions.cpp)
add_library(libnodeordering OBJECT ${NODE_ORDERING_SOURCE_FILES})

# generate targets for each binary
add_executable(kaffpa app/kaffpa.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(kaffpa PRIVATE "-DMODE_KAFFPA")
target_link_libraries(kaffpa ${OpenMP_CXX_LIBRARIES})
install(TARGETS kaffpa DESTINATION bin)
if (USE_TCMALLOC)
  find_library(TCMALLOC_LIB tcmalloc)
  if (TCMALLOC_LIB)
    target_link_libraries(kaffpa ${TCMALLOC_LIB})
    message(STATUS "Using tcmalloc: ${TCMALLOC_LIB}")
  else ()
    message(STATUS "tcmalloc enabled but unavailable on this system")
  endif ()
endif ()

add_executable(global_multisection app/global_multisection.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(global_multisection PRIVATE "-DMODE_KAFFPA -DMODE_GLOBALMS")
target_link_libraries(global_multisection ${OpenMP_CXX_LIBRARIES})
install(TARGETS global_multisection DESTINATION bin)
if (USE_TCMALLOC)
  find_library(TCMALLOC_LIB tcmalloc)
  if (TCMALLOC_LIB)
    target_link_libraries(global_multisection ${TCMALLOC_LIB})
    message(STATUS "Using tcmalloc: ${TCMALLOC_LIB}")
  else ()
    message(STATUS "tcmalloc enabled but unavailable on this system")
  endif ()
endif ()

add_executable(evaluator app/evaluator.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(evaluator PRIVATE "-DMODE_EVALUATOR")
target_link_libraries(evaluator ${OpenMP_CXX_LIBRARIES})
install(TARGETS evaluator DESTINATION bin)

add_executable(node_separator app/node_separator_ml.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(node_separator PRIVATE "-DMODE_NODESEP")
target_link_libraries(node_separator ${OpenMP_CXX_LIBRARIES})
install(TARGETS node_separator DESTINATION bin)

add_executable(label_propagation app/label_propagation.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(label_propagation PRIVATE "-DMODE_LABELPROPAGATION")
target_link_libraries(label_propagation ${OpenMP_CXX_LIBRARIES})
install(TARGETS label_propagation DESTINATION bin)

add_executable(partition_to_vertex_separator app/partition_to_vertex_separator.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(partition_to_vertex_separator PRIVATE "-DMODE_PARTITIONTOVERTEXSEPARATOR")
target_link_libraries(partition_to_vertex_separator ${OpenMP_CXX_LIBRARIES})
install(TARGETS partition_to_vertex_separator DESTINATION bin)

add_executable(interface_test misc/example_library_call/interface_test.cpp interface/kaHIP_interface.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping> $<TARGET_OBJECTS:libnodeordering>)
target_include_directories(interface_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(interface_test PRIVATE "-DMODE_KAFFPA")
target_link_libraries(interface_test ${OpenMP_CXX_LIBRARIES})
if(LIB_METIS)
        target_link_libraries(interface_test ${LIB_METIS})
endif()
install(TARGETS interface_test DESTINATION bin)

if(NOT NOMPI)
  add_executable(kaffpaE app/kaffpaE.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping> $<TARGET_OBJECTS:libkaffpa_parallel>)
  target_compile_definitions(kaffpaE PRIVATE "-DMODE_KAFFPAE")
  target_include_directories(kaffpaE PUBLIC ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(kaffpaE ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_LIBRARIES} OpenMP::OpenMP_CXX )
  install(TARGETS kaffpaE DESTINATION bin)
endif()

add_executable(graphchecker app/graphchecker.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
target_compile_definitions(graphchecker PRIVATE "-DMODE_GRAPHCHECKER")
target_link_libraries(graphchecker ${OpenMP_CXX_LIBRARIES})
install(TARGETS graphchecker DESTINATION bin)

add_executable(edge_partitioning app/spac.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping> $<TARGET_OBJECTS:libspac>)
target_compile_definitions(edge_partitioning PRIVATE "-DMODE_KAFFPA")
target_link_libraries(edge_partitioning ${OpenMP_CXX_LIBRARIES})
install(TARGETS edge_partitioning DESTINATION bin)

add_executable(node_ordering app/node_ordering.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libnodeordering>)
target_compile_definitions(node_ordering PRIVATE "-DMODE_NODESEP -DMODE_NODEORDERING")
target_link_libraries(node_ordering ${OpenMP_CXX_LIBRARIES})
install(TARGETS node_ordering DESTINATION bin)

if(LIB_METIS)
  add_executable(fast_node_ordering app/fast_node_ordering.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libnodeordering>)
  target_compile_definitions(fast_node_ordering PRIVATE "-DMODE_NODESEP -DMODE_NODEORDERING -DFASTORDERING")
  target_link_libraries(fast_node_ordering ${OpenMP_CXX_LIBRARIES} ${LIB_METIS})
  install(TARGETS fast_node_ordering DESTINATION bin)
endif()


# Shared interface library
add_library(kahip SHARED interface/kaHIP_interface.cpp $<TARGET_OBJECTS:libkaffpa>
                                                           $<TARGET_OBJECTS:libmapping>
                                                           $<TARGET_OBJECTS:libnodeordering>)
target_include_directories(kahip PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(kahip PRIVATE "-DMODE_KAFFPA")
target_link_libraries(kahip PUBLIC ${OpenMP_CXX_LIBRARIES})
if(LIB_METIS)
        target_link_libraries(kahip PUBLIC ${LIB_METIS})
endif()
install(TARGETS kahip DESTINATION lib)

# Static interface library
add_library(kahip_static interface/kaHIP_interface.cpp $<TARGET_OBJECTS:libkaffpa>
                                                           $<TARGET_OBJECTS:libmapping>
                                                           $<TARGET_OBJECTS:libnodeordering>)
target_include_directories(kahip_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_compile_definitions(kahip_static PRIVATE "-DMODE_KAFFPA")
target_link_libraries(kahip_static PUBLIC ${OpenMP_CXX_LIBRARIES})
set_target_properties(kahip_static PROPERTIES PUBLIC_HEADER interface/kaHIP_interface.h)

if(LIB_METIS)
        target_link_libraries(kahip_static PUBLIC ${LIB_METIS})
endif()

install(TARGETS kahip_static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

# ParHIP
if(NOT NOMPI AND PARHIP)
  add_subdirectory(parallel/modified_kahip)
  add_subdirectory(parallel/parallel_src)
endif()

if(USE_ILP)
  find_package(Gurobi REQUIRED)
  MESSAGE("Using Gurobi for ILP solver in ilp_improve")
  add_executable(ilp_improve app/ilp_improve.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
  target_include_directories(ilp_improve PUBLIC ${GUROBI_INCLUDE_DIR})
  target_compile_definitions(ilp_improve PRIVATE "-DMODE_ILPIMPROVE")
  target_link_libraries(ilp_improve ${OpenMP_CXX_LIBRARIES} ${GUROBI_LIBRARIES})

  add_executable(ilp_exact app/ilp_exact.cpp $<TARGET_OBJECTS:libkaffpa> $<TARGET_OBJECTS:libmapping>)
  target_include_directories(ilp_exact PUBLIC ${GUROBI_INCLUDE_DIR})
  target_compile_definitions(ilp_exact PRIVATE "-DMODE_ILPIMPROVE -DMODE_ILPEXACT")
  target_link_libraries(ilp_exact ${OpenMP_CXX_LIBRARIES} ${GUROBI_LIBRARIES})

endif()

# pybind11 module
option(BUILDPYTHONMODULE "Build Python binding." OFF)
if (BUILDPYTHONMODULE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/misc/pymodule/pybind11)
  pybind11_add_module(kahip_python_binding ${CMAKE_CURRENT_SOURCE_DIR}/misc/pymodule/kahip.cpp)
  target_link_libraries(kahip_python_binding PUBLIC kahip_static)
  set_target_properties(kahip_python_binding PROPERTIES OUTPUT_NAME "kahip")
endif ()
