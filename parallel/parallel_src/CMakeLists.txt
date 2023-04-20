if(NOT OPTIMIZED_OUTPUT)
  add_definitions("-DNOOUTPUT")
endif()

if(DETERMINISTIC_PARHIP)
  add_definitions("-DDETERMINISTIC_PARHIP")
endif()

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)
include_directories(${MPI_CXX_INCLUDE_PATH})
link_libraries(OpenMP::OpenMP_CXX MPI::MPI_CXX)

set(LIBPARALLEL_SOURCE_FILES
  lib/data_structure/parallel_graph_access.cpp
  lib/data_structure/balance_management.cpp
  lib/data_structure/balance_management_refinement.cpp
  lib/data_structure/balance_management_coarsening.cpp
  lib/parallel_label_compress/node_ordering.cpp
  lib/parallel_contraction_projection/parallel_contraction.cpp
  lib/parallel_contraction_projection/parallel_block_down_propagation.cpp
  lib/parallel_contraction_projection/parallel_projection.cpp
  lib/distributed_partitioning/distributed_partitioner.cpp
  lib/distributed_partitioning/initial_partitioning/initial_partitioning.cpp
  lib/distributed_partitioning/initial_partitioning/distributed_evolutionary_partitioning.cpp
  lib/distributed_partitioning/initial_partitioning/random_initial_partitioning.cpp
  lib/communication/mpi_tools.cpp
  lib/communication/dummy_operations.cpp
  lib/io/parallel_graph_io.cpp
  lib/io/parallel_vector_io.cpp
  lib/tools/random_functions.cpp
  lib/tools/distributed_quality_metrics.cpp
  extern/argtable3-3.2.2/argtable3.c)
add_library(libparallel OBJECT ${LIBPARALLEL_SOURCE_FILES})
target_include_directories(libparallel PUBLIC $<TARGET_PROPERTY:libmodified_kahip_interface,INTERFACE_INCLUDE_DIRECTORIES>)

set(LIBGRAPH2BGF_SOURCE_FILES
  lib/data_structure/parallel_graph_access.cpp
  lib/io/parallel_graph_io.cpp
  lib/data_structure/balance_management.cpp
  lib/data_structure/balance_management_refinement.cpp
  lib/data_structure/balance_management_coarsening.cpp)
add_library(libgraph2bgf OBJECT ${LIBGRAPH2BGF_SOURCE_FILES})

set(LIBEDGELIST_SOURCE_FILES
  lib/data_structure/parallel_graph_access.cpp
  lib/io/parallel_graph_io.cpp
  lib/data_structure/balance_management.cpp
  lib/data_structure/balance_management_refinement.cpp
  lib/data_structure/balance_management_coarsening.cpp
  extern/argtable3-3.2.2/argtable3.c)
add_library(libedgelist OBJECT ${LIBEDGELIST_SOURCE_FILES})

set(LIBDSPAC_SOURCE_FILES
  lib/dspac/dspac.cpp
  lib/dspac/edge_balanced_graph_io.cpp)
add_library(libdspac OBJECT ${LIBDSPAC_SOURCE_FILES})

add_executable(parhip app/parhip.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION")
target_link_libraries(parhip PRIVATE libmodified_kahip_interface)
install(TARGETS parhip DESTINATION bin)

add_executable(toolbox app/toolbox.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(toolbox PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DTOOLBOX")
target_link_libraries(toolbox PRIVATE libmodified_kahip_interface)
install(TARGETS toolbox DESTINATION bin)

add_executable(graph2binary app/graph2binary.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(graph2binary PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION -DGRAPH2DGF")
install(TARGETS graph2binary DESTINATION bin)

add_executable(graph2binary_external app/graph2binary_external.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(graph2binary_external PRIVATE "-DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION -DGRAPH2DGF")
install(TARGETS graph2binary_external DESTINATION bin)

add_executable(readbgf app/readbgf.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(readbgf PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION")
install(TARGETS readbgf DESTINATION bin)

add_executable(edge_list_to_metis_graph app/edge_list_to_metis_graph.cpp $<TARGET_OBJECTS:libedgelist>)
target_compile_definitions(edge_list_to_metis_graph PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DKRONECKER_GENERATOR_PROGRAM")
target_link_libraries(edge_list_to_metis_graph PRIVATE libmodified_kahip_interface)
install(TARGETS edge_list_to_metis_graph DESTINATION bin)

#add_executable(friendster_list_to_metis_graph app/friendster_list_to_metis_graph.cpp $<TARGET_OBJECTS:libedgelist>)
#target_compile_definitions(friendster_list_to_metis_graph PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DKRONECKER_GENERATOR_PROGRAM")
#target_link_libraries(edge_list_to_metis_graph PRIVATE libmodified_kahip_interface)
#install(TARGETS friendster_list_to_metis_graph DESTINATION bin)

add_executable(dspac app/dspac.cpp $<TARGET_OBJECTS:libparallel> $<TARGET_OBJECTS:libdspac>)
target_compile_definitions(dspac PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION")
target_link_libraries(dspac PRIVATE libmodified_kahip_interface)
install(TARGETS dspac DESTINATION bin)

add_library(parhip_interface SHARED interface/parhip_interface.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip_interface PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION")
target_include_directories(parhip_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_link_libraries(parhip_interface PRIVATE libmodified_kahip_interface)
set_target_properties(parhip_interface PROPERTIES PUBLIC_HEADER interface/parhip_interface.h)
install(TARGETS parhip_interface
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

add_library(parhip_interface_static interface/parhip_interface.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip_interface_static PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY -DPARALLEL_LABEL_COMPRESSION")
target_include_directories(parhip_interface_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_link_libraries(parhip_interface_static PRIVATE libmodified_kahip_interface)
install(TARGETS parhip_interface_static DESTINATION lib)
