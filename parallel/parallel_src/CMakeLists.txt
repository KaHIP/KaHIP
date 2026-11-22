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
  lib/io/parallel_graph_io_huge.cpp # last, so the code before it keeps its place
  lib/parallel_label_compress/parallel_label_compress_skip.cpp # last, so the code before it keeps its place
  lib/parallel_contraction_projection/parallel_contraction_send_lists.cpp) # last, so the code before it keeps its place
# argtable3 comes with kahip_parallel_static
add_library(libparallel OBJECT ${LIBPARALLEL_SOURCE_FILES})
target_include_directories(libparallel PUBLIC $<TARGET_PROPERTY:kahip_parallel_static,INTERFACE_INCLUDE_DIRECTORIES>)

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
  ${CMAKE_SOURCE_DIR}/extern/argtable3-3.2.2/argtable3.c)
add_library(libedgelist OBJECT ${LIBEDGELIST_SOURCE_FILES})

set(LIBDSPAC_SOURCE_FILES
  lib/dspac/dspac.cpp
  lib/dspac/edge_balanced_graph_io.cpp)
add_library(libdspac OBJECT ${LIBDSPAC_SOURCE_FILES})

add_executable(parhip app/parhip.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION")
target_link_libraries(parhip PRIVATE kahip_parallel_static)
install(TARGETS parhip DESTINATION bin)

add_executable(toolbox app/toolbox.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(toolbox PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DTOOLBOX")
target_link_libraries(toolbox PRIVATE kahip_parallel_static)
install(TARGETS toolbox DESTINATION bin)

add_executable(graph2binary app/graph2binary.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(graph2binary PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION" "-DGRAPH2DGF")
install(TARGETS graph2binary DESTINATION bin)

add_executable(graph2binary_external app/graph2binary_external.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(graph2binary_external PRIVATE "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION" "-DGRAPH2DGF")
install(TARGETS graph2binary_external DESTINATION bin)

add_executable(readbgf app/readbgf.cpp $<TARGET_OBJECTS:libgraph2bgf>)
target_compile_definitions(readbgf PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION")
install(TARGETS readbgf DESTINATION bin)

add_executable(edge_list_to_metis_graph app/edge_list_to_metis_graph.cpp $<TARGET_OBJECTS:libedgelist>)
target_compile_definitions(edge_list_to_metis_graph PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DKRONECKER_GENERATOR_PROGRAM")
install(TARGETS edge_list_to_metis_graph DESTINATION bin)

# checks that ordered_edge_map iterates exactly like the hashed_graph it replaces
add_executable(ordered_edge_map_test tests/ordered_edge_map_test.cpp)
add_test(NAME ordered_edge_map_test COMMAND ordered_edge_map_test)

#add_executable(friendster_list_to_metis_graph app/friendster_list_to_metis_graph.cpp $<TARGET_OBJECTS:libedgelist>)
#target_compile_definitions(friendster_list_to_metis_graph PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DKRONECKER_GENERATOR_PROGRAM")
#install(TARGETS friendster_list_to_metis_graph DESTINATION bin)

add_executable(dspac app/dspac.cpp $<TARGET_OBJECTS:libparallel> $<TARGET_OBJECTS:libdspac>)
target_compile_definitions(dspac PRIVATE "-DGRAPH_GENERATOR_MPI -DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION")
target_link_libraries(dspac PRIVATE kahip_parallel_static)
install(TARGETS dspac DESTINATION bin)

add_library(parhip_interface SHARED interface/parhip_interface.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip_interface PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION")
target_include_directories(parhip_interface PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_link_libraries(parhip_interface PRIVATE kahip_parallel_static)
set_target_properties(parhip_interface PROPERTIES PUBLIC_HEADER interface/parhip_interface.h)
install(TARGETS parhip_interface
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        )

# pkg-config
# Collect MPI flags for the .pc file
string(REPLACE ";" " -I" PARHIP_MPI_CFLAGS "-I${MPI_CXX_INCLUDE_DIRS}")
set(PARHIP_MPI_LIBS "")
foreach(_lib ${MPI_CXX_LIBRARIES})
  if(IS_ABSOLUTE "${_lib}")
    get_filename_component(_dir "${_lib}" DIRECTORY)
    get_filename_component(_name "${_lib}" NAME_WE)
    string(REGEX REPLACE "^lib" "" _name "${_name}")
    string(APPEND PARHIP_MPI_LIBS " -L${_dir} -l${_name}")
  else()
    string(APPEND PARHIP_MPI_LIBS " ${_lib}")
  endif()
endforeach()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/parhip_interface.pc.in" "${CMAKE_CURRENT_BINARY_DIR}/parhip_interface.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/parhip_interface.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")


add_library(parhip_interface_static interface/parhip_interface.cpp $<TARGET_OBJECTS:libparallel>)
target_compile_definitions(parhip_interface_static PRIVATE "-DGRAPH_GENERATOR_MPI" "-DGRAPHGEN_DISTRIBUTED_MEMORY" "-DPARALLEL_LABEL_COMPRESSION")
target_include_directories(parhip_interface_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_link_libraries(parhip_interface_static PUBLIC kahip_parallel_static)
install(TARGETS parhip_interface_static DESTINATION lib)

# Fortran interface: module parhip plus the communicator conversion it needs
if(BUILDFORTRAN)
  add_library(parhip_fortran STATIC interface/fortran/parhip.F90 interface/fortran/parhip_fortran.cpp)
  set_target_properties(parhip_fortran PROPERTIES Fortran_MODULE_DIRECTORY ${KAHIP_FORTRAN_MODULE_DIR})
  target_include_directories(parhip_fortran PUBLIC $<BUILD_INTERFACE:${KAHIP_FORTRAN_MODULE_DIR}>
                                                   $<INSTALL_INTERFACE:include>)
  target_link_libraries(parhip_fortran PUBLIC parhip_interface_static)
  install(TARGETS parhip_fortran ARCHIVE DESTINATION lib)
  install(FILES ${KAHIP_FORTRAN_MODULE_DIR}/parhip.mod DESTINATION include)

  # the example uses the Fortran MPI module, so it needs a Fortran MPI
  find_package(MPI COMPONENTS Fortran)
  if(MPI_Fortran_FOUND)
    add_executable(parhip_fortran_example ${CMAKE_SOURCE_DIR}/misc/example_library_call/fortran/parhip_example.f90)
    target_link_libraries(parhip_fortran_example PRIVATE parhip_fortran MPI::MPI_Fortran)
    add_test(NAME fortran_interface_parhip
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                     $<TARGET_FILE:parhip_fortran_example> ${MPIEXEC_POSTFLAGS})
  else()
    message(STATUS "No Fortran MPI found: parhip_fortran_example not built")
  endif()
endif()
