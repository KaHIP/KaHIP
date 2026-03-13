cmake_minimum_required(VERSION 3.10)
project(parhip_example CXX)

find_package(MPI REQUIRED)

# Adjust this path to where KaHIP was installed or built
set(KAHIP_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../build" CACHE PATH "Path to KaHIP build directory")

add_executable(parhip_test parhip_test.cpp)
target_include_directories(parhip_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../parallel/parallel_src/interface
    ${MPI_CXX_INCLUDE_DIRS}
)
target_link_libraries(parhip_test PRIVATE
    ${KAHIP_BUILD_DIR}/parallel/parallel_src/libparhip_interface.so
    ${MPI_CXX_LIBRARIES}
)
