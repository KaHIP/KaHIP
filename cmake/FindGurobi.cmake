#### Adapted from http://www.openflipper.org/svnrepo/CoMISo/trunk/CoMISo/cmake/FindGUROBI.cmake

# - Try to find GUROBI
# Once done this will define
#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi

find_path(GUROBI_INCLUDE_DIR 
          NAMES gurobi_c++.h
          PATHS "$ENV{GUROBI_HOME}/include"
          )

find_library( GUROBI_LIBRARY 
              NAMES gurobi90
                    gurobi81
                    gurobi80
                    gurobi75
              PATHS "$ENV{GUROBI_HOME}/lib"
              )


find_library( GUROBI_CXX_LIBRARY 
              NAMES gurobi_g++5.2
              PATHS "$ENV{GUROBI_HOME}/lib" 
              )

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
set(GUROBI_LIBRARIES "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")

if (GUROBI_INCLUDE_DIRS AND GUROBI_LIBRARIES)
  set(GUROBI_FOUND TRUE)
endif()

mark_as_advanced(
  GUROBI_INCLUDE_DIRS
  GUROBI_LIBRARIES
)


