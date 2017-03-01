# Minimum cmake version
cmake_minimum_required (VERSION 3.5.0)

# Name of project and that it is C++ only.
set(code_name "period_elastic")
project (${code_name} )

# ----------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------- 
# Search for zlib as a dependency for SeqAn.
set(CMAKE_MODULE_PATH ~/work/soft/seqan-library-2.2.0/share/cmake/Modules) 
set(CMAKE_PREFIX_PATH ~/work/soft/seqan-library-2.2.0)
find_package (ZLIB)
find_package (Boost REQUIRED COMPONENTS regex)
if (NOT Boost_FOUND )
   message (STATUS "Boost not found, not building.")
   return ()
endif ()
MESSAGE(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
# Would be nice to add an ifdef depending on the machie
execute_process(COMMAND hostname OUTPUT_VARIABLE computer_name)
find_package (SeqAn REQUIRED)
set (OpenMP_C_FLAGS /opt/local/lib/libomp )
set (OpenMP_CXX_FLAGS /opt/local/lib/libomp )
find_package (OpenMP REQUIRED)



# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

#set ( EIGEN3_INCLUDE_DIR "/opt/local/include/eigen3/Eigen" )
#if ( NOT EIGEN3_INCLUDE_DIR )
#    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
#endif()
# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).

# Do use -O3, makes the code way faster
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -I/opt/local/include/libomp -fopenmp ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -g -I/opt/local/include/libomp  -I/opt/local/include/eigen3 -fopenmp -Wall -Wextra -O3 ")

add_executable (${code_name} ${code_name}.cpp do_periodic.hpp nuc_elastic.hpp commandline_parse.hpp )
target_link_libraries (${code_name} ${SEQAN_LIBRARIES} ${OpenPM_LIBRARIES} ${Boost_LIBRARIES} )

