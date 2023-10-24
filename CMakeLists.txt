cmake_minimum_required(VERSION 3.0.0)
project(AssemblyGraph)
set(CMAKE_CXX_FLAGS "-fopenmp -O3 -std=c++0x")
set(CMAKE_BINARY_DIR "bin")
# include directory
include_directories(src)
# handling boost library linking
find_package(Boost 1.46.0 COMPONENTS iostreams filesystem system program_options)
if(NOT Boost_FOUND)
    message(FATAL_ERROR " Cannot find Boost C++ Libraries version 1.46 or higher, abort.")
    return()
endif(NOT Boost_FOUND)

include_directories(${include_directories} ${Boost_INCLUDE_DIRS})
link_directories(${link_directories} ${Boost_LIB_DIRS})

# build Minimizer-Index
file(GLOB SEQINDEX "src/*.cpp")
add_executable(mm_index ${SEQINDEX})
target_link_libraries(mm_index ${Boost_LIBRARIES})
set_target_properties(mm_index PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
