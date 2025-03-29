cmake_minimum_required(VERSION 3.10)
project(mitms)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(LAPACK REQUIRED)

# Add include directories
include_directories(
    /usr/local/include/lapackpp/
    /usr/include/lapackpp/
    ${CMAKE_SOURCE_DIR}/src
)

# Common source files
set(COMMON_SOURCES
    src/configs.cpp
    src/ring.cpp
    src/matrix.cpp
    src/gate.cpp
    src/circuit.cpp
    src/util.cpp
    src/main.cpp
)

# Create object library for common sources
add_library(common_objects OBJECT ${COMMON_SOURCES})

# Original version sources
set(ORIGINAL_SOURCES
    src/database.cpp
    src/search.cpp
)

# Parallel version sources
set(PARALLEL_SOURCES
    src/parallel/database.cpp
    src/parallel/search.cpp
)

# Create executables
add_executable(mitms_original 
    $<TARGET_OBJECTS:common_objects>
    ${ORIGINAL_SOURCES}
)

add_executable(mitms_parallel 
    $<TARGET_OBJECTS:common_objects>
    ${PARALLEL_SOURCES}
)

# Link libraries for both executables
target_link_libraries(mitms_original
    lapackpp
    rt
    pthread
)

target_link_libraries(mitms_parallel
    lapackpp
    rt
    pthread
)
