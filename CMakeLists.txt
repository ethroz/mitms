cmake_minimum_required(VERSION 3.10)
project(mitms)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler flags
add_compile_options(-O0 -g -fpermissive -funroll-loops)

# Find required packages
find_package(LAPACK REQUIRED)

# Add include directories
include_directories(/usr/local/include/lapackpp/ /usr/include/lapackpp/)

# Add source files
set(SOURCES
    src/configs.cpp
    src/ring.cpp
    src/matrix.cpp
    src/gate.cpp
    src/circuit.cpp
    src/util.cpp
    src/database.cpp
    src/search.cpp
    src/main.cpp
)

# Create executable
add_executable(mitms ${SOURCES})

# Link libraries
target_link_libraries(mitms
    lapackpp
    rt
    pthread
)
