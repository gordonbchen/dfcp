cmake_minimum_required(VERSION 3.20)

project(dfcp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

add_executable(dfcp
    src/me.cpp
    src/clusters.cpp
    src/soft_clusters.cpp
    src/max.cpp
    src/expect.cpp
    src/elbo.cpp
    src/tree.cpp
    src/math.cpp
    src/util.cpp
)

target_include_directories(dfcp PRIVATE
    include
)

target_compile_options(dfcp PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -O3
)

target_link_libraries(dfcp PRIVATE
    OpenMP::OpenMP_CXX
)

