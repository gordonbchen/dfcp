cmake_minimum_required(VERSION 3.20)

project(dfcp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

add_executable(impute
    src/me.cpp
    src/seq_array.cpp
    src/clusters.cpp
    src/max.cpp
    src/expect.cpp
    src/elbo.cpp
    src/tree.cpp
    src/math.cpp
    src/util.cpp
    src/pbwt.cpp
    src/fwd_bkwd.cpp
    src/impute_io.cpp
    src/io.cpp
)

target_include_directories(impute PRIVATE
    include
)

target_compile_options(impute PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -O3
)

target_link_libraries(impute PRIVATE
    OpenMP::OpenMP_CXX
)

add_executable(eval_impute
    src/eval_impute.cpp
    src/impute_io.cpp
    src/seq_array.cpp
    src/io.cpp
)

target_include_directories(eval_impute PRIVATE
    include
)

target_compile_options(eval_impute PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -O3
)
