cmake_minimum_required(VERSION 3.18)
project(ALG2D)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 17)

add_executable(ALG2D main.cpp Header/util.cpp Header/util.h Header/vector.cpp Header/vector.h Header/matrix.cpp Header/matrix.h Header/line.cpp Header/line.h Header/point.cpp Header/point.h Header/pointOthers.cpp Header/pointOthers.h Header/lineOthers.cpp Header/lineOthers.h Header/region.cpp Header/region.h Header/section.cpp Header/section.h Header/readFile.cpp Header/readFile.h)

target_compile_options(ALG2D PUBLIC -xCOMMON-AVX512)