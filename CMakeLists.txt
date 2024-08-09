cmake_minimum_required(VERSION 3.20)
project(ALG2D)

set(CMAKE_CXX_STANDARD 17)

set(CL_ALL "en_US.UTF-8")

find_package(OpenMP REQUIRED)

if (OpenMP_C_FOUND)
    message(STATUS "OpenMP found")
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP_C_LIB_NAMES: ${OpenMP_C_LIB_NAMES}")
    message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_CXX_LIB_NAMES: ${OpenMP_CXX_LIB_NAMES}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else ()
    message(FATAL_ERROR "OpenMP not found")
endif ()


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfma")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfma")

# Define the AVX512 option
option(ENABLE_AVX512 "Enable AVX-512 support" ON)

# Check if the AVX512 option is enabled
if (ENABLE_AVX512)
    message(STATUS "AVX-512 support enabled.")
    if (MSVC)
        # For MSVC, use /arch:AVX512 (if supported)
        add_compile_options(/arch:AVX512)
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
            CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR
            CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # For GCC, Clang, or Intel compilers
        add_compile_options(-mavx512f -mavx512dq -mavx512pf -mavx512er)
    endif ()
else ()
    message(STATUS "AVX-512 support not enabled.")
endif ()

add_executable(ALG2D main.cpp Header/line.cpp Header/line.h Header/lineOthers.cpp Header/lineOthers.h Header/matrix.cpp Header/matrix.h Header/point.cpp Header/point.h Header/pointOthers.cpp Header/pointOthers.h Header/readFile.cpp Header/readFile.h Header/region.cpp Header/region.h Header/section.cpp Header/section.h Header/util.cpp Header/util.h Header/vector.cpp Header/vector.h)

target_link_libraries(ALG2D PRIVATE OpenMP::OpenMP_C)
target_link_libraries(ALG2D PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(ALG2D PRIVATE stdc++fs)

get_target_property(SOURCES ALG2D SOURCES)
get_target_property(INCLUDE_DIRS ALG2D INCLUDE_DIRECTORIES)
get_target_property(COMPILE_OPTIONS ALG2D COMPILE_OPTIONS)

message(STATUS "Target: ALG2D")
message(STATUS "  Sources: ${SOURCES}")
message(STATUS "  Include directories: ${INCLUDE_DIRS}")
message(STATUS "  Compile options: ${COMPILE_OPTIONS}")