cmake_minimum_required (VERSION 3.16)

# your project name
project(kiss)

set(CMAKE_CXX_STANDARD 20)
set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++20 -Wall -Weffc++ -Wno-error=effc++ -pedantic \
  -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-newline-eof  \
  -pedantic-errors -Wextra -Waggregate-return -Wcast-align \
  -Wcast-qual -Wconversion \
  -Wdisabled-optimization \
  -Werror -Wfloat-equal -Wformat=2 \
  -Wformat-nonliteral -Wformat-security  \
  -Wformat-y2k \
  -Wimplicit  -Wimport  -Winit-self  -Winline -Winvalid-pch   \
  -Wlong-long \
  -Wmissing-field-initializers -Wmissing-format-attribute   \
  -Wmissing-include-dirs -Wmissing-noreturn \
  -Wpacked  -Wpadded -Wpointer-arith \
  -Wredundant-decls \
  -Wshadow -Wstack-protector \
  -Wstrict-aliasing=2 -Wswitch-default \
  -Wswitch-enum \
  -Wunreachable-code -Wunused \
  -Wunused-parameter \
  -Wvariadic-macros \
  -Wwrite-strings"
)

# Include openmp
find_package(OpenMP 4.5 REQUIRED)

# Include tbb
set(TBB_ROOT "${CMAKE_SOURCE_DIR}/submodules/tbb")
include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbb_preview)

# build
add_definitions(-DLIBSAIS_OPENMP)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/libsais/)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/parallel-divsufsort/)

# link your source code to your project
add_executable(kiss-1 src/kiss-1.cpp)
target_compile_options(kiss-1 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_include_directories(kiss-1 PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(kiss-1 PRIVATE spdlog OpenMP::OpenMP_CXX ${TBB_IMPORTED_TARGETS})
target_include_directories(kiss-1 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(kiss-2 src/kiss-2.cpp)
target_compile_options(kiss-2 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_include_directories(kiss-2 PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(kiss-2 PRIVATE spdlog OpenMP::OpenMP_CXX ${TBB_IMPORTED_TARGETS})
target_include_directories(kiss-2 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(saca-k src/saca-k.cpp src/saca-k-source.cpp)
target_compile_options(saca-k PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_include_directories(saca-k PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(saca-k PRIVATE spdlog OpenMP::OpenMP_CXX)

add_executable(pardss-64 src/pardss.cpp)
target_compile_options(pardss-64 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_include_directories(pardss-64 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(pardss-64 PRIVATE divsufsort spdlog OpenMP::OpenMP_CXX libprange)
target_include_directories(pardss-64 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/submodules/parallel-divsufsort/include")

add_executable(pardss-32 src/pardss_32bit.cpp)
target_compile_options(pardss-32 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_include_directories(pardss-32 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(pardss-32 PRIVATE divsufsort spdlog OpenMP::OpenMP_CXX libprange)
target_include_directories(pardss-32 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/submodules/parallel-divsufsort/include")

add_executable(libsais-test-64 src/libsais.cpp)
target_compile_options(libsais-test-64 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_compile_definitions(libsais-test-64 PRIVATE LIBSAIS_OPENMP)
target_include_directories(libsais-test-64 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(libsais-test-64 PRIVATE libsais spdlog OpenMP::OpenMP_CXX)
target_include_directories(libsais-test-64 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/submodules/libsais/include")
target_include_directories(libsais-test-64 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../libsais_test/Biovoltron/include")

add_executable(libsais-test-32 src/libsais_32bit.cpp)
target_compile_options(libsais-test-32 PRIVATE -fopenmp -Ofast -g3 -std=c++20 -ltcmalloc)
target_compile_definitions(libsais-test-32 PRIVATE LIBSAIS_OPENMP)
target_include_directories(libsais-test-32 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(libsais-test-32 PRIVATE libsais spdlog OpenMP::OpenMP_CXX)
target_include_directories(libsais-test-32 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/submodules/libsais/include")
target_include_directories(libsais-test-32 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../../libsais_test/Biovoltron/include")

