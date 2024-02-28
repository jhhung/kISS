cmake_minimum_required (VERSION 3.16)

# Define project name
project(kISS)

# Define cpp standard
set(CMAKE_CXX_STANDARD 20)

# Include boost
find_package( Boost 1.71.0 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# Include tbb
set(TBB_ROOT "${CMAKE_CURRENT_LIST_DIR}/submodules/tbb")
include(${TBB_ROOT}/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${TBB_ROOT} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1)
find_package(TBB REQUIRED tbb_preview)

# Include Threads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Include openmp
find_package(OpenMP 4.5 REQUIRED)

# Include spdlog
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/submodules/spdlog)

# Add definiton
add_definitions(-mavx -mavx2 -mlzcnt)

# Link your source code to your project
add_executable(kISS src/main.cpp)
target_compile_options(kISS PRIVATE -fopenmp -Ofast -g)

target_link_libraries(kISS PRIVATE ${TBB_IMPORTED_TARGETS})
target_link_libraries(kISS PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(kISS PRIVATE Threads::Threads)
target_link_libraries(kISS PRIVATE spdlog)
target_link_libraries(kISS LINK_PUBLIC ${Boost_LIBRARIES} )
target_include_directories(kISS PRIVATE ${PROJECT_SOURCE_DIR}/include)

include(CTest)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/tests)