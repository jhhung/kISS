cmake_minimum_required (VERSION 3.16)

# Define project name
project(kISS)

# Define cpp standard
set(CMAKE_CXX_STANDARD 20)

# Include Biovoltron's CMakeLists.txt
include(${PROJECT_SOURCE_DIR}/submodules/Biovoltron/CMakeLists.txt)

# include boost
find_package( Boost 1.71.0 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

# Link your source code to your project
add_executable(kISS src/main.cpp)
target_compile_options(kISS PRIVATE -fopenmp -Ofast -g3)

target_link_libraries(kISS LINK_PUBLIC ${Boost_LIBRARIES} )
target_link_libraries(kISS PRIVATE biovoltron)
target_include_directories(kISS PRIVATE ${PROJECT_SOURCE_DIR}/include)
