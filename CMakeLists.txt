cmake_minimum_required (VERSION 3.16)

# Define your project name
project(kISS)

# Include Biovoltron's CMakeLists.txt
include(${PROJECT_SOURCE_DIR}/submodules/Biovoltron/CMakeLists.txt)

# Link your source code to your project
add_executable(kISS src/main.cpp)
target_link_libraries(kISS biovoltron)
