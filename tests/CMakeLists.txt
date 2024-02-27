cmake_minimum_required(VERSION 3.16)
project(kiss-test)

file(GLOB_RECURSE SOURCE_FILES *.cpp)
add_executable(kiss-test ${SOURCE_FILES})

target_link_libraries(kiss-test PRIVATE ${TBB_IMPORTED_TARGETS})
target_link_libraries(kiss-test PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(kiss-test PRIVATE Threads::Threads)
target_link_libraries(kiss-test PRIVATE spdlog)
target_link_libraries(kiss-test LINK_PUBLIC ${Boost_LIBRARIES} )
target_include_directories(kiss-test PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_test(NAME kiss COMMAND kiss-test)