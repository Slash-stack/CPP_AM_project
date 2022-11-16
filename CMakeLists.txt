cmake_minimum_required(VERSION 3.16.3)
project(TP_CPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall")

include_directories(include)
add_subdirectory(src)

add_executable(main demo/main.cpp)
target_link_libraries(main univers)

