cmake_minimum_required(VERSION 3.5)

project(nn_test)
set(CMAKE_CXX_STANDARD 17)

set(SRC
    tensors.cpp
    neural_network.cpp
    siren.cpp
    nn_test.cpp
    tensor_processor.cpp
    )
set(CMAKE_CXX_FLAGS "-g -O3")

add_executable(nn_test ${SRC})
target_include_directories(nn_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})