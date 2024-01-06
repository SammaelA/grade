cmake_minimum_required(VERSION 3.5)

project(nn_test)
set(CMAKE_CXX_STANDARD 17)
include_directories(".")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})

find_package(OpenMP)
find_package(Vulkan)
find_package(FreeImage REQUIRED)
include_directories(${Vulkan_INCLUDE_DIR})
include_directories(${FREEIMAGE_INCLUDE_DIRS})
include_directories("/home/sammael/kernel_slicer/apps/vkutils"
                    "/home/sammael/kernel_slicer/apps/volk"
                    "/home/sammael/kernel_slicer/apps/LiteMath"
                    "/home/sammael/kernel_slicer/apps/utils")
link_directories("/home/sammael/kernel_slicer/apps/volk")

add_compile_definitions(USE_GPU)
add_compile_definitions(USE_VOLK)

set(VKGEN_SOURCE tensor_processor_impl_gpu.cpp 
                 tensor_processor_impl_gpu_ds.cpp
                 tensor_processor_impl_gpu_init.cpp)

set(VKUTILS_SOURCE /home/sammael/kernel_slicer/apps/vkutils/vk_utils.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_copy.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_buffers.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_images.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_context.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_alloc_simple.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_pipeline.cpp
                   /home/sammael/kernel_slicer/apps/vkutils/vk_descriptor_sets.cpp)   

set(SRC
    tensors.cpp
    neural_network.cpp
    siren.cpp
    nn_test.cpp
    tensor_processor.cpp
    tensor_processor_impl.cpp
    tensor_compiler.cpp
    neural_network_2.cpp
    tensor_token.cpp
    siren2.cpp
    )
set(CMAKE_CXX_FLAGS "-g -O3 -Wno-attributes -fopenmp")

add_executable(nn_test ${SRC} ${VKGEN_SOURCE} ${VKUTILS_SOURCE})
target_include_directories(nn_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(nn_test LINK_PUBLIC ${Vulkan_LIBRARY} ${FREEIMAGE_LIBRARIES} OpenMP::OpenMP_CXX volk dl)