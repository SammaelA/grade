cmake_minimum_required(VERSION 3.5)

project(nn_test)
set(CMAKE_CXX_STANDARD 17)
option(MODULE_VULKAN "include GPU support for neural networks" ON)

set(MODULES_SRC "")
set(MODULES_INCLUDE_DIRS "modules")
set(MODULES_LIBRARIES "")
set(MODULES_COMPILE_FLAGS "")

if (MODULE_VULKAN)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  set(SLICER_DIR "/home/sammael/kernel_slicer" CACHE STRING "Custom diff render location")
  
  find_package(OpenMP)
  find_package(Vulkan)
  find_package(FreeImage REQUIRED)

  set(MODULES_INCLUDE_DIRS ${MODULES_INCLUDE_DIRS} ${Vulkan_INCLUDE_DIR}
                                                   ${FREEIMAGE_INCLUDE_DIRS}
                                                   ${SLICER_DIR}/apps/vkutils
                                                   ${SLICER_DIR}/apps/volk
                                                   ${SLICER_DIR}/apps/LiteMath
                                                   ${SLICER_DIR}/apps/utils)

  link_directories(${SLICER_DIR}/apps/volk)

  add_compile_definitions(USE_GPU)
  add_compile_definitions(USE_VOLK)

  set(VKGEN_SOURCE tensor_processor_impl_gpu.cpp 
                  tensor_processor_impl_gpu_ds.cpp
                  tensor_processor_impl_gpu_init.cpp)

  set(VKUTILS_SOURCE ${SLICER_DIR}/apps/vkutils/vk_utils.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_copy.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_buffers.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_images.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_context.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_alloc_simple.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_pipeline.cpp
                    ${SLICER_DIR}/apps/vkutils/vk_descriptor_sets.cpp)   

  set(MODULES_SRC ${MODULES_SRC} ${VKGEN_SOURCE} ${VKUTILS_SOURCE})

  set(MODULES_LIBRARIES ${MODULES_LIBRARIES} ${Vulkan_LIBRARY} ${FREEIMAGE_LIBRARIES} OpenMP::OpenMP_CXX volk dl)
endif()

set(SRC
      direct/tensors.cpp
      direct/neural_network.cpp
      direct/siren.cpp
      direct/nnd_tests.cpp
      
      tensor_processor.cpp
      tensor_processor_impl.cpp
      tensor_compiler.cpp
      neural_network.cpp
      tensor_token.cpp
      siren.cpp
      nn_tests.cpp
      dataset.cpp
      
      standalone_tests.cpp)  

set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -g -O3 -W -Wno-attributes -Wno-sign-compare")
add_executable(nn_test ${MODULES_SRC} ${SRC})
target_include_directories(nn_test PUBLIC ${MODULES_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(nn_test ${MODULES_LIBRARIES})