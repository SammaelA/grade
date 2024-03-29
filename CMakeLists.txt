cmake_minimum_required(VERSION 3.7)

project(main)
message("Building project")

#some examples how to use this file
#cmake CMakeLists.txt -Wno-dev -DTREE_PROJECT=ON -DMODULE_ENZYME=ON  -DMODULE_DR=ON  -DCLANG_VERSION=14
#cmake CMakeLists.txt -Wno-dev -DTREE_PROJECT=ON -DMODULE_ENZYME=OFF -DMODULE_DR=OFF -DMODULE_VULKAN=ON
#cmake CMakeLists.txt -Wno-dev -DDIFF_RECONSTRUCTION=ON -DMODULE_ENZYME=ON -DUSE_CUSTOM_DIFF_RENDER=ON -DDR_SRC=/home/sammael/grade/modules/custom_diff_render
#a list of available options
option(USE_CUSTOM_DIFF_RENDER "Use custom differentiable renderer" OFF)
option(TREE_PROJECT "Build Tree Project" OFF)
option(DIFF_RECONSTRUCTION "Build project for differentiable reconstruction" OFF)
option(CITIES_GENERATOR "Build project for cities generator" OFF)

#list of available modules
option(MODULE_HYDRA "include Hydra module in project" OFF)
option(MODULE_ENZYME "use Enzyme automatic differentiation tool" OFF)
option(MODULE_DR "include Custom differentiable render in project" OFF)
option(MODULE_NN "include neural networks module in project" ON)
option(MODULE_VULKAN "include GPU support for neural networks" OFF)
option(MODULE_LITE_RT "include LiteRT module in project" ON)

#common dependencies for all projects
find_package(SDL2 REQUIRED)
find_library(SDL2_IMAGE_LIBRARY SDL2_image)
find_library(SDL2_TTF_LIBRARY SDL2_ttf)
find_library(SDL2_MIXER_LIBRARY SDL2_mixer)
find_package(Boost REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenMP REQUIRED)

include(FetchContent)
FetchContent_Declare(
  LiteMath
  URL      https://github.com/msu-graphics-group/LiteMath/archive/refs/heads/main.zip
  SOURCE_DIR dependencies/LiteMath
)

FetchContent_MakeAvailable(LiteMath)

#common compiler settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_THREAD_LIBS_INIT "-lpthread") 
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0) 
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -std=c++17 -g -O3 -pthread -Wall -Wno-unused-function \
                     -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result \
                     -Wno-int-to-pointer-cast -Wno-register -Wreturn-type -Wdisabled-optimization \
                     -Wno-narrowing -Wno-attributes -Wno-sign-compare")

#prepare common lists of options to use it when building targets
set(COMMON_LINKER_FLAGS "-lm -lpthread -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lGL -lGLEW -lboost_serialization -lboost_system -lboost_filesystem")
file(GLOB_RECURSE COMMON_SRC RELATIVE ${CMAKE_SOURCE_DIR} "third_party/*.cpp" "common_src/*.cpp" "third_party/*.c" "common_src/*.c")
set(COMMON_INCLUDE_DIRS "common_src", "third_party", "dependencies", "dependencies/SDL2", "dependencies/glm-0.9.9-a2")
set(COMMON_LIBRARIES ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}
                     ${SDL2_TTF_LIBRARY} ${SDL2_MIXER_LIBRARY} 
                     boost_serialization boost_system boost_filesystem
                     OpenMP::OpenMP_CXX)
message("COMMON LIBRARIES: " ${COMMON_LIBRARIES})


#Modules. If module is on, we include its source files, otherwise include stub file
set(MODULES_LINKER_FLAGS "")
set(MODULES_SRC "")
set(MODULES_INCLUDE_DIRS "modules")
set(MODULES_LIBRARIES "")
set(MODULES_COMPILE_FLAGS "")

#####################################################################################################################################################

if (MODULE_HYDRA)
  #integration with Hydra2 is broken and abandoned
  #Soon we will make an integration with Hydra3

  #find_library(HYDRA_LIB1 libcorto.a dependencies/Hydra/HydraAPI/bin)
  #message(HYDRA_LIB1=${HYDRA_LIB1})
  #find_library(HYDRA_LIB2 libhydra_api.a dependencies/Hydra/HydraAPI/bin)
  #message(HYDRA_LIB2=${HYDRA_LIB2})
  #find_library(HYDRA_LIB3 libies_parser.a dependencies/Hydra/HydraAPI/bin)
  #message(HYDRA_LIB3=${HYDRA_LIB3})
  #find_library(HYDRA_LIB4 libmikktspace.a dependencies/Hydra/HydraAPI/bin)
  #message(HYDRA_LIB4=${HYDRA_LIB4})

  #set(MODULES_LIBRARIES ${MODULES_LIBRARIES} glfw freeimage rt stdc++fs OpenCL ${HYDRA_LIB1} ${HYDRA_LIB2} ${HYDRA_LIB3} ${HYDRA_LIB4} dl)
  #set(MODULES_INCLUDE_DIRS ${MODULES_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} "dependencies/Hydra")
  #set(MODULES_SRC ${MODULES_SRC} modules/hydra/hydra_scene_exporter.cpp modules/hydra/export_scene_internal.cpp modules/hydra/init_gl.cpp)
  message("HYDRA MODULE IS TEMPORARY NOT SUPPORTED")
  set(MODULES_SRC ${MODULES_SRC} modules/hydra/stub/stub.cpp)
else()
  message("Module HYDRA                OFF")
  set(MODULES_SRC ${MODULES_SRC} modules/hydra/stub/stub.cpp)
endif()

#####################################################################################################################################################

if (MODULE_ENZYME)
  message("Module ENZYME               ON")
  set(CLANG_VERSION "14" CACHE STRING "clang and LLVM version to use")
  set(CMAKE_CXX_COMPILER "/usr/bin/clang++-${CLANG_VERSION}")
  set(MODULES_SRC ${MODULES_SRC} modules/autodiff/autodiff.cpp)
  set(MODULES_COMPILE_FLAGS ${MODULES_COMPILE_FLAGS} "-fplugin=dependencies/enzyme/ClangEnzyme-${CLANG_VERSION}.so")
else()
  message("Module ENZYME               OFF")
  set(MODULES_SRC ${MODULES_SRC} modules/autodiff/autodiff_stub.cpp)
endif()

#####################################################################################################################################################

if (MODULE_DR)
  message("Module CUSTOM DIFF RENDERER ON")
  if (MODULE_ENZYME)
  else()
    message("TRYING TO BUILD DIFF RENDERER WITHOUT AUTODIFF. PROGRAM MAY NOT WORK AS INTENDED")
  endif()

  set(DR_SRC "dependencies/diff_render" CACHE STRING "Custom diff render location")
  add_definitions(-DUSE_EMBREE)
  add_definitions(-DUSE_CUSTOM_DIFF_RENDER)

  find_library(EMBREE_LIB embree3 dependencies/embree)

  file(GLOB_RECURSE DR_SOURCE RELATIVE ${CMAKE_SOURCE_DIR} "modules/custom_diff_render/raytrace_src/*.cpp" 
                                                           "modules/custom_diff_render/src/*.cpp" 
                                                           "modules/custom_diff_render/common_src/*.cpp")
  set(MODULES_SRC ${MODULES_SRC} ${DR_SOURCE} modules/custom_diff_render/halfgpu_diff_render.cpp)
  set(MODULES_LIBRARIES ${MODULES_LIBRARIES} ${EMBREE_LIB} tbb)
  set(MODULES_INCLUDE_DIRS ${MODULES_INCLUDE_DIRS} "modules/custom_diff_render/src", 
                                                   "modules/custom_diff_render/raytrace_src", 
                                                   "modules/custom_diff_render/common_src", 
                                                   "dependencies/embree/include")
else()
  message("Module CUSTOM DIFF RENDERER OFF")
  set(MODULES_SRC ${MODULES_SRC} modules/custom_diff_render/diff_render_stub.cpp)
endif()

#####################################################################################################################################################

#if (MODULE_NN)
  message("Module NEURALCORE           ON")
  set(MODULES_SRC ${MODULES_SRC} modules/neuralCore/tensor_processor.cpp
                                 modules/neuralCore/tensor_processor_impl.cpp
                                 modules/neuralCore/tensor_compiler.cpp
                                 modules/neuralCore/neural_network.cpp
                                 modules/neuralCore/tensor_token.cpp
                                 modules/neuralCore/siren.cpp
                                 modules/neuralCore/nn_tests.cpp
                                 modules/neuralCore/nn_tests_performance.cpp
                                 modules/neuralCore/dataset.cpp)

  if (MODULE_VULKAN)
    set(MODULES_SRC ${MODULES_SRC} modules/neuralCore/tensor_processor_impl_gpu.cpp 
                                   modules/neuralCore/tensor_processor_impl_gpu_ds.cpp
                                   modules/neuralCore/tensor_processor_impl_gpu_init.cpp)
  endif()
#else()
#endif()

#####################################################################################################################################################

if (MODULE_VULKAN)
  message("Module VULKAN               ON")
  list(APPEND CMAKE_MODULE_PATH "modules/neuralCore/")

  find_package(Vulkan)
  FetchContent_Declare(
    volk
    URL      https://github.com/zeux/volk/archive/refs/heads/master.zip
    SOURCE_DIR dependencies/volk
  )
  FetchContent_Declare(
    vkutils
    URL      https://github.com/msu-graphics-group/vk-utils/archive/refs/heads/master.zip
    SOURCE_DIR dependencies/vkutils
  )
  FetchContent_MakeAvailable(volk vkutils)

  set(MODULES_INCLUDE_DIRS ${MODULES_INCLUDE_DIRS} ${Vulkan_INCLUDE_DIR}
                                                   ${CMAKE_SOURCE_DIR}/dependencies/vkutils
                                                   ${CMAKE_SOURCE_DIR}/dependencies/volk)

  link_directories(${CMAKE_SOURCE_DIR}/dependencies/volk)

  add_compile_definitions(USE_GPU)
  add_compile_definitions(USE_VOLK)

  set(VKUTILS_SOURCE ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_utils.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_copy.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_buffers.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_images.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_context.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_alloc_simple.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_pipeline.cpp
                     ${CMAKE_SOURCE_DIR}/dependencies/vkutils/vk_descriptor_sets.cpp)   

  set(MODULES_SRC ${MODULES_SRC} ${VKUTILS_SOURCE})
  set(MODULES_LIBRARIES ${MODULES_LIBRARIES} ${Vulkan_LIBRARY} OpenMP::OpenMP_CXX volk dl)
else()
  message("Module VULKAN               OFF")
endif()

#####################################################################################################################################################

#if (MODULE_LITE_RT)
message("Module MODULE_LITE_RT       ON")

FetchContent_Declare(
  LiteScene
  URL      https://github.com/msu-graphics-group/LiteScene/archive/refs/heads/main.zip
  SOURCE_DIR dependencies/LiteScene
)
FetchContent_MakeAvailable(LiteScene)

link_directories(${CMAKE_SOURCE_DIR}/modules/LiteRT/dependencies/embree/lib)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/LiteMath)
include_directories(${CMAKE_SOURCE_DIR}/modules/LiteRT/dependencies/embree)

set(SCENE_LOADER_SRC
    ${CMAKE_SOURCE_DIR}/dependencies/LiteScene/pugixml.cpp
    ${CMAKE_SOURCE_DIR}/dependencies/LiteScene/hydraxml.cpp
    ${CMAKE_SOURCE_DIR}/dependencies/LiteScene/cmesh4.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/sdfScene/sdf_scene.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/sdfScene/gpuReady/sdf_scene.cpp
    )

set(CORE_SRC
    ${CMAKE_SOURCE_DIR}/dependencies/LiteMath/Image2d.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/Timer.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BruteForceRT.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/Renderer/eye_ray.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/Renderer/eye_ray_host.cpp
    )

set(BUILDERS_SRC
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BVH/BVH2Common.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BVH/BVH2Common_host.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BVH/cbvh.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BVH/cbvh_fat.cpp
    ${CMAKE_SOURCE_DIR}/modules/LiteRT/BVH/cbvh_embree2.cpp
    )

set(VULKAN_SRC "")

if (MODULE_VULKAN)
set(VULKAN_SRC 
${CMAKE_SOURCE_DIR}/modules/LiteRT/Renderer/eye_ray_gpu.cpp 
${CMAKE_SOURCE_DIR}/modules/LiteRT/Renderer/eye_ray_gpu_ds.cpp 
${CMAKE_SOURCE_DIR}/modules/LiteRT/Renderer/eye_ray_gpu_init.cpp)
endif()

set(MODULES_SRC ${MODULES_SRC} ${SCENE_LOADER_SRC} ${CORE_SRC} ${BUILDERS_SRC} ${VULKAN_SRC})

set(MODULES_LIBRARIES ${MODULES_LIBRARIES} embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tasking)
#else()
#endif()

#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################
                 
if (TREE_PROJECT)
  file(GLOB_RECURSE TP_SRC RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c" "src/*.cpp")  
  add_executable(main ${COMMON_SRC} ${MODULES_SRC} ${TP_SRC})
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COMMON_LINKER_FLAGS} ${MODULES_LINKER_FLAGS}")
  target_include_directories(main PUBLIC ${COMMON_INCLUDE_DIRS} ${MODULES_INCLUDE_DIRS} "src")
  target_link_libraries(main 
                        "${COMMON_LIBRARIES}"
                        "${CMAKE_DL_LIBS}"
                        "${MODULES_LIBRARIES}")
  target_compile_options(main PUBLIC ${MODULES_COMPILE_FLAGS})
endif()

#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################

if (CITIES_GENERATOR)

add_library(libnoise SHARED IMPORTED)
set_target_properties(libnoise PROPERTIES IMPORTED_LOCATION "/usr/lib/libnoise.so")

file(GLOB_RECURSE CG_SRC RELATIVE ${CMAKE_SOURCE_DIR} "cities_generator/*.c" "cities_generator/*.cpp")
add_executable(cg ${COMMON_SRC} ${CG_SRC})
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COMMON_LINKER_FLAGS} -lX11 -lnoise")
target_include_directories(cg PUBLIC ${COMMON_INCLUDE_DIRS} "cities_generator")

target_link_libraries(cg 
                      "${COMMON_LIBRARIES}"
                      glfw
                      libnoise)
endif()

#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################
#####################################################################################################################################################

if (DIFF_RECONSTRUCTION)

  #find_package(OpenCV REQUIRED PATHS dependencies/opencv/build)

  set(RT_EMBREE_LIBS "")
  set(OMP_LIBS "")
  set(DR_SOURCE_RT )
  set(DR_SOURCE )
  if(USE_CUSTOM_DIFF_RENDER)
    set(DR_SRC "/home/sammael/diff_render/diffrender_tutorials" CACHE STRING "Custom diff render location")
    find_package(OpenMP REQUIRED)
    set(AUX_OMP_INC "/usr/lib/llvm-9/lib/clang/9.0.0/include/") # use 'locate omp.h'     on you system to find apropriate dir for clang
    set(AUX_OMP_LIB "/usr/lib/llvm-9/lib/")                     # use 'locate libomp.so' on you system to find apropriate dir for clang

    include_directories("${DR_SRC}/src", "${DR_SRC}/raytrace_src", "${DR_SRC}/common_src", "${DR_SRC}/external/enzyme")

    add_definitions(-DUSE_EMBREE)
    add_definitions(-DUSE_CUSTOM_DIFF_RENDER)
    include_directories(${DR_SRC}/external)
    include_directories(${DR_SRC}/external/embree/include)

    find_package(Threads REQUIRED)
    find_library(EMBREE_LIB embree3 dependencies/embree)
    link_directories(${DR_SRC}/external/embree/lib)
    set(RT_EMBREE_LIBS ${EMBREE_LIB} tbb)

    set(DR_SOURCE_RT
    "${DR_SRC}/raytrace_src/raytrace3d_bf.cpp"
    "${DR_SRC}/raytrace_src/raytrace3d_embree.cpp"
    "${DR_SRC}/raytrace_src/EmbreeRT.cpp"
    "${CMAKE_SOURCE_DIR}/dependencies/LiteMath/Image2d.cpp")

    #set(OMP_LIBS "")
    #if(OpenMP_CXX_FOUND)
      add_definitions(-DUSE_OMP)
      include_directories(${AUX_OMP_INC})
      link_directories(${AUX_OMP_LIB})
      set(OMP_LIBS OpenMP::OpenMP_CXX)
    #endif()

    file(GLOB_RECURSE DR_SOURCE RELATIVE ${CMAKE_SOURCE_DIR} "${DR_SRC}/src/*.cpp" "${DR_SRC}/common_src/*.cpp")

    include_directories("${DR_SRC}/mitsuba_src/stub")
    message(SOURCE ${DR_SRC})
  endif()

  file(GLOB_RECURSE DRec_SRC RELATIVE ${CMAKE_SOURCE_DIR} "diff_reconstruction/*.c" "diff_reconstruction/*.cpp")
  add_executable(dr ${DRec_SRC} ${COMMON_SRC} ${DR_SOURCE} ${DR_SOURCE_RT})
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COMMON_LINKER_FLAGS}")
  target_include_directories(dr PUBLIC ${COMMON_INCLUDE_DIRS} "diff_reconstruction" 
                            "/usr/include/python3.10" "dependencies/python" 
                            ${OpenCV_INCLUDE_DIRS})

  target_link_libraries(dr 
                         "${COMMON_LIBRARIES}"
                          /usr/lib/x86_64-linux-gnu/libpython3.10.so
                          ${OpenCV_LIBS}
                          ${OMP_LIBS} 
                          ${RT_EMBREE_LIBS})

  if(USE_CUSTOM_DIFF_RENDER)
    set(ENZYME_FLAGS "-fplugin=dependencies/enzyme/ClangEnzyme-${CLANG_VERSION}.so")
    message(${ENZYME_FLAGS})
    target_compile_options(dr PUBLIC ${ENZYME_FLAGS})
  else()
  endif()
endif()
