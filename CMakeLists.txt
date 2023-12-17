cmake_minimum_required(VERSION 3.7)

project(main)

#some examples how to use this file
#cmake CMakeLists.txt -Wno-dev -DTREE_PROJECT=ON -DDEBUG=OFF -DMODULE_HYDRA=ON -DMODULE_DR=ON
#cmake CMakeLists.txt -Wno-dev -DDIFF_RECONSTRUCTION=ON -DUSE_CUSTOM_DIFF_RENDER=OFF -DDEBUG=OFF
#cmake CMakeLists.txt -Wno-dev -DDIFF_RECONSTRUCTION=ON -DUSE_CUSTOM_DIFF_RENDER=ON -DTREE_PROJECT=ON -DDEBUG=OFF -DCLANG_VERSION=9
#cmake CMakeLists.txt -Wno-dev -DTREE_PROJECT=ON -DDEBUG=OFF -DMODULE_HYDRA=ON -DMODULE_DR=ON -DDIFF_RECONSTRUCTION=OFF -DAUX_OMP_INC=/usr/include/omp.h -DAUX_OMP_LIB=/usr/lib/x86_64-linux-gnu/libomp.so -DCLANG_VERSION=14

#a list of available options
option(USE_CUSTOM_DIFF_RENDER "Use custom differentiable renderer" OFF)
option(TREE_PROJECT "Build Tree Project" OFF)
option(DIFF_RECONSTRUCTION "Build project for differentiable reconstruction" OFF)
option(CITIES_GENERATOR "Build project for cities generator" OFF)
option(DEBUG "debug mode" OFF)

#list of available modules
option(MODULE_HYDRA "include Hydra module in project" OFF)
option(MODULE_DR "include Custom differentiable render in project" OFF)

#set the compiler and version of C++ standart
set(CMAKE_CXX_STANDARD 17)
if(USE_CUSTOM_DIFF_RENDER AND DIFF_RECONSTRUCTION)
  set(CLANG_VERSION "9" CACHE STRING "clang and LLVM version to use")
  set(CMAKE_CXX_COMPILER "/usr/bin/clang++-${CLANG_VERSION}")
else()
set(CMAKE_CXX_COMPILER)
endif()

#common dependencies for all projects
find_package(SDL2 REQUIRED)
find_library(SDL2_IMAGE_LIBRARY SDL2_image)
message(SDL2_IMAGE_LIBRARY=${SDL2_IMAGE_LIBRARY})
find_library(SDL2_TTF_LIBRARY SDL2_ttf)
message(SDL2_TTF_LIBRARY=${SDL2_TTF_LIBRARY})
find_library(SDL2_MIXER_LIBRARY SDL2_mixer)
message(SDL2_MIXER_LIBRARY=${SDL2_MIXER_LIBRARY})
find_package(Boost REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
message(opgllib = ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}  ${CMAKE_DL_LIBS})
message(include =  ${OPENGL_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})

#common compiler settings
set(CMAKE_THREAD_LIBS_INIT "-lpthread") 
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0) 
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -g -O3 -W -pthread -Wno-sign-compare -Wno-int-to-pointer-cast -Wno-register -Wreturn-type -Wno-narrowing -Wno-unused-result -Wno-unused-parameter")
if (DEBUG)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wshadow")
endif()

#prepare common lists of options to use it when building targets
set(COMMON_LINKER_FLAGS "-lm -lpthread -lSDL2 -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lGL -lGLEW -lboost_serialization -lboost_system -lboost_filesystem")
file(GLOB_RECURSE COMMON_SRC RELATIVE ${CMAKE_SOURCE_DIR} "third_party/*.cpp" "common_src/*.cpp" "third_party/*.c" "common_src/*.c")
set(COMMON_INCLUDE_DIRS "common_src", "third_party", "dependencies/SDL2", "dependencies/glm-0.9.9-a2")
set(COMMON_LIBRARIES ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}
                     ${SDL2_TTF_LIBRARY} ${SDL2_MIXER_LIBRARY} 
                     boost_serialization boost_system boost_filesystem)


#Modules. If module is on, we include its source files, otherwise include stub file
set(MODULES_LINKER_FLAGS "")
set(MODULES_SRC "")
set(MODULES_INCLUDE_DIRS "modules")
set(MODULES_LIBRARIES "")
set(MODULES_COMPILE_FLAGS "")

if (MODULE_HYDRA)
  find_library(HYDRA_LIB1 libcorto.a dependencies/Hydra/HydraAPI/bin)
  message(HYDRA_LIB1=${HYDRA_LIB1})
  find_library(HYDRA_LIB2 libhydra_api.a dependencies/Hydra/HydraAPI/bin)
  message(HYDRA_LIB2=${HYDRA_LIB2})
  find_library(HYDRA_LIB3 libies_parser.a dependencies/Hydra/HydraAPI/bin)
  message(HYDRA_LIB3=${HYDRA_LIB3})
  find_library(HYDRA_LIB4 libmikktspace.a dependencies/Hydra/HydraAPI/bin)
  message(HYDRA_LIB4=${HYDRA_LIB4})

  set(MODULES_LIBRARIES ${MODULES_LIBRARIES} glfw freeimage rt stdc++fs OpenCL ${HYDRA_LIB1} ${HYDRA_LIB2} ${HYDRA_LIB3} ${HYDRA_LIB4} dl)
  set(MODULES_INCLUDE_DIRS ${MODULES_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} "dependencies/Hydra")
  set(MODULES_SRC ${MODULES_SRC} modules/hydra/hydra_scene_exporter.cpp modules/hydra/export_scene_internal.cpp modules/hydra/init_gl.cpp)
else()
  set(MODULES_SRC ${MODULES_SRC} modules/hydra/stub/stub.cpp)
endif()
if (MODULE_DR)
  set(CLANG_VERSION "9" CACHE STRING "clang and LLVM version to use")
  set(CMAKE_CXX_COMPILER "/usr/bin/clang++-${CLANG_VERSION}")

  set(DR_SRC "/home/sammael/diff_render/diffrender_tutorials" CACHE STRING "Custom diff render location")
  find_package(OpenMP REQUIRED)
  set(AUX_OMP_INC "/usr/lib/llvm-9/lib/clang/9.0.0/include/") # use 'locate omp.h'     on you system to find apropriate dir for clang
  set(AUX_OMP_LIB "/usr/lib/llvm-9/lib/")                     # use 'locate libomp.so' on you system to find apropriate dir for clang

  include_directories("${DR_SRC}/src", "${DR_SRC}/LiteMath", "${DR_SRC}/raytrace_src", "${DR_SRC}/common_src", "${DR_SRC}/external/enzyme", "${DR_SRC}/mitsuba_src/stub")

  add_definitions(-DUSE_EMBREE)
  add_definitions(-DUSE_CUSTOM_DIFF_RENDER)
  include_directories(${DR_SRC}/external)
  include_directories(${DR_SRC}/external/embree/include)

  find_package(Threads REQUIRED)
  link_directories(${DR_SRC}/external/embree/lib)
  set(RT_EMBREE_LIBS embree3 tbb)

  set(DR_SOURCE_RT
    "${DR_SRC}/raytrace_src/raytrace3d_bf.cpp"
    "${DR_SRC}/raytrace_src/raytrace3d_embree.cpp"
    "${DR_SRC}/raytrace_src/EmbreeRT.cpp")

  set(OMP_LIBS "")
  if(OpenMP_CXX_FOUND)
    add_definitions(-DUSE_OMP)
    include_directories(${AUX_OMP_INC})
    link_directories(${AUX_OMP_LIB})
    set(OMP_LIBS OpenMP::OpenMP_CXX)
  endif()

  file(GLOB_RECURSE DR_SOURCE RELATIVE ${CMAKE_SOURCE_DIR} "${DR_SRC}/src/*.cpp" "${DR_SRC}/common_src/*.cpp")
  set(MODULES_SRC ${MODULES_SRC} ${DR_SOURCE} ${DR_SOURCE_RT} modules/custom_diff_render/halfgpu_diff_render.cpp modules/custom_diff_render/autodiff.cpp)
  set(MODULES_LIBRARIES ${MODULES_LIBRARIES} ${OMP_LIBS} ${RT_EMBREE_LIBS})
  set(MODULES_COMPILE_FLAGS ${MODULES_COMPILE_FLAGS} "-fplugin=${DR_SRC}/external/enzyme/enzyme/build/Enzyme/ClangEnzyme-${CLANG_VERSION}.so")
  message(${MODULES_COMPILE_FLAGS})
else()
set(MODULES_SRC ${MODULES_SRC} modules/custom_diff_render/diff_render_stub.cpp modules/custom_diff_render/autodiff_stub.cpp)
endif()

#target Tree Project. It is mostly abandoned as of 2023                    
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

if (DIFF_RECONSTRUCTION)

  find_package(OpenCV REQUIRED PATHS dependencies/opencv/build)

  set(RT_EMBREE_LIBS "")
  set(OMP_LIBS "")
  set(DR_SOURCE_RT )
  set(DR_SOURCE )
  if(USE_CUSTOM_DIFF_RENDER)
    set(DR_SRC "/home/sammael/diff_render/diffrender_tutorials" CACHE STRING "Custom diff render location")
    find_package(OpenMP REQUIRED)
    set(AUX_OMP_INC "/usr/lib/llvm-9/lib/clang/9.0.0/include/") # use 'locate omp.h'     on you system to find apropriate dir for clang
    set(AUX_OMP_LIB "/usr/lib/llvm-9/lib/")                     # use 'locate libomp.so' on you system to find apropriate dir for clang

    include_directories("${DR_SRC}/src", "${DR_SRC}/LiteMath", "${DR_SRC}/raytrace_src", "${DR_SRC}/common_src", "${DR_SRC}/external/enzyme")

    add_definitions(-DUSE_EMBREE)
    add_definitions(-DUSE_CUSTOM_DIFF_RENDER)
    include_directories(${DR_SRC}/external)
    include_directories(${DR_SRC}/external/embree/include)

    find_package(Threads REQUIRED)
    link_directories(${DR_SRC}/external/embree/lib)
    set(RT_EMBREE_LIBS embree3 tbb)

    set(DR_SOURCE_RT
    "${DR_SRC}/raytrace_src/raytrace3d_bf.cpp"
    "${DR_SRC}/raytrace_src/raytrace3d_embree.cpp"
    "${DR_SRC}/raytrace_src/EmbreeRT.cpp")

    set(OMP_LIBS "")
    if(OpenMP_CXX_FOUND)
      add_definitions(-DUSE_OMP)
      include_directories(${AUX_OMP_INC})
      link_directories(${AUX_OMP_LIB})
      set(OMP_LIBS OpenMP::OpenMP_CXX)
    endif()

    file(GLOB_RECURSE DR_SOURCE RELATIVE ${CMAKE_SOURCE_DIR} "${DR_SRC}/src/*.cpp" "${DR_SRC}/common_src/*.cpp")

    include_directories("${DR_SRC}/mitsuba_src/stub")
    message(SOURCE ${DR_SRC})
  endif()

  file(GLOB_RECURSE DRec_SRC RELATIVE ${CMAKE_SOURCE_DIR} "diff_reconstruction/*.c" "diff_reconstruction/*.cpp")
  add_executable(dr ${DRec_SRC} ${COMMON_SRC} ${DR_SOURCE} ${DR_SOURCE_RT})
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COMMON_LINKER_FLAGS}")
  target_include_directories(dr PUBLIC ${COMMON_INCLUDE_DIRS} "diff_reconstruction" 
                            "/usr/include/python3.8" "dependencies/python" 
                            ${OpenCV_INCLUDE_DIRS})

  target_link_libraries(dr 
                         "${COMMON_LIBRARIES}"
                          /usr/lib/x86_64-linux-gnu/libpython3.8.so
                          ${OpenCV_LIBS}
                          ${OMP_LIBS} 
                          ${RT_EMBREE_LIBS})

  if(USE_CUSTOM_DIFF_RENDER)
    set(ENZYME_FLAGS "-fplugin=${DR_SRC}/external/enzyme/enzyme/build/Enzyme/ClangEnzyme-${CLANG_VERSION}.so")
    message(${ENZYME_FLAGS})
    target_compile_options(dr PUBLIC ${ENZYME_FLAGS})
  else()
  endif()
endif()
