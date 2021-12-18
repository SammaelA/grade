cmake_minimum_required(VERSION 3.7)

project(main)

find_package(SDL2 REQUIRED)
find_library(SDL2_IMAGE_LIBRARY SDL2_image)
message(SDL2_IMAGE_LIBRARY=${SDL2_IMAGE_LIBRARY})

find_library(SDL2_TTF_LIBRARY SDL2_ttf)
message(SDL2_TTF_LIBRARY=${SDL2_TTF_LIBRARY})

find_library(SDL2_MIXER_LIBRARY SDL2_mixer)
message(SDL2_MIXER_LIBRARY=${SDL2_MIXER_LIBRARY})

find_library(PY_CLUSERING_LIBRARY libpyclustering.so dependencies/pyclustering/ccore)
message(PY_CLUSERING_LIBRARY=${PY_CLUSERING_LIBRARY})

find_library(HYDRA_LIB1 libcorto.a dependencies/Hydra/HydraAPI/bin)
message(HYDRA_LIB1=${HYDRA_LIB1})

find_library(HYDRA_LIB2 libhydra_api.a dependencies/Hydra/HydraAPI/bin)
message(HYDRA_LIB2=${HYDRA_LIB2})

find_library(HYDRA_LIB3 libies_parser.a dependencies/Hydra/HydraAPI/bin)
message(HYDRA_LIB3=${HYDRA_LIB3})

find_library(HYDRA_LIB4 libmikktspace.a dependencies/Hydra/HydraAPI/bin)
message(HYDRA_LIB4=${HYDRA_LIB4})

find_package(Boost REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(PythonLibs REQUIRED)

string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)

message(opgllib = ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}  ${CMAKE_DL_LIBS})
message(include =  ${OPENGL_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_THREAD_LIBS_INIT "-lpthread") 
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0) 
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(CMAKE_CXX_FLAGS "-g -O3 -W -Wno-sign-compare -Wno-int-to-pointer-cast -pthread -pg")
set(TINYLINK "-lX11 -lpthread -lSDL2 -lnoise -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lGL -lGLEW -lboost_serialization -lboost_system -lboost_filesystem")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${TINYLINK}")
file(GLOB IMGUI_SRC src/imgui/*.cpp)
include_directories("dependencies/pyclustering/ccore/include" "dependencies/SDL2" "dependencies/glm-0.9.9-a2" 
                    ${PYTHON_INCLUDE_DIRS} "dependencies/python" "src" "dependencies/Hydra")
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE C_SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c")
add_executable(main ${SOURCES} ${C_SOURCES} ${IMGUI_SRC})
target_link_libraries(main ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}  ${CMAKE_DL_LIBS} 
                        ${SDL2_TTF_LIBRARY} ${SDL2_MIXER_LIBRARY} ${PY_CLUSERING_LIBRARY} 
                        boost_serialization boost_system boost_filesystem
                        ${PYTHON_LIBRARIES} 
                        ${HYDRA_LIB1}
                        ${HYDRA_LIB2}
                        ${HYDRA_LIB3}
                        ${HYDRA_LIB4}
                        freeimage
                        glfw
                        dl
                        rt
                        stdc++fs
                        OpenCL)

