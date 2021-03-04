cmake_minimum_required(VERSION 3.7)

project(main)

find_package(SDL2 REQUIRED)
find_library(SDL2_IMAGE_LIBRARY SDL2_image)

message(SDL2_IMAGE_LIBRARY=${SDL2_IMAGE_LIBRARY})
find_library(SDL2_TTF_LIBRARY SDL2_ttf)

message(SDL2_TTF_LIBRARY=${SDL2_TTF_LIBRARY})
find_library(SDL2_MIXER_LIBRARY SDL2_mixer)

message(SDL2_MIXER_LIBRARY=${SDL2_MIXER_LIBRARY})
find_library(SDL2_MIXER_LIBRARY SDL2_mixer)

message(SDL2_MIXER_LIBRARY=${SDL2_MIXER_LIBRARY})
find_package(Boost REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
message(opgllib = ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}  ${CMAKE_DL_LIBS})
message(include =  ${OPENGL_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")
set(TINYLINK "-lX11 -lpthread -lSDL2 -lnoise -lSDL2_image -lSDL2_mixer -lSDL2_ttf -lGL -lGLEW -lboost_serialization -lboost_system -lboost_filesystem")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${TINYLINK}")
file(GLOB IMGUI_SRC src/imgui/*.cpp)
include_directories(${OPENGL_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")
file(GLOB_RECURSE C_SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.c")
add_executable(main ${SOURCES} ${C_SOURCES} ${IMGUI_SRC})
target_link_libraries(main ${OPENGL_LIBRARY} ${SDL2_LIBRARIES} ${SDL2_IMAGE_LIBRARY} ${GLEW_LIBRARIES}  ${CMAKE_DL_LIBS} 
                        ${SDL2_TTF_LIBRARY} ${SDL2_MIXER_LIBRARY} boost_serialization boost_system boost_filesystem)

