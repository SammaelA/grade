
# Getting started

## Install basic version (no autodiff, diff. render etc.):

Make and Cmake \
`sudo apt-get install make cmake`

Basic graphics libraries \
`sudo apt-get install mesa-common-dev libglu1-mesa-dev mesa-utils`

OpenGL libraries \
`sudo apt-get install libglfw3 libglfw3-dev libglew-dev`

GLM (vectors, matrices) \
`sudo apt-get install libglm-dev`

SDL2 (windows manager) \
`sudo apt-get install libsdl2-dev libsdl2-ttf-dev libsdl2-mixer-dev libsdl2-image-dev`

Boost \
`sudo apt-get install libboost-system-dev libboost-filesystem-dev libboost-serialization-dev`

## Build

Initialize git submodules \
`git submodules init` \
`git submodules update --remote`

Cmake with options what to build \
`cmake CMakeLists.txt -Wno-dev -DTREE_PROJECT=ON`

## List of modules
`MODULE_NN` - neural networks on CPU (no external dependencies) \
`MODULE_SDF_SCENE` - interface to interact with SDF scene (no external dependencies) \
`MODULE_VULKAN` - neural networks on GPU with Vulkan (requires MODULE_NN) \
`MODULE_ENZYME` - Enzyme autodiff \
`MODULE_DR` - Differentiable renderer (requires MODULE_ENZYME) \
`MODULE_HYDRA` - DEPRECATED 
