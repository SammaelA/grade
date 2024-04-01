#!/bin/bash

#empty forlder for shaders
rm -r shaders_gpu
mkdir shaders_gpu

#copy shaders for neuralCore
cd modules/neuralCore/shaders_gpu
find -name "*.spv" | xargs cp --parents -t ../../../shaders_gpu
cd ../../..

#copy shaders for LiteRT
cd modules/LiteRT/Renderer/shaders_gpu
find -name "*.spv" | xargs cp --parents -t ../../../../shaders_gpu
cd ../../../..