#!/bin/bash

#empty forlder for shaders
rm -r shaders_gpu
mkdir shaders_gpu

#build shaders for neuralCore and copy them
cd modules/neuralCore/shaders_gpu
bash build.sh
find -name "*.spv" | xargs cp --parents -t ../../../shaders_gpu
cd ../../..

#build shaders for LiteRT and copy them
cd modules/LiteRT/Renderer/shaders_gpu
bash build.sh
find -name "*.spv" | xargs cp --parents -t ../../../../shaders_gpu
cd ../../../..