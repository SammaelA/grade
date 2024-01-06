#!/bin/bash
cd /home/sammael/kernel_slicer
/home/sammael/kernel_slicer/kslicer /home/sammael/grade/modules/neuralCore/tensor_processor_impl.cpp -mainClass TensorProcessorImpl -stdlibfolder /home/sammael/kernel_slicer/TINYSTL -pattern ipv -reorderLoops YX -I/home/sammael/kernel_slicer/apps/LiteMath ignore -I/home/sammael/kernel_slicer/apps/LiteMathAux ignore -I/home/sammael/kernel_slicer/TINYSTL ignore -shaderCC glsl -suffix _GPU -DKERNEL_SLICER -v 
cd /home/sammael/grade/modules/neuralCore