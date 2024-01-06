#!/bin/sh
glslangValidator -V kernel1D_pow.comp -o kernel1D_pow.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_transpose.comp -o kernel2D_transpose.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_osum.comp -o kernel1D_osum.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_exp.comp -o kernel1D_exp.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_log.comp -o kernel1D_log.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_sum.comp -o kernel1D_sum.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_mul.comp -o kernel2D_mul.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_div.comp -o kernel2D_div.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_sub.comp -o kernel2D_sub.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_sin.comp -o kernel1D_sin.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_add.comp -o kernel2D_add.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_copy.comp -o kernel1D_copy.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_cos.comp -o kernel1D_cos.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_matmul_transposed.comp -o kernel2D_matmul_transposed.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_outer_p_add.comp -o kernel2D_outer_p_add.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel1D_fill.comp -o kernel1D_fill.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V kernel2D_outer_product.comp -o kernel2D_outer_product.comp.spv -DGLSL -I.. -I/home/sammael/kernel_slicer/TINYSTL -I/home/sammael/kernel_slicer/apps/LiteMathAux -I/home/sammael/kernel_slicer/apps/LiteMath 
glslangValidator -V z_memcpy.comp -o z_memcpy.comp.spv
