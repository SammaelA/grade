#pragma once
void __enzyme_fwddiff(void*, ...);
extern int enzyme_dupnoneed;
extern int enzyme_dup;
extern float __enzyme_v[1024];

//expected 
//void F(my_float *in, my_float *out) - free function
//int input size, int output_size, float *input, float *output, float *jac
#define ENZYME_EVALUATE_WITH_DIFF(F, input_size, output_size, input, output, jac) \
{ \
    for(int i = 0; i < input_size-1; i++)\
    {\
      __enzyme_v[i] = 1;\
      __enzyme_fwddiff((void*)(F), input, __enzyme_v, enzyme_dupnoneed, output, jac + i*output_size);\
      __enzyme_v[i] = 0;\
    }\
    __enzyme_v[input_size - 1] = 1;\
    __enzyme_fwddiff((void*)(F), input, __enzyme_v, enzyme_dup, output, jac + (input_size - 1)*output_size);\
    __enzyme_v[input_size - 1] = 0;\
}
