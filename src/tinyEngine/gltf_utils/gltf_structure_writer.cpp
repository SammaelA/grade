#include "gltf_structure_writer.h"

namespace gltf
{
int accessorComponentTypeSizes[6] = {1, 1, 2, 2, 4, 4};
int AccessorTypeComponents[7] = {1,2,3,4,4,9,16};
std::string AccessorTypeNames[7] = {"SCALAR","VEC2","VEC3","VEC4","MAT2","MAT3","MAT4"};
std::string primitiveAttributeTypeNames[8] = {"POSITION","NORMAL","TANGENT","TEXCOORD_0","TEXCOORD_1",
                                              "COLOR_0","JOINTS_0","WEIGHTS_0"};
std::string cameraTypeNames[2]={"perspective","orthographic"};
}