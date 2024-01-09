#include <cstdio>
#include <memory>
#include <algorithm>
#include <cstdint>
#include <cassert>
#include <vector>
#include <functional>
#include <array>
#include <type_traits>
#include <string>
#include <cstring>
#include <chrono>
#include <cmath>

#include "tensor_processor.h"
#include "tensor_compiler.h"
#include "neural_network.h"
#include "siren.h"
#include "nn_tests.h"
#include "direct/nnd_tests.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(int argc, char **argv)
{
  //nnd::perform_tests();
  nn::perform_tests();
  return 0;
}