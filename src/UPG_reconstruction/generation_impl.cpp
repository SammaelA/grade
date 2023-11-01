#include "generation_impl.h"

namespace upg
{
  UniversalGenMesh UniversalGenInstance::generate(std::span <const float> parameters)
  {
    generator.take_params(parameters);
    return generator.generate();
  }
}