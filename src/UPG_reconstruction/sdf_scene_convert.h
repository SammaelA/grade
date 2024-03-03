#pragma once
#include "upg.h"
#include "sdf_scene.h"

namespace upg
{
  SdfScene create_sdf_scene(const UPGStructure &structure, const UPGParametersRaw &params);
}