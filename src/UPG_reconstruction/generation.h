#pragma once

struct Block;
struct ComplexModel;
struct UniversalGenMesh;
namespace upg
{
  bool create_model_from_block(Block &bl, ComplexModel &mod);
  void mesh_to_complex_model(const UniversalGenMesh &mesh, ComplexModel &mod);
}