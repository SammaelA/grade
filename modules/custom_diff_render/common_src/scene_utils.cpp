#include "scene.h"
namespace diff_render
{
void Scene::restore_meshes(bool restore_normals, bool restore_tangents, bool transform_to_unindexed_mesh)
{
  for (auto &mesh : meshes)
  {
    if (restore_normals && mesh.normals.size() != mesh.vertices.size())
    {
      mesh.normals = ::std::vector<float3>(mesh.vertices.size(), float3(0,0,0));
      ::std::vector<int> vert_mul(mesh.vertices.size(), 0);

      for (int i=0;i<mesh.indices.size(); i+=3)
      {
        float3 &v1 = mesh.vertices[mesh.indices[i]];
        float3 &v2 = mesh.vertices[mesh.indices[i+1]];
        float3 &v3 = mesh.vertices[mesh.indices[i+2]];

        float3 l1 = v2-v1;
        float3 l2 = v3-v1;
        float3 n = float3(1,0,0);
        if (length(l1) < 1e-6 || length(l2) < 1e-6)
        {
          logerr("Scene::restore_meshes triangle[%d %d %d] has near-zero size. It may lead to errors",
                 mesh.indices[i], mesh.indices[i+1], mesh.indices[i+2]);
          n = normalize(cross(l1, l2));
        }
        mesh.normals[mesh.indices[i]] += n;
        mesh.normals[mesh.indices[i+1]] += n;
        mesh.normals[mesh.indices[i+2]] += n;

        vert_mul[mesh.indices[i]] += 1;
        vert_mul[mesh.indices[i+1]] += 1;
        vert_mul[mesh.indices[i+2]] += 1;
      }

      for (int i=0; i<mesh.normals.size(); i++)
        mesh.normals[i] /= vert_mul[i];
    }

    if (restore_tangents)
    {

    }

    if (transform_to_unindexed_mesh)
    {
      if (mesh.colors.size() != mesh.vertices.size())
        mesh.colors = ::std::vector<float3>(mesh.vertices.size(), float3(0,0,0));
      if (mesh.tc.size() != mesh.vertices.size())
        mesh.tc = ::std::vector<float2>(mesh.vertices.size(), float2(0,0));
      if (mesh.normals.size() != mesh.vertices.size())
        mesh.normals = ::std::vector<float3>(mesh.vertices.size(), float3(1,0,0));
      if (mesh.tangents.size() != mesh.vertices.size())
        mesh.tangents = ::std::vector<float3>(mesh.vertices.size(), float3(0,1,0));
    
      int v_count = mesh.indices.size();
      auto vertices = ::std::vector<float3>(v_count, float3(0,0,0));
      auto colors = ::std::vector<float3>(v_count, float3(0,0,0));
      auto tc = ::std::vector<float2>(v_count, float2(0,0));
      auto normals = ::std::vector<float3>(v_count, float3(0,0,0));
      auto tangents = ::std::vector<float3>(v_count, float3(0,0,0));
      auto indices = ::std::vector<unsigned int>(v_count, 0);

      int i = 0;
      for (int ind : mesh.indices)
      {
        vertices[i] = mesh.vertices[ind];
        colors[i] = mesh.colors[ind];
        tc[i] = mesh.tc[ind];
        normals[i] = mesh.normals[ind];
        tangents[i] = mesh.tangents[ind];
        indices[i] = i;
        i++;
      }

      mesh.vertices = vertices;
      mesh.colors = colors;
      mesh.tc = tc;
      mesh.normals = normals;
      mesh.tangents = tangents;
      mesh.indices = indices;
    }
  }
}

void transform(TriangleMesh &mesh, const LiteMath::float4x4 &transform)
{
  LiteMath::float4x4 n_tr = LiteMath::transpose(LiteMath::inverse4x4(transform));
  for(auto& v : mesh.vertices)
    v = LiteMath::mul(transform, v);
  for(auto& n : mesh.normals)
    n = LiteMath::to_float3(LiteMath::mul(transform, float4(n.x, n.y, n.z, 0)));
  for(auto& n : mesh.tangents)
    n = LiteMath::to_float3(LiteMath::mul(transform, float4(n.x, n.y, n.z, 0)));
}

void Scene::prepare_for_render() const
{
  if (prepared)
    return;
  prepared = true;
  if (meshes.size() == 0)
    return;

  assert(meshes.size() == restricted_transforms.size());
  assert(meshes.size() == transforms.size());
  assert(meshes.size() == mesh_instancing_types.size());

  transforms_inv.reserve(transforms.size());
  for (int i = 0; i < meshes.size(); i++)
  {
    if (mesh_instancing_types[i] == RESTRICTED_TRANSFORM)
    {
      transforms[i] = ::std::vector<float4x4>(restricted_transforms[i].size());
      for (int j=0;j<restricted_transforms[i].size();j++)
        transforms[i][j] = restricted_transforms[i][j].to_mat();
    }
    transforms_inv.push_back(transforms[i]);
    for (int j=0;j<transforms[i].size();j++)
      transforms_inv.back()[j] = LiteMath::inverse4x4(transforms[i][j]);
  }

  int total_vertices = 0;
  for (int i=0;i<meshes.size();i++)
    total_vertices += transforms[i].size()*meshes[i].vertices.size();

  preparedData.vertices.resize(total_vertices, float3(0,0,0));
  preparedData.orig_vertices.resize(total_vertices, float3(0,0,0));
  preparedData.colors.resize(total_vertices, float3(0,0,0));
  preparedData.tc.resize(total_vertices, float2(0,0));
  preparedData.normals.resize(total_vertices, float3(1,0,0));
  preparedData.tangents.resize(total_vertices, float3(0,1,0));
  
  unsigned offset = 0;
  preparedData.indices.resize(meshes.size());
  for (int i=0;i<meshes.size();i++)
  {
    preparedData.indices[i] = ::std::vector<::std::vector<unsigned>>{transforms[i].size(), ::std::vector<unsigned>()};
    for (int j=0;j<transforms[i].size();j++)
    {
      float4x4 tr = transforms[i][j];
      float4x4 n_tr = LiteMath::transpose(LiteMath::inverse4x4(tr));
      preparedData.indices[i][j] = ::std::vector<unsigned>(meshes[i].indices.size(), (unsigned)-1);
      for (int k = 0; k < meshes[i].indices.size(); k++)
        preparedData.indices[i][j][k] = offset + meshes[i].indices[k];
      
      for (int k = 0; k < meshes[i].vertices.size(); k++)
      {
        preparedData.vertices[offset + k] = tr*meshes[i].vertices[k];
        preparedData.orig_vertices[offset + k] = meshes[i].vertices[k];
      }
      
      assert(meshes[i].colors.size() == meshes[i].vertices.size() || meshes[i].colors.size() == 0);
      for (int k = 0; k < meshes[i].colors.size(); k++)
        preparedData.colors[offset + k] = meshes[i].colors[k];

      assert(meshes[i].tc.size() == meshes[i].vertices.size() || meshes[i].tc.size() == 0);
      for (int k = 0; k < meshes[i].tc.size(); k++)
        preparedData.tc[offset + k] = meshes[i].tc[k];

      assert(meshes[i].normals.size() == meshes[i].vertices.size() || meshes[i].normals.size() == 0);
      for (int k = 0; k < meshes[i].normals.size(); k++)
        preparedData.normals[offset + k] = n_tr*meshes[i].normals[k];

      assert(meshes[i].tangents.size() == meshes[i].vertices.size() || meshes[i].tangents.size() == 0);
      for (int k = 0; k < meshes[i].tangents.size(); k++)
        preparedData.tangents[offset + k] = n_tr*meshes[i].tangents[k];
      
      offset += meshes[i].vertices.size();
    }
  }
}

void Scene::get_prepared_mesh(TriangleMesh &mesh) const
{
  int total_indices = 0;
  for (auto &ind_m : preparedData.indices)
    for (auto &ind_i : ind_m)
      total_indices += ind_i.size();
  mesh.indices.resize(total_indices);
  int ind = 0;
  for (auto &ind_m : preparedData.indices)
    for (auto &ind_i : ind_m)
      for (auto &i : ind_i)
        mesh.indices[ind++] = i;
  mesh.vertices = preparedData.vertices;
  mesh.colors = preparedData.colors;
  mesh.tc = preparedData.tc;
  mesh.normals = preparedData.normals;
  mesh.tangents = preparedData.tangents;
  mesh.textures = meshes[0].textures;
}
} // namespace diff_render
