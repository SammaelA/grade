#include "modeling.h"
#include "tinyEngine/engine.h"
#include "third_party/obj_loader.h"

namespace model_loader
{
Block obj_models_blk;
bool obj_models_blk_loaded = false;
std::string base_path = "resources/models";

Model *create_model_from_block(Block &bl, Texture &tex)
{
    std::string name = bl.get_string("name", "debug_box");
    return create_model_by_name(name, tex);
}

Model *create_model_by_name(std::string name, Texture &tex)
{
    if (name == "debug_box")
    {
        tex = engine::textureManager->get("noise");
        return create_debug_box_model();
    }
    else if (name == "simple_grass")
    {
        return create_simple_grass_model();
    }
    else 
    {
        return load_model_from_obj(name, tex);
    }
}

Model *create_debug_box_model()
{
    Box b = Box(glm::vec3(0,0,0), glm::vec3(1,0,0), glm::vec3(0,1,0), glm::vec3(0,0,1));
    Model *m = new Model;
    visualizer::box_to_model(&b, m);
    return m;
}

Model *create_simple_grass_model()
{
    std::vector<float> vertexes = {-0.5,0,0, -0.5,1,0, 0.5,0,0, 0.5,1,0,   0,0,-0.5, 0,1,-0.5, 0,0,0.5, 0,1,0.5};
    std::vector<float> tc = {0,1,0,0, 0,0,0,0, 1,1,0,0, 1,0,0,0, 0,1,0,0, 0,0,0,0, 1,1,0,0, 1,0,0,0};
    std::vector<float> normals = {0,0,1, 0,0,1, 0,0,1, 0,0,1, 1,0,0, 1,0,0, 1,0,0, 1,0,0};
    std::vector<GLuint> indices = {0, 1, 3, 2, 0, 3, 4,5,7, 6,4,7};

    Model *m = new Model;
    m->positions = vertexes;
    m->normals = normals;
    m->colors = tc;
    m->indices = indices;
    
    return m;
}

void transform_model_to_standart_form(Model *m)
{
    glm::vec3 min_pos = glm::vec3(1e9,1e9,1e9);
    glm::vec3 max_pos = glm::vec3(-1e9,-1e9,-1e9);
    for (int i=0;i<m->positions.size();i+=3)
    {
        min_pos = min(min_pos, glm::vec3(m->positions[i], m->positions[i+1], m->positions[i+2]));
        max_pos = max(max_pos, glm::vec3(m->positions[i], m->positions[i+1], m->positions[i+2]));
    }

    glm::vec3 size = max_pos - min_pos;
    float sz = (MAX(size.x, MAX(size.y,size.z)));
    for (int i=0;i<m->positions.size();i+=3)
    {
        m->positions[i] = (m->positions[i] - min_pos.x)/sz;
        m->positions[i+1] = (m->positions[i+1] - min_pos.y)/sz;
        m->positions[i+2] = (m->positions[i+2] - min_pos.z)/sz;
    }
}

void load_default_blk()
{
    
    load_block_from_file("models.blk", obj_models_blk);
    obj_models_blk_loaded = true;
}

Model *load_model_from_obj(std::string name, Texture &tex)
{
    if (!obj_models_blk_loaded)
    {
        load_default_blk();
    }

    Block *obj = obj_models_blk.get_block(name);
    if (!obj)
    {
        logerr("cannot find model %s. It is not mentioned in models.blk file", name.c_str());
        return nullptr;
    }
    std::string folder_name = obj->get_string("folder_name", name);
    std::string obj_filename = base_path + "/" + folder_name + "/" + obj->get_string("obj", "");
    std::string obj_color_tex = base_path + "/" + folder_name + "/" + obj->get_string("color", "");

    bool success = engine::textureManager->load_tex_to_catalog(name + "_tex", obj_color_tex);
    if (!success)
    {
        logerr("texture manager cannot load file %s", obj_color_tex.c_str());
        return nullptr;
    }
    tex = engine::textureManager->get(name + "_tex");

    objl::Loader loader;
    success = loader.LoadFile(obj_filename);
    if (!success)
    {
        logerr("obj loader cannot load file %s", obj_filename.c_str());
        return nullptr;
    }

    Model *m = new Model;
    int start_index = 0;
    for (auto &mesh : loader.LoadedMeshes)
    {
        logerr("create mesh from %d", loader.LoadedMeshes.size());
        for (auto &lv : mesh.Vertices)
        {
            m->positions.push_back(lv.Position.X);
            m->positions.push_back(lv.Position.Y);
            m->positions.push_back(lv.Position.Z);

            m->normals.push_back(lv.Normal.X);
            m->normals.push_back(lv.Normal.Y);
            m->normals.push_back(lv.Normal.Z);

            m->colors.push_back(lv.TextureCoordinate.X);
            m->colors.push_back(1 - lv.TextureCoordinate.Y);
            m->colors.push_back(0);
            m->colors.push_back(1);
        }

        for (auto &ind : mesh.Indices)
        {
            m->indices.push_back(start_index + ind);
        }
        start_index += mesh.Vertices.size();
    }
    transform_model_to_standart_form(m);
    return m;
}

void transform_model_to_standart_form(ComplexModel &cm)
{
  glm::vec3 min_pos = glm::vec3(1e9,1e9,1e9);
  glm::vec3 max_pos = glm::vec3(-1e9,-1e9,-1e9);
  for (Model *m : cm.models)
  {
    for (int i=0;i<m->positions.size();i+=3)
    {
        min_pos = min(min_pos, glm::vec3(m->positions[i], m->positions[i+1], m->positions[i+2]));
        max_pos = max(max_pos, glm::vec3(m->positions[i], m->positions[i+1], m->positions[i+2]));
    }
  }
  glm::vec3 size = max_pos - min_pos;
  float sz = (MAX(size.x, MAX(size.y,size.z)));
  for (Model *m : cm.models)
  {
    for (int i=0;i<m->positions.size();i+=3)
    {
        m->positions[i] = (m->positions[i] - min_pos.x)/sz;
        m->positions[i+1] = (m->positions[i+1] - min_pos.y)/sz;
        m->positions[i+2] = (m->positions[i+2] - min_pos.z)/sz;
    }
  }
}

bool create_model_from_block(Block &bl, ComplexModel &mod)
{
  std::string name = bl.get_string("name", "");
  if (!obj_models_blk_loaded)
  {
    load_default_blk();
  }

  Block *obj = obj_models_blk.get_block(name);
  if (!obj)
  {
    logerr("cannot find model %s. It is not mentioned in models.blk file", name.c_str());
    return false;
  }
  std::string folder_name = obj->get_string("folder_name", name);
  std::string obj_filename = base_path + "/" + folder_name + "/" + obj->get_string("obj", "");
  std::string obj_color_tex = base_path + "/" + folder_name + "/" + obj->get_string("color", "");

  objl::Loader loader;
  bool success = loader.LoadFile(obj_filename);
  if (!success)
  {
    logerr("obj loader cannot load file %s", obj_filename.c_str());
    return false;
  }

  std::map<std::string, int> meshes_by_material;
  std::vector<int> start_indexes;
  for (auto &mesh : loader.LoadedMeshes)
  {
    int pos;
    auto it = meshes_by_material.find(mesh.MeshMaterial.name);
    if (it == meshes_by_material.end())
    {
      pos = meshes_by_material.size();
      mod.models.push_back(new Model());
      
      std::string base_name_r = "";
      std::string bn_full = mesh.MeshMaterial.map_Ka;
      if (bn_full.empty())
        bn_full = mesh.MeshMaterial.map_Kd;
      for (int ch = bn_full.size()-1;ch>=0;ch--)
      {
        if (bn_full[ch] == '/' || bn_full[ch] == '\\')
          break;
        base_name_r.push_back(bn_full[ch]);
      }
      std::string base_name = base_name_r;
      for (int ch = base_name_r.size()-1;ch>=0;ch--)
        base_name[base_name_r.size()-1-ch] = base_name_r[ch];
      if (base_name == "")
        base_name = obj->get_string("color", "");
      std::string obj_color_tex = base_path + "/" + folder_name + "/" + base_name;
      Texture tex = base_name == "" ? 
                    engine::textureManager->get("white") :
                    engine::textureManager->load_unnamed_tex(obj_color_tex);
      if (!tex.is_valid())
      {
        logerr("texture manager cannot load file %s", obj_color_tex.c_str());
        return false;
      }
      mod.materials.push_back(Material(tex));

      start_indexes.push_back(0);
      meshes_by_material.emplace(mesh.MeshMaterial.name, pos);
    }
    else
    {
      pos = it->second;
    }
    Model *m = mod.models[pos];
    int &start_index = start_indexes[pos];
    //logerr("create mesh from %d", loader.LoadedMeshes.size());
    for (auto &lv : mesh.Vertices)
    {
      m->positions.push_back(lv.Position.X);
      m->positions.push_back(lv.Position.Y);
      m->positions.push_back(lv.Position.Z);

      m->normals.push_back(lv.Normal.X);
      m->normals.push_back(lv.Normal.Y);
      m->normals.push_back(lv.Normal.Z);

      m->colors.push_back(lv.TextureCoordinate.X);
      m->colors.push_back(1 - lv.TextureCoordinate.Y);
      m->colors.push_back(0);
      m->colors.push_back(1);
    }

    for (auto &ind : mesh.Indices)
    {
      m->indices.push_back(start_index + ind);
    }
    start_index += mesh.Vertices.size();
  }
  transform_model_to_standart_form(mod);
  mod.update();

  return true;
}

void save_model_to_obj(const Model *m, const std::string &filename)
{
  std::string header = "#obj file created by custom obj loader\n";
  std::string o_data = "o MainModel\n";
  std::string v_data;
  std::string tc_data;
  std::string n_data;
  std::string s_data = "s off\n";
  std::string f_data;

  int sz = m->positions.size()/3;
  assert(m->normals.size()/3 == sz);
  assert(m->colors.size()/4 == sz);

  for (int i = 0; i < sz; ++i)
  {
    v_data += "v " + std::to_string(m->positions[3*i]) + " " + std::to_string(m->positions[3*i+1]) + " " + std::to_string(m->positions[3*i+2]) + "\n";
    n_data += "vn " + std::to_string(m->normals[3*i]) + " " + std::to_string(m->normals[3*i+1]) + " " + std::to_string(m->normals[3*i+2]) + "\n";
    tc_data += "vt " + std::to_string(m->colors[4*i]) + " " + std::to_string(m->colors[4*i+1]) + "\n";
  }
  for (int i = 0; i < m->indices.size() / 3; ++i)
  {
    f_data += "f " + std::to_string(m->indices[3*i]+1) + "/" + std::to_string(m->indices[3*i]+1) + "/" + std::to_string(m->indices[3*i]+1) + " " +
              std::to_string(m->indices[3*i+1]+1) + "/" + std::to_string(m->indices[3*i+1]+1) + "/" + std::to_string(m->indices[3*i+1]+1) + " " +
              std::to_string(m->indices[3*i+2]+1) + "/" + std::to_string(m->indices[3*i+2]+1) + "/" + std::to_string(m->indices[3*i+2]+1) + "\n";
  }

  std::ofstream out(filename);
  out << header;
  out << o_data;
  out << v_data;
  out << tc_data;
  out << n_data;
  out << s_data;
  out << f_data;
  out.close();
}
}