#include "scene_generator.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "grove_packer.h"
#include "core/tree.h"
#include "save_utils/blk.h"
#include "tree_generators/GE_generator.h"
#include "tree_generators/python_tree_gen.h"
#include "tree_generators/weber_penn_parameters.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/proctree.h"
#include "tree_generators/generated_tree.h"
struct Cell
{
  enum CellStatus
  {
    EMPTY,
    WAITING,
    BORDER,
    FINISHED
  };
  GrovePrototype prototype;
  LightVoxelsCube *voxels_small = nullptr;
  int id = -1;
  CellStatus status;
  std::vector<int> depends;//list of waiting cell (ids) that will use voxels from this cell
  std::vector<int> depends_from;
  AABB influence_bbox;
  explicit Cell(CellStatus _status = CellStatus::EMPTY) {status = _status;}
};

LightVoxelsCube *SceneGenerator::create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                                                     AABB &influence_box)
{
  float min_scale_factor = 1000;
  for (auto &p : prototype.possible_types)
  {
    auto &type = types[p.first];
    min_scale_factor = MIN(min_scale_factor,type.params->get_scale_factor());
  }
  glm::vec3 voxel_sz = 0.5f*(influence_box.max_pos - influence_box.min_pos);
  glm::vec3 voxel_center = influence_box.min_pos + voxel_sz;
  auto *v = new LightVoxelsCube(voxel_center, voxel_sz, 0.5f*min_scale_factor, 1.0f);
  AABB &box = influence_box;
  logerr("created bbox [%f %f %f] - [%f %f %f] for patch [%f %f] - [%f %f]",box.min_pos.x,box.min_pos.y,
  box.min_pos.z, box.max_pos.x,box.max_pos.y,box.max_pos.z,prototype.pos.x - prototype.size.x,
  prototype.pos.y - prototype.size.y, prototype.pos.x + prototype.size.x, prototype.pos.y + prototype.size.y);
  return v;
}
AABB SceneGenerator::get_influence_AABB(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                                        Heightmap &h)
{
  glm::vec3 max_tree_size = glm::vec3(0,0,0);
  for (auto &p : prototype.possible_types)
  {
    auto &type = types[p.first];
    max_tree_size = max(max_tree_size,type.params->get_tree_max_size());
  }
  
  float min_hmap = 0, max_hmap = 0;
  h.get_min_max_imprecise(prototype.pos - prototype.size, prototype.pos + prototype.size, &min_hmap, &max_hmap);
  float br = 5;
  float min_y = min_hmap - br;
  float max_y = max_hmap + max_tree_size.y;
  float y_center = (min_y + max_y)/2;
  float y_sz = (max_y - min_y)/2;
  glm::vec3 voxel_sz = glm::vec3(prototype.size.x + max_tree_size.x, y_sz, prototype.size.y + max_tree_size.z);
  glm::vec3 voxel_center = glm::vec3(prototype.pos.x, y_center, prototype.pos.y);
  return AABB(voxel_center - voxel_sz, voxel_center + voxel_sz);

}
AbstractTreeGenerator *get_gen(std::string &generator_name)
{
  AbstractTreeGenerator *gen;
  if (generator_name == "proctree")
    gen = new Proctree::ProctreeGenerator();
  else if (generator_name == "simple")
    gen = new SimpleTreeGenerator();
  else if (generator_name == "load_from_file")
    gen = new TreeLoaderBlk();
  else if (generator_name == "python_tree_gen")
    gen = new PythonTreeGen();
  else if (generator_name == "ge_gen")
    gen = new GETreeGenerator();
  else
    gen = new mygen::TreeGenerator();
  
  return gen;
}

void SceneGenerator::create_scene_auto()
{
  generate_grove();
}
void SceneGenerator::generate_grove()
{
  logerr("ctx settings size %d",ctx.settings.size());
  for (int i=0;i<ctx.settings.size();i++)
  {
    auto s = ctx.settings.get_name(i);
    logerr("%d set name %s",i,s.c_str());
  }
  GrovePacker packer;
  int max_tc = ctx.settings.get_int("max_trees_per_patch", 1);
  glm::vec2 full_size = ctx.settings.get_vec2("scene_size", glm::vec2(100,100));
  glm::vec2 center = ctx.settings.get_vec2("scene_center", glm::vec2(0,0));
  glm::vec2 start_pos = center - 0.5f*full_size;
  glm::vec2 cell_size = ctx.settings.get_vec2("cell_size", glm::vec2(50,50));
  glm::vec2 mask_pos = center;
  glm::vec2 heightmap_size = ctx.settings.get_vec2("heightmap_size", glm::vec2(1000,1000));
  glm::vec3 center3 = glm::vec3(center.x,0,center.y);
  float hmap_cell_size = ctx.settings.get_double("heightmap_cell_size", 10.0f);


  ctx.scene->heightmap = new Heightmap(center3, 0.5f*heightmap_size,hmap_cell_size);
  ctx.scene->heightmap->random_generate(0,0,10);
  logerr("heightmap size %f %f", ctx.scene->heightmap->get_size().x,ctx.scene->heightmap->get_size().y);
  ctx.scene->grove.center = center3;
  ctx.scene->grove.ggd_name = "blank";
  
  std::vector<TreeTypeData> types;
  Block *types_bl = ctx.settings.get_block("types");
  if (!types_bl)
  {
    logerr("error: scene generation settings should have \"types\" block");
    return;
  }  
  else
  {
    for (int i=0;i<types_bl->size();i++)
    {
      std::string type_name = types_bl->get_name(i);
      auto it = ctx.tree_types.find(type_name);
      if (it == ctx.tree_types.end())
      {
        logerr("error: using unknown tree type %s",type_name.c_str());
      }
      else
      {
        types.push_back(it->second);
      }
    }
    if (types.empty())
    {
      logerr("error: \"types\" block should contain at least one valid tree type");
      return;
    }
  }
  
  GroveGenerationData ggd;
  ggd.types = types; 


  GroveMask mask = GroveMask(glm::vec3(mask_pos.x,0,mask_pos.y), 0.5f*full_size, 3);
  mask.set_round(MIN(full_size.x, full_size.y));

  int cells_x = ceil(full_size.x/cell_size.x);
  int cells_y = ceil(full_size.y/cell_size.y);
  std::vector<Cell> cells = std::vector<Cell>(cells_x*cells_y,Cell(Cell::CellStatus::WAITING));
  std::list<int> waiting_cells;
  std::list<int> border_cells;

  for (int i=0;i<cells_x;i++)
  {
    for (int j=0;j<cells_y;j++)
    {
      int id = i*cells_y + j;
      cells[id].id = id;
      //TODO: do we need a cell here?
      cells[id].status = (i % 2 && j % 2) ? Cell::CellStatus::WAITING : Cell::CellStatus::EMPTY;
      //cells[id].status = Cell::CellStatus::WAITING;
      if (cells[id].status == Cell::CellStatus::WAITING)
      {
        glm::vec2 center = start_pos + cell_size*glm::vec2(i+0.5,j+0.5);
        cells[id].prototype.pos = center;
        cells[id].prototype.size = cell_size;
        cells[id].prototype.possible_types = {std::pair<int, float>(0,1)};
        cells[id].prototype.trees_count = MAX(urand()*max_tc,1);
        cells[id].influence_bbox = get_influence_AABB(cells[id].prototype, ggd.types, *(ctx.scene->heightmap));
        waiting_cells.push_back(id);
      }
    }
  }

  for (int c_id : waiting_cells)
  {
    //find dependencies
    int j0 = c_id % cells_y;
    int i0 = c_id / cells_y;
    bool search = true;
    int d = 3;
    int d_prev = 0;
    while (search)
    {
      search = false;
      auto func = [&](int i1, int j1)
      {
          int i = i0 + i1;
          int j = j0 + j1;

          logerr("test %d %d",i,j);
          int ncid = i*cells_y + j;
          if (i >= 0 && j >= 0 && i < cells_x && j <= cells_y && ncid > c_id)
          {
            auto &c = cells[ncid];
            if (c.status == Cell::CellStatus::WAITING && c.influence_bbox.intersects(cells[c_id].influence_bbox))
            {
              cells[c_id].depends.push_back(ncid);
              c.depends_from.push_back(c_id);
              search = true;
            }
          }
      };
      for (int i1=-d;i1<=d;i1++)
      {
        for (int j1=-d;j1<=d;j1++)
        {
          if (i1 <= -d_prev || j1 <= -d_prev || i1 >= d_prev || j1 >= d_prev)
            func(i1,j1);
        }
      }
      d++;
      d_prev = d;
    }
    debug("depends of cell %d: ",c_id);
    for (auto &d : cells[c_id].depends)
      debug("%d ",d);
    debugnl();
  }

  for (int c_id : waiting_cells)
  {
    auto &c = cells[c_id];

    //temp stuff
    ggd.pos.x = c.prototype.pos.x;
    ggd.pos.z = c.prototype.pos.y;
    ggd.pos.y = ctx.scene->heightmap->get_height(ggd.pos);
    ggd.size.x = cell_size.x;
    ggd.size.z = cell_size.y;
    ggd.trees_count = c.prototype.trees_count;

    ::Tree *trees = new ::Tree[ggd.trees_count];
    GroveGenerator grove_gen;
    GrovePrototype prototype;
    prototype.pos = glm::vec2(ggd.pos.x, ggd.pos.z);
    prototype.size = glm::vec2(ggd.size.x, ggd.size.z);
    prototype.trees_count = ggd.trees_count;
    prototype.possible_types = {std::pair<int, float>(0, 1)};
    LightVoxelsCube *voxels = create_grove_voxels(prototype, ggd.types, c.influence_bbox);
    voxels->add_heightmap(*(ctx.scene->heightmap));
    for (int i = 0; i < ggd.obstacles.size(); i++)
    {
      voxels->add_body(ggd.obstacles[i],10000);
    }
    for (auto &dep_cid : c.depends_from)
    {
      voxels->add_voxels_cube(cells[dep_cid].voxels_small);
    }
    //debug_voxels = voxels;

    grove_gen.prepare_patch(prototype, ggd.types, *(ctx.scene->heightmap), mask, *voxels, trees);

    packer.add_trees_to_grove(ggd, ctx.scene->grove, trees, ctx.scene->heightmap);

    if (!c.depends.empty())
    {
      c.status = Cell::CellStatus::BORDER;
      c.voxels_small = new LightVoxelsCube(voxels,glm::ivec3(0,0,0), voxels->get_vox_sizes(), 4,
                                           glm::vec2(0,1e8));
      border_cells.push_back(c_id);
    }
    else
    {
      c.status = Cell::CellStatus::FINISHED;
    }
    delete voxels;
    delete[] trees;
    //TODO: can be done not every iteration

    auto it = border_cells.begin();

    while (it != border_cells.end())
    {
      bool have_deps = false;
      for (int &dep : cells[*it].depends)
      {
        if (cells[dep].status == Cell::CellStatus::WAITING)
        {
          have_deps = true;
          break;
        }
      }
      if (!have_deps)
      {
        logerr("removed dependency %d",*it);
        cells[*it].status = Cell::CellStatus::FINISHED;
        delete cells[*it].voxels_small;
        it = border_cells.erase(it);
      }
      else
      {
        it++;
      }   
    }
  }
}
