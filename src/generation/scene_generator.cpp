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
#include "grass_generator.h"
#include "graphics_utils/debug_transfer.h"
#include "scene_generator_helper.h"
#include <algorithm>
#include <thread>
#include <mutex>

template <typename T>
struct atomwrapper
{
  /*Operations on the vector itself (e.g. adding or removing elements) must not be performed concurrently.*/
  std::atomic<T> _a;
  atomwrapper() : _a() {}
  atomwrapper(const std::atomic<T> &a) : _a(a.load()) {}
  atomwrapper(const atomwrapper &other) : _a(other._a.load()) {}
  atomwrapper &operator=(const atomwrapper &other) { _a.store(other._a.load()); }
};

std::mutex ctx_lock;
std::vector<atomwrapper<bool>> thread_finished;
LightVoxelsCube *create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
                                     AABB &influence_box)
{
  float min_scale_factor = 1000;
  for (auto &p : prototype.possible_types)
  {
    auto &type = types[p.first];
    min_scale_factor = MIN(min_scale_factor,type.params->get_scale_factor());
  }
  float vox_scale = 0.5/0.8;
  glm::vec3 voxel_sz = vox_scale*(influence_box.max_pos - influence_box.min_pos);
  glm::vec3 voxel_center = influence_box.min_pos + voxel_sz;
  auto *v = new LightVoxelsCube(voxel_center, voxel_sz, vox_scale*min_scale_factor, 1.0f);
  AABB &box = influence_box;
  float Mvoxels = 1e-6*v->get_size_cnt();
  debugl(1, "created voxels array [%.1f %.1f %.1f] - [%.1f %.1f %.1f] for patch [%.1f %.1f] - [%.1f %.1f] with %.2f Mvoxels\n",
  box.min_pos.x,box.min_pos.y,
  box.min_pos.z, box.max_pos.x,box.max_pos.y,box.max_pos.z,prototype.pos.x - prototype.size.x,
  prototype.pos.y - prototype.size.y, prototype.pos.x + prototype.size.x, prototype.pos.y + prototype.size.y,
  Mvoxels);
  return v;
}

LightVoxelsCube *create_cell_small_voxels(Cell &c, SceneGenerator::SceneGenerationContext &ctx)
{
  float vox_scale = 4;
  glm::vec3 voxel_sz = vox_scale*(c.influence_bbox.max_pos - c.influence_bbox.min_pos);
  glm::vec3 voxel_center = c.influence_bbox.min_pos + voxel_sz;
  
  auto *voxels = new LightVoxelsCube(voxel_center, voxel_sz, vox_scale, 1.0f);

  auto func = [&](const std::pair<AABB, uint64_t> &p)
  {
    logerr("added bbox");
    voxels->add_AABB(p.first,false, 10000);
  };
  logerr("added bbox 1");
  ctx.objects_bvh.iterate_over_intersected_bboxes(voxels->get_bbox(), func);

  return voxels;
}

AABB get_influence_AABB(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
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

class RawTreesDatabase
{
public:
  struct TreeToken
  {
    Tree *trees;
    int count;
    int array_n;
    int start_pos;
  };
  RawTreesDatabase() {};
  RawTreesDatabase(const RawTreesDatabase&) = delete;
  RawTreesDatabase( RawTreesDatabase&) = delete;
  TreeToken get_empty_trees(int cnt, GroveGenerationData &ggd);
  void generation_finished(TreeToken &token);
  bool pack_ready(GrovePacker &packer, SceneGenerator::SceneGenerationContext &ctx, bool forced);

  std::mutex database_lock;
private:
  static constexpr int max_trees_in_array = 32;
  struct TreeArray
  {
    Tree *data = nullptr;
    int cnt_max = max_trees_in_array;
    int cnt_real = 0;
    std::vector<bool> finished;//length = cnt;
    bool all_finished = false;
    bool filled = false;
    bool closed = false;
    GroveGenerationData ggd;
  };
  std::vector<TreeArray> arrays;
};

RawTreesDatabase::TreeToken RawTreesDatabase::get_empty_trees(int cnt, GroveGenerationData &_ggd)
{
  for (int i=0;i<arrays.size();i++)
  {
    if (arrays[i].closed || arrays[i].filled)
      continue;
    if (arrays[i].cnt_real + cnt <= arrays[i].cnt_max)
    {
      TreeToken tok = TreeToken{arrays[i].data + arrays[i].cnt_real, cnt, i, arrays[i].cnt_real};
      arrays[i].cnt_real += cnt;
      for (int j=0;j<cnt;j++)
        arrays[i].finished.push_back(false);
      return tok;
    }
    else
    {
      logerr("%d do not fin in %d/%d", cnt, arrays[i].cnt_real, arrays[i].cnt_max);
      //cannot place more trees here, array is filled
      arrays[i].filled = true;
    }
  }

  //haven't found place in existing arrays, create new one

  arrays.emplace_back();
  arrays.back().cnt_max = MAX(max_trees_in_array, cnt);
  arrays.back().data = new Tree[arrays.back().cnt_max];
  arrays.back().cnt_real = cnt;
  arrays.back().filled = (cnt == arrays.back().cnt_max);
  arrays.back().ggd = _ggd;
  arrays.back().finished = std::vector<bool>(cnt, false);

  return TreeToken{arrays.back().data, cnt, (int)arrays.size() - 1, 0};
};

void RawTreesDatabase::generation_finished(TreeToken &token)
{
  auto &arr = arrays[token.array_n];
  for (int i=token.start_pos;i<token.start_pos + token.count;i++)
    arr.finished[i] = true;
  
  if (arr.filled)
  {
    arr.all_finished = true;
    for (auto fin : arr.finished)
      arr.all_finished = arr.all_finished && fin;
  }
}

bool RawTreesDatabase::pack_ready(GrovePacker &packer, SceneGenerator::SceneGenerationContext &ctx, bool forced)
{
  if (arrays.empty())
    return false;
  bool job_finished = true;
  for (auto &arr : arrays)
  {
    if (!arr.closed && (arr.all_finished || forced))
    {
      arr.ggd.trees_count = arr.cnt_real;
      packer.add_trees_to_grove(arr.ggd, ctx.scene->grove, arr.data, ctx.scene->heightmap);
      delete[] arr.data;
      arr.closed = true;
    }
    job_finished = job_finished && arr.closed;
  }

  return job_finished;
}

struct GenerationJob
{
  SceneGenerator::SceneGenerationContext &ctx;

  std::vector<Cell> &cells;

  std::list<int> waiting_cells;
  std::list<int> border_cells;
  GroveGenerationData ggd;
  RawTreesDatabase &rawTreesDatabase;

  GrassGenerator &grassGenerator;
  std::mutex grassGenerator_lock;

  GroveMask &mask;

  int cells_x;
  int cells_y;
  bool grass_needed;
  int id;
  GenerationJob(SceneGenerator::SceneGenerationContext &_ctx, std::vector<Cell> &_cells,
                std::vector<TreeTypeData> &_types, RawTreesDatabase &_database, GrassGenerator &_grassGenerator,
                GroveMask &_mask, int _cells_x, int _cells_y, bool _grass_needed, int _id) : ctx(_ctx),
                                                                cells(_cells),
                                                                rawTreesDatabase(_database),
                                                                grassGenerator(_grassGenerator),
                                                                mask(_mask)
  {
    ggd.types = _types;
    cells_x = _cells_x;
    cells_y = _cells_y;
    grass_needed = _grass_needed;
    id = _id;
  };
  GenerationJob(const GenerationJob &) = delete;
  GenerationJob(GenerationJob &) = delete;
  void prepare_dependencies(); //done in main thread;
  void generate();//can be done in separate thread     
  bool remove_unused_borders();      
};

void SceneGenerator::create_scene_auto()
{
  generate_grove();
}
void SceneGenerator::generate_grove()
{
  GrovePacker packer;
  RawTreesDatabase rawTreesDatabase;

  int max_tc = ctx.settings.get_int("max_trees_per_patch", 1);
  int fixed_patches_count = ctx.settings.get_int("fixed_patches_count", 0);
  float patches_density = ctx.settings.get_double("patches_density", 1);
  glm::vec2 heightmap_size = ctx.settings.get_vec2("heightmap_size", glm::vec2(1000,1000));
  glm::vec2 full_size = ctx.settings.get_vec2("scene_size", glm::vec2(100,100));
  glm::vec2 grass_field_size = ctx.settings.get_vec2("grass_field_size", glm::vec2(500,500));
  grass_field_size = max(min(grass_field_size, heightmap_size), full_size);
  glm::vec2 center = ctx.settings.get_vec2("scene_center", glm::vec2(0,0));
  glm::vec2 cell_size = ctx.settings.get_vec2("cell_size", glm::vec2(50,50));
  glm::vec2 mask_pos = center;
  glm::vec3 center3 = glm::vec3(center.x,0,center.y);
  float hmap_cell_size = ctx.settings.get_double("heightmap_cell_size", 10.0f);

  Block *grass_blk = ctx.settings.get_block("grass");
  bool grass_needed = (grass_blk != nullptr);


  ctx.scene->heightmap = new Heightmap(center3, heightmap_size,hmap_cell_size);
  ctx.scene->heightmap->random_generate(0,0,100);
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


  int cells_x = ceil(grass_field_size.x/cell_size.x);
  int cells_y = ceil(grass_field_size.y/cell_size.y);
  glm::vec2 start_pos = center - 0.5f*cell_size*glm::vec2(cells_x, cells_y);
  std::vector<Cell> cells = std::vector<Cell>(cells_x*cells_y,Cell(Cell::CellStatus::EMPTY));
  
  GrassGenerator grassGenerator;
  glm::vec2 mask_size = full_size;
  if (cells_x <= 2 || cells_y <= 2)
    mask_size *= 2.0f;
  GroveMask mask = GroveMask(glm::vec3(mask_pos.x,0,mask_pos.y), mask_size, 3);
  mask.set_square(mask_size.x, mask_size.y);

  //std::list<int> waiting_cells;
  //std::list<int> border_cells;

  std::vector<GenerationJob *> generationJobs;
  const int max_jobs_cnt = 8;
  int job_size = MAX(ceil((float)cells_x/max_jobs_cnt), 1);
  int jobs_cnt = MAX(ceil((float)cells_x/job_size), 1);
  logerr("starting generation %dx%d cells with %d jobs (job size = %d)", cells_x, cells_y, jobs_cnt, job_size);

  for (int i=0;i<jobs_cnt;i++)
  {
    generationJobs.push_back(new GenerationJob(ctx, cells, types, rawTreesDatabase, grassGenerator, mask, cells_x, cells_y, 
                                               grass_needed, i));
    std::atomic<bool> ab(false);
    thread_finished.push_back(ab);                                      
  }

  for (int i=0;i<cells_x;i++)
  {
    for (int j=0;j<cells_y;j++)
    {
      int id = i*cells_y + j;
      cells[id].id = id;
      cells[id].bbox = AABB2D(start_pos + cell_size*glm::vec2(i,j), start_pos + cell_size*glm::vec2(i+1,j+1));
    }
  }

  if (fixed_patches_count > 0)
  {
    std::vector<int> cells_n = {};
    for (int i=0;i<cells.size();i++)
    {
      glm::vec2 ct = 0.5f*(cells[i].bbox.min_pos + cells[i].bbox.max_pos);
      if (mask.get_bilinear(glm::vec3(ct.x,0,ct.y)) > 0.1)
        cells_n.push_back(i);
    }
    std::random_shuffle(cells_n.begin(), cells_n.end());
    for (int i=0;i<MIN(fixed_patches_count, cells_n.size());i++)
    {
      cells[cells_n[i]].status = Cell::CellStatus::WAITING;
    }
  }

  for (int i=0;i<cells_x;i++)
  {
    for (int j=0;j<cells_y;j++)
    {
      int id = i*cells_y + j;

      int trees_count = 0;
      if (urand() <= patches_density || (cells[id].status == Cell::CellStatus::WAITING))
        trees_count = MAX(urand(0.4,0.6)*max_tc,1);
      if (grass_needed)
        cells[id].status = Cell::CellStatus::WAITING;
      if (cells[id].status == Cell::CellStatus::WAITING)
      {
        glm::vec2 center = start_pos + cell_size*glm::vec2(i+0.5,j+0.5);
        cells[id].prototype.pos = center;
        cells[id].prototype.size = 0.5f*cell_size;
        cells[id].prototype.possible_types = {};
        for (int i=0;i<types.size();i++)
        {
          cells[id].prototype.possible_types.push_back(std::pair<int, float>(i,1.0f/types.size()));
        }
        cells[id].prototype.trees_count = trees_count;
        cells[id].influence_bbox = get_influence_AABB(cells[id].prototype, ggd.types, *(ctx.scene->heightmap));
        generationJobs[i/job_size]->waiting_cells.push_back(id);
      }
    }
  }

  if (grass_needed)
  {
    grassGenerator.set_grass_types(ctx.grass_types, *grass_blk);
    grassGenerator.prepare_grass_patches(cells, cells_x, cells_y);
  }

  std::vector<std::thread> threads;
  for (auto *j : generationJobs)
  {
    j->prepare_dependencies();
  }
  for (auto *j : generationJobs)
  {
    threads.push_back(std::thread(&GenerationJob::generate, j));
  }

  for (auto &t : threads)
  {
    logerr("launching thread");
    t.detach();
  }

  bool threads_working = true;
  while (threads_working)
  {
    std::lock(rawTreesDatabase.database_lock, ctx_lock);
    rawTreesDatabase.pack_ready(packer, ctx, false);
    rawTreesDatabase.database_lock.unlock();
    ctx_lock.unlock();

    threads_working = false;
    for (auto fin : thread_finished)
    {
      threads_working = threads_working || !fin._a.load();
    }
  }
  logerr("finished packing");
  for (auto *j : generationJobs)
  {
    delete j;
  }

  rawTreesDatabase.pack_ready(packer, ctx, true);

  if (grass_needed)
  {
    for (auto &c : cells)
    {
      if ((c.status != Cell::CellStatus::EMPTY) && !c.grass_patches.empty())
        grassGenerator.generate_grass_in_cell(c, c.planar_occlusion);
    }
    grassGenerator.pack_all_grass(ctx.scene->grass, *(ctx.scene->heightmap));
  }

  for (auto &c : cells)
  {
    if (c.planar_occlusion)
      delete c.planar_occlusion;
    if (c.voxels_small)
    {
      logerr("missed dependencies. Cell %d", c.id);
      delete c.voxels_small;
    }
  }
}

void GenerationJob::prepare_dependencies()
{
  for (int c_id : waiting_cells)
  {
    if (cells[c_id].prototype.trees_count == 0)
      continue;
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

          int ncid = i*cells_y + j;
          if (i >= 0 && j >= 0 && i < cells_x && j < cells_y && ncid > c_id)
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
}
void GenerationJob::generate()
{
  for (int c_id : waiting_cells)
  {
    //generation trees and grass
    logerr("generating cell %d", c_id);
    bool short_lived_voxels = false;
    auto &c = cells[c_id];
    
    //c.cell_lock.lock();
    GrovePrototype prototype = c.prototype;
    auto influence_bb = c.influence_bbox;
    //c.cell_lock.unlock();

    if (prototype.trees_count > 0)
    {
      //temp stuff
      ggd.pos.x = prototype.pos.x;
      ggd.pos.z = prototype.pos.y;

      ctx_lock.lock();
      ggd.pos.y = ctx.scene->heightmap->get_height(ggd.pos);
      ctx_lock.unlock();

      ggd.size.x = prototype.size.x;
      ggd.size.z = prototype.size.y;
      ggd.trees_count = prototype.trees_count;

      //::Tree *trees = new ::Tree[ggd.trees_count];
      GroveGenerator grove_gen;
      LightVoxelsCube *voxels = create_grove_voxels(prototype, ggd.types, influence_bb);
      
      ctx_lock.lock();
      auto func = [&](const std::pair<AABB, uint64_t> &p)
      {
        voxels->add_AABB(p.first,false, 10000);
      };
      ctx.objects_bvh.iterate_over_intersected_bboxes(voxels->get_bbox(), func);
      ctx_lock.unlock();

      c.cell_lock.lock();
      std::vector<int> deps = c.depends_from;
      c.cell_lock.unlock();
      for (auto &dep_cid : deps)
      {
        cells[dep_cid].cell_lock.lock();
        voxels->add_voxels_cube(cells[dep_cid].voxels_small);
        cells[dep_cid].cell_lock.unlock();
      }

      ctx_lock.lock();
      voxels->add_heightmap(*(ctx.scene->heightmap));
      ctx_lock.unlock();

      rawTreesDatabase.database_lock.lock();
      auto token = rawTreesDatabase.get_empty_trees(ggd.trees_count, ggd);
      rawTreesDatabase.database_lock.unlock();

      grove_gen.prepare_patch(prototype, ggd.types, *(ctx.scene->heightmap), mask, *voxels, token.trees);
      debugl(1, "created patch with %d trees\n", ggd.trees_count);

      rawTreesDatabase.database_lock.lock();
      rawTreesDatabase.generation_finished(token);
      rawTreesDatabase.database_lock.unlock();

      //std::lock(ctx_lock, packer_lock);
      //packer.add_trees_to_grove(ggd, ctx.scene->grove, trees, ctx.scene->heightmap);
      //ctx_lock.unlock();
      //packer_lock.unlock();

      c.cell_lock.lock();
      if (!c.depends.empty())
      {
        c.status = Cell::CellStatus::BORDER;
        c.voxels_small = new LightVoxelsCube(voxels,glm::ivec3(0,0,0), voxels->get_vox_sizes(), 4,
                                            glm::vec2(0,1e8));
        border_cells.push_back(c_id);
      }
      else
      {
        c.status = Cell::CellStatus::FINISHED_PLANTS;
      }
      c.cell_lock.unlock();
      if (debugTransferSettings.save_detailed_voxels_count > debugTransferData.debug_voxels.size())
        debugTransferData.debug_voxels.push_back(voxels);
      else
      delete voxels;
      //delete[] trees;
    } 
    else 
    {
      c.cell_lock.lock();
      if (!c.depends.empty())
      {
        c.status = Cell::CellStatus::BORDER;
        c.voxels_small = create_cell_small_voxels(c, ctx);
        border_cells.push_back(c_id);
      }
      else
      {
        c.status = Cell::CellStatus::FINISHED_PLANTS;
      }
      c.cell_lock.unlock();
    }

    //TODO: can be done not every iteration
    remove_unused_borders();
    logerr("finished generating cell %d", c_id);
  }

  while (!remove_unused_borders())
    std::this_thread::yield();
  if (grass_needed)
  {
    auto it = waiting_cells.begin();
    while (it != waiting_cells.end())
    {
      if (!cells[*it].planar_occlusion)
      {
        AABB2D bbox = cells[*it].bbox;
        glm::vec3 occ_center = 0.5f * glm::vec3(bbox.max_pos.x + bbox.min_pos.x, 0, bbox.max_pos.y + bbox.min_pos.y);
        glm::vec2 bord = glm::vec2(bbox.max_pos.x - bbox.min_pos.x, bbox.max_pos.y - bbox.min_pos.y);
        float vox_size = MIN(bord.x, bord.y) / 64;
        cells[*it].planar_occlusion = new Field_2d(occ_center, 0.5f * bord, vox_size);
        std::function<float(glm::vec2 &)> func = [&](glm::vec2 &p) -> float
        {
          glm::vec3 pos = glm::vec3(p.x, 0, p.y);
          pos.y = ctx.scene->heightmap->get_height(pos) + 1;
          if (ctx.objects_bvh.contains(pos))
            return 1e9;
          else
            return 0;
        };
        cells[*it].planar_occlusion->fill_func(func);
      }
      it++;
    }
  }
  thread_finished[id]._a.store(true);
}

bool GenerationJob::remove_unused_borders()
{
  auto it = border_cells.begin();
  while (it != border_cells.end())
  {
    bool have_deps = false;
    cells[*it].cell_lock.lock();
    std::vector<int> deps = cells[*it].depends;
    cells[*it].cell_lock.unlock();
    for (int &dep : deps)
    {
      if (cells[dep].status.load() == Cell::CellStatus::WAITING)
      {
        have_deps = true;
        break;
      }
    }
    if (!have_deps)
    {
      cells[*it].cell_lock.lock();
      //logerr("removed dependency %d",*it);
      if (grass_needed)
      {
        if (cells[*it].voxels_small)
        {
          glm::vec3 occ_center = cells[*it].voxels_small->get_center();
          AABB bbox = cells[*it].voxels_small->get_bbox();
          glm::vec2 bord = glm::vec2(bbox.max_pos.x - bbox.min_pos.x, bbox.max_pos.z - bbox.min_pos.z);
          float vox_size = cells[*it].voxels_small->get_voxel_size();
          cells[*it].planar_occlusion = new Field_2d(occ_center, 0.5f * bord, vox_size);
          std::function<float(glm::vec2 &)> func = [&](glm::vec2 &p) -> float
          {
            return cells[*it].voxels_small->get_occlusion_projection(glm::vec3(p.x, 0, p.y));
          };
          cells[*it].planar_occlusion->fill_func(func);
        }
      }
      cells[*it].status = Cell::CellStatus::FINISHED_PLANTS;

      if (debugTransferSettings.save_small_voxels_count > debugTransferData.debug_voxels.size())
        debugTransferData.debug_voxels.push_back(cells[*it].voxels_small);
      else
        delete cells[*it].voxels_small;

      cells[*it].voxels_small = nullptr;
      it = border_cells.erase(it);
      cells[*it].cell_lock.unlock();
    }
    else
    {
      it++;
    }
  }

  return border_cells.empty();
}

uint64_t SceneGenerator::add_object_blk(Block &b)
{
  std::string name = b.get_string("name", "debug_box");
  glm::mat4 transform = b.get_mat4("transform");
  glm::vec4 from = transform*glm::vec4(0,0,0,1);
  glm::vec4 to = transform*glm::vec4(1,1,1,1);
  logerr("adding object %f %f %f - %f %f %f",from.x, from.y, from.z, to.x, to.y, to.z);
  bool new_model = true;
  unsigned model_num = 0;
  unsigned pos = 0;
  for (auto &im : ctx.scene->instanced_models)
  {
    if (im.name == name)
    {
      im.instances.push_back(transform);
      new_model = false;
      pos = im.instances.size() - 1;
      return SceneGenHelper::pack_id(0, (int)Scene::ERROR, 0, 0);
    }
    model_num++;
  }

  if (new_model)
  {
    ModelLoader loader;
    ctx.scene->instanced_models.emplace_back();
    ctx.scene->instanced_models.back().model = loader.create_model_from_block(b, ctx.scene->instanced_models.back().tex);
    ctx.scene->instanced_models.back().model->update();
    ctx.scene->instanced_models.back().instances.push_back(transform);
    pos = 0;
  }
  
  logerr("model num %d %d", model_num, ctx.scene->instanced_models.size());
  auto &im = ctx.scene->instanced_models[model_num];
  std::vector<AABB> boxes;
  SceneGenHelper::get_AABB_list_from_instance(im.model, transform, boxes, 4, 1.1);
  uint64_t id = SceneGenHelper::pack_id(0,(int)Scene::SIMPLE_OBJECT,model_num,pos);
  ctx.objects_bvh.add_bboxes(boxes, id);
  return id;
}

bool SceneGenerator::remove_object(uint64_t packed_id)
{
  unsigned _empty;
  unsigned category;
  unsigned type;
  unsigned id;

  SceneGenHelper::unpack_id(packed_id, _empty, category, type, id);
  if (category == (int)Scene::SIMPLE_OBJECT && type < ctx.scene->instanced_models.size() && 
      id < ctx.scene->instanced_models[type].instances.size())
  {
    auto beg = ctx.scene->instanced_models[type].instances.begin();
    ctx.scene->instanced_models[type].instances.erase(std::next(beg, id));
    ctx.objects_bvh.remove_bboxes(packed_id);
    return true;
  }
  else
  {
    return false;
  }
}