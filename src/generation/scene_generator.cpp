#include "scene_generator.h"
#include "graphics_utils/volumetric_occlusion.h"
#include "grove_packer.h"
#include "core/tree.h"
#include "save_utils/blk.h"
#include "tree_generators/all_generators.h"
#include "grass_generator.h"
#include "graphics_utils/debug_transfer.h"
#include "scene_generator_helper.h"
#include "metainfo_manager.h"
#include <algorithm>
#include <thread>
#include <mutex>
#include <chrono>

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
    min_scale_factor = MIN(min_scale_factor,type.get_params()->get_scale_factor());
  }
  float vox_scale = 0.5/0.8;
  glm::vec3 voxel_sz = 0.5f*(influence_box.max_pos - influence_box.min_pos);
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
  float vox_scale = 4*ctx.biome_map_pixel_size;
  glm::vec3 voxel_sz = 0.5f*(c.influence_bbox.max_pos - c.influence_bbox.min_pos);
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

AABB get_influence_AABB(std::vector<GrovePrototype> &prototypes, std::vector<TreeTypeData> &types,
                        Heightmap &h)
{
  if (prototypes.empty())
    return AABB();
  glm::vec2 pos = prototypes[0].pos;
  glm::vec2 size = prototypes[0].size;

  glm::vec3 max_tree_size = glm::vec3(0,0,0);
  for (auto &prototype : prototypes)
  {
    for (auto &p : prototype.possible_types)
    {
      auto &type = types[p.first];
      max_tree_size = max(max_tree_size,type.get_params()->get_tree_max_size());
    }
  }

  float min_hmap = 0, max_hmap = 0;
  h.get_min_max_imprecise(pos - size, pos + size, &min_hmap, &max_hmap);
  float br = 5;
  float min_y = min_hmap - br;
  float max_y = max_hmap + max_tree_size.y;
  float y_center = (min_y + max_y)/2;
  float y_sz = (max_y - min_y)/2;

  glm::vec3 voxel_sz = glm::vec3(size.x + max_tree_size.x, y_sz, size.y + max_tree_size.z);
  glm::vec3 voxel_center = glm::vec3(pos.x, y_center, pos.y);
  return AABB(voxel_center - voxel_sz, voxel_center + voxel_sz);

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
      //logerr("%d do not fin in %d/%d", cnt, arrays[i].cnt_real, arrays[i].cnt_max);
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
    if (!arr.closed && (arr.cnt_real > 0) && (arr.all_finished || forced))
    {
      logerr("adding %d trees to grove",arr.cnt_real);
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

void align(float &from, float by)
{
  from = ceil(from/by)*by;
}
void align(glm::vec2 &from, glm::vec2 by)
{
  from.x = ceil(from.x/by.x)*by.x;
  from.y = ceil(from.y/by.y)*by.y;
}
void align(glm::vec2 &from, float by)
{
  from.x = ceil(from.x/by)*by;
  from.y = ceil(from.y/by)*by;
}

void SceneGenerator::init_scene(Block &_settings)
{
  ctx.settings = _settings;

  ctx.heightmap_size = ctx.settings.get_vec2("heightmap_size", glm::vec2(1000,1000));
  ctx.hmap_pixel_size = ctx.settings.get_double("heightmap_cell_size", 10.0f);
  ctx.biome_map_pixel_size = ctx.settings.get_double("biome_map_pixel_size", 1.0f);
  ctx.full_size = ctx.settings.get_vec2("scene_size", glm::vec2(100,100));
  ctx.grass_field_size = ctx.settings.get_vec2("grass_field_size", glm::vec2(1750,1750));
  ctx.grass_field_size = max(min(ctx.grass_field_size, ctx.heightmap_size), ctx.full_size);
  ctx.center = ctx.settings.get_vec2("scene_center", glm::vec2(0,0));
  ctx.cell_size = ctx.settings.get_vec2("cell_size", glm::vec2(150,150));
  ctx.center3 = glm::vec3(ctx.center.x,0,ctx.center.y);

  align(ctx.hmap_pixel_size, ctx.biome_map_pixel_size);
  align(ctx.heightmap_size, ctx.biome_map_pixel_size);
  align(ctx.grass_field_size, ctx.biome_map_pixel_size);
  align(ctx.full_size, ctx.biome_map_pixel_size);
  align(ctx.cell_size, ctx.biome_map_pixel_size);
  align(ctx.grass_field_size, ctx.cell_size);
  
  ctx.cells_x = ceil(ctx.grass_field_size.x/ctx.cell_size.x);
  ctx.cells_y = ceil(ctx.grass_field_size.y/ctx.cell_size.y);
  ctx.start_pos = ctx.center - 0.5f*ctx.cell_size*glm::vec2(ctx.cells_x, ctx.cells_y);
  
  ctx.cells = std::vector<Cell>(ctx.cells_x*ctx.cells_y,Cell(Cell::CellStatus::EMPTY));
  for (int i=0;i<ctx.cells_x;i++)
  {
    for (int j=0;j<ctx.cells_y;j++)
    {
      int id = i*ctx.cells_y + j;
      ctx.cells[id].id = id;
      ctx.cells[id].bbox = AABB2D(ctx.start_pos + ctx.cell_size*glm::vec2(i,j), ctx.start_pos + ctx.cell_size*glm::vec2(i+1,j+1));
    }
  }

  ctx.biome_map.create(AABB2D(ctx.center - 0.5f*ctx.grass_field_size, ctx.center + 0.5f*ctx.grass_field_size), 1);

  debug("Initialized scene\n");
  debug("Heightmap size %.1fx%.1f\n", ctx.heightmap_size.x, ctx.heightmap_size.y);
  debug("Grass field size %.1fx%.1f\n", ctx.grass_field_size.x, ctx.grass_field_size.y);
  debug("Vegetation scene size %.1fx%.1f\n", ctx.full_size.x, ctx.full_size.y);
  debug("created %dx%d cells with %.1fx%.1f size each\n", ctx.cells_x, ctx.cells_y, ctx.cell_size.x, ctx.cell_size.y);
  debug("created biome map %dx%d pixels\n", ctx.biome_map.pixels_w(), ctx.biome_map.pixels_h());
}

void SceneGenerator::create_heightmap_simple_auto()
{
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  ctx.scene->heightmap = new Heightmap(ctx.center3, ctx.heightmap_size,ctx.hmap_pixel_size);
  ctx.scene->heightmap->load_from_image(0,0,25,"heightmap1.jpg");
  std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
  float ms = 1e-4*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  debug("created heightmap. Took %.2f ms\n", ms);
}

void SceneGenerator::create_global_grove_mask(GroveMask &mask)
{
  glm::vec2 mask_size = ctx.grass_field_size;
  std::function<float(glm::vec2 &)> func = [&](glm::vec2 &po) -> float
  {
    glm::vec3 p0 = glm::vec3(po.x + 0.5*mask.get_cell_size(), 0, po.y + 0.5*mask.get_cell_size());
    p0.y = ctx.scene->heightmap->get_height(p0) + ctx.biome_map_pixel_size;
    int hits = 0;
    float kernel[9] = {1, 3, 1,
                       3, 6, 3,
                       1, 3, 1};
    for (int i=-1;i<=1;i++)
    {
      for (int j=-1;j<=1;j++)
      {
        glm::vec3 p = p0 + glm::vec3(i*0.5*mask.get_cell_size(), 0, j*0.5*mask.get_cell_size());
        hits += kernel[3*(i+1) + j + 1]*((int)ctx.objects_bvh.contains(p));
      }
    }
    return 1 - hits / 22.0; 
  };
  mask.fill_func(func);
  mask.save_as_image("grass_mask",0,1);
}

void SceneGenerator::generate_grove()
{
  GrovePacker packer;
  RawTreesDatabase rawTreesDatabase;
  glm::vec2 mask_pos = ctx.center;
  int max_tc = ctx.settings.get_int("max_trees_per_patch", 1);
  int fixed_patches_count = ctx.settings.get_int("fixed_patches_count", 0);
  float patches_density = ctx.settings.get_double("patches_density", 1);

  Block *grass_blk = ctx.settings.get_block("grass");
  bool grass_needed = (grass_blk != nullptr);

  ctx.scene->grove.center = ctx.center3;
  ctx.scene->grove.ggd_name = "blank";
  
  std::vector<TreeTypeData> types = metainfoManager.get_all_tree_types();
  GroveGenerationData ggd;
  ggd.types = types; 

  std::vector<Cell> &cells = ctx.cells;
  int &cells_x = ctx.cells_x;
  int &cells_y = ctx.cells_y;

/*
  for (auto &c : cells)
  {
    ctx.biome_map.get_stat(c.biome_stat, c.bbox);
    
    debug("cell %d stat:", c.id);
    for (auto &p : c.biome_stat)
    {
      debug(" (%d, %d)",p.first,p.second);
    }
    debugnl();
    
  }
*/
  GrassGenerator grassGenerator;
  glm::vec2 mask_size = ctx.grass_field_size;
  GroveMask global_mask = GroveMask(glm::vec3(mask_pos.x,0,mask_pos.y), mask_size, ctx.biome_map_pixel_size);
  create_global_grove_mask(global_mask);

  std::vector<GenerationJob *> generationJobs;
  const int max_jobs_cnt = 8;
  int job_size = MAX(ceil((float)cells_x/max_jobs_cnt), 1);
  int jobs_cnt = MAX(ceil((float)cells_x/job_size), 1);
  logerr("starting generation %dx%d cells with %d jobs (job size = %d)", cells_x, cells_y, jobs_cnt, job_size);

  for (int i=0;i<jobs_cnt;i++)
  {
    generationJobs.push_back(new GenerationJob(ctx, cells, types, rawTreesDatabase, grassGenerator, global_mask, cells_x, cells_y, 
                                               grass_needed, i));
    std::atomic<bool> ab(false);
    thread_finished.push_back(ab);                                      
  }

  if (fixed_patches_count > 0)
  {
    std::vector<int> cells_n = {};
    for (int i=0;i<cells.size();i++)
    {
      glm::vec2 ct = 0.5f*(cells[i].bbox.min_pos + cells[i].bbox.max_pos);
      if (global_mask.get_bilinear(glm::vec3(ct.x,0,ct.y)) > 0.1)
        cells_n.push_back(i);
    }
    std::random_shuffle(cells_n.begin(), cells_n.end());
    for (int i=0;i<MIN(fixed_patches_count, cells_n.size());i++)
    {
      cells[cells_n[i]].status = Cell::CellStatus::WAITING;
    }
  }

  prepare_patches(0);
  prepare_patches(1);

  for (int i=0;i<cells_x;i++)
  {
    for (int j=0;j<cells_y;j++)
    {
      glm::vec2 cell_center = ctx.start_pos + ctx.cell_size*glm::vec2(i+0.5,j+0.5);
      int id = i*cells_y + j;
      ctx.biome_map.get_stat(cells[id].biome_stat, cells[id].bbox);
      int cnt_all = 0;
      for (auto &p : cells[id].biome_stat)
      {
        cnt_all += p.second;
      }
      
      for (auto &p : cells[id].biome_stat)
      {
        //prepare main prototype for biome
        float fract = p.second/(float)cnt_all;

        if (fract < 0.01)
          continue;

        Biome &biome = metainfoManager.get_biome(p.first);

        for (auto &pn : cells[id].trees_patches)
        {
          if (ctx.trees_patches[pn].biome_id == p.first)
          {
            auto &patch = ctx.trees_patches[pn];

            GroveMask *mask = new GroveMask(glm::vec3(cell_center.x,0,cell_center.y), 
                                                      0.5f*ctx.cell_size, ctx.biome_map_pixel_size);

            std::function<float(glm::vec2 &)> func = [&](glm::vec2 &po) -> float
            {
              return patch.border.contains(po) ? 1 : 0; 
            };
            mask->fill_func(func);
            if (cells[id].prototypes.size() > 0 && cells[id].prototypes.back().biome_mask)
              mask->mul(*(cells[id].prototypes.back().biome_mask));
            else
            {
              GroveMask *biome_mask = new GroveMask(glm::vec3(cell_center.x,0,cell_center.y), 
                                                              0.5f*ctx.cell_size, ctx.biome_map_pixel_size);
              ctx.biome_map.set_mask(*biome_mask, p.first);
              mask->mul(*(biome_mask));
            }

            auto patch_prototype = GrovePrototype();
            cells[id].prototypes.back();
            
            int cells_cnt = 0;
            std::function<void(glm::vec2 &, float)> reader = [&](glm::vec2 &po, float val)
            {
              if (val > 0)
                cells_cnt++;
            };
            mask->read_func(reader);
            float tk = 1e-4*patch.density*cells_cnt;
            float tk_frac = tk - (int)tk;

            patch_prototype.pos = cell_center;
            patch_prototype.size = 0.5f*ctx.cell_size;
            patch_prototype.possible_types = patch.types;
            patch_prototype.trees_count = (int)tk + (tk_frac > urand());
            patch_prototype.biome_mask = mask;
            //logerr("cell %d prepared biome %d patch (%d cells) Expected %d trees",id,p.first,cells_cnt,
            //patch_prototype.trees_count);
            if (patch_prototype.trees_count > 0)
            {
              cells[id].prototypes.push_back(patch_prototype);
            }
            else if (patch_prototype.biome_mask)
            {
              delete patch_prototype.biome_mask;
            }
          }
        } 

        if (!biome.trees.main_types.empty())
        {
          auto prototype = GrovePrototype();
          
          prototype.pos = cell_center;
          prototype.size = 0.5f*ctx.cell_size;
          prototype.possible_types = biome.trees.main_types;
          float tk = 1e-4*biome.trees.main_density*p.second;
          float tk_frac = tk - (int)tk;
          prototype.trees_count = (int)tk + (tk_frac > urand());

          if (fract < 0.95)
          {
            prototype.biome_mask = new GroveMask(glm::vec3(prototype.pos.x,0,prototype.pos.y), prototype.size, ctx.biome_map_pixel_size);
            ctx.biome_map.set_mask(*(prototype.biome_mask), p.first);
          }
          if (prototype.trees_count > 0)
          {
            cells[id].prototypes.push_back(prototype);
          }
          else if (prototype.biome_mask)
          {
            delete prototype.biome_mask;
          }
          //logerr("cell %d prepared biome %d Expected %d trees",id,p.first,prototype.trees_count);
        }
      }
      
      if (cells[id].prototypes.size() > 0 || grass_needed)
        cells[id].status = Cell::CellStatus::WAITING;

      if (cells[id].status == Cell::CellStatus::WAITING)
      {
        cells[id].influence_bbox = get_influence_AABB(cells[id].prototypes, ggd.types, *(ctx.scene->heightmap));
        generationJobs[i/job_size]->waiting_cells.push_back(id);
      }
    }
  }

  if (grass_needed)
  {
    //grassGenerator.set_grass_types(*grass_blk);
    //grassGenerator.prepare_grass_patches(cells, cells_x, cells_y);
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
  for (auto *j : generationJobs)
  {
    delete j;
  }

  rawTreesDatabase.pack_ready(packer, ctx, true);

  if (grass_needed)
  {
    for (auto &c : cells)
    {
      //if ((c.status != Cell::CellStatus::EMPTY) && !c.grass_patches.empty())
      //  grassGenerator.generate_grass_in_cell(c, c.planar_occlusion);
      
      if (c.status != Cell::CellStatus::EMPTY)
      {
        glm::vec2 cell_center = 0.5f*(c.bbox.max_pos + c.bbox.min_pos);
        int cnt_all = 0;
        for (auto &p : c.biome_stat)
        {
          cnt_all += p.second;
        }
        for (auto &p : c.biome_stat)
        {
          //prepare grass for biome
          float fract = p.second/(float)cnt_all;

          if (fract < 0.01)
            continue;

          Biome &biome = metainfoManager.get_biome(p.first);
          GroveMask *biome_mask = new GroveMask(glm::vec3(cell_center.x,0,cell_center.y), 0.5f*ctx.cell_size, 
                                               ctx.biome_map_pixel_size);
          ctx.biome_map.set_mask(*(biome_mask), p.first);

          for (int patch_id : c.grass_patches)
          {
            if (ctx.grass_patches[patch_id].biome_id == p.first)
            {
              auto &patch = ctx.grass_patches[patch_id];
              std::function<float(glm::vec2 &)> patch_func = [&](glm::vec2 &po) -> float
              {
                return patch.border.contains(po) ? 1 : 0; 
              };
              GroveMask *patch_mask = new GroveMask(glm::vec3(cell_center.x,0,cell_center.y), 0.5f*ctx.cell_size, 
                                                              ctx.biome_map_pixel_size);
              patch_mask->fill_func(patch_func);
              patch_mask->mul(*biome_mask);
              grassGenerator.generate_grass_in_cell(c, c.planar_occlusion,patch_mask, &global_mask, 
                                                    patch.density, patch.types);

              delete patch_mask;
            }
          }
          
          std::function<float(glm::vec2 &, float)> func = [&](glm::vec2 &po, float occ) -> float
          {
            float res = occ;
            for (int patch_id : c.grass_patches)
            {
              if (ctx.grass_patches[patch_id].biome_id == p.first && ctx.grass_patches[patch_id].border.contains(po))
              {
                auto &patch_desc = biome.grass.patch_descs[ctx.grass_patches[patch_id].patch_type_id];
                res = CLAMP(MIN(res,1 - patch_desc.push),0,1);
              }
            }
            return res;
          };

          biome_mask->fill_func(func);
          grassGenerator.generate_grass_in_cell(c, c.planar_occlusion,biome_mask,  &global_mask, 
                                                biome.grass.main_density, biome.grass.main_types);

          delete biome_mask;        
        }
      }
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

  packer.prepare_grove_atlas(ctx.scene->grove, 512, 512, true, true, true);
}

void GenerationJob::prepare_dependencies()
{
  for (int c_id : waiting_cells)
  {
    int tc = 0;
    for (auto &prot : cells[c_id].prototypes)
      tc += prot.trees_count;
    if (tc == 0)
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
    /*
    debug("depends of cell %d: ",c_id);
    for (auto &d : cells[c_id].depends)
      debug("%d ",d);
    debugnl();
    */
  }
}
void GenerationJob::generate()
{
  for (int c_id : waiting_cells)
  {
    //generation trees and grass
    bool short_lived_voxels = false;
    auto &c = cells[c_id];
    auto influence_bb = c.influence_bbox;

    int tc = 0;
    for (auto &prot : c.prototypes)
      tc += prot.trees_count;
    
    if (tc > 0)
    {
      //temp stuff
      GrovePrototype base_prototype = c.prototypes[0];//all prototypes should have the same size
      ggd.pos.x = base_prototype.pos.x;
      ggd.pos.z = base_prototype.pos.y;

      ctx_lock.lock();
      ggd.pos.y = ctx.scene->heightmap->get_height(ggd.pos);
      ctx_lock.unlock();

      ggd.size.x = base_prototype.size.x;
      ggd.size.z = base_prototype.size.y;
      ggd.trees_count = base_prototype.trees_count;

      GroveGenerator grove_gen;
      LightVoxelsCube *voxels = create_grove_voxels(base_prototype, ggd.types, influence_bb);
      
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

      for (auto &prototype : c.prototypes)
      {
        ggd.trees_count = prototype.trees_count;
        rawTreesDatabase.database_lock.lock();
        auto token = rawTreesDatabase.get_empty_trees(ggd.trees_count, ggd);
        rawTreesDatabase.database_lock.unlock();

        grove_gen.prepare_patch(prototype, ggd.types, *(ctx.scene->heightmap), mask, *voxels, token.trees);
        debugl(1, "created patch with %d trees\n", ggd.trees_count);

        rawTreesDatabase.database_lock.lock();
        rawTreesDatabase.generation_finished(token);
        rawTreesDatabase.database_lock.unlock();
      }
      c.cell_lock.lock();
      if (!c.depends.empty())
      {
        c.status = Cell::CellStatus::BORDER;
        c.voxels_small = new LightVoxelsCube(voxels,glm::ivec3(0,0,0), voxels->get_vox_sizes(),
                                            voxels->get_block_size() > 1 ? voxels->get_block_size() : 5,
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
    //logerr("finished generating cell %d", c_id);
  }
  int tries = 0;
  int max_tries = 100;
  while (!remove_unused_borders())
  {
    tries = MIN(tries + 1, max_tries);
    for (int t=0;t<tries;t++)
      std::this_thread::yield();
  }
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
      cells[*it].cell_lock.unlock();
      it = border_cells.erase(it);
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
  /*
  glm::vec4 from = transform*glm::vec4(0,0,0,1);
  glm::vec4 to = transform*glm::vec4(1,1,1,1);
  logerr("adding object %f %f %f - %f %f %f",from.x, from.y, from.z, to.x, to.y, to.z);*/

  bool place_on_terrain = b.get_bool("on_terrain", false);

  if (place_on_terrain)
  {
    glm::vec4 from = transform*glm::vec4(0,0,0,1);
    glm::vec3 pos = glm::vec3(from);
    float min_h, float_max_h;
    pos.y = ctx.scene->heightmap->get_height(pos);
    transform[3][1] += (pos.y - from.y);
  }
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
      break;
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
    ctx.scene->instanced_models.back().name = name;
    pos = 0;
  }
  
  //logerr("model num %d %d", model_num, ctx.scene->instanced_models.size());
  auto &im = ctx.scene->instanced_models[model_num];
  std::vector<AABB> boxes;
  SceneGenHelper::get_AABB_list_from_instance(im.model, transform, boxes, 2*ctx.biome_map_pixel_size, 1.05);
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

SceneGenerator::SceneGenerator(SceneGenerationContext &_ctx):
ctx(_ctx)
{
  if (!metainfoManager.reload_all())
  {
    logerr("failed to load meta information for scene generator");
  }
}

void SceneGenerator::set_default_biome(std::string biome_name)
{
  int id = metainfoManager.get_biome_id_by_name(biome_name);
  ctx.base_biome_id = id;
  if (id>=0)
    ctx.biome_map.set_rect(ctx.biome_map.borders(),id);
}

void SceneGenerator::set_biome_round(glm::vec2 pos, float r, std::string biome_name)
{
  int id = metainfoManager.get_biome_id_by_name(biome_name);
  float outer_r = 1.5*r;
  if (id>=0)
    ctx.biome_map.set_round(pos,r,outer_r,id);
}

void SceneGenerator::plant_tree(glm::vec2 pos, int type)
{
  if (type < 0 || type >= metainfoManager.get_all_tree_types().size())
  {
    logerr("Failed to place tree manually. Invalid type = %d",type);
    return;
  }
  glm::ivec2 c_ij = (pos - ctx.start_pos)/ctx.cell_size;
  if (c_ij.x >= 0 && c_ij.x < ctx.cells_x && c_ij.y >= 0 && c_ij.y < ctx.cells_y)
  {
    glm::vec3 p = glm::vec3(pos.x, 0, pos.y);
    p.y = ctx.scene->heightmap->get_bilinear(p);
    Cell &c = ctx.cells[c_ij.x*ctx.cells_y + c_ij.y];
    if (c.prototypes.empty())
    {
      c.prototypes.emplace_back();
      c.prototypes.back().pos = 0.5f*(c.bbox.max_pos + c.bbox.min_pos);
      c.prototypes.back().size = 0.5f*(c.bbox.max_pos - c.bbox.min_pos);
      c.prototypes.back().trees_count = 0;
    }
    bool have_type = false;
    for (auto &tp : c.prototypes.back().possible_types)
    {
      if (tp.first == type)
      {
        have_type = true;
        break;
      }
    }
    if (!have_type)
    {
      c.prototypes.back().possible_types.push_back(std::pair<int,float>(type,0));
    }
    c.prototypes.back().preplanted_trees.push_back(std::pair<int,glm::vec3>(type,p));
    c.prototypes.back().trees_count++;
  } 
}

SceneGenerator::Patch::Patch(glm::vec2 pos, Biome::PatchDesc &patchDesc, int _biome_id, int _patch_type_id)
{
  Normal norm = Normal(0, 1);
  float r = patchDesc.size + patchDesc.size_std_dev*norm.get();
  border = Sphere2D(pos, r);
  density = patchDesc.density + patchDesc.density_std_dev*norm.get();
  types = patchDesc.types;
  patch_type_id = _patch_type_id;
  biome_id = _biome_id;
}

void SceneGenerator::prepare_patches(int patch_type)
{
  auto biomes = metainfoManager.get_all_biomes();
  std::vector<std::vector<std::pair<Biome::PatchDesc, float>>> patch_descs;
  
  for (auto &b : biomes)
  {
    patch_descs.emplace_back();
    auto &descs = ((patch_type == 0) ? b.trees.patch_descs : b.grass.patch_descs);
    for (auto &pd : descs)
    {
      float chance = pd.coverage_part*SQR(ctx.biome_map_pixel_size)/(PI*SQR(pd.size));
      if (chance > 1e-6)
        patch_descs.back().push_back(std::pair<Biome::PatchDesc, float>(pd,chance));
    }
  }

  auto &patches = (patch_type == 0) ? ctx.trees_patches : ctx.grass_patches;
  for (int i=0;i<ctx.biome_map.pixels_h();i++)
  {
    for (int j=0;j<ctx.biome_map.pixels_w();j++)
    {
      int type = ctx.biome_map.get(j,i);
      if (type >= 0 && type < patch_descs.size() && patch_descs[type].size() > 0)
      {
        glm::vec2 pos = ctx.biome_map.borders().min_pos + ctx.biome_map_pixel_size*glm::vec2(j,i);
        float rnd = urand();
        int desc_n = 0;
        for (auto &p : patch_descs[type])
        {
          if (p.second > rnd)
          {
            patches.push_back(Patch(pos, p.first, type, desc_n));
            for (auto &c : ctx.cells)//TODO: optimize me
            {
              if (patches.back().border.intersects(c.bbox))
              {
                if (patch_type == 0)
                  c.trees_patches.push_back(patches.size() - 1);
                else
                {
                  c.grass_patches.push_back(patches.size() - 1);
                }
              }
            }
            break;
          }
          else
          {
            rnd -= p.second;
          }
          desc_n++;
        }
      }
    }
  }
}