#include "scene_generator_plants.h"
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

namespace scene_gen
{
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
  void generate_plants_cells(SceneGenerator::SceneGenerationContext &ctx, std::vector<int> cell_ids)
  {
    GrovePacker packer;
    RawTreesDatabase rawTreesDatabase;
    std::vector<TreeTypeData> types = metainfoManager.get_all_tree_types();
    GroveGenerationData ggd;
    ggd.types = types; 
    ggd.task = GENERATE | CLUSTERIZE | MODELS;
    std::vector<Cell> &cells = ctx.cells;
    int &cells_x = ctx.cells_x;
    int &cells_y = ctx.cells_y;

    std::vector<GenerationJob *> generationJobs;
    const int max_jobs_cnt = 8;
    int job_size = MAX(ceil((float)cells_x/max_jobs_cnt), 1);
    int jobs_cnt = MAX(ceil((float)cells_x/job_size), 1);
    logerr("starting generation %dx%d cells with %d jobs (job size = %d)", cells_x, cells_y, jobs_cnt, job_size);

    GrassGenerator grassGenerator;//tmp
    bool grass_needed = false;//tmp

    glm::vec2 mask_pos = ctx.center;
    glm::vec2 mask_size = ctx.grass_field_size;
    GroveMask global_mask = GroveMask(glm::vec3(mask_pos.x,0,mask_pos.y), mask_size, ctx.biome_map_pixel_size);
    global_mask.fill_const(1);

    thread_finished = {};
    for (int i=0;i<jobs_cnt;i++)
    {
      generationJobs.push_back(new GenerationJob(ctx, cells, types, rawTreesDatabase, grassGenerator, global_mask, cells_x, cells_y, 
                                                grass_needed, i));
      std::atomic<bool> ab(false);
      thread_finished.push_back(ab);                                      
    }

    for (int i = 0; i < cells_x; i++)
    {
      for (int j = 0; j < cells_y; j++)
      {
        int id = i * cells_y + j;
        if (cells[id].prototypes.size() > 0 || grass_needed)
          cells[id].status = Cell::CellStatus::WAITING;

        if (cells[id].status == Cell::CellStatus::WAITING)
        {
          cells[id].influence_bbox = get_influence_AABB(cells[id].prototypes, ggd.types, *(ctx.scene->heightmap));
          generationJobs[i / job_size]->waiting_cells.push_back(id);
        }
      }
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
      std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
      rawTreesDatabase.pack_ready(packer, ctx, false);
      std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
      float ms = 1e-4*std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
      if (ms > 0.1)
        logerr("clustering took %.1f ms", ms);
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

    for (auto &c : cells)
    {
      if (c.planar_occlusion)
        delete c.planar_occlusion;
      if (c.voxels_small)
      {
        logerr("missed dependencies. Cell %d", c.id);
        delete c.voxels_small;
      }
      c.prototypes.clear();
    }

    packer.prepare_grove_atlas(ctx.scene->grove, 512, 512, true, true, true);
  }
};