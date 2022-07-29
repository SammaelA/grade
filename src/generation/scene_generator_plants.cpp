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
  
  float get_small_voxels_size(Cell &c)
  {
    float vox_scale = 0.5/0.8;
    int block_sz = LightVoxelsCube::get_default_block_size();
    float small_voxel_size = block_sz*vox_scale;
    float cell_size = (c.bbox.max_pos.x - c.bbox.min_pos.x);
    int vox_cnt = ceil(0.5*(cell_size/small_voxel_size - 1));
    vox_cnt = block_sz*ceil(vox_cnt/(block_sz + 0.0));
    small_voxel_size = cell_size/(2*vox_cnt + 1);

    return small_voxel_size;
  }
  
  LightVoxelsCube *create_grove_voxels(Cell &c)
  {
    glm::vec3 voxel_sz = 0.5f*(c.influence_bbox.max_pos - c.influence_bbox.min_pos);
    glm::vec3 voxel_center = 0.5f*(c.influence_bbox.max_pos + c.influence_bbox.min_pos);

    //we need to make small light voxels cube (made from this one with voxels 5 times bigger) have exactly the same size
    //as the cell it belongs to
    float small_voxel_size = get_small_voxels_size(c);
    int vox_cnt_y = ceil(0.5*(2*voxel_sz.y/small_voxel_size - 1));
    voxel_sz.y = 0.5*((2*vox_cnt_y + 1)*small_voxel_size);

    auto *v = new LightVoxelsCube(voxel_center, voxel_sz, small_voxel_size/LightVoxelsCube::get_default_block_size(), 1.0f);
    //AABB box = v->get_bbox();
    //debug("created voxels array [%.1f %.1f %.1f] - [%.1f %.1f %.1f]\n",
    //box.min_pos.x,box.min_pos.y,
    //box.min_pos.z, box.max_pos.x,box.max_pos.y,box.max_pos.z);
    return v;
  }

  void create_cell_small_voxels(Cell &c, SceneGenerator::SceneGenerationContext &ctx)
  {
    float vox_scale = get_small_voxels_size(c);
    glm::vec3 center = glm::vec3(0.5f*(c.bbox.max_pos.x + c.bbox.min_pos.x), 
                                 0.5f*(c.influence_bbox.max_pos.y + c.influence_bbox.min_pos.y),
                                 0.5f*(c.bbox.max_pos.y + c.bbox.min_pos.y));
    glm::vec3 size = glm::vec3(0.5f*(c.bbox.max_pos.x - c.bbox.min_pos.x), 
                               0.5f*(c.influence_bbox.max_pos.y - c.influence_bbox.min_pos.y),
                               0.5f*(c.bbox.max_pos.y - c.bbox.min_pos.y));
    c.voxels_small = new LightVoxelsCube(center, size, vox_scale, 1.0, 1, 2);
    auto func = [&](const std::pair<AABB, uint64_t> &p)
    {
      c.voxels_small->add_AABB(p.first,true, 10000);
    };
    ctx.objects_bvh.iterate_over_intersected_bboxes(c.voxels_small->get_bbox(), func);
  }

  enum {BRANCH_OCC_MUL = 100};
  
  void add_trees_to_voxels_array(int cnt, Tree *trees, LightVoxelsCube *voxels)
  {
    float blm = BRANCH_OCC_MUL/pow(voxels->get_voxel_size(),3);
    for (int i=0;i<cnt;i++)
    {
      for (auto bh : trees[i].branchHeaps)
      {
        for (auto &b : bh->branches)
        {
          for (auto &seg : b.segments)
          {
            float v = (1.0/3)*PI*length(seg.begin - seg.end)*
                      (seg.rel_r_begin*seg.rel_r_begin + seg.rel_r_begin*seg.rel_r_end + seg.rel_r_end*seg.rel_r_end);
            voxels->set_occluder_simple(seg.end, blm*v);
          }
        }
      }
    }
  }

  void remove_compressed_tree_from_voxels_array(SceneGenerator::SceneGenerationContext &ctx, CompressedTree &t, LightVoxelsCube *voxels)
  {
    if (!voxels)
      return;
    
    float blm = BRANCH_OCC_MUL/pow(voxels->get_voxel_size(),3);
    std::vector<CompressedTree::Node> nodes = t.LOD_roots;
    while (!nodes.empty())
    {
      std::vector<CompressedTree::Node> new_nodes;
      for (auto &node : nodes)
      {
        if (node.type == CompressedTree::MODEL)
        {
          auto &tr = ctx.scene->grove.instancedBranchesDirect[node.model_num]->IDA.transforms[node.instance_num];
          float size_mul_sq = 1;
          for (auto &br_id : ctx.scene->grove.instancedBranchesDirect[node.model_num]->branches)
          {
            auto &b = ctx.scene->grove.instancedCatalogue.get(br_id);
            if (b.joints.size() < 2)
              continue;
            glm::vec3 prev_pos = tr * glm::vec4(b.joints[0].pos, 1);
            for (int i=1;i<b.joints.size();i++)
            {
              glm::vec3 pos = tr * glm::vec4(b.joints[i].pos, 1);
              float v = (1.0/3)*PI*length(pos - prev_pos)*size_mul_sq*
                      (b.joints[i-1].r*b.joints[i-1].r + b.joints[i-1].r*b.joints[i].r + b.joints[i].r*b.joints[i].r);
              voxels->set_occluder_simple(pos, -blm*v);
              prev_pos = pos;
            }
          }
        }
        new_nodes.insert(new_nodes.end(), node.children.begin(), node.children.end());
      }
      nodes = std::move(new_nodes);
    }
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
    GroveGenerationData ggd;
    RawTreesDatabase &rawTreesDatabase;
    GroveMask &mask;

    int cells_x;
    int cells_y;
    int id;
    GenerationJob(SceneGenerator::SceneGenerationContext &_ctx, std::vector<Cell> &_cells,
                  std::vector<TreeTypeData> &_types, RawTreesDatabase &_database,
                  GroveMask &_mask, int _cells_x, int _cells_y, int _id) : ctx(_ctx),
                                                                  cells(_cells),
                                                                  rawTreesDatabase(_database),
                                                                  mask(_mask)
    {
      ggd.types = _types;
      cells_x = _cells_x;
      cells_y = _cells_y;
      id = _id;
    };
    GenerationJob(const GenerationJob &) = delete;
    GenerationJob(GenerationJob &) = delete;
    void prepare_dependencies(); //done in main thread;
    void generate();//can be done in separate thread           
  };
  void GenerationJob::prepare_dependencies()
  {
    for (int c_id : waiting_cells)
    {
      if (!cells[c_id].depends.empty())
        continue;
      int tc = 0;
      for (auto &prot : cells[c_id].prototypes)
        tc += prot.trees_count;
      if (tc == 0)
        continue;
      //find dependencies
      int j0 = c_id % cells_y;
      int i0 = c_id / cells_y;
      bool search = true;
      int d = 1;
      int d_prev = 0;
      while (search)
      {
        search = false;
        auto func = [&](int i1, int j1)
        {
            int i = i0 + i1;
            int j = j0 + j1;

            int ncid = i*cells_y + j;
            if (i >= 0 && j >= 0 && i < cells_x && j < cells_y && c_id != ncid)
            {
              auto &c = cells[ncid];
              if (c.influence_bbox.intersects(cells[c_id].influence_bbox))
              {
                cells[c_id].depends.push_back(ncid);
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
      //generation trees
      auto &c = cells[c_id];

      int tc = 0;
      for (auto &prot : c.prototypes)
        tc += prot.trees_count;
      
      if (tc > 0)
      {
        logerr("cell %d creating %d trees in %d tree prototypes", c.id, tc, c.prototypes.size());

        ggd.pos.x = 0.5f * (c.bbox.max_pos.x + c.bbox.min_pos.x);
        ggd.pos.z = 0.5f * (c.bbox.max_pos.y + c.bbox.min_pos.y);

        ctx_lock.lock();
        ggd.pos.y = ctx.scene->heightmap->get_height(ggd.pos);
        ctx_lock.unlock();

        ggd.size.x = 0.5f * (c.bbox.max_pos.x - c.bbox.min_pos.x);
        ggd.size.z = 0.5f * (c.bbox.max_pos.y - c.bbox.min_pos.y);
        ggd.trees_count = tc;

        GroveGenerator grove_gen;
        LightVoxelsCube *voxels = create_grove_voxels(c);
        c.cell_lock.lock();
        std::vector<int> deps = c.depends;
        if (!c.voxels_small)
          create_cell_small_voxels(c, ctx);
        voxels->add_voxels_cube(c.voxels_small);
        c.cell_lock.unlock();
        for (auto &dep_cid : deps)
        {
          cells[dep_cid].cell_lock.lock();
          if (!cells[dep_cid].voxels_small)
            create_cell_small_voxels(cells[dep_cid], ctx);
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

          if (c.voxels_small)
          {
            c.cell_lock.lock();
            add_trees_to_voxels_array(ggd.trees_count, token.trees, c.voxels_small);
            c.cell_lock.unlock();
          }
          for (auto &dep_cid : deps)
          {
            cells[dep_cid].cell_lock.lock();
            add_trees_to_voxels_array(ggd.trees_count, token.trees, cells[dep_cid].voxels_small);
            cells[dep_cid].cell_lock.unlock();
          }
        }
        delete voxels;
      } 

      c.cell_lock.lock();
      c.prototypes.clear();
      c.status = Cell::CellStatus::FINISHED_PLANTS;
      c.cell_lock.unlock();
    }

    thread_finished[id]._a.store(true);
  }

  void generate_plants_cells(SceneGenerator::SceneGenerationContext &ctx, std::vector<int> all_cell_ids)
  {
    std::vector<int> cell_ids;//only valid cells where we have something to generate
    for (auto &id : all_cell_ids)
    {
      if (id < 0 || id >= ctx.cells.size())
        logerr("trying to make plants in invalid cell id = %d",id);
      else if (ctx.cells[id].prototypes.size() > 0)
      {
        cell_ids.push_back(id);
        ctx.cells[id].status = Cell::CellStatus::WAITING;
      }
    }
    if (cell_ids.empty())
      return;
    
    GrovePacker packer;
    RawTreesDatabase rawTreesDatabase;
    std::vector<TreeTypeData> types = metainfoManager.get_all_tree_types();
    GroveGenerationData ggd;
    ggd.types = types; 
    ggd.task = MINIMUM_FOR_RENDER;
    std::vector<Cell> &cells = ctx.cells;

    std::vector<GenerationJob *> generationJobs;
    const int max_jobs_cnt = 8;
    int job_size = MAX(ceil((float)cell_ids.size()/max_jobs_cnt), 1);
    int jobs_cnt = MAX(ceil((float)cell_ids.size()/job_size), 1);
    logerr("starting generation %d cells with %d jobs (job size = %d)", cell_ids.size(), jobs_cnt, job_size);

    thread_finished = {};
    for (int i=0;i<jobs_cnt;i++)
    {
      generationJobs.push_back(new GenerationJob(ctx, cells, types, rawTreesDatabase, ctx.global_mask, ctx.cells_x, ctx.cells_y, i));
      std::atomic<bool> ab(false);
      thread_finished.push_back(ab);                                                                         
    }

    int cur_job = 0;
    for (auto &id : cell_ids)
    {
      if (cells[id].status == Cell::CellStatus::WAITING)
      {
        generationJobs[cur_job]->waiting_cells.push_back(id);
        if (jobs_cnt > 1)
          cur_job = (cur_job + 1)%jobs_cnt;
      }
    }

    for (auto *j : generationJobs)
    {
      j->prepare_dependencies();
    }
    if (jobs_cnt > 1)
    {
      std::vector<std::thread> threads;
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
    }
    else if (jobs_cnt == 1)
    {
      generationJobs[0]->generate();
    }
    for (auto *j : generationJobs)
      delete j;

    rawTreesDatabase.pack_ready(packer, ctx, true);
    packer.prepare_grove_atlas(ctx.scene->grove, 512, 512, true, true, true);
  }

  void remove_trees_from_scene(SceneGenerator::SceneGenerationContext &ctx, std::vector<int> &ids)
  {
    for (int &id : ids)
    {
      auto it = ctx.scene->grove.trees_by_global_id.find(id);
      logerr("%d %d", id, ctx.scene->grove.trees_by_global_id.size());
      if (it != ctx.scene->grove.trees_by_global_id.end())
      {
        auto &t = ctx.scene->grove.compressedTrees[it->second];
        glm::vec2 pos_xz = glm::vec2(t.pos.x, t.pos.z);
        glm::ivec2 c_ij = (pos_xz - ctx.start_pos)/ctx.cell_size;
        int cell_id = c_ij.x*ctx.cells_y + c_ij.y;
        if (cell_id >= 0 && cell_id < ctx.cells.size())
        {
          logerr("remove tree cell %d %d", cell_id, t.global_id);
          remove_compressed_tree_from_voxels_array(ctx, t, ctx.cells[cell_id].voxels_small);
          for (int &cid : ctx.cells[cell_id].depends)
          {
            if (ctx.cells[cid].bbox.intersectsXZ(t.bbox))
              logerr("remove tree cell %d %d", cid, t.global_id);
              remove_compressed_tree_from_voxels_array(ctx, t, ctx.cells[cid].voxels_small);
          }
        }
      }
    }
    GrovePacker::remove_trees_from_grove(ctx.scene->grove, ids);
  }
}