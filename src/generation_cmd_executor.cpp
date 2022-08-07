#include "cmd_executors.h"
#include "common_utils/utility.h"
#include "generation/scene_generator_helper.h"
#include "generation/metainfo_manager.h"
#include "generation/scene_generator_plants.h"
#include "generation/grove_packer.h"
#include "parameter_selection/parameter_selection.h"
namespace scene_gen
{
  void align(float &from, float by)
  {
    from = ceil(from / by) * by;
  }
  void align(glm::vec2 &from, glm::vec2 by)
  {
    from.x = ceil(from.x / by.x) * by.x;
    from.y = ceil(from.y / by.y) * by.y;
  }
  void align(glm::vec2 &from, float by)
  {
    from.x = ceil(from.x / by) * by;
    from.y = ceil(from.y / by) * by;
  }

  void init_scene(Block &_settings, SceneGenerationContext &ctx)
  {
    metainfoManager.reload_all();
    ctx.settings = _settings;

    ctx.heightmap_size = ctx.settings.get_vec2("heightmap_size", glm::vec2(1000, 1000));
    ctx.hmap_pixel_size = ctx.settings.get_double("heightmap_cell_size", 8.0f);
    ctx.biome_map_pixel_size = ctx.settings.get_double("biome_map_pixel_size", 4.0f);
    ctx.full_size = ctx.settings.get_vec2("scene_size", glm::vec2(100, 100));
    ctx.grass_field_size = ctx.settings.get_vec2("grass_field_size", glm::vec2(1750, 1750));
    ctx.grass_field_size = max(min(ctx.grass_field_size, ctx.heightmap_size), ctx.full_size);
    ctx.center = ctx.settings.get_vec2("scene_center", glm::vec2(0, 0));
    ctx.cell_size = ctx.settings.get_vec2("cell_size", glm::vec2(150, 150));
    ctx.center3 = glm::vec3(ctx.center.x, 0, ctx.center.y);

    align(ctx.hmap_pixel_size, ctx.biome_map_pixel_size);
    align(ctx.heightmap_size, ctx.biome_map_pixel_size);
    align(ctx.grass_field_size, ctx.biome_map_pixel_size);
    align(ctx.full_size, ctx.biome_map_pixel_size);
    align(ctx.cell_size, ctx.biome_map_pixel_size);
    align(ctx.grass_field_size, ctx.cell_size);

    ctx.cells_x = ceil(ctx.grass_field_size.x / ctx.cell_size.x);
    ctx.cells_y = ceil(ctx.grass_field_size.y / ctx.cell_size.y);
    ctx.start_pos = ctx.center - 0.5f * ctx.cell_size * glm::vec2(ctx.cells_x, ctx.cells_y);

    ctx.cells = std::vector<Cell>(ctx.cells_x * ctx.cells_y, Cell(Cell::CellStatus::EMPTY));
    for (int i = 0; i < ctx.cells_x; i++)
    {
      for (int j = 0; j < ctx.cells_y; j++)
      {
        int id = i * ctx.cells_y + j;
        ctx.cells[id].id = id;
        ctx.cells[id].bbox = AABB2D(ctx.start_pos + ctx.cell_size * glm::vec2(i, j), ctx.start_pos + ctx.cell_size * glm::vec2(i + 1, j + 1));
        ctx.cells[id].influence_bbox = AABB(glm::vec3(ctx.cells[id].bbox.min_pos.x, -10, ctx.cells[id].bbox.min_pos.y),
                                            glm::vec3(ctx.cells[id].bbox.max_pos.x, 300, ctx.cells[id].bbox.max_pos.y));
      }
    }

    ctx.biome_map.create(AABB2D(ctx.center - ctx.grass_field_size, ctx.center + ctx.grass_field_size), ctx.biome_map_pixel_size);
    ctx.global_mask.create(glm::vec3(ctx.center.x,0,ctx.center.y), ctx.grass_field_size, ctx.biome_map_pixel_size);
    ctx.global_mask.fill_const(1);

    debug("Initialized scene\n");
    debug("Heightmap size %.1fx%.1f\n", ctx.heightmap_size.x, ctx.heightmap_size.y);
    debug("Grass field size %.1fx%.1f\n", ctx.grass_field_size.x, ctx.grass_field_size.y);
    debug("Vegetation scene size %.1fx%.1f\n", ctx.full_size.x, ctx.full_size.y);
    debug("created %dx%d cells with %.1fx%.1f size each\n", ctx.cells_x, ctx.cells_y, ctx.cell_size.x, ctx.cell_size.y);
    debug("created biome map %dx%d pixels\n", ctx.biome_map.pixels_w(), ctx.biome_map.pixels_h());

    ctx.scene.grove.center = ctx.center3;
    ctx.scene.grove.ggd_name = "blank";
  }


  AABB get_influence_AABB(Cell &c, const Heightmap &h)
  {
    glm::vec3 max_tree_size = glm::vec3(0.33f*(c.bbox.max_pos.x - c.bbox.min_pos.x),400,0.33f*(c.bbox.max_pos.y - c.bbox.min_pos.y));
    float min_hmap = 10000, max_hmap = 0;
    h.get_min_max_imprecise(c.bbox.min_pos, c.bbox.max_pos, &min_hmap, &max_hmap);
    float br = 5;
    float min_y = min_hmap - br;
    float max_y = max_hmap + max_tree_size.y;
    //TODO: handle different y of neighbour cells
    min_y = 0;
    max_y = max_tree_size.y;
    return AABB(glm::vec3(c.bbox.min_pos.x - max_tree_size.x, min_y, c.bbox.min_pos.y - max_tree_size.z),
                glm::vec3(c.bbox.max_pos.x + max_tree_size.x, max_y, c.bbox.max_pos.y + max_tree_size.z));
  }

  void create_heightmap(Block &settings, SceneGenerationContext &ctx)
  {
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    ctx.scene.heightmap = new Heightmap(ctx.center3, ctx.heightmap_size, ctx.hmap_pixel_size);
    bool load_from_image = settings.get_bool("load_from_image", true);
    float base = settings.get_double("base_height", 0);
    float mn = settings.get_double("min_height", 0);
    float mx = settings.get_double("max_height", 100);
    if (load_from_image)
    {
      std::string hmap_tex_name = settings.get_string("tex_name", "heightmap1.jpg");
      ctx.scene.heightmap->load_from_image(base, mn, mx, hmap_tex_name);
    }
    else
    {
      ctx.scene.heightmap->random_generate(base, mn, mx);
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float ms = 1e-4 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    debug("created heightmap. Took %.2f ms\n", ms);
  }

  void create_global_grove_mask(GroveMask &mask, SceneGenerationContext &ctx)
  {
    glm::vec2 mask_size = ctx.grass_field_size;
    std::function<float(glm::vec2 &)> func = [&](glm::vec2 &po) -> float
    {
      glm::vec3 p0 = glm::vec3(po.x + 0.5*mask.get_cell_size(), 0, po.y + 0.5*mask.get_cell_size());
      p0.y = ctx.scene.heightmap->get_height(p0) + ctx.biome_map_pixel_size;
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
  }

  uint64_t add_object_blk(Block &b, SceneGenerationContext &ctx)
  {
    std::string name = b.get_string("name", "debug_box");
    glm::mat4 transform = b.get_mat4("transform");

    bool place_on_terrain = b.get_bool("on_terrain", false);

    if (place_on_terrain)
    {
      glm::vec4 from = transform * glm::vec4(0, 0, 0, 1);
      glm::vec3 pos = glm::vec3(from);
      float min_h, float_max_h;
      pos.y = ctx.scene.heightmap->get_height(pos);
      transform[3][1] += (pos.y - from.y);
    }
    bool new_model = true;
    unsigned model_num = 0;
    unsigned pos = 0;
    for (auto &im : ctx.scene.instanced_models)
    {
      if (im.name == name)
      {
        new_model = false;
        pos = im.instances.size();
        break;
      }
      model_num++;
    }

    if (new_model)
    {
      ctx.scene.instanced_models.emplace_back();
      ctx.scene.instanced_models.back().model = model_loader::create_model_from_block(b, ctx.scene.instanced_models.back().tex);
      ctx.scene.instanced_models.back().model->update();
      ctx.scene.instanced_models.back().name = name;
      pos = 0;
    }

    //add object to scene BVH
    auto &im = ctx.scene.instanced_models[model_num];
    std::vector<AABB> boxes;
    SceneGenHelper::get_AABB_list_from_instance(im.model, transform, boxes, 4 * ctx.biome_map_pixel_size, 1.1);
    uint64_t id = SceneGenHelper::pack_id(0, (int)Scene::SIMPLE_OBJECT, model_num, pos);
    ctx.objects_bvh.add_bboxes(boxes, id);
    ctx.objects_bvh.rebuild();

    // get instance AABB
    glm::vec3 mn_pos = boxes[0].min_pos;
    glm::vec3 mx_pos = boxes[0].max_pos;
    for (AABB &box : boxes)
    {
      mn_pos = min(mn_pos, box.min_pos);
      mx_pos = max(mx_pos, box.max_pos);
    }
    
    im.bboxes.push_back(AABB(mn_pos, mx_pos));
    im.instances.push_back(transform);

    //update voxel arrays in cells
    glm::ivec2 c_ij_min = (glm::vec2(mn_pos.x, mn_pos.z) - ctx.start_pos) / ctx.cell_size;
    c_ij_min = max(c_ij_min, glm::ivec2(0,0));
    c_ij_min = min(c_ij_min, glm::ivec2(ctx.cells_x-1, ctx.cells_y-1));
    glm::ivec2 c_ij_max = (glm::vec2(mx_pos.x, mx_pos.z) - ctx.start_pos) / ctx.cell_size;
    c_ij_max = max(c_ij_max, glm::ivec2(0,0));
    c_ij_max = min(c_ij_max, glm::ivec2(ctx.cells_x-1, ctx.cells_y-1));
    for (int i=c_ij_min.x; i<=c_ij_max.x;i++)
    {
      for (int j=c_ij_min.y; j<=c_ij_max.y;j++)
      {
        Cell &c = ctx.cells[i * ctx.cells_y + j];
        if (c.voxels_small)
        {
          for (auto &box : boxes)
          {
            if (c.voxels_small->get_bbox().intersects(box))
              c.voxels_small->add_AABB(box,true, 10000);
          }
        }
      }
    }
    return id;
  }

  void plant_tree(glm::vec2 pos, int type, SceneGenerationContext &ctx)
  {
    if (type < 0 || type >= metainfoManager.get_all_tree_types().size())
    {
      logerr("Failed to place tree manually. Invalid type = %d", type);
      return;
    }
    glm::ivec2 c_ij = (pos - ctx.start_pos) / ctx.cell_size;
    if (c_ij.x >= 0 && c_ij.x < ctx.cells_x && c_ij.y >= 0 && c_ij.y < ctx.cells_y)
    {
      glm::vec3 p = glm::vec3(pos.x, 0, pos.y);
      p.y = ctx.scene.heightmap->get_bilinear(p);
      Cell &c = ctx.cells[c_ij.x * ctx.cells_y + c_ij.y];
      if (c.prototypes.empty())
      {
        c.prototypes.emplace_back();
        c.prototypes.back().pos = 0.5f * (c.bbox.max_pos + c.bbox.min_pos);
        c.prototypes.back().size = 0.5f * (c.bbox.max_pos - c.bbox.min_pos);
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
        c.prototypes.back().possible_types.push_back(std::pair<int, float>(type, 0));
      }
      c.prototypes.back().preplanted_trees.push_back(std::pair<int, glm::vec3>(type, p));
      c.prototypes.back().trees_count++;
    }
  }
}
void GenerationCmdExecutor::execute(int max_cmd_count)
{
  int cmd_left = max_cmd_count;
  while (!genCmdBuffer.empty() && cmd_left != 0)
  {
    auto &cmd = genCmdBuffer.front();
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    switch (cmd.type)
    {
    case GC_GEN_HMAP:
      if (genCtx.scene.heightmap)
        delete genCtx.scene.heightmap;
      scene_gen::create_heightmap(cmd.args, genCtx);
      if (genCtx.scene.heightmap)
      {
        for (auto &c : genCtx.cells)
          c.influence_bbox = scene_gen::get_influence_AABB(c, *(genCtx.scene.heightmap));
      }
      break;
    case GC_ADD_OBJECT:
      if (genCtx.scene.heightmap)
      {
        uint64_t id = scene_gen::add_object_blk(cmd.args, genCtx);
        debug("added object id = %lu \n", id);
      }
      break;
    case GC_CLEAR_SCENE:
      if (genCtx.inited)
      {
        genCtx = SceneGenerationContext();
        genCtx.inited = false;
      }
      break;
    case GC_INIT_SCENE:
      if (!genCtx.inited)
      {
        scene_gen::init_scene(cmd.args, genCtx);
        genCtx.inited = true;
      }
      break;
    case GC_REMOVE_BY_ID:
      for (int i = 0; i < cmd.args.size(); i++)
      {
        glm::ivec4 unpacked_id = cmd.args.get_ivec4(i, glm::ivec4(-1, 0, 0, 0));
        if (unpacked_id.x < 0)
          continue;
        else if (unpacked_id.y < 0)
        {
          // remove everything
          genCtx.scene.instanced_models.clear();
          genCtx.objects_bvh.clear();
        }
        else if (unpacked_id.z < 0)
        {
          // remove everything in this category
          if (unpacked_id.z == Scene::ObjCategories::SIMPLE_OBJECT)
          {
            genCtx.scene.instanced_models.clear();
            genCtx.objects_bvh.clear();
          }
        }
        else if (unpacked_id.w < 0)
        {
          // TODO
        }
        else
        {
          uint64_t id = SceneGenHelper::pack_id(unpacked_id.x, unpacked_id.y, unpacked_id.z, unpacked_id.w);
          auto func = [&](const std::pair<AABB, uint64_t> &p)
          {
            auto &ctx = genCtx;
            glm::ivec2 c_ij_min = (glm::vec2(p.first.min_pos.x, p.first.min_pos.z) - ctx.start_pos) / ctx.cell_size;
            c_ij_min = max(c_ij_min, glm::ivec2(0,0));
            c_ij_min = min(c_ij_min, glm::ivec2(ctx.cells_x-1, ctx.cells_y-1));
            glm::ivec2 c_ij_max = (glm::vec2(p.first.max_pos.x, p.first.max_pos.z) - ctx.start_pos) / ctx.cell_size;
            c_ij_max = max(c_ij_max, glm::ivec2(0,0));
            c_ij_max = min(c_ij_max, glm::ivec2(ctx.cells_x-1, ctx.cells_y-1));
            for (int i=c_ij_min.x; i<=c_ij_max.x;i++)
            {
              for (int j=c_ij_min.y; j<=c_ij_max.y;j++)
              {
                Cell &c = ctx.cells[i * ctx.cells_y + j];
                if (c.voxels_small)
                {
                  c.voxels_small->add_AABB(p.first,true, -10000);
                }
              }
            }
          };
          genCtx.objects_bvh.remove_bboxes_iterate(id, func);
          genCtx.objects_bvh.rebuild();
          if (unpacked_id.y == Scene::ObjCategories::SIMPLE_OBJECT)
          {
            int m_num = unpacked_id.z;
            int inst_num = unpacked_id.w;
            if (m_num >= 0 && m_num < genCtx.scene.instanced_models.size())
            {
              auto &im = genCtx.scene.instanced_models[m_num];
              if (inst_num >= 0 && inst_num < im.instances.size())
              {
                if (im.instances.size() > 1)
                {
                  im.instances.erase(im.instances.begin() + inst_num);
                  im.bboxes.erase(im.bboxes.begin() + inst_num);
                }
                else
                  genCtx.scene.instanced_models.erase(genCtx.scene.instanced_models.begin() + m_num);
              }
            }
          }
        }
      }
      break;
    case GC_PLANT_TREE:
      {
      std::string type_name = cmd.args.get_string("type_name");
      int type_id = metainfoManager.get_tree_type_id_by_name(type_name);
      if (type_id >= 0)
      {
        scene_gen::plant_tree(cmd.args.get_vec2("pos"), type_id, genCtx);
      }
      else
      {
        logerr("trying to create tree with unknown type %s", type_name.c_str());
      }
      }
      break;
    case GC_GEN_TREES_CELL:
    {
      std::vector<int> ids = {};
      for (int i=0;i<cmd.args.size();i++)
      {
        int id = cmd.args.get_int(i,-1);
        if (id >= 0)
          ids.push_back(id);
      }
      scene_gen::generate_plants_cells(genCtx, ids);
    }
      break;
    case GC_UPDATE_GLOBAL_MASK:
      scene_gen::create_global_grove_mask(genCtx.global_mask, genCtx);
      break;
    case GC_REMOVE_TREES:
    {
      std::vector<int> trees;
      cmd.args.get_arr("ids",trees);
      logerr("removing %d trees", trees.size());
      scene_gen::remove_trees_from_scene(genCtx, trees);
    }
      break;
    case GC_TREE_GEN_PARAMETER_SELECTION:
    {
      Block *set_bl = cmd.args.get_block("settings");
      Block *ref_bl = cmd.args.get_block("reference");
      if (set_bl && ref_bl)
      {
        ParameterSelector sel;
        auto result = sel.parameter_selection(*ref_bl, *set_bl, nullptr);

        int i = 0;
        for (auto &res : result.best_candidates)
        {
          metainfoManager.add_tree_type(res,std::string("__tmp_")+std::to_string(i));
          i++;
        }
      }
      else
      {
        logerr("cannot perform parameter selection settings or/and reference blocks are not found");
      }
    }
      break;
    case GC_GEN_BIOME_MAP:
      genCtx.biome_map.set_rect(genCtx.biome_map.borders(), metainfoManager.get_biome_id_by_name("empty_biome"));
      genCtx.biome_map.set_default_biome(metainfoManager.get_biome_id_by_name("empty_biome"));
      genCtx.biome_map.set_round(glm::vec3(0,0,0), 100, 150, metainfoManager.get_biome_id_by_name("simple_forest"));
      break;
    case GC_SET_DEFAULT_BIOME:
    {
      int id = cmd.args.get_int("id",-1);
      if (id != -1)
        genCtx.biome_map.set_default_biome(id);
    }
      break;
    case GC_SET_BIOME_ROUND:
    {
      int id = cmd.args.get_int("id",-1);
      if (id == -1)//we should remove biome i.e. set it to default
        id = genCtx.biome_map.get_default_biome();
      if (id != -1)
      {
        float outer_radius = cmd.args.get_double("outer_radius",-1);
        float inner_radius = cmd.args.get_double("inner_radius",-1);
        glm::vec3 pos = cmd.args.get_vec3("pos");
        if (inner_radius > 0 && outer_radius > inner_radius)
          genCtx.biome_map.set_round(glm::vec2(pos.x, pos.z), inner_radius, outer_radius, id);
      }
    }
      break;
    case GC_PREPARE_PLANT_PROTOTYPES:
    {
      scene_gen::prepare_tree_prototypes(genCtx);
      int t_cnt = 0;
      int c_cnt = 0;
      int p_cnt = 0;
      for (auto &c : genCtx.cells)
      {
        if (!c.prototypes.empty())
          c_cnt++;
        p_cnt += c.prototypes.size();
        for (auto &p : c.prototypes)
        {
          t_cnt += p.trees_count;
        }
      }
      logerr("created %d prototypes in %d cells. %d trees to generate",p_cnt, c_cnt, t_cnt);
    }
      break;
    default:
      logerr("GenerationCmdExecutor: command %d is not implemented yet", (int)(cmd.type));
      break;
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    float ms = 1e-4 * std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    logerr("%s took %.3f ms", ToString(cmd.type), ms);
    genCmdBuffer.pop();
    cmd_left--;
  }
}