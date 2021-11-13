
#define GLEW_EXPERIMENTAL
#include "tinyEngine/TinyEngine.h"
#include "glm/trigonometric.hpp"
#include "tinyEngine/image.h"
#include "tinyEngine/color.h"
#include "tinyEngine/helper.h"
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include <string>
#include "tinyEngine/camera.h"
#include "graphics_utils/modeling.h"
#include "graphics_utils/texture_manager.h"
#include "common_utils/utility.h"
#include "core/grove.h"
#include "render/grove_renderer.h"
#include "save_utils/config.h"
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include "graphics_utils/terrain.h"
#include "render/shadow.h"
#include "generation/grove_packer.h"
#include "app.h"
#include "render/grass_renderer.h"
#include "tinyEngine/deffered_target.h"
#include "render/ambient_occlusion.h"
#include "tinyEngine/cubemap.h"
#include "common_utils/python_interaction.h"
#include "gltf_utils/general_gltf_writer.h"
#include <thread>
#include "generation/parameter_selection.h"
#include "save_utils/blk.h"
#include "clustering/clustering_benchmark.h"
#include "tree_generators/load_tree_structure.h"
#include "clustering/clustering_debug_status.h"
#include "generation/grove_generator.h"
#include "tree_generators/GE_generator.h"
#include "tree_generators/python_tree_gen.h"
#include "tree_generators/weber_penn_parameters.h"
#include "tree_generators/simple_generator.h"
#include "tree_generators/proctree.h"
#include "tree_generators/generated_tree.h"
#include "render/billboard_cloud_renderer.h"
#include "render/visualizer.h"

View Tiny::view;   //Window and Interface  (Requires Initialization)
Event Tiny::event; //Event Handler
Audio Tiny::audio; //Audio Processor       (Requires Initialization)
TextureManager textureManager;
Config config;
AppContext appContext;
ShadowMap shadowMap;
DefferedTarget defferedTarget;
GrovePacker packer;
AbstractTreeGenerator *gen = nullptr;
DebugVisualizer *debugVisualizer = nullptr;
GrovePacked grove;
GroveGenerationData ggd;
GroveRenderer *GR = nullptr;
ClusteringBenchmark benchmark;
LightVoxelsCube *debug_voxels = nullptr;
struct RenderData
{
  bool regenerate_shadows = true;

  HBAORenderer *hbaoRenderer;
  Cubemap *cubemap;
  PostFx *defferedLight;
  PostFx *startScreenShader;
  Shader *defaultShader;
  Shader *debugShader;

  Heightmap *heightmap;
  GroveRenderer *groveRenderer;
  HeightmapTex *heightmapTex;
  GrassRenderer *grassRenderer;
  TerrainRenderer *terrainRenderer;

  void clear();
} data;

bool generation_needed = false;
bool saving_needed = false;
bool loading_needed = false;
bool print_perf = false;
bool only_gen = false;
bool visualize_voxels = false;
bool visualize_initial_voxels = false;
bool statistics_run = false;
bool gltf_export = false;
bool need_initialization = true;
bool parameter_selection = false;
bool clustering_benchmark = false;

std::string generator_name = "default";
std::string generator_fixed_preset_name = "";
std::string parameter_selector_name = "default_selection";
std::string clustering_benchmark_path = "benchmark.blk";
struct StatRunLaunchParams
{
  int trees = 1;
  float max_ind_dist = 0.7;
} statRunLaunchParams;
GroveRenderer::Precision pres = GroveRenderer::MEDIUM;
std::string grove_type_name = "default";
std::string save_path = ".";
std::string load_path = ".";

int parse_arguments(int argc, char *argv[])
{
  int k = 1;
  while (k < argc)
  {
    if (std::string(argv[k]) == "-g")
    {
      generation_needed = true;
      if (argc == k + 1)
      {
        logerr("write grove type name after -g");
        return 1;
      }
      else
        grove_type_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-perf")
    {
      print_perf = true;
      k++;
    }
    else if (std::string(argv[k]) == "-parameter_selection")
    {
      parameter_selection = true;
      if (argc > k+1)
      {
        parameter_selector_name = std::string(argv[k+1]);
        k++;
      }
      k++;
    }
    else if (std::string(argv[k]) == "-benchmark")
    {
      clustering_benchmark = true;
      if (argc > k+1)
      {
        clustering_benchmark_path = std::string(argv[k+1]);
        k++;
      }
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_voxels")
    {
      visualize_voxels = true;
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_clusters")
    {
      clusteringDebugInfo.visualize_clusters = true;
      k++;
    }
    else if (std::string(argv[k]) == "-visualize_initial_voxels")
    {
      visualize_initial_voxels = true;
      k++;
    }
    else if (std::string(argv[k]) == "-only_gen")
    {
      only_gen = true;
      k++;
    }
    else if (std::string(argv[k]) == "-export_gltf")
    {
      gltf_export = true;
      k++;
    }
    else if (std::string(argv[k]) == "-s")
    {
      saving_needed = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -s");
        return 1;
      }
      else
        save_path = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-l")
    {
      loading_needed = true;
      if (argc == k + 1)
      {
        logerr("write path to load from after -l");
        return 1;
      }
      else
        load_path = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-low_precision")
    {
      pres = GroveRenderer::LOW;
      k++;
    }
    else if (std::string(argv[k]) == "-high_precision")
    {
      pres = GroveRenderer::DEBUG;
      k++;
    }
    else if (std::string(argv[k]) == "-h")
    {
      logerr("-g <grove type name> -s <save path> -l <load path>");
      return 1;
    }
    else if (std::string(argv[k]) == "-stat_run")
    {
      if (argc <= k + 2)
      {
        logerr("-stat_run <num_tree> <max_individual distance>");
      }
      else
      {
        statRunLaunchParams.trees = atoi(argv[k + 1]);
        statRunLaunchParams.max_ind_dist = atof(argv[k + 2]);
        statistics_run = true;
        //only_gen = true;
        generation_needed = true;
        debug("stat run %d trees %f MID\n", statRunLaunchParams.trees, statRunLaunchParams.max_ind_dist);
        k += 3;
      }
    }
    else if (std::string(argv[k]) == "-generator")
    {
      bool ok = argc >= k+2;
      if (ok)
      {
        if (std::string(argv[k+1]) == "=")
          generator_name = std::string(argv[k+2]);
        else
          ok = false;
      }
      if (!ok)
      {
        logerr("-generator = <generator_name>");
        logerr("possible names : default, proctree, simple, python_tree_gen");
      }
      k += 3;
      if (generator_name == "python_tree_gen")
      {
        bool ok = argc >= k+2;
        if (ok && std::string(argv[k]) == "-parameters")
        {
          if (std::string(argv[k+1]) == "=")
            generator_fixed_preset_name = std::string(argv[k+2]);
          else
            ok = false;
        }
        else 
          ok = false;

        if (ok)
          k += 3;
      }
    }
    else if (std::string(argv[k]) == "-prepare_dataset")
    {
      clusteringDebugInfo.prepare_dataset = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -prepare_dataset");
        return 1;
      }
      else
        clusteringDebugInfo.dataset_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-save_csv")
    {
      clusteringDebugInfo.save_csv = true;
      if (argc == k + 1)
      {
        logerr("write path to save after -save_csv");
        return 1;
      }
      else
        clusteringDebugInfo.csv_file_name = argv[k + 1];
      k += 2;
    }
    else if (std::string(argv[k]) == "-progress_bar")
    {
      show_progress = 1;
      k++;
    }
    else
    {
      logerr("unknown command \"%s\". Write -h for help", argv[k]);
      return 1;
    }
  }
  if (saving_needed && !generation_needed)
  {
    logerr("Error: -s should be used only with -g flag");
    return 1;
  }
  if (!generation_needed && !loading_needed)
  {
    need_initialization = false;
    //logerr("You choosed to do nothing. Program will be closed");
    //return 0;
  }
  return -1;
}

void clear_current_grove()
{
  grove = GrovePacked();
  if (GR)
  {
    delete GR;
    GR = nullptr;
    data.groveRenderer = nullptr;
  }
  packer = GrovePacker();
}


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

LightVoxelsCube *create_grove_voxels(GrovePrototype &prototype, std::vector<TreeTypeData> &types,
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
void generate_grove()
{
  ggd.types[0].generator_name = generator_name;
  int max_tc = ggd.trees_count;
  glm::vec2 full_size = glm::vec2(200,200);
  glm::vec2 start_pos = glm::vec2(-100, -100);
  glm::vec2 cell_size = glm::vec2(ggd.size.x,ggd.size.z);
  glm::vec2 mask_pos = start_pos + 0.5f*full_size;
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
        cells[id].influence_bbox = get_influence_AABB(cells[id].prototype, ggd.types, *data.heightmap);
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
    ggd.pos.y = data.heightmap->get_height(ggd.pos);
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
    voxels->add_heightmap(*data.heightmap);
    for (auto &dep_cid : c.depends_from)
    {
      voxels->add_voxels_cube(cells[dep_cid].voxels_small);
    }
    //debug_voxels = voxels;

    grove_gen.prepare_patch(prototype, ggd.types, *data.heightmap, mask, *voxels, trees);

    packer.add_trees_to_grove(ggd, grove, trees, data.heightmap);

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
    if (debug_voxels) 
      delete debug_voxels;
    debug_voxels = voxels;
    //delete voxels;
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
  if (debug_voxels)
  {
    LightVoxelsCube *vox = new LightVoxelsCube(debug_voxels,glm::ivec3(0,0,0), debug_voxels->get_vox_sizes(), 4,
                                               glm::vec2(0,1e8));
    debug_voxels = vox;
  }
  /*
  for (int i = 0;i<2;i++)
  {
    for (int j=0;j<2;j++)
    {
      ggd.pos.x = (i + urand())*3*ggd.size.x;
      ggd.pos.z = (j + urand())*3*ggd.size.z;
      ggd.trees_count = MAX(urand()*max_tc,1);
      ggd.trees_count = 1;
      ::Tree *trees = new ::Tree[ggd.trees_count];

      GroveGenerator grove_gen;
      GrovePrototype prototype;
      prototype.pos = glm::vec2(ggd.pos.x, ggd.pos.z);
      prototype.size = glm::vec2(ggd.size.x, ggd.size.z);
      prototype.trees_count = ggd.trees_count;
      prototype.possible_types = {std::pair<int, float>(0,1)};
      LightVoxelsCube *voxels = create_grove_voxels(prototype, ggd.types,*data.heightmap);
      voxels->add_heightmap(*data.heightmap);
      //debug_voxels = voxels;

      GroveMask mask = GroveMask(ggd.pos, prototype.size,3);
      mask.set_round(MIN(prototype.size.x,prototype.size.y));
      grove_gen.prepare_patch(prototype, ggd.types, *data.heightmap, mask, *voxels, trees);

      packer.add_trees_to_grove(ggd, grove, trees, data.heightmap);
    
      delete[] trees;
      if (i == 1 && j== 1)
        debug_voxels = voxels;
      else
        delete voxels;
    }
  }
  */
}
void generate_single_tree(ParametersSet *par, GrovePacked &res)
{
    GrovePacker packer;
    GroveGenerationData tree_ggd;
    DistributionData dd;
    distibutionGenerator.d = &dd;

    tree_ggd.trees_count = 1;
    TreeTypeData type = ggd.types[0];
    type.params = par;
    tree_ggd.types = {type};
    tree_ggd.synts_count = 0;
    tree_ggd.synts_precision = 1;
    tree_ggd.pos = glm::vec3(0,0,0);
    tree_ggd.size = glm::vec3(150,150,250);
    tree_ggd.obstacles = {};
    tree_ggd.clustering_max_individual_distance = 0.0;
    tree_ggd.name = "single_tree";
    tree_ggd.task = GenerationTask::IMPOSTORS;
    Tree single_tree;
    gen->create_grove(tree_ggd, &single_tree, *data.heightmap);
    packer.add_trees_to_grove(tree_ggd, res, &single_tree, data.heightmap, false);
    print_alloc_info();
    distibutionGenerator.d = nullptr;
    dd.clear();
}
void generate_grove_renderer(GrovePacked *source, GroveGenerationData *gen_data)
{
  std::vector<float> LODs_dists = {15000, 1500, 500, 200, 30};
  if (pres == GroveRenderer::Precision::LOW)
    LODs_dists.back() = -10;
  if (data.groveRenderer)
  {
    delete data.groveRenderer;
  }
  data.groveRenderer = new GroveRenderer(source, gen_data, 5, LODs_dists, print_perf, pres);
  GR = data.groveRenderer;
}
int base_initialization()
{
  //base initialization
  glewInit();
  Tiny::view.lineWidth = 1.0f;
  Tiny::window("Procedural Tree", appContext.WIDTH, appContext.HEIGHT);
  Tiny::event.handler = [&]()
  { eventHandler(appContext, Tiny::event); };
  textureManager = TextureManager("./resources/textures/");
  data.startScreenShader = new PostFx("simple_render.fs");

  //PythonHelper ph;
  //ph.test();
  return -1;
}
int full_initialization()
{
  config.load_config();
  config.load_ggds();

  appContext.camera.pos = glm::vec3(-300, 70, 0);
  appContext.light.dir = glm::normalize(glm::vec3(-0.2, 0.5, -0));
  appContext.light.color = glm::vec3(0.99, 0.9, 0.7);
  appContext.light.intensity = 1;
  appContext.light.ambient_q = 0.5;
  appContext.light.diffuse_q = 0.4;
  appContext.light.specular_q = 0.1;
  appContext.light.has_shadow_map = true;
  appContext.light.shadow_map_size = glm::vec2(4096, 4096);
  shadowMap.create(4096, 4096);
  defferedTarget.create(2 * 1920, 2 * 1080);
  defferedTarget.set_clear_color(glm::vec4(0.0, 0.0, 0.0, 0.0));

  data.hbaoRenderer = new HBAORenderer();
  data.hbaoRenderer->create(800, 600);
  data.cubemap = new Cubemap(1920, 1080);
  data.defferedLight = new PostFx("deffered_light.fs");
  data.debugShader = new Shader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  data.defaultShader = new Shader({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"});
  debugVisualizer = new DebugVisualizer(textureManager.get("wood"), data.defaultShader);

  srand(1);

  data.heightmap = new Heightmap(glm::vec3(0, 0, 0), glm::vec2(2000, 2000), 5);
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
  data.heightmap->random_generate(0, 0, 10);

  if (generation_needed)
  {
    SDL_GL_SwapWindow(Tiny::view.gWindow);
    ggd = config.get_ggd(grove_type_name);
    if (generator_name == "proctree")
    {
      ggd.types[0].params = new Proctree::Properties();
    }
    if (generator_name == "python_tree_gen")
    {
      auto *p = new WeberPennParameters();
      if (!generator_fixed_preset_name.empty())
      {
        p->name = generator_fixed_preset_name;
        p->settings_already_in_file = true;
      }
      else
      {
        p->name = "default";
        p->settings_already_in_file = false;
      }
      ggd.types[0].params = p;
    }
    if (parameter_selection)
    {
      std::function<void(ParametersSet *, GrovePacked &)> _generate = generate_single_tree;
      ParameterSelector sel(_generate);
      ParametersSet *start = ggd.types[0].params;
      Quality imp_qual = ggd.impostor_quality;
      ggd.impostor_quality = Quality::ULTRALOW;
      //ggd.bill_1_quality = Quality::ULTRALOW;
      //ggd.bill_2_quality = Quality::ULTRALOW;
      sel.select(start,parameter_selector_name);
      ggd.impostor_quality = imp_qual;
    }
    
    if (statistics_run)
    {
      ggd.clustering_max_individual_distance = statRunLaunchParams.max_ind_dist;
      ggd.trees_count = statRunLaunchParams.trees;
    }
    if (clustering_benchmark)
    {
      benchmark.perform_benchmark(clustering_benchmark_path, gen, ggd, data.heightmap);
      if (benchmark.grove_count() > 0)
      {
        appContext.benchmark_grove_needed = 0;
      }
    }
    else
    {
      generate_grove();
    }
    if (visualize_initial_voxels)
    {
      mygen::TreeGenerator *mygen_gen = dynamic_cast<mygen::TreeGenerator *>(gen);
      if (mygen_gen)
        debug_voxels = mygen_gen->voxels;
    }
    if (debug_voxels)
    {
        debugVisualizer->visualize_light_voxels(debug_voxels);
    }
    /*LightVoxelsCube *voxx = new LightVoxelsCube(glm::vec3(0,0,0), glm::ivec3(8,8,8),1,3,2);
    for (int i=-8;i<=8;i++)
    for (int j=-8;j<=8;j++)
    for (int k=-8;k<=8;k++)
    voxx->set_occluder_trilinear(glm::vec3(i,j,k),abs(10*i));
    voxx->prepare_mip_levels();
    for (int i=0;i<0;i++)
    {
      logerr("visualize voxels");
      if (test_voxels_cube[0])
      {
        auto *voxels = test_voxels_cube[0];
        debugVisualizer->visualize_light_voxels(voxels,
                            voxels->get_center() - voxels->get_voxel_size()*glm::vec3(voxels->get_vox_sizes()),
                            voxels->get_voxel_size()*(2.0f*glm::vec3(voxels->get_vox_sizes()) + glm::vec3(1)),
                            powf(2,i)*glm::vec3(voxels->get_voxel_size()),
                            powf(2,i)*0.85f*voxels->get_voxel_size(),
                            0.01,glm::vec3(0,100*(i+1),0),glm::vec3(1,1,1),i);
      }
    }*/
  }
  if (saving_needed)
  {
    bool status = prepare_directory(save_path);
    if (status)
    {
      std::string f_path = save_path + std::string("/grove.dat");
      FILE *f = fopen(f_path.c_str(), "wb");
      saver::set_textures_path(save_path);
      saver::save(f, grove);
      fclose(f);
    }
    else
    {
      logerr("error occured while saving to path %s", save_path.c_str());
    }
  }
  if (loading_needed)
  {
    try
    {
      if (boost::filesystem::exists(load_path))
      {
        if (boost::filesystem::is_directory(load_path))
        {
          grove = GrovePacked();
          struct stat sb;
          if (stat(load_path.c_str(), &sb) != 0)
          {
            mkdir(load_path.c_str(), 0777);
          }
          if (!S_ISDIR(sb.st_mode))
          {
            logerr("given load path \"%s\" is not a directory", load_path.c_str());
          }
          std::string f_path = load_path + std::string("/grove.dat");
          FILE *f = fopen(f_path.c_str(), "rb");
          saver::set_textures_path(load_path);
          saver::load(f, grove);
          ggd = config.get_ggd(grove.ggd_name);
          logerr("grove %s", grove.ggd_name.c_str());
          fclose(f);
        }
        else
        {
          logerr("given load path %s is not a directory", load_path.c_str());
        }
      }
      else
      {
        logerr("given load path %s does not exist", load_path.c_str());
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
    }
  }
  generate_grove_renderer(&grove, &ggd);

  for (int i = 0; i < ggd.obstacles.size(); i++)
    debugVisualizer->add_bodies(ggd.obstacles[i], 1);
  data.terrainRenderer = new TerrainRenderer(*data.heightmap, glm::vec3(0, 0, 0), glm::vec2(2500, 2500), glm::vec2(25, 25));

  data.heightmapTex = new HeightmapTex(*data.heightmap, 2048, 2048);
  data.grassRenderer = new GrassRenderer();

  if (gltf_export)
  {
    gltf::GeneralGltfWriter ggw;
    ggw.add_packed_grove(grove, ggd);
    ggw.add_model(data.terrainRenderer->flat_terrain);
    ggw.convert_to_gltf("scene");
    ggw.clear();
    debug("Scene exported to .gltf format\n");
  }
  if (only_gen)
  {
    debug("Grove successfully generated. Exiting.");
    return 0;
  }

  appContext.renderMode = RenderMode::Rendering;
  Texture tex = textureManager.get("reference_tree_test");
  return -1;
}

void simple_render_pipeline()
{
  auto &ctx = appContext;

  ctx.fpsCounter.tick();
  if (ctx.fpsCounter.get_frame_n() % 100 == 0 && print_perf)
  {
    fprintf(stderr, "FPS: %4.1f\n", ctx.fpsCounter.get_average_fps());
  }
  if (data.regenerate_shadows)
  {
    //shadows pass
    data.regenerate_shadows = false;
    shadowMap.use(ctx.light);
    glm::mat4 sh_viewproj = shadowMap.get_transform();
    if (data.groveRenderer)
    {
      data.groveRenderer->render(data.groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), ctx.camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                ctx.light, ctx.groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (data.terrainRenderer)
    {
      data.terrainRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(),
                                  0, ctx.camera.pos, ctx.light, true);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                ctx.camera.pos, *data.heightmapTex, ctx.light, true);
    }
    shadowMap.start_trans_pass();
    if (data.groveRenderer)
    {
      data.groveRenderer->render(data.groveRenderer->get_max_LOD(), shadowMap.get_projection(),
                                shadowMap.get_view(), ctx.camera,
                                glm::vec2(shadowMap.SHADOW_WIDTH, shadowMap.SHADOW_HEIGHT),
                                ctx.light, ctx.groveRendererDebugParams, sh_viewproj, 0, true);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(shadowMap.get_projection(), shadowMap.get_view(), shadowMap.get_transform(), 0,
                                ctx.camera.pos, *data.heightmapTex, ctx.light, true);
    }
    shadowMap.finish_trans_pass();

    shadowMap.blur();
  }
  //Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  defferedTarget.target();

  if (ctx.render_mode <= ctx.DEBUG_RENDER_MODE)
  {
    data.debugShader->use();
    data.debugShader->uniform("projectionCamera", ctx.projection * ctx.camera.camera());
    if (ctx.render_mode == ctx.DEBUG_RENDER_MODE)
    {
      data.debugShader->texture("tex", textureManager.get(ctx.debug_tex));
      data.debugShader->texture("tex_arr", textureManager.get(ctx.debug_tex));
      data.debugShader->uniform("need_tex", true);
      data.debugShader->uniform("need_arr_tex", false);
      data.debugShader->uniform("need_coord", false);
      data.debugShader->uniform("slice", 0);
    }
    else if (ctx.render_mode == ctx.ARRAY_TEX_DEBUG_RENDER_MODE)
    {
      data.debugShader->texture("tex", textureManager.get_arr(ctx.debug_tex));
      data.debugShader->texture("tex_arr", textureManager.get_arr(ctx.debug_tex));
      data.debugShader->uniform("need_tex", false);
      data.debugShader->uniform("need_arr_tex", true);
      data.debugShader->uniform("need_coord", false);
      data.debugShader->uniform("slice", ctx.debug_layer);
    }
  }
  else
  {
    //depth prepass
    /*tr.render(ctx.projection * ctx.camera.camera(),shadowMap.get_transform(),shadowMap.getTex(),
                ctx.camera.pos,ctx.light, true);

      glClearColor(clearcolor.x, clearcolor.y, clearcolor.z, 1.0f);*/
    //color pass
    if (data.terrainRenderer)
    {
      data.terrainRenderer->render(ctx.projection, ctx.camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                  ctx.camera.pos, ctx.light);
    }
    if (data.grassRenderer)
    {
      data.grassRenderer->render(ctx.projection, ctx.camera.camera(), shadowMap.get_transform(), 0 * shadowMap.getTex(),
                                ctx.camera.pos, *data.heightmapTex, ctx.light);
    }
    if (ctx.render_mode != 2 && data.groveRenderer)
    {
      data.groveRenderer->render(ctx.forced_LOD, ctx.projection, ctx.camera.camera(), ctx.camera,
                                 glm::vec2(Tiny::view.WIDTH, Tiny::view.HEIGHT), ctx.light,
                                 ctx.groveRendererDebugParams, shadowMap.get_transform(), 0 * shadowMap.getTex());
    }
    if (debugVisualizer)
      debugVisualizer->render(ctx.projection * ctx.camera.camera(), ctx.render_mode);
  }

  //postfx
  /*uniform vec3 dir_to_sun;
uniform vec3 camera_pos;
uniform vec3 ambient_diffuse_specular;
uniform vec3 light_color;
uniform vec2 sts_inv;
uniform sampler2D shadowMap;
uniform bool need_shadow;
uniform mat4 shadow_mat;*/
  //hbaoRenderer.render(ctx,defferedTarget.get_view_pos());

  data.cubemap->render(ctx.projection, ctx.camera.camera(), ctx.camera);
  glm::vec3 ads = glm::vec3(ctx.light.ambient_q, ctx.light.diffuse_q, ctx.light.specular_q);
  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  data.defferedLight->use();
  data.defferedLight->get_shader().texture("colorTex", defferedTarget.get_color());
  data.defferedLight->get_shader().texture("normalsTex", defferedTarget.get_normals());
  data.defferedLight->get_shader().texture("viewPosTex", defferedTarget.get_view_pos());
  data.defferedLight->get_shader().texture("worldPosTex", defferedTarget.get_world_pos());
  data.defferedLight->get_shader().texture("aoTex", data.hbaoRenderer->get_tex());
  data.defferedLight->get_shader().texture("shadowMap", shadowMap.getTex());
  data.defferedLight->get_shader().texture("cubeTex", data.cubemap->get_tex());
  data.defferedLight->get_shader().uniform("dir_to_sun", ctx.light.dir);
  data.defferedLight->get_shader().uniform("camera_pos", ctx.camera.pos);
  data.defferedLight->get_shader().uniform("ambient_diffuse_specular", ads);
  data.defferedLight->get_shader().uniform("light_color", ctx.light.color);
  data.defferedLight->get_shader().uniform("need_shadow", shadowMap.getTex() != 0);
  data.defferedLight->get_shader().uniform("shadow_mat", shadowMap.get_transform());
  data.defferedLight->get_shader().uniform("sts_inv", 1.0f / ctx.light.shadow_map_size);
  data.defferedLight->render();
}

void start_screen_pipeline()
{
  Texture start_screen = textureManager.get("start_screen");
  Tiny::view.target(glm::vec3(0.6, 0.7, 1));
  data.startScreenShader->use();
  data.startScreenShader->get_shader().texture("tex", start_screen);
  data.startScreenShader->get_shader().uniform("tex_transform", glm::vec4(0, 1, 1, -1));
  data.startScreenShader->render();
}

int main(int argc, char *argv[])
{
  int status = parse_arguments(argc, argv);
  if (status >= 0)
    return status;

  status = base_initialization();
  if (status >= 0)
    return status;

  if (need_initialization && appContext.renderMode != RenderMode::Rendering)
  {
    status = full_initialization();
    if (status >= 0)
      return status;
  }

  Tiny::view.pipeline = [&]()
  {
    if (need_initialization && appContext.renderMode != RenderMode::Rendering)
    {
      full_initialization();
      //std::thread t1([&]{full_initialization();});
      //t1.join();
    }
    if (appContext.regeneration_needed)
    {
      clear_current_grove();
      generate_grove();
      generate_grove_renderer(&grove, &ggd);
      //print_alloc_info();
      appContext.regeneration_needed = false;
    }
    else if (appContext.add_generation_needed)
    {
      //print_alloc_info();
      generate_grove();
      //print_alloc_info();
      generate_grove_renderer(&grove, &ggd);
      //print_alloc_info();
      appContext.add_generation_needed = false;
    }
    if (appContext.benchmark_grove_current != appContext.benchmark_grove_needed &&
        appContext.benchmark_grove_needed >= 0)
    {
      int g = appContext.benchmark_grove_needed % benchmark.grove_count();
      appContext.benchmark_grove_current = appContext.benchmark_grove_needed;
      generate_grove_renderer(&(benchmark.get_grove(g)), &ggd);
    }
    if (appContext.renderMode == RenderMode::StartingScreen)
      start_screen_pipeline();
    else if (appContext.renderMode == RenderMode::Rendering)
      simple_render_pipeline();
  };

  Tiny::view.interface = [&]()
  {
    bool p = true;
    //ImGui::ShowDemoWindow(&p);
  };

  Tiny::loop([&]() {});
  {
    for (int i = 0; i < ggd.obstacles.size(); i++)
    {
      delete ggd.obstacles[i];
      ggd.obstacles[i] = nullptr;
    }
    ggd.obstacles.clear();
  }

  data.clear();
  distibutionGenerator.clear();
  Tiny::quit();

  return 0;
}

void RenderData::clear()
{
  regenerate_shadows = true;
  if (hbaoRenderer)
    delete hbaoRenderer;
  if (cubemap)
    delete cubemap;
  if (defferedLight)
    delete defferedLight;
  if (startScreenShader)
    delete startScreenShader;
  if (defaultShader)
    delete defaultShader;
  if (debugShader)
    delete debugShader;
  if (groveRenderer)
    delete groveRenderer;
  if (heightmapTex)
    delete heightmapTex;
  if (grassRenderer)
    delete grassRenderer;
  if (terrainRenderer)
    delete terrainRenderer;
  if (heightmap)
    delete heightmap;
}
