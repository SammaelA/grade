#include "gui.h"
#include "third_party/imgui/imgui.h"
#include "cmd_buffers.h"
#include "common_utils/utility.h"
#include "generation/metainfo_manager.h"
#include "tree_utils/tree_modeling.h"
#include "third_party/icons.h"
#include "hydra_utils/hydra_scene_exporter.h"
#include <map>

char __input_buf[4096];

void GUI::render()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplSDL2_NewFrame(engine::view->gWindow);
  ImGui::NewFrame();

  render_main_toolbar();
  read_from_console_nonblock();

  ImGui::Render();
  glViewport(0, 0, (int)engine::view->io.DisplaySize.x, (int)engine::view->io.DisplaySize.y);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  while(glGetError() != GL_NO_ERROR){}//ImGUI produces invalid value error every frame, I dont know why
}

bool GUI::input_text(const char *label, std::string &text)
{
  int len = text.size();
  strncpy(__input_buf, text.c_str(), text.size() + 1);
  bool get = ImGui::InputText(label, __input_buf, 4096, ImGuiInputTextFlags_EnterReturnsTrue);
  text = std::string(__input_buf);
  if (len != strlen(__input_buf))
    input_active = true;

  return get;
}

void GUI::blk_modification_interface(Block *b, const std::string &title)
{
  if (ImGui::CollapsingHeader(title.c_str()))
  {
    if (b)
    {
      for (int i=0;i<b->size();i++)
      {
        Block::Value &val = b->values[i];
        const char *name = b->names[i].c_str();
        switch (val.type)
        {
        case Block::ValueType::BOOL:
          ImGui::Checkbox(name, &(val.b));
        break;
        case Block::ValueType::INT:
        {
          int v = val.i;
          ImGui::InputInt(name, &v);
          val.i = v;
        }
        break;
        case Block::ValueType::UINT64:
        {
          int v = val.u;
          ImGui::InputInt(name, &v);
          val.u = v;
        }
        break;
        case Block::ValueType::DOUBLE:
          ImGui::InputDouble(name, &(val.d));
        break;
        case Block::ValueType::VEC2:
        {
          float v[2];
          v[0] = val.v2.x;
          v[1] = val.v2.y;
          ImGui::InputFloat2(name, v);
          val.v2.x = v[0];
          val.v2.y = v[1];
        }
        break;
        case Block::ValueType::VEC3:
        {
          float v[3];
          v[0] = val.v3.x;
          v[1] = val.v3.y;
          v[2] = val.v3.z;
          ImGui::InputFloat3(name, v);
          val.v3.x = v[0];
          val.v3.y = v[1];
          val.v3.z = v[2];
        }
        break;
        case Block::ValueType::VEC4:
        {
          float v[4];
          v[0] = val.v4.x;
          v[1] = val.v4.y;
          v[2] = val.v4.z;
          v[3] = val.v4.w;
          ImGui::InputFloat4(name, v);
          val.v4.x = v[0];
          val.v4.y = v[1];
          val.v4.z = v[2];
          val.v4.w = v[3];
        }
        break;
        case Block::ValueType::IVEC2:
        {
          int v[2];
          v[0] = val.iv2.x;
          v[1] = val.iv2.y;
          ImGui::InputInt2(name, v);
          val.iv2.x = v[0];
          val.iv2.y = v[1];
        }
        break;
        case Block::ValueType::IVEC3:
        {
          int v[3];
          v[0] = val.iv3.x;
          v[1] = val.iv3.y;
          v[2] = val.iv3.z;
          ImGui::InputInt3(name, v);
          val.iv3.x = v[0];
          val.iv3.y = v[1];
          val.iv3.z = v[2];
        }
        break;
        case Block::ValueType::IVEC4:
        {
          int v[4];
          v[0] = val.iv4.x;
          v[1] = val.iv4.y;
          v[2] = val.iv4.z;
          v[3] = val.iv4.w;
          ImGui::InputInt4(name, v);
          val.iv4.x = v[0];
          val.iv4.y = v[1];
          val.iv4.z = v[2];
          val.iv4.w = v[3];
        }
        break;
        case Block::ValueType::MAT4:
        {
          std::string n0 = std::string(name)+"[0]";
          std::string n1 = std::string(name)+"[1]";
          std::string n2 = std::string(name)+"[2]";
          std::string n3 = std::string(name)+"[3]";
          ImGui::InputFloat4(n0.c_str(), &(val.m4[0].x));
          ImGui::InputFloat4(n1.c_str(), &(val.m4[1].x));
          ImGui::InputFloat4(n2.c_str(), &(val.m4[2].x));
          ImGui::InputFloat4(n3.c_str(), &(val.m4[3].x));
        }
        break;
        case Block::ValueType::STRING:
        {
          input_text(name, *(val.s));
        }
        break;
        case Block::ValueType::BLOCK:
          blk_modification_interface(val.bl, b->names[i]);
          break;
        default:
          ImGui::Text("not supported");
          break;
        }
      }
    }
    else
    {
      ImGui::Text("empty");
    }
  }
}

GUI::GUI(AppContext &app_ctx, const SceneGenerationContext &gen_ctx) : appCtx(app_ctx),
                                                                       genCtx(gen_ctx)
{
  ImGuiIO &io = ImGui::GetIO();
  io.Fonts->AddFontDefault();

  // merge in icons from Font Awesome
  static const ImWchar icons_ranges[] = {ICON_MIN_FA, ICON_MAX_16_FA, 0};
  ImFontConfig icons_config;
  icons_config.MergeMode = true;
  icons_config.PixelSnapH = true;
  io.Fonts->AddFontFromFileTTF(FONT_ICON_FILE_NAME_FAS, 16.0f, &icons_config, icons_ranges);

  for (int i=0;i<InputCommands::IC_COMMANDS_COUNT;i++)
  {
    command_names.emplace(ToString((InputCommands)i), (InputCommands)i);
  }
};

void GUI::render_parameter_selection_menu()
{
  static bool selection_finished = false;

  if (selection_finished)
  {
    static std::string selected_name = "";
    if (selected_name != "")
    {
      static std::string new_name = "";
      bool get = input_text("New name", new_name);
      if (get)
      {
        TreeTypeData type = metainfoManager->get_tree_type(selected_name);
        metainfoManager->add_tree_type(type, new_name);
        selected_name = "";
        new_name = "";
      }
    }
    else
    {
      static int cur_item = -1;
      const char * types[16] = {nullptr};
      auto &types_map = metainfoManager->see_all_tree_type_names();
      int i=0;
      for (auto &p : types_map)
      {
        if (i>255)
          break;
        if (strncmp(p.first.c_str(),"__tmp", 5) == 0)
        {
          types[i] = p.first.c_str();
          i++;
        }
      }
      ImGui::ListBox("Selected types", &cur_item, types, i);
      bool add = ImGui::Button("Add this type");
      bool finish = ImGui::Button("Finish selection");
      if (add && cur_item >= 0)
        selected_name = std::string(types[cur_item]);
      if (finish)
      {
        metainfoManager->save_all();
        selection_finished = false;
      }
    }
  }
  else
  {
    static bool settings_loaded = false;
    static Block set_info, ref_info;
    if (!settings_loaded)
    {
      settings_loaded = true;
      
      set_info.clear();
      ref_info.clear();
      load_block_from_file("parameter_selection_settings.blk", set_info);
      load_block_from_file("parameter_selection_reference.blk", ref_info);
    }
    ImGui::Begin("Create new tree type by parameter selection");
    blk_modification_interface(&set_info, "Parameter selection settings");
    blk_modification_interface(&ref_info, "Parameter selection reference");

    if (ImGui::Button("Load settings"))
      settings_loaded = false;
    if (ImGui::Button("Save settings"))
    {
      
      save_block_to_file("parameter_selection_settings.blk", set_info);
      save_block_to_file("parameter_selection_reference.blk", ref_info);
    }
    bool start = ImGui::Button("Start");
    ImGui::End();

    if (start)
    {
      Block b;
      b.add_block("settings", &set_info);
      b.add_block("reference", &ref_info);
      inputCmdBuffer->push(IC_TREE_GEN_PARAMETER_SELECTION, b);
      selection_finished = true;
    }
  }
}

void GUI::render_main_toolbar()
{
  input_active = false;

  ImGuiWindowFlags window_flags = 0;
  window_flags |= ImGuiWindowFlags_NoMove;
  window_flags |= ImGuiWindowFlags_NoResize;
  window_flags |= ImGuiWindowFlags_NoCollapse;
  window_flags |= ImGuiWindowFlags_NoTitleBar;
  window_flags |= ImGuiWindowFlags_NoScrollbar;
  ImGui::SetNextWindowPos(ImVec2(0, -10), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImVec2(appCtx.window_width, 60), ImGuiCond_Always);
  ImGui::SetNextWindowContentSize(ImVec2(appCtx.window_width, 60));
  ImGui::Begin(" ", nullptr, window_flags);
  ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(15, 15));
  ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(10, 10));

  bool press_1 = ImGui::Button(ICON_FA_GEAR);ImGui::SameLine();
  bool press_2 = ImGui::Button(ICON_FA_HOUSE);ImGui::SameLine();
  bool press_3 = ImGui::Button(ICON_FA_TREE);ImGui::SameLine();
  bool press_4 = ImGui::Button(ICON_FA_COMPUTER);ImGui::SameLine();
  bool press_5 = ImGui::Button(ICON_FA_DIAGRAM_PROJECT);ImGui::SameLine();
  bool press_6 = ImGui::Button(ICON_FA_MOUNTAIN);ImGui::SameLine();
  bool press_7 = ImGui::Button(ICON_FA_IMAGE);ImGui::SameLine();
  bool press_8 = ImGui::Button(ICON_FA_FILE);
  ImGui::PopStyleVar(2);
  ImGui::End();
  {
    static bool show = false;
    if (press_1)
      show = !show;
    if (show)
      render_debug_settings();
  }
  {
    static bool show = false;
    if (press_2)
      show = !show;
    if (show)
      render_model_creation_info();
  }
  {
    static bool show = false;
    if (press_3)
      show = !show;
    if (show)
      render_tree_plant_info();
  }
  {
    static bool show = false;
    if (press_4)
      show = !show;
    if (show)
      text_input();
  }
  {
    static bool show = false;
    if (press_5)
      show = !show;
    if (show)
      render_parameter_selection_menu();
  }
  {
    static bool show = false;
    if (press_6)
      show = !show;
    if (show)
      render_biome_toolbar();    
  }
  {
    static bool show = false;
    if (press_7)
      show = !show;
    if (show)
      render_hydra_toolbar();    
  }
  {
    static bool show = false;
    if (press_8)
      show = !show;
    if (show)
      show = !render_save_settings();    
  }
  if (appCtx.active_cell_id >= 0)
    render_cell_info();
  
  if (input_active)
    appCtx.frames_from_last_input = 0;
  else 
    appCtx.frames_from_last_input++;
}

void GUI::render_biome_toolbar()
{
  static bool update = true;
  static std::vector<std::string> names;
  static const char *types[256];
  static int types_cnt = 0;
  static int cur_item = -1;
  static bool biome_brush = false;
  bool start_brush = false, end_brush = false, create_new_map = false, b_update = false;
  ImGui::Begin("Biome Map");

  if (ImGui::CollapsingHeader("BiomesList"))
  {
    ImGui::ListBox("Biomes", &cur_item, types, types_cnt);
  }
  else
    cur_item = -1;
  if (cur_item != -1)
  {
    if (ImGui::Button("Set as default biome"))
    {}
  }
  if (biome_brush && appCtx.biome_brush == -1)
  {
    biome_brush = false;
    end_brush = true;
  }
  if (biome_brush)
  {
    ImGui::SliderFloat("Brush scale", &appCtx.biome_brush_size, 10, 500);
    if (ImGui::Button("Stop modification [B]"))
    {
      end_brush = true;
    }
  }
  else 
  {
    if (ImGui::Button("Modify map"))
    {
      Block b;
      b.add_bool("render_biome_mask_debug", true);
      inputCmdBuffer->push(InputCommands::IC_UPDATE_RENDER_DEBUG_PARAMS, b);
      biome_brush = true;
      end_brush = false;
    }
    create_new_map = ImGui::Button("Create biome map");
    if (ImGui::Button("Prepare for trees generation"))
      inputCmdBuffer->push(InputCommands::IC_PREPARE_PLANT_PROTOTYPES);
    if (ImGui::Button("Regenerate all trees"))
    {
      inputCmdBuffer->push(InputCommands::IC_REMOVE_ALL_PLANTS);
      inputCmdBuffer->push(InputCommands::IC_GEN_ALL_PLANTED_TREES);
    }
    b_update = ImGui::Button("Update biome list");
  }
  ImGui::End();
  if (end_brush)
  {
    Block b;
    b.add_bool("render_biome_mask_debug", false);
    inputCmdBuffer->push(InputCommands::IC_UPDATE_RENDER_DEBUG_PARAMS, b);
    biome_brush = false;
  }

  if (b_update)
  {
    metainfoManager->reload_all();
  }
  if (update || b_update)
  {
    names.clear();
    for (int i=0;i<256;i++)
      types[i] = nullptr;
    auto &b_names = metainfoManager->see_all_biome_names();

    types_cnt = b_names.size();
    names.reserve(types_cnt);
    int j = 0;
    for (auto &p : b_names)
    {
      if (j>255)
        break;
      names.push_back(p.first);
      types[j] = names.back().c_str();
      j++;
    }
    update = false;
  }
  if (create_new_map)
    inputCmdBuffer->push(IC_GEN_BIOME_MAP);
  if (biome_brush && cur_item >= 0)
    appCtx.biome_brush = metainfoManager->get_biome_id_by_name(names[cur_item]);
  else  
    appCtx.biome_brush = -1;
}

void GUI::render_debug_settings()
{
    constexpr int options_cnt = 4;
    static bool checks[options_cnt] = {false, false, false, false};
    static bool prev_checks[options_cnt] = {false, false, false, false};
    static std::string check_names[options_cnt] = {"Show grid", "Show global mask", "Show_biome_mask", "Show bvh"};
    static std::string options_names[options_cnt] = {"render_grid_debug", "render_grove_mask_debug", "render_biome_mask_debug", 
                                                     "render_bvh_debug"};
    ImGui::Begin("Render debug params");
    for (int i=0;i<options_cnt;i++)
        ImGui::Checkbox(check_names[i].c_str(), &(checks[i]));
    ImGui::End();

    bool change = false;
    Block b;
    for (int i=0;i<options_cnt;i++)
    {
        if (checks[i] != prev_checks[i])
        {
            change = true;
            b.add_bool(options_names[i], checks[i]);
        }
        prev_checks[i] = checks[i];
    }

    if (change)
    {
        inputCmdBuffer->push(InputCommands::IC_UPDATE_RENDER_DEBUG_PARAMS, b);
    }
}

void GUI::render_cell_info()
{
  if (appCtx.active_cell_id < 0)
    return;
  auto &cell = appCtx.cells.at(appCtx.active_cell_id);
  bool visualize_prev = cell.visualize_voxels_array;
  ImGui::Begin("Cell info"); 
  ImGui::Text("Cell id: %d", appCtx.active_cell_id);
  ImGui::Checkbox("Show voxels array", &(cell.visualize_voxels_array));
  bool remove_trees = ImGui::Button("Remove trees");
  ImGui::End();

  if (visualize_prev && !cell.visualize_voxels_array)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer->push(InputCommands::IC_REMOVE_VOXELS_DEBUG, b);
  }
  else if (!visualize_prev && cell.visualize_voxels_array)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer->push(InputCommands::IC_VISUALIZE_VOXELS_DEBUG, b);    
  }
  if (remove_trees)
  {
    Block b;
    b.add_int("cell_id",appCtx.active_cell_id);
    inputCmdBuffer->push(InputCommands::IC_CLEAR_CELL, b);  
  }
}

void GUI::render_model_creation_info()
{
  static int cur_item = -1;
  int prev_cur_item = cur_item;
  static std::vector<std::string> names;
  static const char *types[256];
  static int types_cnt = 0;
  static bool update = true;

  static float angle_x = 0, angle_y = 0, angle_z = 0;
  static float scale = 1;
  static bool on_terrain = false;

  ImGui::Begin("Create object"); 
  ImGui::ListBox("Available models", &cur_item, types, types_cnt);
  ImGui::SliderAngle("Rotate X", &angle_x);
  ImGui::SliderAngle("Rotate Y", &angle_y);
  ImGui::SliderAngle("Rotate Z", &angle_z);
  ImGui::SliderFloat("Scale", &scale, 1, 1000);
  ImGui::Checkbox("On terrain", &on_terrain);
  appCtx.cur_obj_angles = glm::vec3(angle_x, angle_y, angle_z);
  appCtx.cur_obj_scale = scale;
  appCtx.cur_object_on_terrain = on_terrain;
  bool b_update = ImGui::Button("Update objects list");
  ImGui::End();

  if ((update || b_update) && model_loader::obj_models_blk_loaded)
  {
    names.clear();
    for (int i=0;i<256;i++)
      types[i] = nullptr;
    types_cnt = model_loader::obj_models_blk.size();
    names.reserve(types_cnt);
    for (int j=0;j<types_cnt;j++)
    {
      if (j>255)
        break;
      names.push_back(model_loader::obj_models_blk.get_name(j));
      types[j] = names.back().c_str();
    }
    update = false;
  }

  if (cur_item != prev_cur_item)
  {
    prev_cur_item = cur_item;
    if (cur_item >= 0)
      appCtx.active_object_name = types[cur_item];
  }
}

void GUI::render_tree_plant_info()
{
  static int cur_item = -1;
  int prev_cur_item = cur_item;
  const char * types[256] = {nullptr};
  auto &types_map = metainfoManager->see_all_tree_type_names();
  int i=0;
  for (auto &p : types_map)
  {
    if (i>255)
      break;
    if (strncmp(p.first.c_str(),"__", 2) != 0)
    {
      types[i] = p.first.c_str();
      i++;
    }
  }
  ImGui::Begin("Plant tree"); 
  ImGui::ListBox("Available Types", &cur_item, types, i);
  ImGui::End();
  if (cur_item != prev_cur_item)
  {
    prev_cur_item = cur_item;
    if (cur_item >= 0)
      appCtx.active_tree_type = types[cur_item];
  }
}

void GUI::render_hydra_toolbar()
{
  static bool inited = false;
  static Block settings;
  if (!inited)
  {
    hydra::get_default_settings(settings);
    inited = true;
  }

  settings.set_vec3("camera_pos", appCtx.camera.pos);
  settings.set_vec3("camera_look_at",appCtx.camera.pos + appCtx.camera.front);
  settings.set_vec3("camera_up", appCtx.camera.up);

  blk_modification_interface(&settings, "hydra exporter settings");
  bool exp = ImGui::Button("Export");
  if (exp)
  {
    inputCmdBuffer->push(IC_EXPORT_SCENE_TO_HYDRA, settings);
  }
  {
    appCtx.camera.pos = settings.get_vec3("camera_pos");
    appCtx.camera.front = glm::normalize(settings.get_vec3("camera_look_at") - appCtx.camera.pos);
    appCtx.camera.up = glm::normalize(settings.get_vec3("camera_up"));
  }
}

bool GUI::render_save_settings()
{
  static std::string save_name = "scene";
  static std::string load_name = "scene";

  ImGui::Begin("Save/Load scene");

  input_text("Save as", save_name);
  bool sv = ImGui::Button("Save");

  input_text("Load from file", load_name);
  bool ld = ImGui::Button("Load");
  bool clear = ImGui::Button("Clear scene");
  ImGui::End();

  if (sv)
  {
    Block b;
    b.add_string("name", save_name);
    inputCmdBuffer->push(IC_SAVE_SCENE, b);
  }
  if (ld)
  {
    Block b;
    b.add_string("name", load_name);
    inputCmdBuffer->push(IC_LOAD_SCENE, b);
  }
  if (clear)
  {
    inputCmdBuffer->push(IC_CLEAR_SCENE);
  }
  return sv || ld || clear;
}