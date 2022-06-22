#include "visualizer.h"
#include "graphics_utils/texture_manager.h"

DebugVisualizer::DebugVisualizer():
Visualizer(),
debugShader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
bodyShader({"debug.vs", "body.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{

}
DebugVisualizer::DebugVisualizer(Texture _tree_tex, Shader *_tree_shader) : 
Visualizer(_tree_tex, textureManager.empty(), _tree_shader),
debugShader({"debug.vs", "debug.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
bodyShader({"debug.vs", "body.fs"}, {"in_Position", "in_Normal", "in_Tex"})
{
}
DebugVisualizer::~DebugVisualizer()
{
    for (int i = 0; i < debugModels.size(); i++)
    {
        delete (debugModels[i]);
    }
}
void DebugVisualizer::render(glm::mat4 view_proj, int mode)
{
    if (!tree_shader)
    {
        logerr("empty tree shader in debug visualizer\n");
    }
    else if ((mode == -1) || (mode == 1))
    {
        tree_shader->use();
        tree_shader->texture("tex", tree_tex);
        tree_shader->uniform("projectionCamera", view_proj);

        for (int i = 0; i < debugModels.size(); i++)
        {
            if (currentModes[i] != 1)
                continue;
            tree_shader->uniform("model", debugModels[i]->model);
            debugModels[i]->update();
            debugModels[i]->render(GL_TRIANGLES);
        }
    }
    if ((mode == -1) || (mode == 2))
    {
        bodyShader.use();
        bodyShader.uniform("projectionCamera", view_proj);
        bodyShader.uniform("need_coord",true);
        for (int i = 0; i < debugModels.size(); i++)
        {
            if (currentModes[i] != 2)
                continue;
            bodyShader.uniform("main_color",glm::vec3((i % 3)/3.0,(i % 5)/5.0,(i % 7)/7.0));
            bodyShader.uniform("model", debugModels[i]->model);
            debugModels[i]->update();
            debugModels[i]->render(GL_TRIANGLES);
        }
    }
        bodyShader.use();
        bodyShader.uniform("projectionCamera", view_proj);
        bodyShader.uniform("need_coord",false);
        for (int i = 0; i < debugModels.size(); i++)
        {
            if (currentModes[i] != 3)
                continue;
            bodyShader.uniform("main_color",glm::vec3((i % 3)/3.0,(i % 5)/5.0,(i % 7)/7.0));
            bodyShader.uniform("model", debugModels[i]->model);
            debugModels[i]->update();
            debugModels[i]->render(GL_TRIANGLES);
        }
}
void DebugVisualizer::branch_to_model_debug(Branch *b, int level, Model &m)
{
    if (b->level < 0)
        return;
    if (b->level == level || level < 0)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, MAX(3 - b->level,0));
        for (auto &segment : b->segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3, 1);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b->segments.back(), vets.back(), ringsize, (float)(i % 3) / 3, 1);
        }
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, &m);
        }
    }
    for (auto &joint : b->joints)
    {
        for (auto branch : joint.childBranches)
            branch_to_model_debug(branch, level, m);
    }
}
void DebugVisualizer::add_branch_debug(Branch *b, glm::vec3 scale, glm::vec3 shift, int level)
{
    Model *m = new Model();
    debugModels.push_back(m);
    currentModes.push_back(1);
    branch_to_model_debug(b, level, *m);
    m->shift(shift);
    m->scale = scale;
}
void DebugVisualizer::enable_all()
{
    for (int i = 0; i < currentModes.size(); i++)
        currentModes[i] = 1;
}
void DebugVisualizer::disable_all()
{
    for (int i = 0; i < currentModes.size(); i++)
        currentModes[i] = 0;
}
void DebugVisualizer::add_bodies(Body *b_ptr, int count)
{
    Model *m = new Model();
    debugModels.push_back(m);
    currentModes.push_back(3);
    for (int i=0;i<count;i++)
    {
        body_to_model(b_ptr + i,m);
    }
}
DebugVisualizer& DebugVisualizer::operator=(const DebugVisualizer& dv)
{
    tree_tex = dv.tree_tex;
    leaves_tex = dv.leaves_tex;
    tree_shader = dv.tree_shader;
}
void DebugVisualizer::visualize_light_voxels(LightVoxelsCube *voxels, glm::vec3 shift, glm::vec3 scale)
{
    visualize_light_voxels(voxels,
                           voxels->get_center() - voxels->get_voxel_size()*glm::vec3(voxels->get_vox_sizes()),
                           voxels->get_voxel_size()*(2.0f*glm::vec3(voxels->get_vox_sizes()) + glm::vec3(1)),
                           2.0f*glm::vec3(voxels->get_voxel_size()),
                           0.5f*voxels->get_voxel_size(),
                           0.1,
                           shift,
                           scale);


}
void DebugVisualizer::visualize_light_voxels(LightVoxelsCube *voxels,glm::vec3 pos, glm::vec3 size, glm::vec3 step,
                                             float dot_size, float threshold, glm::vec3 shift, glm::vec3 scale, int mip)
{
    int count = ((int)(size.x/step.x)) * ((int)(size.y/step.y)) * ((int)(size.z/step.z));
    Model *m = new Model();
    debugModels.push_back(m);
    currentModes.push_back(2);

    for (float x = pos.x; x < pos.x + size.x; x += step.x)
    {
        for (float y = pos.y; y < pos.y + size.y; y += step.y)
        {
            for (float z = pos.z; z < pos.z + size.z; z += step.z)
            {
                float occ = voxels->get_occlusion_simple_mip(glm::vec3(x,y,z),mip);
                if (occ < threshold || occ > 1e8)
                    continue;
                glm::vec4 tex;
                tex.w = 1;
                tex.z = MIN(1,occ/(10*threshold));
                tex.y = MIN(1,occ/(100*threshold));
                tex.x = MIN(1,occ/(1000*threshold));
                Box b = Box(shift + glm::vec3(scale.x*x,scale.y*y,scale.z*z),glm::vec3(dot_size,0,0),glm::vec3(0,dot_size,0),glm::vec3(0,0,dot_size));
                body_to_model(&b,m,true,tex);
            }
        }
    }
}

Model *visualizer::visualize_light_voxels(LightVoxelsCube *voxels, glm::vec3 shift, glm::vec3 scale)
{
  return visualize_light_voxels(voxels,
                                voxels->get_center() - voxels->get_voxel_size() * glm::vec3(voxels->get_vox_sizes()),
                                voxels->get_voxel_size() * (2.0f * glm::vec3(voxels->get_vox_sizes()) + glm::vec3(1)),
                                2.0f * glm::vec3(voxels->get_voxel_size()),
                                0.5f * voxels->get_voxel_size(),
                                0.1,
                                shift,
                                scale);
}
Model *visualizer::visualize_light_voxels(LightVoxelsCube *voxels, glm::vec3 pos, glm::vec3 size, glm::vec3 step,
                                          float dot_size, float threshold, glm::vec3 shift, glm::vec3 scale, int mip)
{
  int count = ((int)(size.x / step.x)) * ((int)(size.y / step.y)) * ((int)(size.z / step.z));
  Model *m = new Model();
  Visualizer vis;
  for (float x = pos.x; x < pos.x + size.x; x += step.x)
  {
    for (float y = pos.y; y < pos.y + size.y; y += step.y)
    {
      for (float z = pos.z; z < pos.z + size.z; z += step.z)
      {
        float occ = voxels->get_occlusion_simple_mip(glm::vec3(x, y, z), mip);
        if (occ < threshold || occ > 1e8)
          continue;
        glm::vec4 tex;
        tex.w = 1;
        tex.z = MIN(1, occ / (10 * threshold));
        tex.y = MIN(1, occ / (100 * threshold));
        tex.x = MIN(1, occ / (1000 * threshold));
        Box b = Box(shift + glm::vec3(scale.x * x, scale.y * y, scale.z * z), glm::vec3(dot_size, 0, 0), glm::vec3(0, dot_size, 0), glm::vec3(0, 0, dot_size));
        vis.body_to_model(&b, m, true, tex);
      }
    }
  }
  logerr("created model %d", m->positions.size());
  m->update();
  return m;
}

Model *visualizer::visualize_aabb(AABB &box, glm::vec3 &color)
{
  ::std::vector<AABB> boxes = {box};
  ::std::vector<glm::vec3> colors = {color};
  return visualize_aabb(boxes, colors);
}

Model *visualizer::visualize_aabb(::std::vector<AABB> &boxes, ::std::vector<glm::vec3> &colors)
{
  if (colors.size() == boxes.size())
  {
    Model *m = new Model();
    Visualizer vis;
    for (int i=0;i<colors.size();i++)
    {
      glm::vec3 sz = boxes[i].max_pos - boxes[i].min_pos;
      Box *box_ptr = new Box(boxes[i].min_pos, glm::vec3(sz.x,0,0), glm::vec3(0,sz.y,0), glm::vec3(0,0,sz.z));
      vis.body_to_model(box_ptr, m, true, glm::vec4(colors[i], 1));
      delete box_ptr;
    }
    m->update();
    return m;
  }
  else
  {
    logerr("visualize_aabb colors and boxes arrays should have the same size");
  }
}