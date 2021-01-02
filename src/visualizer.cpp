#include "generated_tree.h"
#include "visualizer.h"
#include "tinyEngine/utility.h"
#define PI 3.14159265f

void Visualizer::leaf_to_model(Leaf &l, Model *m)
{
  if (l.dead)
    return;
  if (l.edges.size()<4)
  {
    return;
  }
  glm::vec3 a = l.edges[0];
  glm::vec3 b = l.edges[1];
  glm::vec3 c = l.edges[2];
  glm::vec3 n = glm::normalize(glm::cross(a-b,c-b));
  std::vector<float> tex_c{1,0,0,0,0,1,1,1};
  int _b = m->positions.size()/3;
  for (int i=0;i<4;i++)
  {
    glm::vec3 v = l.edges[i];
    m->positions.push_back(v.x);
    m->positions.push_back(v.y);
    m->positions.push_back(v.z);
    m->normals.push_back(n.x);
    m->normals.push_back(n.y);
    m->normals.push_back(n.z);
    m->colors.push_back(tex_c[2*i]);
    m->colors.push_back(tex_c[2*i+1]);
    m->colors.push_back(0);
    m->colors.push_back(1);
  }

  m->indices.push_back(_b);
  m->indices.push_back(_b+1);
  m->indices.push_back(_b+2);
  m->indices.push_back(_b+2);
  m->indices.push_back(_b+3);
  m->indices.push_back(_b);
}
void Visualizer::get_base_ring(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos)
{
  sv.ringsize = ring_size;
  glm::vec3 start = s.begin;
  glm::vec3 end   = s.end;
  glm::vec3 dir = end - start;
  
  //Get Some Normal Vector
  glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
  glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

  //Add the Correct Number of Indices
  glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)(2*PI/ring_size), dir);

  for (int i = 0; i < ring_size; i++)
  {
    VertexData vd;
    vd.pos = start + s.rel_r_begin*glm::vec3(n.x,n.y,n.z);
    vd.normal = n;
    vd.tex_coord.x = ((float)i)/ring_size;
    vd.tex_coord.y = rel_ring_pos;
    sv.bigRing.push_back(vd);
    n = r*n;
  }
}
void Visualizer::get_last_seg_vertexes(Segment &s, SegmentVertexes &sv, int ring_size, float rel_ring_pos)
{
  sv.ringsize = ring_size;
  glm::vec3 start = s.begin;
  glm::vec3 end   = s.end;
  glm::vec3 dir = end - start;
  
  //Get Some Normal Vector
  glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
  glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

  //Add the Correct Number of Indices
  glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)(2*PI/ring_size), dir);

  for (int i = 0; i < ring_size; i++)
  {
    VertexData vd;
    vd.pos = end + s.rel_r_end*glm::vec3(n.x,n.y,n.z);
    vd.normal = n;
    vd.tex_coord.x = ((float)i)/ring_size;
    vd.tex_coord.y = rel_ring_pos;
    sv.smallRing.push_back(vd);
    n = r*n;
  }
}
void Visualizer::seg_vertexes_to_model(SegmentVertexes &sv, Model *m)
{
  Model *h = m;
  int _b = h->positions.size()/3;
  if (sv.smallRing.size()<sv.ringsize || sv.bigRing.size()<sv.ringsize)
    return;
  for (auto pos: sv.smallRing)
      {
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        h->colors.push_back(pos.tex_coord.x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(0.0);
        h->colors.push_back(1.0);
      }
      for (auto pos: sv.bigRing)
      {
        h->positions.push_back(pos.pos.x);
        h->positions.push_back(pos.pos.y);
        h->positions.push_back(pos.pos.z);
        h->normals.push_back(pos.normal.x);
        h->normals.push_back(pos.normal.y);
        h->normals.push_back(pos.normal.z);
        h->colors.push_back(pos.tex_coord.x);
        h->colors.push_back(pos.tex_coord.y);
        h->colors.push_back(0.0);
        h->colors.push_back(1.0);
      }
  int ringsize = sv.ringsize;
  for(int i = 0; i < ringsize; i++)
  {
    //Bottom Triangle
    h->indices.push_back(_b+i);
    h->indices.push_back(_b+(i + 1)%(ringsize));
    h->indices.push_back(_b+i+ringsize);
    //Upper Triangle
    h->indices.push_back(_b+i+ringsize);
    h->indices.push_back(_b+(i + 1)%(ringsize));
    h->indices.push_back(_b+(i + 1)%(ringsize) + ringsize);
  }
}
void Visualizer::joint_to_model(Joint &j, Model *m, bool leaves)
{

}
void Visualizer::recursive_branch_to_model(Branch &b, Model *m, bool leaves)
{
    if (b.level<0)
        return;
    if (!leaves)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, curParams.max_depth() - b.level);
        for (auto &segment : b.segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b.segments.back(),vets.back(),ringsize,(float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m);
        }
    }
    else
    {
        for (auto &joint : b.joints)
        {
            if (joint.leaf)
                leaf_to_model(*(joint.leaf), m);
        }
    }
    for (auto &joint : b.joints)
    {
        for (auto branch : joint.childBranches)
            recursive_branch_to_model(*branch, m, leaves);
    }
}
void Visualizer::branch_to_model(Branch &b, Model *m, bool leaves)
{
    if (!leaves)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, curParams.max_depth() - b.level);
        for (auto &segment : b.segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b.segments.back(), vets.back(), ringsize, (float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            seg_vertexes_to_model(vt, m);
        }
    }
  else
  {
    for (auto &joint : b.joints)
        joint_to_model(joint,m,leaves);
  }
}
void Visualizer::segment_to_model(Segment &s, Model *m, bool leaves)
{
    
    glm::vec3 start = s.begin;
    glm::vec3 end   = s.end;
    glm::vec3 dir = end - start;
    int ringsize = 8;
    Model *h = m;
    //Get Some Normal Vector
    glm::vec3 x = glm::normalize(dir + glm::vec3(1.0, 1.0, 1.0));
    glm::vec4 n = glm::vec4(glm::normalize(glm::cross(dir, x)), 1.0);

    //Add the Correct Number of Indices
    glm::mat4 r = glm::rotate(glm::mat4(1.0), (float)3.141/ringsize, dir);

    //Index Buffer
    int _b = h->positions.size()/3;

    //GL TRIANGLES
    for(int i = 0; i < ringsize; i++){
      //Bottom Triangle
      h->indices.push_back(_b+i*2+0);
      h->indices.push_back(_b+(i*2+2)%(2*ringsize));
      h->indices.push_back(_b+i*2+1);
      //Upper Triangle
      h->indices.push_back(_b+(i*2+2)%(2*ringsize));
      h->indices.push_back(_b+(i*2+3)%(2*ringsize));
      h->indices.push_back(_b+i*2+1);
    }

    for(int i = 0; i < ringsize; i++){

      h->positions.push_back(start.x + s.rel_r_begin*n.x);
      h->positions.push_back(start.y + s.rel_r_begin*n.y);
      h->positions.push_back(start.z + s.rel_r_begin*n.z);
      h->normals.push_back(n.x);
      h->normals.push_back(n.y);
      h->normals.push_back(n.z);
      n = r*n;

      h->positions.push_back(end.x + s.rel_r_end*n.x);
      h->positions.push_back(end.y + s.rel_r_end*n.y);
      h->positions.push_back(end.z + s.rel_r_end*n.z);
      h->normals.push_back(n.x);
      h->normals.push_back(n.y);
      h->normals.push_back(n.z);
      n = r*n;

    }
}
void Visualizer::add_branch_layer(Tree &t, int layer, Model *m)
{
    if (layer >= 0 && layer < t.branchHeaps.size())
    {
        for (auto &branch : t.branchHeaps[layer]->branches)
        {
            if (!branch.dead)
                branch_to_model(branch, m, false);
        }
    }
}
Visualizer::Visualizer(Texture *_tree_tex, Texture *_leaves_tex, Shader *_tree_shader)
{
    tree_tex = _tree_tex;
    tree_shader = _tree_shader;
    leaves_tex = _leaves_tex;
}
DebugVisualizer::DebugVisualizer(Texture *_tree_tex, Shader *_tree_shader):
Visualizer(_tree_tex,nullptr,_tree_shader)
{
}
DebugVisualizer::~DebugVisualizer()
{
    for (int i=0;i<debugModels.size();i++)
    {
        delete(debugModels[i]);
    }
}
void DebugVisualizer::render(glm::mat4 view_proj)
{
    if (!tree_shader)
    {
        logerr("empty debug shader\n");
        return;
    }
    tree_shader->use();
    if (tree_tex)
        tree_shader->texture("tex",*tree_tex);
    tree_shader->uniform("projectionCamera", view_proj);
    
    for (int i=0;i<debugModels.size();i++)
    {
        if (currentModes[i]==0)
            continue;
        tree_shader->uniform("model", debugModels[i]->model);
        debugModels[i]->update();
        debugModels[i]->render(GL_TRIANGLES);
    }
}
void DebugVisualizer::branch_to_model_debug(Branch *b, int level, Model &m)
{
    if (b->level<0)
        return;
    if (b->level == level || level < 0)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, curParams.max_depth() - b->level);
        for (auto &segment : b->segments)
        {
            SegmentVertexes vt;
            get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            get_last_seg_vertexes(b->segments.back(),vets.back(),ringsize,(float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
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
    branch_to_model_debug(b,level,*m);
    m->shift(shift);
    m->scale = scale;
}
void DebugVisualizer::enable_all()
{
    for (int i=0;i<currentModes.size();i++)
        currentModes[i] = 1;
}
void DebugVisualizer::disable_all()
{
    for (int i=0;i<currentModes.size();i++)
        currentModes[i] = 0;
}