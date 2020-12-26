#include "debug_visualizer.h"
#include "generated_tree.h"
DebugVisualizer::DebugVisualizer(Texture *_tree_tex, Shader *_tree_shader)
{
    tree_tex = _tree_tex;
    tree_shader = _tree_shader;
}
void DebugVisualizer::render(glm::mat4 view_proj)
{
    if (!tree_shader)
        return;
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
void branch_to_model(Branch *b, int level, Model &m)
{
    Tree t;
    TreeGenerator gen = TreeGenerator(t);
    if (b->level<0)
        return;
    if (b->level == level || level < 0)
    {
        std::vector<SegmentVertexes> vets;
        int i = 0;
        int ringsize = 3 * pow(2, gen.curParams.max_depth() - b->level);
        for (auto &segment : b->segments)
        {
            SegmentVertexes vt;
            gen.get_base_ring(segment, vt, ringsize, (float)(i % 3) / 3);
            if (!vets.empty())
                vets.back().smallRing = vt.bigRing;
            //segment_to_model(segment,m,leaves);
            vets.push_back(vt);
            i++;
        }
        if (!vets.empty())
        {
            gen.get_last_seg_vertexes(b->segments.back(),vets.back(),ringsize,(float)(i % 3) / 3);
        }

        //seg_vertexes_to_model(vets.front(),m);
        for (auto &vt : vets)
        {
            gen.seg_vertexes_to_model(vt, &m);
        }
    }
    for (auto &joint : b->joints)
    {
        for (auto branch : joint.childBranches)
            branch_to_model(branch, level, m);
    }

    
}
void DebugVisualizer::add_branch(Branch *b, glm::vec3 scale, glm::vec3 shift, int level)
{
    Model *m = new Model();
    debugModels.push_back(m);
    currentModes.push_back(1);
    branch_to_model(b,level,*m);
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