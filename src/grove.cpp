#include "grove.h"
#include "visualizer.h"
#include "texture_manager.h"
#include "billboard_cloud.h"
GroveRenderer::GroveRenderer(GrovePacked *_source, int LODs_count):
renderer({"default.vs", "default.fs"}, {"in_Position", "in_Normal", "in_Tex"}),
leaf(textureManager.get("leaf")),
wood(textureManager.get("wood"))
{
    Visualizer v = Visualizer();
    source = _source;
    debug("creating grove renderer with %d LODs\n",_source->clouds.size());
    for (int i=0;i<_source->clouds.size();i++)
    {
        LODs.emplace_back();
        LODs.back().cloud = _source->clouds[i];
        Model *m = new Model();
        for (int j=0;j<i;j++)
        {
            auto packed_branches = _source->branchCatalogue.get_level(j);
            for (PackedBranch &pb : packed_branches)
            {
                v.packed_branch_to_model(pb,m,false);
            }
        }
        LODs.back().m = m;
    }
}
void GroveRenderer::render(int lod, glm::mat4 prc)
{
    if (lod < 0 || lod >= LODs.size())
    {
        logerr("trying to render grove with wrong LOD number %d. Grove has %d LODs",lod, LODs.size());
        return;
    }
    Model *m = LODs[lod].m;
    renderer.use();
	renderer.uniform("projectionCamera", prc);
    renderer.texture("tex",wood);
	renderer.uniform("wireframe", false);
    renderer.uniform("model", m->model);
	m->update();
	m->render(GL_TRIANGLES);

    LODs[lod].cloud->render(prc);
}