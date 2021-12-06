#include "modeling.h"
#include "texture_manager.h"

Model *ModelLoader::create_model_from_block(Block &bl, Texture &tex)
{
    std::string name = bl.get_string("name", "debug_box");
    if (name == "debug_box")
    {
        tex = textureManager.get("noise");
        return create_debug_box_model();
    }
    else 
    {
        logerr("Unknows model name %s", name.c_str());
        tex = textureManager.get("noise");
        return create_debug_box_model();
    }
}

Model *ModelLoader::create_debug_box_model()
{
    Box b = Box(glm::vec3(0,0,0), glm::vec3(1,0,0), glm::vec3(0,1,0), glm::vec3(0,0,1));
    //Ellipsoid cyl = Ellipsoid(glm::vec3(0,0,0), glm::vec3(1,0,0), glm::vec3(0,1,0), glm::vec3(0,0,1));
    Visualizer v;
    Model *m = new Model;
    //v.ellipsoid_to_model(&cyl, m, 16, 16);
    v.box_to_model(&b, m);
    return m;
}