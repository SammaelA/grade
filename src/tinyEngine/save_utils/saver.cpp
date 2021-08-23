#include "saver.h"
#include "../../texture_manager.h"
#include "../../grove.h"
namespace saver
{
std::string textures_path = ".";  
void set_textures_path (std::string &tp) {textures_path = tp;} 
bool save(FILE *f, std::string &t)
{
            try
        {
            short k = t.size();
            int c = fwrite(&k, sizeof(short), 1, f);
            bool ok = (c == 1);
            if (k > 0)
            {
                c = fwrite(t.c_str(), sizeof(char), k, f);
                ok = ok && (c == k);
            }
            return ok;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
}
bool load(FILE *f, std::string &t)
{
        try
        {
            short k;
            int c = fread(&k, sizeof(short), 1, f);
            bool ok = (c == 1);
            if (k > 0 && ok)
            {
                char *ptr = new char[k];
                c = fread(ptr, sizeof(char), k, f);
                ok = ok && (c == k);
                t.clear();
                t = std::string(ptr);
                delete[] ptr;
            }
            return ok;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return false;
        }
}
bool save(FILE *f, std::vector<GLfloat> &t)
{
    return savesv(f,t);
}
bool load(FILE *f, std::vector<GLfloat> &t)
{
    return loadsv(f,t);
}
    bool save(FILE *f, std::vector<PackedBranch> &t)
    {
        return savev(f,t);
    }
    bool load(FILE *f, std::vector<PackedBranch> &t)
    {
        return loadv(f,t);
    }
bool save(FILE *f, PackedLeaf &t)
{
    bool status = true;
    status = status && savesv(f,t.edges); 
    return status;
}
bool load(FILE *f, PackedLeaf &t)
{
    bool status = true;
    status = status && loadsv(f,t.edges); 
    return status;
}

bool save(FILE *f, PackedJoint &t)
{
    bool status = true;
    status = status && save(f,t.pos) && save(f,t.r); 
    return status;
}
bool load(FILE *f, PackedJoint &t)
{
    bool status = true;
    status = status && load(f,t.pos) && load(f,t.r); 
    return status;
}

bool save(FILE *f, PackedBranch &t)
{
    bool status = true;
    status = status && savev(f,t.joints) && savev(f,t.leaves) &&
             save(f,t.level) && save(f,t.plane_coef) && savev(f,t.r_mults)
             && save(f,t.type_id); 
    return status;
}
bool load(FILE *f, PackedBranch &t)
{
    bool status = true;
    status = status && loadv(f,t.joints) && loadv(f,t.leaves) &&
             load(f,t.level) && load(f,t.plane_coef) && loadv(f,t.r_mults)
             && load(f,t.type_id); 
    return status;
}

bool save(FILE *f, BranchCatalogue &t)
{
    bool status = true;
    status = status && savev(f,t.branches);
    return status;
}
bool load(FILE *f, BranchCatalogue &t)
{
    bool status = true;
    status = status && loadv(f,t.branches);
    return status;
}

bool save(FILE *f, BranchStructure &t)
{
    bool status = true;
    status = status && savev(f,t.childBranchesInstanced) && savev(f,t.childBranches) && save(f,t.pos);
    return status;
}
bool load(FILE *f, BranchStructure &t)
{
    bool status = true;
    status = status && loadv(f,t.childBranchesInstanced) && loadv(f,t.childBranches) && load(f,t.pos);
    return status;
}

bool save(FILE *f, InstancedBranch &t)
{
    bool status = true;
    status = status && savesv(f,t.branches) && save(f,t.IDA);
    return status;
}
bool load(FILE *f, InstancedBranch &t)
{
    bool status = true;
    status = status && loadsv(f,t.branches) && load(f,t.IDA);
    return status;
}

bool save(FILE *f, Texture &t)
{
    //it is not actual saving of texture - only its id
    unsigned char *pixels = new unsigned char[4*t.W*t.H*t.layers];
    glBindTexture(t.type,t.texture);
    glGetTexImage(t.type,0,GL_RGBA,GL_UNSIGNED_BYTE,pixels);
    std::string name = textures_path + "/atlas" + std::to_string(t.texture) + ".raw";
    FILE *tex = fopen(name.c_str(), "wb"); 
    fwrite(pixels,sizeof(unsigned char),4*t.W*t.H*t.layers,tex);
    fclose(tex);
    delete[] pixels;
    bool status = true;
    status = status && save(f,t.texture) && save(f,t.type) && save(f,t.W) &&
             save(f,t.H) && save(f,t.layers);
    return status;
}
bool load(FILE *f, Texture &t)
{
    bool status = true;
    status = status && load(f,t.texture) && load(f,t.type) && load(f,t.W) &&
             load(f,t.H) && load(f,t.layers);

    unsigned char *pixels = new unsigned char[4*t.W*t.H*t.layers];
    std::string name = textures_path + "/atlas" + std::to_string(t.texture) + ".raw";
    FILE *tex = fopen(name.c_str(), "rb"); 
    int c = fread(pixels,sizeof(unsigned char),4*t.W*t.H*t.layers,tex);
    status = status && 4*t.W*t.H*t.layers;

    if (status && t.type == GL_TEXTURE_2D)
        t = textureManager.load_unnamed(t,pixels);
    else if (status && t.type == GL_TEXTURE_2D_ARRAY)
    {
        t = textureManager.load_unnamed_arr(t,pixels);
    }
    else
    {
        status = false;
        t = textureManager.empty();
    }
    fclose(tex);
    delete[] pixels;
    logerr("loaded texture %d",t.texture);
    return status;
}

bool save(FILE *f, TextureAtlas &t)
{
    bool status = true;
    status = status && save(f,t.clearColor) && save(f,t.colorTex) && save(f,t.curNum) &&
             save(f,t.height) && save(f,t.width) && save(f,t.layers) && save(f,t.gridHN) &&
             save(f,t.gridWN) && save(f,t.isGrid);

    return status;
}
bool load(FILE *f, TextureAtlas &t)
{
    bool status = true;
    status = status && load(f,t.clearColor) && load(f,t.colorTex) && load(f,t.curNum) &&
             load(f,t.height) && load(f,t.width) && load(f,t.layers) && load(f,t.gridHN) &&
             load(f,t.gridWN) && load(f,t.isGrid);

    return status;   
}

bool save(FILE *f, InstanceDataArrays &t)
{
    bool status = true;
    status = status && savesv(f,t.centers_par) && savesv(f,t.centers_self) && savev(f,t.transforms) && 
                       savesv(f,t.type_ids) && savesv(f, t.tree_ids);
    return status;
}
bool load(FILE *f, InstanceDataArrays &t)
{
        bool status = true;
    status = status && loadsv(f,t.centers_par) && loadsv(f,t.centers_self) && loadv(f,t.transforms) &&
                       loadsv(f,t.type_ids) && loadsv(f, t.tree_ids);
    return status;
}

bool save(FILE *f, Billboard &t)
{
    bool status = true;
    status = status && save(f,t.branch_id) && save(f,t.id) && save(f,t.instancing) &&
            save(f,t.planeCoef) && savev(f,t.positions);
    return status;
}
bool load(FILE *f, Billboard &t)
{
    bool status = true;
    status = status && load(f,t.branch_id) && load(f,t.id) && load(f,t.instancing) &&
            load(f,t.planeCoef) && loadv(f,t.positions);
    return status;
}

bool save(FILE *f, BillboardData &t)
{
    bool status = true;
    status = status && save(f,t.IDA) && savev(f,t.billboards) && save(f,t.id) && save(f, t.base_position);
    return status;
}
bool load(FILE *f, BillboardData &t)
{
    bool status = true;
    status = status && load(f,t.IDA) && loadv(f,t.billboards) && load(f,t.id) && load(f, t.base_position);
    return status;
}

bool save(FILE *f, BillboardCloudData &t)
{
    bool status = true;
    status = status && save(f,t.atlas) && savelst(f,t.billboards) && save(f,t.level) && save(f,t.valid);
    return status;
}
bool load(FILE *f, BillboardCloudData &t)
{
    bool status = true;
    status = status && load(f,t.atlas) && loadlst(f,t.billboards) && load(f,t.level) && load(f,t.valid);
    return status;
}

bool save(FILE *f, GrovePacked &t)
{
    bool status = true;
    status = status && save(f,t.center) && savev(f,t.clouds) && savev(f,t.impostors) && savelst(f,t.instancedBranches) &&
             save(f,t.instancedCatalogue) && save(f,t.ggd_name);
    return status;
}
bool load(FILE *f, GrovePacked &t)
{
    bool status = true;
    status = status && load(f,t.center) && loadv(f,t.clouds) &&loadv(f,t.impostors) && loadlst(f,t.instancedBranches) &&
             load(f,t.instancedCatalogue) && load(f,t.ggd_name);
    return status; 
}
bool save(FILE *f, ImpostorsData &t)
{
    bool status = true;
    status = save(f, t.atlas) && save(f, t.level) && save(f, t.valid) && savev(f, t.impostors);
    return status;
}
bool load(FILE *f, ImpostorsData &t)
{
    bool status = true;
    status = load(f, t.atlas) && load(f, t.level) && load(f, t.valid) && loadv(f, t.impostors);
    return status;
}

bool save(FILE *f, Impostor &t)
{
    bool status = true;
    status = save(f, t.id) && save(f, t.bcyl) && save(f, t.IDA) && savev(f, t.slices) && save(f, t.top_slice);
    return status;
}
bool load(FILE *f, Impostor &t)
{
    bool status = true;
    status = load(f, t.id) && load(f, t.bcyl) && load(f, t.IDA) && loadv(f, t.slices) && load(f, t.top_slice);
    return status;
}
}