#include "ggd_parser.h"
#include "parser.h"
#include "config.h"
#include <malloc.h>
#include "../../tree.h"
#include "../../grove.h"
extern "C" int parse_file(char *file);
bool Config::load_ggds()
{

    GroveGenerationData def;
    def.size = glm::vec3(150,200,150);
    def.pos = glm::vec3(0,0,0);
    def.trees_count = 1;
    def.synts_count = 0;
    def.synts_precision = 1;
    def.types = {TreeTypeData(0,TreeStructureParameters(),std::string("wood"),std::string("leaf"))};
    def.name = "default";
    ggds.emplace("default",def);
    pd.presets_n = 0;
    dat.ggds_c = 0;
    parse_file("groves.txt");
    for (int i=0; i< dat.ggds_c;i++)
    {
        #define GETV(a) glm::vec3(a[0],a[1],a[2])
        ggd &cur_g = dat.ggds[i];
        GroveGenerationData gen;
        gen.trees_count = cur_g.count;
        gen.synts_count = cur_g.synts_count;
        gen.synts_precision = cur_g.synts_precision;
        gen.size = GETV(cur_g.size);
        gen.pos = GETV(cur_g.pos);
        for (int j=0;j<cur_g.obsts_c;j++)
        {
            obst &cur_o = cur_g.obsts[j];
            glm::vec3 pos = GETV(cur_o.pos);
            glm::vec3 a = GETV(cur_o.a);
            glm::vec3 b = GETV(cur_o.b);
            glm::vec3 c = GETV(cur_o.c);
            if (cur_o.type == 1)
                gen.obstacles.push_back(new Box(pos,a,b,c));
            else if (cur_o.type == 2) 
                gen.obstacles.push_back(new Ellipsoid(pos,a,b,c));
            else if (cur_o.type == 3) 
                gen.obstacles.push_back(new Cylinder(pos,a,b,c));
        }
        for (int j=0;j<cur_g.ttds_c;j++)
        {
            ttd &cur_tt = cur_g.ttds[j];
            TreeStructureParameters p = get(cur_tt.name);
            TreeTypeData ttd(cur_tt.id,p,cur_tt.wood,cur_tt.leaf);
            gen.types.push_back(ttd);
        }
        gen.name = cur_g.name;
        ggds.emplace(cur_g.name,gen);
    }
    debugl(10,"loaded %d grove data's",ggds.size());
}