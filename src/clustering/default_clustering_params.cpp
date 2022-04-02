#include "default_clustering_params.h"
#include "common_utils/utility.h"

Block trunks_default_params;
Block branches_default_params;
Block trees_default_params;
Block dummy;
std::string default_params_name = "clustering_defaults.blk";
bool params_loaded = false;
ClusteringStep current_clustering_step;
Block &get_default_block()
{
    load_default_blocks();

    if (current_clustering_step == ClusteringStep::TRUNKS)
        return trunks_default_params;
    else if (current_clustering_step == ClusteringStep::BRANCHES)
        return branches_default_params;
    else if (current_clustering_step == ClusteringStep::TREES)
        return trees_default_params;
    else 
        return dummy;
}
void load_default_blocks()
{
    if (params_loaded)
        return;
    params_loaded = true;
    BlkManager man;
    Block base;
    man.load_block_from_file(default_params_name,base);
    Block *b = base.get_block("trunks_default_params");
    if (b)
    {
        trunks_default_params.copy(b);
    }
    else
    {
        logerr("trunks clustering params block was not found in %s",default_params_name.c_str());
    }

    b = base.get_block("branches_default_params");
    if (b)
    {
        branches_default_params.copy(b);
    }
    else
    {
        logerr("branches clustering params block was not found in %s",default_params_name.c_str());
    }

    b = base.get_block("trees_default_params");
    if (b)
    {
        trees_default_params.copy(b);
    }
    else
    {
        logerr("trees clustering params block was not found in %s",default_params_name.c_str());
    }
}