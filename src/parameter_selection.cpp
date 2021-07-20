#include "parameter_selection.h"
void ParameterSelector::select(TreeStructureParameters &param, SelectionType sel_type, Metric metric)
{
    if (sel_type == BruteForce && metric == CompressionRatio)
    {
        std::vector<ParameterDesc> mask;
        std::vector<double> data;
        std::vector<double> data_max;
        param.get_mask_and_data(mask, data);
        GrovePacked tree;
        float max_compression = 0;
        for (int i = 0;i<10;i++)
        {
            GrovePacked tree = GrovePacked();
            data[0]--;
            param.load_from_mask_and_data(mask,data);
            generate(param, tree);
            int origins = 0;
            int instances = tree.instancedBranches.size();
            for (int j = 0;j < tree.instancedCatalogue.levels(); j++)
            {
                origins += tree.instancedCatalogue.get_level(j).size();
            }
            float comp = (float)instances/origins;
            logerr("instanced branches %d/%d = %f",instances,origins,comp);
            if (comp > max_compression)
            {
                max_compression = comp;
                data_max = data;
            }
        }
        param.load_from_mask_and_data(mask,data);
        logerr("parameter selection finished max_compression %f with values %f",max_compression,data[0]);
    } 
}