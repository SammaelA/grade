#include "sandbox.h"
#include "generation/scene_generator.h"
#include "generation/grove_packer.h"
#include "generation/metainfo_manager.h"
#include "graphics_utils/texture_manager.h"
#include "tinyEngine/image.h"
#include "parameter_selection/impostor_similarity.h"
#include "parameter_selection/genetic_algorithm.h"
#include "tree_generators/GE_generator.h"
#include "parameter_selection/parameter_selection.h"
#include <thread>
#include <chrono>
struct BS_Grid
{
    std::vector<std::pair<int, int>> bins;
    ParameterList base_set;
};

ParameterList grid_params(std::vector<std::pair<int, int>> bins, ParameterList &base_set)
{
    ParameterList params = base_set;
                        int k=0;
                    for (auto &p : params.categorialParameters)
                    {
                        if (!p.second.fixed())
                        {
                            p.second.val = p.second.possible_values[bins[k].first];
                            p.second.possible_values = {p.second.val};
                        }
                        k++;
                    }
                    for (auto &p : params.ordinalParameters)
                    {
                        if (!p.second.fixed())
                        {
                            p.second.val = p.second.min_val + bins[k].first; 
                            p.second.min_val = p.second.val;
                            p.second.max_val = p.second.val;
                        } 
                        k++;
                    }
                    for (auto &p : params.continuousParameters)
                    {
                        if (!p.second.fixed())
                        {
                            float step = (p.second.max_val - p.second.min_val)/bins[k].second;
                            float min = p.second.min_val;
                            p.second.val = min + (bins[k].first + urand())*step; 
                            p.second.min_val = min + (bins[k].first)*step;
                            p.second.max_val = min + (bins[k].first + 1)*step;
                        }
                        k++;
                    }
    return params;
}

void fill_bs_grid_bins(BS_Grid &grid, ParameterList &params, int num_bins)
{
    for (auto &p : params.categorialParameters)
    {
        if (!p.second.fixed())
            grid.bins.push_back(std::pair<int, int>(0, p.second.possible_values.size()));
        else
            grid.bins.push_back(std::pair<int, int>(0, 0));
    }
    for (auto &p : params.ordinalParameters)
    {
        if (!p.second.fixed())
            grid.bins.push_back(std::pair<int, int>(0, p.second.max_val - p.second.min_val));  
        else
            grid.bins.push_back(std::pair<int, int>(0, 0));
    }
    for (auto &p : params.continuousParameters)
    {
        if (!p.second.fixed())
            grid.bins.push_back(std::pair<int, int>(0, num_bins));  
        else
            grid.bins.push_back(std::pair<int, int>(0, 0));
    }
}

void bruteforce_selection(const std::function<std::vector<float>(std::vector<ParameterList> &)> &f, int num_bins, 
                          int detalization_count, int detalization_depth, int num_samples, float &best_val, 
                          ParameterList &bestParams)
{
    std::vector<BS_Grid> progress_bars;
    progress_bars.emplace_back();
    progress_bars[0].base_set = bestParams;
    fill_bs_grid_bins(progress_bars[0], bestParams, num_bins);

    for (int layer = 0;layer<detalization_depth;layer++)
    {
        if (layer == detalization_depth - 1)
            detalization_count = 1;

        std::vector<BS_Grid> new_progress_bars;
        for (auto &grid : progress_bars)
        {
            std::vector<std::pair<float, std::vector<std::pair<int, int>>>> cur_best;
            int i=0;
            while (i < grid.bins.size())
            {
                float sum_metric = 0;
                std::vector<ParameterList> params; 
                for (int sample = 0;sample<num_samples;sample++)
                {
                    params.push_back(grid_params(grid.bins, grid.base_set));
                }
                auto vec = f(params);
                for (auto &v : vec)
                {
                    sum_metric += v;
                    logerr("v = %f", v);
                }
                logerr("");
                sum_metric /= num_samples;
                if (cur_best.size() < detalization_count)
                {
                    cur_best.push_back(std::pair<float, std::vector<std::pair<int, int>>>(sum_metric, grid.bins));
                }
                else
                {
                    float min_val = 1e9;
                    int min_pos = 0;
                    for (int j = 0;j<cur_best.size();j++)
                    {
                        if (cur_best[j].first < min_val)
                        {
                            min_val = cur_best[j].first;
                            min_pos = j;
                        }
                    }
                    if (min_val < sum_metric)
                        cur_best[min_pos] = std::pair<float, std::vector<std::pair<int, int>>>(sum_metric, grid.bins);
                }

                if (grid.bins[i].first < grid.bins[i].second - 1)
                {
                    grid.bins[i].first++;
                }
                else
                {
                    i++;
                    while (i < grid.bins.size() && grid.bins[i].second <= MAX(1, grid.bins[i].first + 1))
                        i++;
                    if (i < grid.bins.size())
                    {
                        grid.bins[i].first++;
                        int t = -1;
                        for (int j=0;j<i;j++)
                        {
                            grid.bins[j].first = 0;
                            if (grid.bins[j].second > 1 && t<0)
                                t = j;
                        };
                        if (t >= 0)
                        i = t;
                    }
                }
                /*
                for (int j=0;j<grid.bins.size();j++)
                {
                    if (grid.bins[j].second > 1)
                        debug("%d/%d ", grid.bins[j].first, grid.bins[j].second);
                }
                debug("  i = %d\n",i);
                */
            }
            if (layer < detalization_depth - 1)
            {
                for (auto &best : cur_best)
                {
                    new_progress_bars.emplace_back();
                    new_progress_bars.back().base_set = grid_params(best.second, grid.base_set);
                    //new_progress_bars.back().base_set.print();
                    fill_bs_grid_bins(new_progress_bars.back(), new_progress_bars.back().base_set, num_bins);
                }
            }
            else
            {
                best_val = cur_best[0].first;
                bestParams = grid_params(cur_best[0].second, grid.base_set);
            }
        }
        progress_bars = std::move(new_progress_bars);
    }
}

float dot_metric(Tree &single_tree, float dst_dot)
{
    double sum_dot = 0;
    int dot_cnt = 0;

    for (auto &bh : single_tree.branchHeaps)
    {
        for (auto &b : bh->branches)
        {
            glm::vec3 dir = normalize(b.segments.front().begin - b.segments.front().end);
            for (auto &j : b.joints)
            {
                for (auto *chb : j.childBranches)
                {
                    glm::vec3 ch_dir = normalize(chb->segments.front().begin - chb->segments.front().end);
                    float w = pow(10, 3 - b.level);
                    sum_dot += w * SQR(0.5 - CLAMP(glm::dot(dir, ch_dir), 0, 1));
                    dot_cnt += w;
                }
            }
        }
    }
    float dt = sum_dot / dot_cnt;
    float metric = 1 - dt;
    logerr("dot metric %f %f", dt, metric);
    return metric;
}

void sandbox_main(int argc, char **argv, Scene &scene)
{
    metainfoManager.reload_all();
    TreeTypeData type = metainfoManager.get_tree_type("simpliest_tree_default");
    Block b, ref_info;
    BlkManager man;
    man.load_block_from_file("parameter_selection_settings.blk", b);
    man.load_block_from_file("parameter_selection_reference.blk", ref_info);
    ParameterSelector sel;
    //auto res = sel.parameter_selection(type, b, &scene);
    auto res = sel.parameter_selection(ref_info, b, &scene);
}