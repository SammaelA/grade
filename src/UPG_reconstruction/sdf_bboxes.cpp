#include "sdf_bboxes.h"


namespace upg
{
  void bbox_borders_search(unsigned int hi, unsigned int hj, unsigned int hk,
                           float& max_i, float& max_j, float& max_k,
                           bool& ind_i, bool& ind_j, bool& ind_k,
                           const std::vector<std::vector<std::vector<unsigned int>>> &is_marked_3D);

  /*
    the argument 'cuts_quant' in the 'get_bbox_list' (see below) function shows
    how many time we want to cut each bbox in half along some axis
  */
  std::vector<AABB> get_bbox_list(std::function<float(const float3 &)> sdf,
                                  const AABB &sdf_bbox, int bbox_count, int cuts_quant){
    std::vector<AABB> bbox_list;
    
    int bbox_segment_quant = pow(2, cuts_quant);
    float3 size_bbox = sdf_bbox.size();
    float3 segment_bbox = size_bbox / bbox_segment_quant;

    // Thanks to the point_xyz variable, we move between the "bboxes"
    float3 point_xyz = segment_bbox / 2;
    
    int semidiag_len = sqrt(pow(point_xyz.x, 2) + pow(point_xyz.y, 2) +
                            pow(point_xyz.z, 2));
    
    std::vector<AABB> bbox_list_1D;
    std::vector<std::vector<AABB>> bbox_list_2D;
    std::vector<std::vector<std::vector<AABB>>> bbox_list_3D;

    std::vector<float3> max_bbox_1D;
    std::vector<std::vector<float3>> max_bbox_2D;
    std::vector<std::vector<std::vector<float3>>> max_bbox_3D;

    std::vector<float3> min_bbox_1D;
    std::vector<std::vector<float3>> min_bbox_2D;
    std::vector<std::vector<std::vector<float3>>> min_bbox_3D;

    std::vector<unsigned int> is_marked_1D;
    std::vector<std::vector<unsigned int>> is_marked_2D;
    std::vector<std::vector<std::vector<unsigned int>>> is_marked_3D;
    /*
      values in is_marked_nD vectors:
        0 - an empty bbox

        3 - the bbox is in contact with three others
        2 - the bbox is in contact with other two
        1 - the bbox is in contact with another

        4 - unvisited bbox
        5 - visited bbox
    */

    for(unsigned int axis_x = 0; axis_x < bbox_segment_quant; axis_x += 1){
      for(unsigned int axis_y = 0; axis_y < bbox_segment_quant; axis_y += 1){
        int ind = 0;
        for(unsigned int axis_z = 0; axis_z < bbox_segment_quant; axis_z += 1){
          bbox_list_1D.push_back(AABB(point_xyz - segment_bbox / 2, point_xyz +
                                      segment_bbox / 2));
          max_bbox_1D.push_back(point_xyz + segment_bbox / 2);
          min_bbox_1D.push_back(point_xyz - segment_bbox / 2);

          if (sdf(point_xyz) >= semidiag_len & ind == 0)
            is_marked_1D.push_back(4);
          else
            is_marked_1D.push_back(0);

          point_xyz.z += segment_bbox.z;
        }

        is_marked_1D.push_back(0);
        is_marked_1D.insert(is_marked_1D.begin(), 0);

        bbox_list_2D.push_back(bbox_list_1D);
        max_bbox_2D.push_back(max_bbox_1D);
        min_bbox_2D.push_back(min_bbox_1D);
        is_marked_2D.push_back(is_marked_1D);

        bbox_list_1D.erase(bbox_list_1D.begin(), bbox_list_1D.end());
        max_bbox_1D.erase(max_bbox_1D.begin(), max_bbox_1D.end());
        min_bbox_1D.erase(min_bbox_1D.begin(), min_bbox_1D.end());
        is_marked_1D.erase(is_marked_1D.begin(), is_marked_1D.end());

        point_xyz.z = 0;
        point_xyz.y += segment_bbox.y;
      }

      is_marked_2D.push_back(is_marked_2D[0]);
      is_marked_2D.insert(is_marked_2D.begin(), is_marked_2D[0]);
      is_marked_2D[0].assign(is_marked_2D[0].size(), 0);
      is_marked_2D[is_marked_2D.size() - 1].assign(is_marked_2D[0].size(), 0);


      bbox_list_3D.push_back(bbox_list_2D);
      max_bbox_3D.push_back(max_bbox_2D);
      min_bbox_3D.push_back(min_bbox_2D);
      is_marked_3D.push_back(is_marked_2D);

      bbox_list_2D.erase(bbox_list_2D.begin(), bbox_list_2D.end());
      max_bbox_2D.erase(max_bbox_2D.begin(), max_bbox_2D.end());
      min_bbox_2D.erase(min_bbox_2D.begin(), min_bbox_2D.end());
      is_marked_2D.erase(is_marked_2D.begin(), is_marked_2D.end());

      point_xyz.y = 0;
      point_xyz.x += segment_bbox.x;
    }

    std::vector<unsigned int> row_of_zeros(is_marked_3D[0][0].size());
    is_marked_3D.push_back(is_marked_3D[0]);
    is_marked_3D.insert(is_marked_3D.begin(), is_marked_3D[0]);
    is_marked_3D[0].assign(is_marked_3D[0].size(), row_of_zeros);
    is_marked_3D[is_marked_3D.size() - 1].assign(is_marked_3D[0].size(), row_of_zeros);

    
    // Looking for corners
    unsigned int touches_quantity;
    for(unsigned int i = 0; i < bbox_list_3D.size(); i++){
      for(unsigned int j = 0; j < bbox_list_3D[0].size(); j++){
        for(unsigned int k = 0; k < bbox_list_3D[0][0].size(); k++){
          if(is_marked_3D[i][j][k] == 4){
            touches_quantity = is_marked_3D[i + 1][j][k] + is_marked_3D[i - 1][j][k] + \
                               is_marked_3D[i][j + 1][k] + is_marked_3D[i][j - 1][k] + \
                               is_marked_3D[i][j][k + 1] + is_marked_3D[i][j][k - 1];
            touches_quantity /= 4;

            if(touches_quantity > 3)
                continue;
            else if(touches_quantity == 1)
                is_marked_3D[i][j][k] = 1;

            else{
              if(is_marked_3D[i][j + 1][k] == 4 && is_marked_3D[i][j - 1][k] == 4)
                continue;
              else if(is_marked_3D[i][j][k + 1] == 4 && is_marked_3D[i][j][k - 1] == 4)
                continue;
              else if(is_marked_3D[i + 1][j][k] == 4 && is_marked_3D[i - 1][j][k] == 4)
                continue;

              else{
                if(touches_quantity == 3)
                  is_marked_3D[i][j][k] = 3;
                else
                  is_marked_3D[i][j][k] = 2;
              }
            }

          }
        }
      }
    }

    bool ind_i, ind_j, ind_k; // 1 - it is possible to expand the box along the x axis
                              // 0 - it is impossible
    bool ind; // 0 - it is possible to combine the boxes
              // 1 - it is not yet possible
    float max_i, max_j, max_k;
    float max_size = std::max(std::max(is_marked_3D.size(), is_marked_3D[0].size()),
                         is_marked_3D[0][0].size());
    
    for(unsigned int i = 1; i < is_marked_3D.size() - 1; i++){
      for(unsigned int j = 1; j < is_marked_3D[0].size() - 1; j++){
        for(unsigned int k = 1; k < is_marked_3D[0][0].size() - 1; k++){
          if(is_marked_3D[i][j][k] == 3){
            ind_i = 1, ind_j = 1, ind_k = 1;
            ind = 1;
            for(unsigned int h = 1; h < max_size; h++){
              if(ind_i == 1)
                max_i = i + h;
              if(ind_j == 1)
                max_j = j + h;
              if(ind_k == 1)
                max_k = k + h;

              if(ind == 1){
                for(unsigned int hj = max_j; hj >= 0; hj--){
                  for(unsigned int hk = max_k; hk >= 0; hk--)
                    bbox_borders_search(0, hj, hk, max_i, max_j, max_k,
                                        ind_i, ind_j, ind_k, is_marked_3D);
                  if(ind == 0)
                    break;
                }
                for(unsigned int hi = max_i; hi >= 0; hi--){
                  for(unsigned int hk = max_k; hk >= 0; hk--)
                    bbox_borders_search(hi, 0, hk, max_i, max_j, max_k,
                                        ind_i, ind_j, ind_k, is_marked_3D);
                  if(ind == 0)
                    break;
                }
                for(unsigned int hi = max_i; hi >= 0; hi--){
                  for(unsigned int hj = max_k; hj >= 0; hj--)
                    bbox_borders_search(hi, hj, 0, max_i, max_j, max_k,
                                        ind_i, ind_j, ind_k, is_marked_3D);
                  if(ind == 0)
                    break;
                }
              }
              else{
                for(unsigned int hi = i; hi <= max_i; hi++)
                  for(unsigned int hj = j; hj <= max_j; hj++)
                    for(unsigned int hk = k; hk <= max_k; hk++)
                      is_marked_3D[hi][hj][hk] = 5;
                
                bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                         max_bbox_3D[max_i - 1][max_j - 1][max_k - 1]));

              }
            }
          }

          else if(is_marked_3D[i][j][k] == 2){
            // fixing the i axis
            if((is_marked_3D[i][j + 1][k] != 0 && is_marked_3D[i][j + 1][k] != 1 && is_marked_3D[i][j + 1][k] != 5) ||
               (is_marked_3D[i][j][k + 1] != 0 && is_marked_3D[i][j][k + 1] != 1 && is_marked_3D[i][j][k + 1] != 5) ||
               (is_marked_3D[i][j + 1][k + 1] != 0 && is_marked_3D[i][j + 1][k + 1] != 1 && is_marked_3D[i][j + 1][k + 1] != 5)){
              ind_j = 0, ind_k = 0;
              ind = 1;

              for(unsigned int h = 0; h < std::max(is_marked_3D[0].size(),
                  is_marked_3D[0][0].size()); h++){
                if(ind_j == 1)
                  max_j = j + h;
                if(ind_k == 1)
                  max_k = k + h;

                if(ind == 1){
                  // fixing max_j
                  for(unsigned int hk = max_k; hk >= 0; hk--){
                    if(is_marked_3D[i][max_j][max_k - hk] == 1 || 
                       is_marked_3D[i][max_j][max_k - hk] == 0 ||
                       is_marked_3D[i][max_j][max_k - hk] == 4){
                      if(ind_j == 1){
                        max_j = max_j - 1;
                        ind_j = 0;
                      }
                      else if(ind_k == 1){
                        max_k = max_k - 1;
                        ind_k = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }

                  // fixing max_k
                  for(unsigned int hj = max_j; hj >= 0; hj--){
                    if(is_marked_3D[i][max_j - hj][max_k] == 1 || 
                       is_marked_3D[i][max_j - hj][max_k] == 0 ||
                       is_marked_3D[i][max_j - hj][max_k] == 4){
                      if(ind_j == 1){
                        max_j = max_j - 1;
                        ind_j = 0;
                      }
                      else if(ind_k == 1){
                        max_k = max_k - 1;
                        ind_k = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }
                }
                else{
                  for(unsigned int hj = 0; hj < max_j; hj++)
                    for(unsigned int hk = 0; hk < max_k; hk++)
                      is_marked_3D[i][hj][hk] = 5;
                
                  bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                          max_bbox_3D[i - 1][max_j - 1][max_k - 1]));
                }
              }
            }

            // fixing the j axis
            else if((is_marked_3D[i + 1][j][k] != 0 && is_marked_3D[i + 1][j][k] != 1 && is_marked_3D[i + 1][j][k] != 5) ||
                    (is_marked_3D[i][j][k + 1] != 0 && is_marked_3D[i][j][k + 1] != 1 && is_marked_3D[i][j][k + 1] != 5) ||
                    (is_marked_3D[i + 1][j][k + 1] != 0 && is_marked_3D[i + 1][j][k + 1] != 1 && is_marked_3D[i + 1][j][k + 1] != 5)){
              ind_i = 0, ind_k = 0;
              ind = 1;

              for(unsigned int h = 0;
                  h < std::max(is_marked_3D.size(), is_marked_3D[0][0].size()); h++){
                if(ind_i == 1)
                  max_i = i + h;
                if(ind_k == 1)
                  max_k = k + h;

                if(ind == 1){
                  // fixing max_i
                  for(unsigned int hk = max_k; hk >= 0; hk--){
                    if(is_marked_3D[max_i][j][max_k - hk] == 1 || 
                       is_marked_3D[max_i][j][max_k - hk] == 0 ||
                       is_marked_3D[max_i][j][max_k - hk] == 4){
                      if(ind_i == 1){
                        max_i = max_i - 1;
                        ind_i = 0;
                      }
                      else if(ind_k == 1){
                        max_k = max_k - 1;
                        ind_k = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }

                  // fixing max_k
                  for(unsigned int hi = max_i; hi >= 0; hi--){
                    if(is_marked_3D[max_i - hi][j][max_k] == 1 || 
                       is_marked_3D[max_i - hi][j][max_k] == 0 ||
                       is_marked_3D[max_i - hi][j][max_k] == 4){
                      if(ind_i == 1){
                        max_i = max_i - 1;
                        ind_i = 0;
                      }
                      else if(ind_k == 1){
                        max_k = max_k - 1;
                        ind_k = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }
                }
                else{
                  for(unsigned int hi = 0; hi < max_i; hi++)
                    for(unsigned int hk = 0; hk < max_k; hk++)
                      is_marked_3D[hi][j][hk] = 5;
                
                  bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                          max_bbox_3D[max_i - 1][j - 1][max_k - 1]));
                }
              }
            }

            // fixing the k axis
            else if((is_marked_3D[i + 1][j][k] != 0 && is_marked_3D[i + 1][j][k] != 1 && is_marked_3D[i + 1][j][k] != 5) ||
                    (is_marked_3D[i][j + 1][k] != 0 && is_marked_3D[i][j + 1][k] != 1 && is_marked_3D[i][j + 1][k] != 5) ||
                    (is_marked_3D[i + 1][j + 1][k] != 0 && is_marked_3D[i + 1][j + 1][k] != 1 && is_marked_3D[i + 1][j + 1][k] != 5)){
              ind_i = 0, ind_k = 0;
              ind = 1;
  
              for(unsigned int h = 0;
                  h < std::max(is_marked_3D.size(), is_marked_3D[0].size()); h++){
                if(ind_i == 1)
                  max_i = i + h;
                if(ind_j == 1)
                  max_j = j + h;

                if(ind == 1){
                  // fixing max_i
                  for(unsigned int hj = max_j; hj >= 0; hj--){
                    if(is_marked_3D[max_i][j - hj][k] == 1 || 
                       is_marked_3D[max_i][j - hj][k] == 0 ||
                       is_marked_3D[max_i][j - hj][k] == 4){
                      if(ind_i == 1){
                        max_i = max_i - 1;
                        ind_i = 0;
                      }
                      else if(ind_j == 1){
                        max_j = max_j - 1;
                        ind_j = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }

                  // fixing max_j
                  for(unsigned int hi = max_i; hi >= 0; hi--){
                    if(is_marked_3D[max_i - hi][max_j][k] == 1 || 
                       is_marked_3D[max_i - hi][max_j][k] == 0 ||
                       is_marked_3D[max_i - hi][max_j][k] == 4){
                      if(ind_i == 1){
                        max_i = max_i - 1;
                        ind_i = 0;
                      }
                      else if(ind_j == 1){
                        max_j = max_j - 1;
                        ind_j = 0;
                      }
                      else{
                        ind = 0;
                        break;
                      }
                    }
                  }
                }
                else{
                  for(unsigned int hi = 0; hi < max_i; hi++)
                    for(unsigned int hj = 0; hj < max_k; hj++)
                      is_marked_3D[hi][hj][k] = 5;
                
                  bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                          max_bbox_3D[max_i - 1][max_j - 1][k - 1]));
                }
              }
            }
          }

          else if(is_marked_3D[i][j][k] == 1){
            unsigned int coeff;

            if(is_marked_3D[i][j][k + 1] == 1){
              is_marked_3D[i][j][k] = 5;
              is_marked_3D[i][j][k + 1] = 5;

              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i - 1][j - 1][k]));
            }

            else if(is_marked_3D[i][j][k + 1] != 0 && is_marked_3D[i][j][k + 1] != 5){
              coeff = 0;
              while(is_marked_3D[i][j][k + coeff] != 0 && is_marked_3D[i][j][k + coeff] != 5){
                is_marked_3D[i][j][k + coeff] = 5;
                coeff++;
              }
              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i - 1][j - 1][k + coeff - 1]));
            }

            else if(is_marked_3D[i][j + 1][k] == 1){
              is_marked_3D[i][j][k] = 5;
              is_marked_3D[i][j + 1][k] = 5;

              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i - 1][j][k - 1]));
            }

            else if(is_marked_3D[i][j + 1][k] != 0 && is_marked_3D[i][j + 1][k] != 5){
              coeff = 0;
              while(is_marked_3D[i][j + coeff][k] != 0 && is_marked_3D[i][j + coeff][k] != 5){
                is_marked_3D[i][j + coeff][k] = 5;
                coeff++;
              }
              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i - 1][j + coeff - 1][k - 1]));
            }

            else if(is_marked_3D[i + 1][j][k] == 1){
              is_marked_3D[i][j][k] = 5;
              is_marked_3D[i + 1][j][k] = 5;

              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i][j - 1][k - 1]));
            }

            else if(is_marked_3D[i + 1][j][k] != 0 && is_marked_3D[i + 1][j][k] != 5){
              coeff = 0;
              while(is_marked_3D[i + coeff][j][k] != 0 && is_marked_3D[i + coeff][j][k] != 5){
                is_marked_3D[i + coeff][j][k] = 5;
                coeff++;
              }
              bbox_list.push_back(AABB(min_bbox_3D[i - 1][j - 1][k - 1],
                                       max_bbox_3D[i + coeff - 1][j - 1][k - 1]));
            }

          }
        }
      }
    }

    // Checking for the remaining boxes
    for(unsigned int i = 1; i < is_marked_3D.size() - 1; i++)
      for(unsigned int j = 1; j < is_marked_3D[0].size() - 1; j++)
        for(unsigned int k = 1; k < is_marked_3D[0][0].size() - 1; k++)
          if(is_marked_3D[i][j][k] != 0 && is_marked_3D[i][j][k] != 5){
            is_marked_3D[i][j][k] = 5;
            bbox_list.push_back(AABB(min_bbox_3D[i][j][k], max_bbox_3D[i][j][k]));
          }

    return bbox_list;
  }

  void bbox_borders_search(unsigned int hi, unsigned int hj, unsigned int hk,
                           float& max_i, float& max_j, float& max_k,
                           bool& ind_i, bool& ind_j, bool& ind_k,
                           const std::vector<std::vector<std::vector<unsigned int>>> &is_marked_3D){
    if(is_marked_3D[max_i - hi][max_j - hj][max_k - hk] != 3 &&
       is_marked_3D[max_i - hi][max_j - hj][max_k - hk] != 4){
      if(ind_i == 1){
        max_i -= 1;
        ind_i = 0;
      }
      else if(ind_j == 1){
        max_j -= 1;
        ind_j = 0;
      }
      else if(ind_k == 1){
        max_k -= 1;
        ind_k = 0;
      }
    }
  }
}