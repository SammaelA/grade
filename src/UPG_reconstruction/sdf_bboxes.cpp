#include "sdf_bboxes.h"
namespace upg
{
  /*
    the argument 'cuts_quant' in the 'get_bbox_list' (see below) function shows
    how many time we want to cut each bbox in half along some axis
  */
  std::vector<AABB> get_bbox_list(std::function<float(const float3 &)> sdf,
                                  const AABB &sdf_bbox, int bbox_count, int cuts_quant){
    std::vector<AABB> bbox_list;

    int bbox_segment_quant = pow(2, cuts_quant);
    float3 size_bbox = sdf_bbox.size();
    float3 segment_len = size_bbox / bbox_segment_quant;
    float3 semi_segment_len = segment_len / 2;

    float3 point_xyz = sdf_bbox.min_pos + semi_segment_len;

    float semidiag_len = sqrt(pow(semi_segment_len.x, 2) + pow(semi_segment_len.y, 2) +
                            pow(semi_segment_len.y, 2));

    // for(float axis_x = point_xyz.x; axis_x < sdf_bbox.max_pos.x; axis_x += segment_len.x)
    //   for(float axis_y = point_xyz.y; axis_y < sdf_bbox.max_pos.y; axis_y += segment_len.y)
    //     for(float axis_z = point_xyz.z; axis_z < sdf_bbox.max_pos.z; axis_z += segment_len.z)
    //       if(sdf(float3(axis_x, axis_y, axis_z)) <= semidiag_len){
    //         bbox_list.push_back(AABB(float3(axis_x, axis_y, axis_z) - semi_segment_len,
    //                                  float3(axis_x, axis_y, axis_z) + semi_segment_len));
    //       }

    std::vector<float3> max_bbox;
    std::vector<float3> min_bbox;
    std::vector<unsigned int> is_marked;
    /*
      values in is_marked vectors:
        0 - an empty bbox

        3 - the bbox is in contact with three others
        2 - the bbox is in contact with other two
        1 - the bbox is in contact with another

        4 - unvisited bbox
        5 - visited bbox
    */


    // Adding "empty" boxes
    for (unsigned int n = 0; n < pow(bbox_segment_quant + 2, 2); n++){
      bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
      max_bbox.push_back(float3(0, 0, 0));
      min_bbox.push_back(float3(0, 0, 0));
      is_marked.push_back(0);
    }
    
    for(float axis_x = point_xyz.x; axis_x < sdf_bbox.max_pos.x; axis_x += segment_len.x){
      for(unsigned int n = 0; n < bbox_segment_quant + 2; n++){
        bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
        max_bbox.push_back(float3(0, 0, 0));
        min_bbox.push_back(float3(0, 0, 0));
        is_marked.push_back(0);
      }

      for(float axis_y = point_xyz.y; axis_y < sdf_bbox.max_pos.y; axis_y += segment_len.y){
        bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
        max_bbox.push_back(float3(0, 0, 0));
        min_bbox.push_back(float3(0, 0, 0));
        is_marked.push_back(0);

        for(float axis_z = point_xyz.z; axis_z < sdf_bbox.max_pos.z; axis_z += segment_len.z){
          bbox_list.push_back(AABB(float3(axis_x, axis_y, axis_z) - semi_segment_len,
                                   float3(axis_x, axis_y, axis_z) + semi_segment_len));
          min_bbox.push_back(float3(axis_x, axis_y, axis_z) - semi_segment_len);
          max_bbox.push_back(float3(axis_x, axis_y, axis_z) + semi_segment_len);
          if(sdf(float3(axis_x, axis_y, axis_z)) <= semidiag_len)
            is_marked.push_back(4);
          else
            is_marked.push_back(0);
        }

        bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
        max_bbox.push_back(float3(0, 0, 0));
        min_bbox.push_back(float3(0, 0, 0));
        is_marked.push_back(0);
      }

      for(unsigned int n = 0; n < bbox_segment_quant + 2; n++){
        bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
        max_bbox.push_back(float3(0, 0, 0));
        min_bbox.push_back(float3(0, 0, 0));
        is_marked.push_back(0);
      }

    }

    for (unsigned int n = 0; n < pow(bbox_segment_quant + 2, 2); n++){
      bbox_list.push_back(AABB(float3(0, 0, 0), float3(1, 1, 1)));
      max_bbox.push_back(float3(0, 0, 0));
      min_bbox.push_back(float3(0, 0, 0));
      is_marked.push_back(0);
    }

    unsigned int ratio_i = pow(bbox_segment_quant + 2, 2);
    unsigned int ratio_j = (bbox_segment_quant + 2);

    float touches_quantity;
    for(unsigned int i = 0; i < bbox_segment_quant + 2; i++){
      for(unsigned int j = 0; j < bbox_segment_quant + 2; j++){
        for(unsigned int k = 0; k < bbox_segment_quant + 2; k++){
          if(is_marked[i * ratio_i + j * ratio_j + k] == 4){
            touches_quantity = is_marked[(i + 1) * ratio_i + j * ratio_j + k] + \
                               is_marked[(i - 1) * ratio_i + j * ratio_j + k] + \
                               is_marked[i * ratio_i + (j + 1) * ratio_j + k] + \
                               is_marked[i * ratio_i + (j - 1) * ratio_j + k] + \
                               is_marked[i * ratio_i + j * ratio_j + (k + 1)] + \
                               is_marked[i * ratio_i + j * ratio_j + (k - 1)];
            touches_quantity /= 4;

            if(touches_quantity > 3)
              continue;
            else if(touches_quantity == 1)
              is_marked[i * ratio_i + j * ratio_j + k] = 1;
            else{
              if((is_marked[i * ratio_i + j * ratio_j + k] == 4 &&
                  is_marked[i * ratio_i + (j - 1) * ratio_j + k] == 4) ||

                 (is_marked[i * ratio_i + j * ratio_j + (k + 1)] == 4 &&
                  is_marked[i * ratio_i + j * ratio_j + (k - 1)] == 4) ||

                 (is_marked[(i + 1) * ratio_i + j * ratio_j + k] == 4 &&
                  is_marked[(i - 1) * ratio_i + j * ratio_j + k] == 4)){
                continue;
              }
              else{
                if(touches_quantity == 3)
                  is_marked[i * ratio_i + j * ratio_j + k] = 3;
                else 
                  is_marked[i * ratio_i + j * ratio_j + k] = 2;
              }
            }
          }
        }
      }
    }


    // Joining bboxes
    std::vector<AABB> bbox_join;

    bool ind_i, ind_j, ind_k; // 1 - it is possible to expand the bbox along the x axis
                              // 0 - it is impossible
    bool ind; // 0 - it is possible to combine the bboxes
              // 1 - it is not yet possible
    float max_i, max_j, max_k;

    for(unsigned int i = 0; i < bbox_segment_quant + 2; i++){
      for(unsigned int j = 0; j < bbox_segment_quant + 2; j++){
        for(unsigned int k = 0; k < bbox_segment_quant + 2; k++){
          if(is_marked[i * ratio_i + j * ratio_j + k] == 3){
            ind_i = ind_j = ind_k = 1;
            ind = 1;

            for(unsigned int h = 1; h < bbox_segment_quant + 2; h++){
              if(ind == 1){
                if(ind_i == 1)
                  max_i = i + h;
                if(ind_j == 1)
                  max_j = j + h;
                if(ind_k == 1)
                  max_k = k + h;

                int ax = 3;
                for(int hi = max_i; hi >= 0; hi--){
                  for(int hj = max_j; hj >= 0; hj--){
                    for(int hk = max_k; hk >= 0; hk--){
                      int is_marked_max = is_marked[(max_i - hi) * ratio_i +
                                                    (max_j - hj) * ratio_j +
                                                    (max_k - hk)];
                      if(is_marked_max != 3 && is_marked_max != 4){
                        if(ind_i == 1){
                          max_i -= 1;
                          ind_i = 0;
                          ax--;
                          break;
                        }
                        else if(ind_j == 1 && ax == 2){
                          max_j -= 1;
                          ind_j = 0;
                          break;
                        }
                        else if(ind_k == 1 && ax == 1){
                          max_k -= 1;
                          ind_k = 0;
                          ind = 0;
                        }
                      }
                    }
                    if(ax == 2){
                      ax--;
                      break;
                    }
                    else if(ax == 1)
                      break;
                  }
                  if(ax == 1){
                    ax--;
                    break;
                  }
                }
              }
              else if(ind == 0){
                bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                         max_bbox[max_i * ratio_i + max_j * ratio_j
                                                  + max_k]));

                for(unsigned int hi = 0; hi < max_i; hi++)
                  for(unsigned int hj = 0; hj < max_j; hj++)
                    for(unsigned int hk = 0; hk < max_k; hk++)
                      is_marked[hi * ratio_i + hj * ratio_j + hk] = 5;
                
                break;
              }
            }
          }
        }
      }
    }

    for(unsigned int i = 0; i < bbox_segment_quant + 2; i++){
      for(unsigned int j = 0; j < bbox_segment_quant + 2; j++){
        for(unsigned int k = 0; k < bbox_segment_quant + 2; k++){
          if(is_marked[i * ratio_i + j * ratio_j + k] == 2 ||
             is_marked[i * ratio_i + j * ratio_j + k] == 3){
            // fixing the i axis
            ind_j = ind_k = 1;
            ind = 1;

            for(unsigned int h = 1; h < bbox_segment_quant + 2; h++){
              if(ind == 1){
                if(ind_j == 1)
                  max_j = j + h;
                if(ind_k == 1)
                  max_k = k + h;
                
                int ax = 2;
                for(int hj = max_j; hj >= 0; hj--){
                  for(int hk = max_k; hk >= 0; hk--){
                    int is_marked_max = is_marked[i * ratio_i +
                                                  (max_j - hj) * ratio_j +
                                                  (max_k - hk)];
                    if(is_marked_max == 1 || is_marked_max == 5 || is_marked_max == 0){
                      if(ind_j == 1){
                        max_j -= 1;
                        ind_j = 0;
                        ax--;
                        break;
                      }
                      else if(ind_k == 1 && ax == 1){
                        max_k -= 1;
                        ind_k = 0;
                        ind = 0;
                      }
                    }
                  }
                  if(ax == 1){
                    ax--;
                    break;
                  }
                }
              }
              else if(ind == 0){
                bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                         max_bbox[i * ratio_i + max_j * ratio_j
                                                  + max_k]));
                
                for(unsigned int hj = 0; hj < max_j; hj++)
                  for(unsigned int hk = 0; hk < max_k; hk++)
                    is_marked[i * ratio_i + hj * ratio_j + hk] = 5;
              }
            }


            // fixing the j axis
            ind_i = ind_k = 1;
            ind = 1;
            for(unsigned int h = 1; h < bbox_segment_quant + 2; h++){
              if(ind == 1){
                if(ind_i == 1)
                  max_i = i + h;
                if(ind_k == 1)
                  max_k = k + h;
                
                int ax = 2;
                for(int hi = max_i; hi >= 0; hi--){
                  for(int hk = max_k; hk >= 0; hk--){
                    int is_marked_max = is_marked[(max_i - hi) * ratio_i +
                                                  j * ratio_j +
                                                  (max_k - hk)];
                    if(is_marked_max == 1 || is_marked_max == 5){
                      if(ind_i == 1){
                        max_i -= 1;
                        ind_i = 0;
                        ax--;
                      }
                      else if(ind_k == 1 && ax == 1){
                        max_k -= 1;
                        ind_k = 0;
                        ind = 0;
                      }
                    }
                  }
                  if(ax == 1){
                    ax--;
                    break;
                  }
                }
              }
              else if(ind == 0){
                bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                         max_bbox[max_i * ratio_i + j * ratio_j
                                                  + max_k]));
                
                for(unsigned int hi = 0; hi < max_j; hi++)
                  for(unsigned int hk = 0; hk < max_k; hk++)
                    is_marked[hi * ratio_i + j * ratio_j + hk] = 5;
              }
            }

            // fixing the k axis
            ind_i = ind_j = 1;
            ind = 1;
            for(unsigned int h = 1; h < bbox_segment_quant + 2; h++){
              if(ind == 1){
                if(ind_i == 1)
                  max_i = i + h;
                if(ind_j == 1)
                  max_j = j + h;
                
                int ax = 2;
                for(int hi = max_i; hi >= 0; hi--){
                  for(int hj = max_j; hj >= 0; hj--){
                    int is_marked_max = is_marked[(max_i - hi) * ratio_i +
                                                  (max_j - hj) * ratio_j + k];
                    if(is_marked_max == 1 || is_marked_max == 5 || is_marked_max == 0){
                      if(ind_i == 1){
                        max_i -= 1;
                        ind_i = 0;
                        ax--;
                        break;
                      }
                      else if(ind_j == 1 && ax == 1){
                        max_j -= 1;
                        ind_j = 0;
                        ind = 0;
                      }
                    }
                  }
                  if(ax == 1){
                    ax--;
                    break;
                  }
                }
              }
              else if(ind == 0){
                bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                         max_bbox[max_i * ratio_i + max_j * ratio_j
                                                  + k]));
                
                for(unsigned int hi = 0; hi < max_i; hi++)
                  for(unsigned int hj = 0; hj < max_j; hj++)
                    is_marked[hi * ratio_i + hj * ratio_j + k] = 5;
              }
            }
          }
        }
      }
    }

    for(unsigned int i = 0; i < bbox_segment_quant + 2; i++){
      for(unsigned int j = 0; j < bbox_segment_quant + 2; j++){
        for(unsigned int k = 0; k < bbox_segment_quant + 2; k++){
          unsigned int val_1_cont;
          if(is_marked[i * ratio_i + j * ratio_j + k] == 1 || 
             is_marked[i * ratio_i + j * ratio_j + k] == 2 ||
             is_marked[i * ratio_i + j * ratio_j + k] == 3){
            int coeff = 1;
            if(is_marked[(i + coeff) * ratio_i + j * ratio_j + k] != 5 &&
               is_marked[(i + coeff) * ratio_i + j * ratio_j + k] != 0){

              val_1_cont = is_marked[(i + coeff) * ratio_i + j * ratio_j + k];
              is_marked[(i + coeff) * ratio_i + j * ratio_j + k] = 5;

              while(true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[(i + coeff) * ratio_i + j *
                                           ratio_j + k]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[(i + coeff) * ratio_i + j * ratio_j + k];
                  is_marked[(i + coeff) * ratio_i + j * ratio_j + k] = 5;
                }
                  
              }
            }

            else if(is_marked[(i - coeff) * ratio_i + j * ratio_j + k] != 5 &&
                    is_marked[(i - coeff) * ratio_i + j * ratio_j + k] != 0){
              val_1_cont = is_marked[(i - coeff) * ratio_i + j * ratio_j + k];
              is_marked[(i - coeff) * ratio_i + j * ratio_j + k] = 5;
              
              while (true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[(i - coeff) * ratio_i + j *
                                           ratio_j + k]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[(i - coeff) * ratio_i + j * ratio_j + k];
                }
              }
            }


            else if(is_marked[i * ratio_i + (j + coeff) * ratio_j + k] != 5 &&
                    is_marked[i * ratio_i + (j + coeff) * ratio_j + k] != 0){
              val_1_cont = is_marked[i * ratio_i + (j + coeff) * ratio_j + k];
              is_marked[i * ratio_i + (j + coeff) * ratio_j + k] = 5;

              while(true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[i * ratio_i + (j + coeff) *
                                           ratio_j + k]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[i * ratio_i + (j + coeff) * ratio_j + k];
                }
              }
            }

            else if(is_marked[i * ratio_i + (j - coeff) * ratio_j + k] != 5 &&
                    is_marked[i * ratio_i + (j - coeff) * ratio_j + k] != 0){
              val_1_cont = is_marked[i * ratio_i + (j - coeff) * ratio_j + k];
              is_marked[i * ratio_i + (j - coeff) * ratio_j + k] = 5;

              while(true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[i * ratio_i + (j - coeff) *
                                           ratio_j + k]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[i * ratio_i + (j - coeff) * ratio_j + k];
                }
              }
            }


            else if(is_marked[i * ratio_i + j * ratio_j + (k + coeff)] != 5 &&
                    is_marked[i * ratio_i + j * ratio_j + (k + coeff)] != 0){
              val_1_cont = is_marked[i * ratio_i + j * ratio_j + (k + coeff)];
              is_marked[i * ratio_i + j * ratio_j + (k + coeff)] = 5;

              while(true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[i * ratio_i + j * ratio_j +
                                           (k + coeff)]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[i * ratio_i + j * ratio_j + (k + coeff)];
                }
              }
            }

            else if(is_marked[i * ratio_i + j * ratio_j + (k - coeff)] != 5 &&
                    is_marked[i * ratio_i + j * ratio_j + (k - coeff)] != 0){
              val_1_cont = is_marked[i * ratio_i + j * ratio_j + (k - coeff)];
              is_marked[i * ratio_i + j * ratio_j + (k - coeff)] = 5;

              while(true){
                if(val_1_cont == 0 || val_1_cont == 5){
                  bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                           max_bbox[i * ratio_i + j * ratio_j +
                                           (k - coeff)]));
                  break;
                }
                else{
                  coeff++;
                  val_1_cont = is_marked[i * ratio_i + j * ratio_j + (k - coeff)];
                }
              }
            }
          }
        }
      }
    }

    // Checking for the remaining boxes
    for(unsigned int i = 1; i < bbox_segment_quant + 1; i++){
      for(unsigned int j = 1; j < bbox_segment_quant + 1; j++){
        for(unsigned int k = 1; k < bbox_segment_quant + 1; k++){
          if(is_marked[i * ratio_i + j * ratio_j + k] != 0 &&
            is_marked[i * ratio_i + j * ratio_j + k] != 5){
            is_marked[i * ratio_i + j * ratio_j + k] = 5;
            bbox_join.push_back(AABB(min_bbox[i * ratio_i + j * ratio_j + k],
                                     max_bbox[i * ratio_i + j * ratio_j + k]));
          }
        }
      }
    }

    return bbox_join;
  }
}