#pragma once

void compare_sandbox(int argc, char **argv);
void render_mygen_cup(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                      std::string texture_path, std::string save_dir);