#pragma once

void compare_sandbox(int argc, char **argv);
void render_mygen_cup_not_diff(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                           std::string texture_path, std::string save_dir);
void render_mygen_cup(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                      std::string texture_path, std::string save_dir);
void render_mygen_cup_demo(MitsubaInterface &mi, CameraSettings &camera, std::vector<float> params,
                           std::string texture_path, std::string save_dir);