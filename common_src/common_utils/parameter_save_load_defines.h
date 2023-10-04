#pragma once

#define P_INT(par, type)                        \
  if (type == SaveLoadMode::BLK_SAVE)           \
    b.set_double(#par, par);                    \
  else if (type == SaveLoadMode::BLK_LOAD)      \
    par = b.get_double(#par, par);              \
  else if (type == SaveLoadMode::PAR_LIST_SAVE) \
    par = list.ordinalParameters.at(#par);      \
  else                                          \
    list.ordinalParameters.emplace(#par, par);

#define P_CAT(par, type)                        \
  if (type == SaveLoadMode::BLK_SAVE)           \
    b.set_double(#par, par);                    \
  else if (type == SaveLoadMode::BLK_LOAD)      \
    par = b.get_double(#par, par);              \
  else if (type == SaveLoadMode::PAR_LIST_SAVE) \
    par = list.categorialParameters.at(#par);   \
  else                                          \
    list.categorialParameters.emplace(#par, par);

#define P_FLOAT(par, type)                      \
  if (type == SaveLoadMode::BLK_SAVE)           \
    b.set_double(#par, par);                    \
  else if (type == SaveLoadMode::BLK_LOAD)      \
    par = b.get_double(#par, par);              \
  else if (type == SaveLoadMode::PAR_LIST_SAVE) \
    par = list.continuousParameters.at(#par);   \
  else                                          \
    list.continuousParameters.emplace(#par, par);

#define P_VEC2(par, type)                                               \
  if (type == SaveLoadMode::BLK_SAVE)                                   \
    b.set_vec2(#par, par);                                              \
  else if (type == SaveLoadMode::BLK_LOAD)                              \
    par = b.get_vec2(#par, par);                                        \
  else if (type == SaveLoadMode::PAR_LIST_SAVE)                         \
  {                                                                     \
    par.x = list.continuousParameters.at(std::string(#par) + "_x");     \
    par.y = list.continuousParameters.at(std::string(#par) + "_y");     \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    list.continuousParameters.emplace(std::string(#par) + "_x", par.x); \
    list.continuousParameters.emplace(std::string(#par) + "_y", par.y); \
  }

#define P_VEC3(par, type)                                               \
  if (type == SaveLoadMode::BLK_SAVE)                                   \
    b.set_vec3(#par, par);                                              \
  else if (type == SaveLoadMode::BLK_LOAD)                              \
    par = b.get_vec3(#par, par);                                        \
  else if (type == SaveLoadMode::PAR_LIST_SAVE)                         \
  {                                                                     \
    par.x = list.continuousParameters.at(std::string(#par) + "_x");     \
    par.y = list.continuousParameters.at(std::string(#par) + "_y");     \
    par.z = list.continuousParameters.at(std::string(#par) + "_z");     \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    list.continuousParameters.emplace(std::string(#par) + "_x", par.x); \
    list.continuousParameters.emplace(std::string(#par) + "_y", par.y); \
    list.continuousParameters.emplace(std::string(#par) + "_z", par.z); \
  }

#define P_VEC4(par, type)                                               \
  if (type == SaveLoadMode::BLK_SAVE)                                   \
    b.set_vec4(#par, par);                                              \
  else if (type == SaveLoadMode::BLK_LOAD)                              \
    par = b.get_vec4(#par, par);                                        \
  else if (type == SaveLoadMode::PAR_LIST_SAVE)                         \
  {                                                                     \
    par.x = list.continuousParameters.at(std::string(#par) + "_x");     \
    par.y = list.continuousParameters.at(std::string(#par) + "_y");     \
    par.z = list.continuousParameters.at(std::string(#par) + "_z");     \
    par.w = list.continuousParameters.at(std::string(#par) + "_w");     \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    list.continuousParameters.emplace(std::string(#par) + "_x", par.x); \
    list.continuousParameters.emplace(std::string(#par) + "_y", par.y); \
    list.continuousParameters.emplace(std::string(#par) + "_z", par.z); \
    list.continuousParameters.emplace(std::string(#par) + "_w", par.w); \
  }
