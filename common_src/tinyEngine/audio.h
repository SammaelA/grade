#pragma once
#include <deque>
#include <unordered_map>
#include <string>
#include <initializer_list>
#include <SDL2/SDL.h>
#include <SDL2/SDL_mixer.h>
class Audio
{
public:
  bool enabled = false;

  std::unordered_map<std::string, Mix_Chunk *> sounds;
  std::deque<Mix_Chunk *> unprocessed;

  bool init();
  bool quit();

  void load(std::initializer_list<std::string> in);
  void play(std::string sound);
  void process();
};