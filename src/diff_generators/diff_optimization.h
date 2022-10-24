#pragma once

struct Block;
class MitsubaInterface;
namespace dopt
{
  void test();
  float image_based_optimization(Block &settings_blk, MitsubaInterface &mi);
}