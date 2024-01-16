#pragma once
#include <cstring>
#include <cstdio>
#include <errno.h>
#include <vector>

namespace nsdf
{
  struct PrimitiveSDFScene
  {
    std::vector<float> scene_data; //pos_0.x, pos_0.y, pos_0.z, radius_0, pos_1.x, ...

    bool to_file(const char *filename)
    {
      FILE *f = fopen(filename, "w");
      if (!f)
      {
        fprintf(stderr, "failed to open/create file %s. Errno %d\n",filename, (int)errno);
        return false;
      }
      fprintf(f, "spheres count = %d\n", (int)scene_data.size()/4);
      for (int i=0;i<scene_data.size();i+=4)
        fprintf(f, "pos = (%f, %f, %f), radius = %f\n", scene_data[i],scene_data[i+1],scene_data[i+2],scene_data[i+3]);

      int res = fclose(f);
      if (res != 0)
      {
        fprintf(stderr, "failed to close file %s. fclose returned %d\n",filename, res);
        return false;
      }
      return true;
    }

    bool from_file(const char *filename)
    {
      FILE *f = fopen(filename, "r");
      if (!f)
      {
        fprintf(stderr, "failed to open file %s. Errno %d\n",filename, (int)errno);
        return false;
      }
      int spheres_count = 0;
      fscanf(f, "spheres count = %d\n", &spheres_count);
      scene_data.resize(4*spheres_count, 0);
      for (int i=0;i<scene_data.size();i+=4)
        fscanf(f, "pos = (%f, %f, %f), radius = %f\n", &(scene_data[i]),&(scene_data[i+1]),&(scene_data[i+2]),&(scene_data[i+3]));

      int res = fclose(f);
      if (res != 0)
      {
        fprintf(stderr, "failed to close file %s. fclose returned %d\n",filename, res);
        return false;
      }
      return true;
    }
  };
}