#include "utils.h"
#include "dmesh.h"
#include <iostream>

#ifdef WIN32
  #include <direct.h>     // for windows mkdir
#else
  #include <sys/stat.h>   // for linux mkdir
  #include <sys/types.h>
  #include <ftw.h>
  #include <errno.h>
#endif
namespace diff_render
{
float LossAndDiffLoss(const Img& b, const Img& a, Img& a_outDiff)
{
  assert(a.width()*a.height() == b.width()*b.height());
  double accumMSE = 0.0f;
  const size_t imgSize = a.width()*a.height();
  for(size_t i=0;i<imgSize;i++)
  {
    const float3 diffVec = b.data()[i] - a.data()[i];
    a_outDiff.data()[i] = 2.0f*diffVec;                    // (I[x,y] - I_target[x,y])    // dirrerential of the loss function 
    accumMSE += double(dot(diffVec, diffVec));             // (I[x,y] - I_target[x,y])^2  // the loss function itself
  }
  return float(accumMSE);
}

void CHECK_NaN(float f)
{
  if (f != f)
  {
    printf("exception!!!\n");
    float *t = nullptr;
    printf("%f\n",*t);
  }
}

void CHECK_NaN(float3 f)
{
  if (f.x != f.x || f.y != f.y || f.z != f.z)
  {
    printf("exception!!!\n");
    float *t = nullptr;
    printf("%f\n",*t);
  }
}

#ifndef WIN32
int unlink_cb_fun(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    int rv;

    if (ftwbuf->level == 0)
        return 0;

    rv = remove(fpath);

    if (rv)
        perror(fpath);

    return rv;
}
#endif

void prepare_and_clear_directory(const ::std::string &dir)
{
  #ifdef WIN32
  mkdir(dir.c_str());
  //TODO: clear direactory
  #else
  int res = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (res == 0 || errno == EEXIST)
  {
    res = nftw(dir.c_str(), unlink_cb_fun, 64, FTW_DEPTH | FTW_PHYS);
    if (res != 0)
      logerr("failed to clear directory %s : %s", dir.c_str(), strerror(errno));
  }
  else
  {
    logerr("failed to create directory %s : %s", dir.c_str(), strerror(errno));
  }
  #endif
}
}