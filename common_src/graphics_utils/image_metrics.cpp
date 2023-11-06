#include "image_metrics.h"
#include "tinyEngine/texture.h"

float ImageMetric::get(const Texture &t1, const Texture &t2, Metric metric)
{
  if (t1.get_W() != t2.get_W() || t1.get_H() != t2.get_H() || t1.get_layers() != 1 || t2.get_layers() != 1)
  {
    logerr("ImageMetric requires textures with 1 layer and same sizes");
    return 0;
  }
  int sz = 4*t1.get_W()*t1.get_H();
  float *data_1 = new float[sz];
  float *data_2 = new float[sz];

  glBindTexture(GL_TEXTURE_2D, t1.texture);
  glGetTexImage(GL_TEXTURE_2D,
                0,
                GL_RGBA,
                GL_FLOAT,
                data_1);
  glBindTexture(GL_TEXTURE_2D, t2.texture);
  glGetTexImage(GL_TEXTURE_2D,
                0,
                GL_RGBA,
                GL_FLOAT,
                data_2);
  glBindTexture(t1.type, 0);

  float res = 0;
  switch (metric)
  {
  case Metric::MAE :
  {
    double d = 0;
    for (int i=0;i<sz;i++)
      d += abs(data_1[i] - data_2[i]);
    res = d/sz;
  }
  break;
  case Metric::MSE :
  {
    double d = 0;
    for (int i=0;i<sz;i++)
      d += SQR(data_1[i] - data_2[i]);
    res = d/sz;
  }
  break;
  case Metric::PSNR :
  {
    double d = 0;
    for (int i=0;i<sz;i++)
      d += SQR(data_1[i] - data_2[i]);
    res = -10*log10(MAX(1e-9f,d/sz));
  }
  break;
  case Metric::IOU :
  {
    int un = 0;
    int in = 0;
    for (int i=0;i<sz;i++)
    {
      un += (data_1[i] > 0 || data_2[i] > 0);
      in += (data_1[i] > 0 && data_2[i] > 0);
    }
    res = ((float)in)/un;
  }
  break;
  default:
    break;
  }


  delete[] data_1;
  delete[] data_2;

  return res;
}