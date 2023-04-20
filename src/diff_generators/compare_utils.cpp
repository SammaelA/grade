#include <string>
#include <glm/glm.hpp>
#include "third_party/stb_image.h"
#include "cstdio"

namespace compare_utils
{
  #define CU_MSE(x, y) ((x)/255.0 - (y)/255.0)*((x)/255.0 - (y)/255.0)
  float loss_psnr(const std::string &img1, const std::string &img2)
  {
    int w1=0,h1=0,c1=0;
    auto *t1 = stbi_load(img1.c_str(),&w1,&h1,&c1,0);
    int w2=0,h2=0,c2=0;
    auto *t2 = stbi_load(img2.c_str(),&w2,&h2,&c2,0);

    //printf("%s %d %d %d %d %d %d\n", img1.c_str(), w1, w2, h1, h2, c1, c2);

    assert(w1 == w2 && h1 == h2 && c1 == c2 && c1>=3);

    double mse = 0;
    double cnt = 0;
    for (int i=0;i<w1*h1*c1;i+=c1)
    {
      mse += CU_MSE(t1[i], t2[i]) + CU_MSE(t1[i+1], t2[i+1]) + CU_MSE(t1[i+2], t2[i+2]);
      cnt += 3;
    }
    mse /= cnt;

    float psnr = 10*(log(1.0/mse)/log(10.0));
    //printf("psnr %f\n", psnr);
    return psnr;
  }

  float loss_silhouette_psnr(const std::string &img1, const std::string &img2, glm::vec3 background_color)
  {
    unsigned char r0 = 255*background_color.r;
    unsigned char g0 = 255*background_color.g;
    unsigned char b0 = 255*background_color.b;

    int w1=0,h1=0,c1=0;
    auto *t1 = stbi_load(img1.c_str(),&w1,&h1,&c1,0);
    int w2=0,h2=0,c2=0;
    auto *t2 = stbi_load(img2.c_str(),&w2,&h2,&c2,0);

    //printf("%s %d %d %d %d %d %d\n", img1.c_str(), w1, w2, h1, h2, c1, c2);

    assert(w1 == w2 && h1 == h2 && c1 == c2 && c1>=3);

    double mse = 0;
    double cnt = 0;
    for (int i=0;i<w1*h1*c1;i+=c1)
    {
      if (!(t1[i] == r0 && t1[i+1] == g0 && t1[i+2] == b0) || 
          !(t2[i] == r0 && t2[i+1] == g0 && t2[i+2] == b0))
      {
        mse += CU_MSE(t1[i], t2[i]) + CU_MSE(t1[i+1], t2[i+1]) + CU_MSE(t1[i+2], t2[i+2]);
        cnt += 3;
      }
    }
    mse /= cnt;

    float psnr = 10*(log(1.0/mse)/log(10.0));
    //printf("psnr %f\n", psnr);
    return psnr;
  }
  float loss_silhouette_iou(const std::string &img1, const std::string &img2, glm::vec3 background_color)
  {
    unsigned char r0 = 255*background_color.r;
    unsigned char g0 = 255*background_color.g;
    unsigned char b0 = 255*background_color.b;

    int w1=0,h1=0,c1=0;
    auto *t1 = stbi_load(img1.c_str(),&w1,&h1,&c1,0);
    int w2=0,h2=0,c2=0;
    auto *t2 = stbi_load(img2.c_str(),&w2,&h2,&c2,0);

    //printf("%s %d %d %d %d %d %d\n", img1.c_str(), w1, w2, h1, h2, c1, c2);

    assert(w1 == w2 && h1 == h2 && c1 == c2 && c1>=3);

    double i_cnt = 0;
    double u_cnt = 0;
    for (int i=0;i<w1*h1*c1;i+=c1)
    {
      bool in_1 = !(t1[i] == r0 && t1[i+1] == g0 && t1[i+2] == b0);
      bool in_2 = !(t2[i] == r0 && t2[i+1] == g0 && t2[i+2] == b0);
      i_cnt += in_1 && in_2;
      u_cnt += in_1 || in_2;
    }

    //printf("IoU %f\n", i_cnt/u_cnt);
    return i_cnt/u_cnt;
  }

  void turntable_loss(const std::string &path1, const std::string &path2, int image_count)
  {
    float l1=0,l2=0,l3=0;
    for (int i=0;i<image_count;i++)
    {
      char path[1024];
      sprintf(path, "%s/frame-%04d.png", path1.c_str(), i);
      std::string img1(path);
      sprintf(path, "%s/frame-%04d.png", path2.c_str(), i);
      std::string img2(path);

      l1 += loss_psnr(img1, img2) / image_count;
      l2 += loss_silhouette_psnr(img1, img2, glm::vec3(1,1,1)) / image_count;
      l3 += loss_silhouette_iou(img1, img2, glm::vec3(1,1,1)) / image_count;
    }
    printf("Turntable Loss for %d images\n", image_count);
    printf("Textured PSNR = %.2f\n", l1);
    printf("With silhouette PSNR = %.2f\n", l2);
    printf("Silhouette IoU = %.4f\n", l3);
  }
} // namespace compare_utils