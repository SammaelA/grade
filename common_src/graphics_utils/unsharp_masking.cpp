#include "unsharp_masking.h"
#include "tinyEngine/engine.h"
#include "gauss_blur_precise.h"
#include "image_arithmetic.h"

Texture UnsharpMasking::perform(Texture &t, float sigma, float strength)
{
  assert(t.type == GL_TEXTURE_2D);

  float w = t.get_W();
  float h = t.get_H();

  Texture res_tex = engine::textureManager->create_texture(w, h);
  GaussFilter gauss(sigma);
  Texture blurred_tex = gauss.perform_gauss_blur(t);
  ImageArithmetics::add(res_tex, t, blurred_tex, 1 + strength, -strength);

  return res_tex;
}