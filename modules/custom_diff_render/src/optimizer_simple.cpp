#include "optimizer.h"
#include "utils.h"
#include "scene.h"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace LiteMath;
namespace diff_render
{
void OptimizerParameters::set_default()
{
  if (alg == OPT_ALGORITHM::GD_Adam)
  {
    decayPeriod = 100;
    decay_mult = 0.75;
    base_lr = 0.1;
    position_lr = 0.1;
    textures_lr = 0.1;
  }
}

struct OptSimple : public IOptimizer
{
  OptSimple(){}

  void         Init(const Scene& a_scene, ::std::shared_ptr<IDiffRender> a_pDRImpl, 
                    const CamInfo* a_cams, const Img* a_images, int a_numViews, OptimizerParameters a_params) override;

  Scene Run (size_t a_numIters, float &final_error, ::std::vector<Scene> *iterations_dump = nullptr) override;

protected:
  
  typedef ::std::vector<::std::vector<::std::pair<int, float>>> IntervalLearningRate; //{offset, learning rate} for every diff. mesh

  float EvalFunction(const Scene& mesh, DScene& gradScene);
  IntervalLearningRate GetLR(DScene& gradScene);
  void  OptStep(DScene &gradScene, const IntervalLearningRate &lr);
  void  OptUpdateScene(DScene &gradScene, Scene* mesh);
  void  StepDecay(int a_iterId, IntervalLearningRate &lr) const;

  Scene        m_scene; ///<! global mesh optimized mesh
  size_t       m_iter = 0;
  OptimizerParameters m_params;

  ::std::vector<GradReal> m_GSquare; ///<! m_GSquare is a vector of the sum of the squared gradients at or before iteration 'i'
  ::std::vector<GradReal> m_vec; 

  ::std::shared_ptr<IDiffRender> m_pDR = nullptr;
  const Img*     m_targets  = nullptr; 
  const CamInfo* m_cams     = nullptr; 
  int            m_numViews = 0;
};

IOptimizer* CreateSimpleOptimizer() { return new OptSimple; };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

OptSimple::IntervalLearningRate OptSimple::GetLR(DScene& gradScene)
{
  IntervalLearningRate lrr;
  for (auto &dmesh : gradScene.get_dmeshes())
  {
    GradReal *ptr = dmesh.full_data();
    lrr.emplace_back();
    auto &lr = lrr.back();
    lr.push_back({0, m_params.position_lr});
    if (dmesh.color(0))
      lr.push_back({dmesh.color(0) - ptr, m_params.base_lr});
    if (dmesh.tex_count() > 0)
      lr.push_back({dmesh.tex(0,0) - ptr, m_params.textures_lr});
    if (dmesh.transform_mat(0))
      lr.push_back({dmesh.transform_mat(0) - ptr, m_params.transforms_lr});
    if (dmesh.restricted_transform(0))
      lr.push_back({dmesh.restricted_transform(0) - ptr, m_params.transforms_lr});
  }
  return lrr;
}

void  OptSimple::StepDecay(int a_iterId, IntervalLearningRate &lrr) const
{
  if(a_iterId > 0 && a_iterId % m_params.decayPeriod == 0)
  {
    for (auto &lr : lrr)
      for (auto &p : lr)
        p.second *= m_params.decay_mult;
  }
}

void OptSimple::OptUpdateScene(DScene &gradScene, Scene* scene)
{
  for (int mesh_id = 0; mesh_id < scene->get_meshes().size(); mesh_id++)
  {
    auto &mesh = scene->get_mesh_modify(mesh_id);
    auto *dmesh_ptr = gradScene.get_dmesh(mesh_id);
    if (!dmesh_ptr)
      continue;
    auto &dmesh = *dmesh_ptr;
    assert(dmesh.vertex_count() == mesh.vertex_count());
    
    GradReal *d_pos = dmesh.pos(0);
    if (d_pos)
    {
      for(int vertId=0; vertId< mesh.vertex_count(); vertId++)
      {
        mesh.vertices[vertId].x -= d_pos[3*vertId+0];
        mesh.vertices[vertId].y -= d_pos[3*vertId+1];
        mesh.vertices[vertId].z -= d_pos[3*vertId+2];
      }
    }
    
    GradReal *d_tr = dmesh.transform_mat(0);
    if (d_tr)
    {
      for (int tr_n=0; tr_n<scene->get_transform_modify(mesh_id).size(); tr_n++)
      {
        auto &tr = scene->get_transform_modify(mesh_id)[tr_n];
        d_tr = dmesh.transform_mat(tr_n);

        tr[0][0] -= d_tr[0];
        tr[0][1] -= d_tr[1];
        tr[0][2] -= d_tr[2];
        tr[0][3] -= d_tr[3];

        tr[1][0] -= d_tr[4];
        tr[1][1] -= d_tr[5];
        tr[1][2] -= d_tr[6];
        tr[1][3] -= d_tr[7];

        tr[2][0] -= d_tr[8];
        tr[2][1] -= d_tr[9];
        tr[2][2] -= d_tr[10];
        tr[2][3] -= d_tr[11];
      }
    }
    else
    {
      d_tr = dmesh.restricted_transform(0);
      if (d_tr)
      {
        for (int tr_n=0; tr_n<scene->get_restricted_transform_modify(mesh_id).size(); tr_n++)
        {
          auto &tr = scene->get_restricted_transform_modify(mesh_id)[tr_n];
          d_tr = dmesh.restricted_transform(tr_n);

          tr.data.s.translate.x -= d_tr[0];
          tr.data.s.translate.y -= d_tr[1];
          tr.data.s.translate.z -= d_tr[2];

          tr.data.s.rotate.x -= d_tr[3];
          tr.data.s.rotate.y -= d_tr[4];
          tr.data.s.rotate.z -= d_tr[5];

          tr.data.s.scale -= d_tr[6];
        }
      }
    }

    GradReal *d_col = dmesh.color(0);
    if (d_col)
    {
      assert(dmesh.vertex_count() == mesh.colors.size());
      for(int faceId=0; faceId < mesh.colors.size(); faceId++)
      {
        mesh.colors[faceId].x -= d_col[3*faceId+0];
        mesh.colors[faceId].y -= d_col[3*faceId+1];
        mesh.colors[faceId].z -= d_col[3*faceId+2];
      }
    }
    
    if (gradScene.get_shading_model() != SHADING_MODEL::SILHOUETTE &&
        gradScene.get_shading_model() != SHADING_MODEL::VERTEX_COLOR)
    {
      for (int i=0;i<dmesh.tex_count();i++)
      {
        int sz = mesh.textures[i].data.size();
        GradReal *d_tex = dmesh.tex(i, 0);
        for (int j=0;j<sz;j++)
          mesh.textures[i].data[j] -= d_tex[j];
      }
    }
  }
}

void OptSimple::OptStep(DScene &gradScene, const IntervalLearningRate &lr)
{
  int off = 0;
  int lr_n = 0;
  for (auto &gradMesh : gradScene.get_dmeshes())
  {
    GradReal *grad = gradMesh.full_data();
    if(m_params.alg >= OptimizerParameters::GD_AdaGrad)
    {
      if(m_params.alg == OptimizerParameters::GD_AdaGrad)  // ==> GSquare[i] = gradF[i]*gradF[i]
      {
        for(size_t i=0;i<gradMesh.full_size();i++)
          m_GSquare[off+i] += (grad[i]*grad[i]);
      }
      else if(m_params.alg == OptimizerParameters::GD_RMSProp) // ==> GSquare[i] = GSquarePrev[i]*a + (1.0f-a)*gradF[i]*gradF[i]
      {
        const float alpha = 0.5f;
        for(size_t i=0;i<gradMesh.full_size();i++)
          m_GSquare[off+i] = 2.0f*(m_GSquare[off+i]*alpha + (grad[i]*grad[i])*(1.0f-alpha)); // does not works without 2.0f
      }
      else if(m_params.alg == OptimizerParameters::GD_Adam) // ==> Adam(m[i] = b*mPrev[i] + (1-b)*gradF[i], GSquare[i] = GSquarePrev[i]*a + (1.0f-a)*gradF[i]*gradF[i])
      {
        const float alpha = 0.5f;
        const float beta  = 0.25f;
        for(size_t i=0;i<gradMesh.full_size();i++)
          m_vec[off+i] = m_vec[off+i]*beta + grad[i]*(1.0f-beta);

        for(size_t i=0;i<gradMesh.full_size();i++)
          m_GSquare[off+i] = 2.0f*(m_GSquare[off+i]*alpha + (grad[i]*grad[i])*(1.0f-alpha)); // does not works without 2.0f

        for(size_t i=0;i<gradMesh.full_size();i++)
          grad[i] = m_vec[off+i];
      }
      
      //xNext[i] = x[i] - gamma/(sqrt(GSquare[i] + epsilon));
      for (int i=0;i<gradMesh.full_size();i++)
        grad[i] = grad[i]/( ::std::sqrt(m_GSquare[i] + GradReal(1e-8f)));
    }
    
    for (int i=0;i<lr[lr_n].size();i++)
    {
      int next_offset = (i == lr[lr_n].size() - 1) ? gradMesh.full_size() : lr[lr_n][i+1].first;
      for (int j=lr[lr_n][i].first; j<next_offset; j++)
        grad[j] *= lr[lr_n][i].second;
    }
    off += gradMesh.full_size();
    lr_n++;
  }
}

float OptSimple::EvalFunction(const Scene& scene, DScene& gradScene)
{
  ::std::vector<Img> images(m_numViews);
  float mse = m_pDR->d_render_and_compare(scene, m_cams, m_targets, m_numViews, m_targets[0].width()*m_targets[0].height(), gradScene, 
                                          m_params.verbose ? images.data() : nullptr);
  if (m_params.verbose)
  {
    for(int i=0;i<m_numViews;i++) 
    {
      ::std::stringstream strOut;
      strOut  << "output/rendered_opt" << i << "/render_" << ::std::setfill('0') << ::std::setw(4) << m_iter << ".bmp";
      auto temp = strOut.str();
      LiteImage::SaveImage(temp.c_str(), images[i]);
    }
  }
  m_iter++;
  return mse;
}

void OptSimple::Init(const Scene& a_scene, ::std::shared_ptr<IDiffRender> a_pDRImpl, 
                     const CamInfo* a_cams, const Img* a_images, int a_numViews, OptimizerParameters a_params) 
{ 
  assert(a_numViews > 0);

  m_scene       = a_scene; 
  m_iter        = 0; 
  m_pDR         = a_pDRImpl;

  m_targets     = a_images;
  m_cams        = a_cams;
  m_numViews    = a_numViews;
  m_params      = a_params;
}

Scene OptSimple::Run(size_t a_numIters, float &final_error, ::std::vector<Scene> *iterations_dump) 
{ 
  DScene gradScene;
  gradScene.reset(m_scene, m_pDR->mode, m_params.differentiable_mesh_ids);

  if(m_params.alg >= OptimizerParameters::GD_AdaGrad) {
    m_GSquare.resize(gradScene.full_size());
    memset(m_GSquare.data(), 0, sizeof(GradReal)*m_GSquare.size());
  }

  if(m_params.alg == OptimizerParameters::GD_Adam) {
    m_vec.resize(gradScene.full_size());
    memset(m_vec.data(), 0, sizeof(GradReal)*m_vec.size());
  }

  auto lr = GetLR(gradScene);
  
  m_iter = 0;
  float error = 0;
  float psnr = 0;

  for(size_t iter=0, trueIter = 0; iter < a_numIters; iter++, trueIter++)
  {
    error = EvalFunction(m_scene, gradScene);
    psnr = -10*log10(max(1e-9f,error));
    if (m_params.verbose)
      ::std::cout << "iter " << trueIter << ", PSNR = " << psnr << ::std::endl;
    OptStep(gradScene, lr);
    OptUpdateScene(gradScene, &m_scene);
    StepDecay(iter, lr);

    if (iterations_dump)
      iterations_dump->push_back(m_scene);
  }

  final_error = error;
  return m_scene;
}
}