#pragma once

#include "common_utils/LiteMath_ext.h"
#include "dmesh.h"
#include "dmodels.h"
#include "../raytrace_src/raytrace.h"
#include "virtual_drender.h"

#if DEBUG_RENDER
constexpr static int  MAXTHREADS    = 1;
#else
constexpr static int  MAXTHREADS    = 14;
#endif

#include "myomp.h"
#include "../common_src/qmc.h"
#include <vector>
#include <set>
#include <memory>
#include <iostream>

void dmat_dtransforms_jac(float const mat[7], float jac[7][12]);
namespace diff_render
{

struct Edge 
{
  int v0, v1; // vertex ID, v0 < v1
  int vn0, vn1;
  int mesh_n;
  int instance_n;
  Edge(int a_v0, int a_v1, int a_vn0, int a_vn1, int a_mesh_n, int a_instance_n)
  {
    mesh_n = a_mesh_n;
    instance_n = a_instance_n;
    if (a_v0 < a_v1)
    {
      v0 = a_v0;
      vn0 = a_vn0;
      v1 = a_v1;
      vn1 = a_vn1;
    }
    else
    {
      v0 = a_v1;
      vn0 = a_vn1;
      v1 = a_v0;
      vn1 = a_vn0;
    }  
  }
  bool operator<(const Edge &e) const { return this->v0 != e.v0 ? this->v0 < e.v0 : this->v1 < e.v1; } // for sorting edges
};

// for sampling edges with inverse transform sampling
struct Sampler 
{
  ::std::vector<float> pmf; // probability mass function
  ::std::vector<float> cdf;
};

// build a discrete CDF using edge length
Sampler build_edge_sampler(const Scene &scene, const ::std::vector<Edge> &edges);

// binary search for inverting the CDF in the sampler
int sample(const Sampler &sampler, const float u);

inline void edge_grad(const Scene &scene, const Edge &e, const float2 d_v0, const float2 d_v1, const AuxData aux,
                      ::std::vector<::std::vector<GradReal>> &d_pos, ::std::vector<::std::vector<GradReal>> &d_tr);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<SHADING_MODEL material>
struct DiffRender : public IDiffRender
{
  DiffRender() { };
  virtual ~DiffRender() override { };
  virtual void init(const DiffRenderSettings &settings) override
  {
    m_samples_per_pixel = settings.spp;
    m_pTracer = MakeRayTracer3D("");
    m_hammSamples.resize(2*settings.spp);
    mode = material;

    qmc::init(m_table);
    qmc::planeHammersley(m_hammSamples.data(), settings.spp);
  }
  
  
  void commit(const Scene &scene) override
  {
    scene.prepare_for_render();
    m_pTracer->Init(&scene); // Build Acceleration structurres and e.t.c. if needed
    m_pLastPreparedScene = &scene;
  }

  void render(const Scene &scene, const CamInfo* cams, Img *imgames, int a_viewsNum) override // TODO: add BSPImage rendering
  {
    auto sqrt_num_samples  = (int)sqrt((float)m_samples_per_pixel);
    auto samples_per_pixel = sqrt_num_samples * sqrt_num_samples;

    if(&scene != m_pLastPreparedScene)
    {
      ::std::cout << "[DiffRender::render]: error, renderer was not prepared for this scene!" << ::std::endl;
      return;
    }
    
    for(int camId =0; camId<a_viewsNum; camId++) { // TODO: also can make parallel if m_pTracer->clone() is implemented

      const CamInfo& cam = cams[camId];
      Img&           img = imgames[camId];

      m_pTracer->SetCamera(cam);
      
      #if (DEBUG_RENDER==0)
      #pragma omp parallel for collapse (2) num_threads(MAXTHREADS) 
      #endif 
      for (int y = 0; y < img.height(); y++) { // for each pixel 
        for (int x = 0; x < img.width(); x++) {
          
          float3 pixelColor = float3(0,0,0);
  
          for (int dy = 0; dy < sqrt_num_samples; dy++) { // for each subpixel
            for (int dx = 0; dx < sqrt_num_samples; dx++) {
  
              auto xoff = (dx + 0.5f) / float(sqrt_num_samples);
              auto yoff = (dy + 0.5f) / float(sqrt_num_samples);
              auto screen_pos = float2{x + xoff, y + yoff};
      
              auto color = shade<material>(scene, m_pTracer.get(), screen_pos);
  
              pixelColor += (color / samples_per_pixel);
            }
          }
  
          img[int2(x,y)] = pixelColor;
        }
      }
    }
  }

  
  void d_render(const Scene &scene, const CamInfo* cams, const Img *adjoints, int a_viewsNum, const int edge_samples_in_total,
                DScene &d_mesh,
                Img* debugImages = nullptr, int debugImageNum = 0) override
  {  
    if(&scene != m_pLastPreparedScene)
    {
      ::std::cout << "[DiffRender::render]: error, renderer was not prepared for this scene!" << ::std::endl;
      return;
    }

    for(int camId=0; camId<a_viewsNum; camId++) {
      
      m_pTracer->SetCamera(cams[camId]);
  
      m_aux.pCamInfo      = cams + camId;
      m_aux.debugImages   = debugImages;
      m_aux.debugImageNum = debugImageNum;

      if (material != SHADING_MODEL::SILHOUETTE) //constexpr
        interior_derivatives(scene, adjoints[camId], d_mesh);
  
      edge_derivatives(scene, adjoints[camId], edge_samples_in_total, d_mesh);
    }
  }

  virtual float d_render_and_compare(const Scene &scene, const CamInfo* cams, const Img *target_images, int a_viewsNum, 
                                     const int edge_samples_in_total, DScene &d_mesh, Img* outImages = nullptr) override
  {
    ::std::vector<Img> local_images(a_viewsNum);
    Img *images = outImages ? outImages : local_images.data();
    for(int i=0;i<a_viewsNum;i++)
      images[i].resize(cams[i].width,cams[i].height);

    commit(scene);
    render(scene, cams, images, a_viewsNum);

    ::std::vector<Img> adjoints(a_viewsNum);
    for(auto& im : adjoints)
      im = Img(images[0].width(), images[0].height(), float3{0, 0, 0});

    float mse = 0.0f;
    #pragma omp parallel for num_threads(a_viewsNum) reduction(+:mse)
    for(int i=0;i<a_viewsNum;i++)
      mse += LossAndDiffLoss(images[i], target_images[i], adjoints[i]);
    d_mesh.clear();
    d_render(scene, cams, adjoints.data(), a_viewsNum, edge_samples_in_total, d_mesh);
    return mse/float(a_viewsNum*images[0].width()*images[0].height());
  }

private:

  const Scene* m_pLastPreparedScene = nullptr;

  void interior_derivatives(const Scene &scene, const Img &adjoint,
                            DScene &d_mesh) 
  {
    auto sqrt_num_samples = (int)sqrt((float)m_samples_per_pixel);
    auto samples_per_pixel = sqrt_num_samples * sqrt_num_samples;
  
    DScene grads[MAXTHREADS];
    for(int i=0;i<MAXTHREADS;i++) {
      grads[i] = d_mesh;
      grads[i].clear(); // TODO: make this more effitiient
    }
    
    #if (DEBUG_RENDER==0)
    #pragma omp parallel for collapse (2) num_threads(MAXTHREADS)
    #endif
    for (int y = 0; y < adjoint.height(); y++) { // for each pixel  
      for (int x = 0; x < adjoint.width(); x++)  {
  
        for (int samId = 0; samId < samples_per_pixel; samId++) // for each subpixel
        {         
          float xoff = m_hammSamples[2*samId+0];
          float yoff = m_hammSamples[2*samId+1];
          
          const auto val = adjoint[int2(x,y)] / samples_per_pixel;
          shade_grad<material>(scene, m_pTracer.get(), float2(x + xoff, y + yoff), val, m_aux, grads[omp_get_thread_num()]);     
               
        } // for (int samId = 0; samId < samples_per_pixel; samId++)
      }
    }
  
    // accumulate gradient from different threads (parallel reduction/hist)
    //
    for(int i=0;i<MAXTHREADS;i++) 
      for (int m_n=0; m_n<d_mesh.get_dmeshes().size(); m_n++)
        for(size_t j=0;j<d_mesh.get_dmeshes()[m_n].full_size(); j++)
          d_mesh.get_dmeshes()[m_n].full_data()[j] += grads[i].get_dmeshes()[m_n].full_data()[j];
  }

  void edge_derivatives_was(const Scene &scene, const Img &adjoint,
                            DScene &d_mesh)
  {
    DScene grads[MAXTHREADS];
    for(int i=0;i<MAXTHREADS;i++) 
    {
      grads[i] = d_mesh;
      grads[i].clear(); // TODO: make this more effitiient
    }
    
    //#if (DEBUG_RENDER==0)
    //#pragma omp parallel for collapse (2) num_threads(MAXTHREADS)
    //#endif
    for (int y = 0; y < adjoint.height(); y++) 
    {
      for (int x = 0; x < adjoint.width(); x++)  
      {
        bool is_edge_pixel = true;
        if (!is_edge_pixel)
          continue;
        /*
        for (int samId = 0; samId < samples_per_pixel; samId++) // for each subpixel
        {         
          float xoff = m_hammSamples[2*samId+0];
          float yoff = m_hammSamples[2*samId+1];
          
          const auto val = adjoint[int2(x,y)] / samples_per_pixel;
          shade_grad<material>(scene, m_pTracer.get(), float2(x + xoff, y + yoff), val, m_aux, grads[omp_get_thread_num()]);     
               
        } // for (int samId = 0; samId < samples_per_pixel; samId++)
        */
      }
    }
  
    // accumulate gradient from different threads (parallel reduction/hist)
    //
    for(int i=0;i<MAXTHREADS;i++) 
      for (int m_n=0; m_n<d_mesh.get_dmeshes().size(); m_n++)
        for(size_t j=0;j<d_mesh.get_dmeshes()[m_n].full_size(); j++)
          d_mesh.get_dmeshes()[m_n].full_data()[j] += grads[i].get_dmeshes()[m_n].full_data()[j];
  }

  void edge_derivatives(
        const Scene &scene,
        const Img &adjoint,
        const int num_edge_samples,
        DScene &d_scene) 
  {
    // (1) We need to project 3d mesh to screen for correct edje sampling  
    //TODO: scene is prepared, we can project only the prepared arrays
    Scene scene2d_diff, scene2d_full;
    ::std::set<Edge> edges_s;
    {
      ::std::vector<int> diff_mesh_n = {0};
      ::std::vector<bool> diff_mesh_b(scene.get_meshes().size(), false);
      for (int i = 0; i< scene.get_meshes().size(); i++)
      {
        if (d_scene.get_dmesh(i))
          diff_mesh_b[i] = true;
      }
      for (int i=0; i<scene.get_meshes().size(); i++)
      {
        //TODO: support instancing for diff render
        //if (diff_mesh_b[i])
        //  assert(scene.get_transform(i).size() == 1);
        for (int j=0;j<scene.get_transform(i).size();j++)
        {
          TriangleMesh m = scene.get_mesh(i);
          transform(m, scene.get_transform(i)[j]);
          for (auto &v : m.vertices)
          {
            auto vCopy = v;
            VertexShader(*(m_aux.pCamInfo), vCopy.x, vCopy.y, vCopy.z, v.M);
          }

          scene2d_full.add_mesh(m);
          if (d_scene.get_dmesh(i))
          {
            //this mesh is differentiable, we shoud add it to diff scene and collect edges for it
            scene2d_diff.add_mesh(m);

            for (size_t k=0; k<scene.get_mesh(i).indices.size();k+=3) 
            {
              auto A = scene.get_index(i, j, k);
              auto B = scene.get_index(i, j, k+1);
              auto C = scene.get_index(i, j, k+2); 
              auto An = scene.get_vertex_n(i, k);
              auto Bn = scene.get_vertex_n(i, k+1);
              auto Cn = scene.get_vertex_n(i, k+2); 
              edges_s.insert(Edge(A, B, An, Bn, i, j));
              edges_s.insert(Edge(B, C, Bn, Cn, i, j));
              edges_s.insert(Edge(C, A, Cn, An, i, j));
              //logerr("tri (%d %d)(%d %d)(%d %d)", A, An, B, Bn, C, Cn);
            }
          }
        }
      }
    }
    ::std::vector<int> map = scene.get_instance_id_mapping();
    ::std::vector<int> map_2d(map.size(), 0);
    scene2d_full.set_instance_id_mapping(map_2d);
    scene2d_full.prepare_for_render();
  
    scene2d_diff.set_instance_id_mapping(map_2d);
    scene2d_diff.prepare_for_render();

    ::std::vector<Edge> edges = ::std::vector<Edge>(edges_s.begin(), edges_s.end());

    // (2) prepare edge sampler
    auto edge_sampler = build_edge_sampler(scene2d_diff, edges);
  
    // (3) do edje sampling
    // 
    prng::RandomGen gens[MAXTHREADS];
    ::std::vector<::std::vector<GradReal>> d_pos[MAXTHREADS];
    ::std::vector<::std::vector<GradReal>> d_tr[MAXTHREADS];
  
    for(int i=0;i<MAXTHREADS;i++)
    {
      gens [i] = prng::RandomGenInit(7777 + i*i + 1);
      for (auto &dmesh : d_scene.get_dmeshes())
      {
        d_pos[i].emplace_back();
        d_pos[i].back().resize(3*dmesh.vertex_count(), 0);
        d_tr[i].emplace_back();
        d_tr[i].back().resize(DMesh::TRANSFORM_SIZE*dmesh.instance_count(), 0);
      }
    }
  
    //float maxRelativeError = 0.0f;
    #if (DEBUG_RENDER==0)
    #pragma omp parallel for num_threads(MAXTHREADS)
    #endif
    for (int i = 0; i < num_edge_samples; i++) 
    { 
      auto& gen = gens[omp_get_thread_num()];

      const float rnd0 = prng::rndFloat(&gen);
      const float rnd1 = prng::rndFloat(&gen);
      // pick an edge
      auto edge_id = sample(edge_sampler, rnd0);
      auto edge    = edges[edge_id];
      auto pmf     = edge_sampler.pmf[edge_id];
      
      // pick a point p on the edge
      auto v0 = LiteMath::to_float2(scene2d_diff.get_pos(edge.v0));
      auto v1 = LiteMath::to_float2(scene2d_diff.get_pos(edge.v1));
      auto t = rnd1;
      auto p = v0 + t * (v1 - v0);
      int xi = int(p.x); 
      int yi = int(p.y); // integer coordinates
      if (xi < 0 || yi < 0 || xi >= adjoint.width() || yi >= adjoint.height() || length(adjoint[int2(xi,yi)]) < 1e-3)
        continue;

      // sample the two sides of the edge
      auto n = normal2D((v1 - v0) / length(v1 - v0));
      
      const float2 coordIn  = p - 1e-3f * n;
      const float2 coordOut = p + 1e-3f * n;

      const auto color_in  = shade<material>(scene2d_full, m_pTracer.get(), coordIn);
      const auto color_out = shade<material>(scene2d_full, m_pTracer.get(), coordOut);

      // get corresponding adjoint from the adjoint image,
      // multiply with the color difference and divide by the pdf & number of samples.
      float pdf    = pmf  / (length(v1 - v0));
      float weight = 1.0f / (pdf * float(num_edge_samples));
      float adj    = dot(color_in - color_out, adjoint[int2(xi,yi)]);
      
      if(adj*adj > 0.0f)
      {
        // the boundary point is p = v0 + t * (v1 - v0)
        // according to Reynolds transport theorem, the derivatives w.r.t. q is color_diff * dot(n, dp/dq)
        // dp/dv0.x = (1 - t, 0), dp/dv0.y = (0, 1 - t)
        // dp/dv1.x = (    t, 0), dp/dv1.y = (0,     t)
        
        auto d_v0 = float2{(1 - t) * n.x, (1 - t) * n.y} * adj * weight; // v0: (dF/dx_proj, dF/dy_proj)
        auto d_v1 = float2{     t  * n.x,      t  * n.y} * adj * weight; // v1: (dF/dx_proj, dF/dy_proj)
        
        edge_grad(scene, edge, d_v0, d_v1, m_aux, d_pos[omp_get_thread_num()], d_tr[omp_get_thread_num()]);
      }
    }    
  
    //::std::cout << " (VS_X_grad/VS_Y_grad) maxError = " << maxRelativeError*100.0f << "%" << ::std::endl;
  
    // accumulate gradient from different threads (parallel reduction/hist)
    //, logerr("dpos[%d %d] = %f",i,j,d_pos[i][j])
    int mesh_pos = 0;
    for (auto &dmesh : d_scene.get_dmeshes())
    {
      GradReal *accum_pos = dmesh.pos(0);
      if (accum_pos)
      {
        for(int i=0;i<MAXTHREADS;i++) 
          for(size_t j=0;j<3*dmesh.vertex_count(); j++)
            accum_pos[j] += d_pos[i][mesh_pos][j];
      }
      GradReal *accum_tr = dmesh.transform_mat(0);
      if (accum_tr)
      {
        //optimize transform matrix directly
        for(int i=0;i<MAXTHREADS;i++) 
          for(size_t j=0;j<DMesh::TRANSFORM_SIZE*dmesh.instance_count(); j++)
            accum_tr[j] += d_tr[i][mesh_pos][j];
      }
      else
      {
        //optimize simple transformation parameters (translation + rotation + scale)
        accum_tr = dmesh.restricted_transform(0);
        assert(accum_tr);
        ::std::vector<GradReal> accum_tr_mat(DMesh::TRANSFORM_SIZE*dmesh.instance_count(), 0);
        GradReal jac[DMesh::RESTRICTED_TRANSFORM_SIZE][DMesh::TRANSFORM_SIZE];
      
        for(int i=0;i<MAXTHREADS;i++) 
          for(size_t j=0;j<DMesh::TRANSFORM_SIZE*dmesh.instance_count(); j++)
            accum_tr_mat[j] += d_tr[i][mesh_pos][j];

        for(size_t j=0;j<dmesh.instance_count(); j++)
        {
          //calculate jacobian
          for (int k=0;k<DMesh::RESTRICTED_TRANSFORM_SIZE; k++)
            for (int l=0;l<DMesh::TRANSFORM_SIZE; l++)
              jac[k][l] = 0;
          dmat_dtransforms_jac(scene.get_restricted_transform(dmesh.get_mesh_id())[j].data.M, jac);

          for (int k=0;k<DMesh::RESTRICTED_TRANSFORM_SIZE; k++)
            for (int l=0;l<DMesh::TRANSFORM_SIZE; l++)
              accum_tr[DMesh::RESTRICTED_TRANSFORM_SIZE*j + k] += jac[k][l]*accum_tr_mat[DMesh::TRANSFORM_SIZE*j + l];
        }
      }
      mesh_pos++;
    }
  }


  ::std::shared_ptr<IRayTracer> m_pTracer = nullptr;
  int m_samples_per_pixel;
  AuxData m_aux;

  unsigned int m_table[qmc::QRNG_DIMENSIONS][qmc::QRNG_RESOLUTION];
  ::std::vector<float> m_hammSamples;
};


::std::shared_ptr<IDiffRender> MakeDifferentialRenderer(const DiffRenderSettings &settings);
}