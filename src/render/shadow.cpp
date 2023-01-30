#include "shadow.h"
#include "tinyEngine/engine.h"
    void ShadowMap::use(DirectedLight &light)
    {
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp.texture, 0);
        glClearColor(1, 0, 0, 1);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);    
        shadow_camera.pos = 1500.f*light.dir;
        shadow_camera.front =  -shadow_camera.pos;
        shadow_camera.up = glm::vec3( 0.0f, 1.0f,  0.0f);
        view = shadow_camera.camera();
        float near_plane = 100, far_plane = 2000;
        projection = glm::ortho(-1000.0f, 1000.0f, -1000.0f, 1000.0f, near_plane, far_plane); 
        //projection = glm::perspective(glm::radians(90.0f), (float)SHADOW_WIDTH / SHADOW_HEIGHT, 1.0f, 3000.0f);
        viewproj = projection * view;

    }
    glm::mat4 &ShadowMap::get_transform()
    {   
        return viewproj;
    }
    void ShadowMap::start_trans_pass()
    {
        glEnable(GL_BLEND); 
        glDisable(GL_DEPTH_TEST); 
        glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBlendEquationSeparate(GL_MIN, GL_FUNC_ADD);
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, srcDepthTex.texture, 0);
        glClearColor(1, 0, 0, 1);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);  
    }
    void ShadowMap::finish_trans_pass()
    {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBlendFuncSeparate(GL_ONE, GL_ZERO, GL_ONE, GL_ZERO);
        glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
        glEnable(GL_DEPTH_TEST); 
        glDisable(GL_BLEND); 
    }
    void ShadowMap::create(int w, int h)
    {
        postFx = new PostFx("gaussian_blur.fs");
        SHADOW_WIDTH = w;
        SHADOW_HEIGHT = h;
        depthMapFBO = create_framebuffer();

        depthMap = engine::textureManager->create_texture(SHADOW_WIDTH, SHADOW_HEIGHT, GL_DEPTH_COMPONENT16, 1, nullptr, GL_DEPTH_COMPONENT, GL_FLOAT);
        glBindTexture(GL_TEXTURE_2D, depthMap.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);  
        
        VSMdepthTex = engine::textureManager->create_texture(SHADOW_WIDTH, SHADOW_HEIGHT, GL_RGBA32F, 1, nullptr, GL_RGBA, GL_UNSIGNED_BYTE);
        glBindTexture(GL_TEXTURE_2D, VSMdepthTex.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor); 
       
        srcDepthTex = engine::textureManager->create_texture(SHADOW_WIDTH, SHADOW_HEIGHT, GL_RGBA32F, 1, nullptr, GL_RGBA, GL_UNSIGNED_BYTE);
        glBindTexture(GL_TEXTURE_2D, srcDepthTex.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        VSMdepthTexTemp = engine::textureManager->create_texture(SHADOW_WIDTH, SHADOW_HEIGHT, GL_RGBA32F, 1, nullptr, GL_RGBA, GL_UNSIGNED_BYTE);
        glBindTexture(GL_TEXTURE_2D, VSMdepthTexTemp.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        int prev_FBO = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap.texture, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp.texture, 0);
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
        }
        else
        {
            debugl(10,"Shadow map created %d %d",SHADOW_HEIGHT,SHADOW_WIDTH);
        }
        glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);  
    }
    ShadowMap::~ShadowMap() 
    {
        delete postFx;
        engine::textureManager->delete_tex(depthMap);
        engine::textureManager->delete_tex(VSMdepthTex);
    }
    void ShadowMap::blur()
    {
        PostFx merge = PostFx("depth_merge.fs");
        int prev_FBO = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prev_FBO);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glDisable(GL_DEPTH_TEST);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTex.texture, 0);
        merge.use();
        merge.get_shader().texture("trans",srcDepthTex);
        merge.get_shader().texture("depth",VSMdepthTexTemp);
        merge.render();

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp.texture, 0);
        postFx->use();
        postFx->get_shader().texture("tex",VSMdepthTex);
        postFx->get_shader().uniform("pass",0);
        postFx->get_shader().uniform("tex_size_inv",glm::vec2(0.1/SHADOW_WIDTH, 0.1/SHADOW_HEIGHT));
        postFx->render();

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTex.texture, 0);
        postFx->use();
        postFx->get_shader().texture("tex",VSMdepthTexTemp);
        postFx->get_shader().uniform("pass",1);
        postFx->get_shader().uniform("tex_size_inv",glm::vec2(0.1/SHADOW_WIDTH, 0.1/SHADOW_HEIGHT));
        postFx->render();

        glEnable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_FRAMEBUFFER, prev_FBO);

        //glDeleteTextures(1,&VSMdepthTexTemp);
        //glDeleteTextures(1,&srcDepthTex);
    }