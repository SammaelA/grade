#include "shadow.h"
#include "graphics_utils/texture_manager.h"
    void ShadowMap::use(DirectedLight &light)
    {
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp, 0);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);    
        shadow_camera.pos = 1500.f*light.dir;
        shadow_camera.front =  -shadow_camera.pos;
        shadow_camera.up = glm::vec3( 0.0f, 1.0f,  0.0f);
        view = shadow_camera.camera();
        float near_plane = 1.0f, far_plane = 7.5f;
        projection = glm::ortho(-10.0f, 10.0f, -10.0f, 10.0f, near_plane, far_plane); 
        projection = glm::perspective(glm::radians(90.0f), (float)SHADOW_WIDTH / SHADOW_HEIGHT, 1.0f, 3000.0f);
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
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, srcDepthTex, 0);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);  
    }
    void ShadowMap::finish_trans_pass()
    {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glEnable(GL_DEPTH_TEST); 
        glDisable(GL_BLEND); 
    }
    void ShadowMap::create(int w, int h)
    {
        postFx = new PostFx("gaussian_blur.fs");
        SHADOW_WIDTH = w;
        SHADOW_HEIGHT = h;
        glGenFramebuffers(1, &depthMapFBO);  

        glGenTextures(1, &depthMap);
        glBindTexture(GL_TEXTURE_2D, depthMap);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);  
        
        glGenTextures(1, &VSMdepthTex);
        glBindTexture(GL_TEXTURE_2D, VSMdepthTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor); 
       
        glGenTextures(1, &srcDepthTex);
        glBindTexture(GL_TEXTURE_2D, srcDepthTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glGenTextures(1, &VSMdepthTexTemp);
        glBindTexture(GL_TEXTURE_2D, VSMdepthTexTemp);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp, 0);
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            print_FB_status(glCheckFramebufferStatus(GL_FRAMEBUFFER));
        }
        else
        {
            debugl(10,"Shadow map created %d %d",SHADOW_HEIGHT,SHADOW_WIDTH);
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);  
    }
    ShadowMap::~ShadowMap() 
    {
        delete postFx;
        glDeleteTextures(1,&depthMap);
        glDeleteTextures(1,&VSMdepthTex);
    }
    void ShadowMap::blur()
    {
        PostFx merge = PostFx("depth_merge.fs");
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glDisable(GL_DEPTH_TEST);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTex, 0);
        merge.use();
        merge.get_shader().texture("trans",srcDepthTex);
        merge.get_shader().texture("depth",VSMdepthTexTemp);
        merge.render();

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTexTemp, 0);
        postFx->use();
        postFx->get_shader().texture("tex",VSMdepthTex);
        postFx->get_shader().uniform("pass",0);
        postFx->get_shader().uniform("tex_size_inv",glm::vec2(1.0/SHADOW_WIDTH, 1.0/SHADOW_HEIGHT));
        postFx->render();

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, VSMdepthTex, 0);
        postFx->use();
        postFx->get_shader().texture("tex",VSMdepthTexTemp);
        postFx->get_shader().uniform("pass",1);
        postFx->get_shader().uniform("tex_size_inv",glm::vec2(1.0/SHADOW_WIDTH, 1.0/SHADOW_HEIGHT));
        postFx->render();

        glEnable(GL_DEPTH_TEST);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        //glDeleteTextures(1,&VSMdepthTexTemp);
        //glDeleteTextures(1,&srcDepthTex);
    }