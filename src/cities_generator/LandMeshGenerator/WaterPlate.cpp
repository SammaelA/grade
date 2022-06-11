#include "WaterPlate.h"
#include "cities_generator/ImageLoader.h"
#include "cities_generator/GameManager.h"

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

WaterPlate::WaterPlate(glm::vec3 pos, glm::vec2 size)
{
    position = glm::vec3(pos.x + size.x / 2, pos.y, pos.z + size.y / 2);
    scale = glm::vec3(size.x, 1, size.y);
    timer = 0;

    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("waterPlateShader");
    
    GLfloat vertexes[] = {
        -0.5f,  0.0f, -0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
         0.5f,  0.0f, -0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
         0.5f,  0.0f,  0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,

         0.5f,  0.0f,  0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        -0.5f,  0.0f,  0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        -0.5f,  0.0f, -0.5f,  // 0.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
    };

    vertexesNum = sizeof(vertexes) / sizeof(vertexes[0]);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertexes), vertexes, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
//     glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
//     glEnableVertexAttribArray(1);
//     glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
//     glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    //CREATE FB TEXTURE
	glGenTextures(1, &depthTexture);
	glBindTexture(GL_TEXTURE_2D, depthTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, SCREEN_WIDTH, SCREEN_HEIGHT, 0,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
    glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &aboveWaterTexture);
    glBindTexture(GL_TEXTURE_2D, aboveWaterTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCREEN_WIDTH, SCREEN_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &bloomTexture);
    glBindTexture(GL_TEXTURE_2D, bloomTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCREEN_WIDTH, SCREEN_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenTextures(1, &underWaterTexture);
    glBindTexture(GL_TEXTURE_2D, underWaterTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCREEN_WIDTH, SCREEN_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    //CREATE FBO AND ATTACH TEXTURE TO IT
	glGenFramebuffers(1, &FBO);
	glBindFramebuffer(GL_FRAMEBUFFER, FBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,GL_TEXTURE_2D, depthTexture,0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, (GLuint)aboveWaterTexture, 0); 

    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE)
    {
		debug("FBO NO OK!!");
        throw std::exception{};
    }

    //DUDV and normal maps
    DUDVTexture = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/Water/waterDUDV.png"
    );
    glBindTexture(GL_TEXTURE_2D, DUDVTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glBindTexture(GL_TEXTURE_2D, 0);

    normalMapTexture = ImageLoader::StbLoadImage(
        std::string(ASSETS_FOLDER) + "Textures/Water/waterNormal.png"
    );
    glBindTexture(GL_TEXTURE_2D, normalMapTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glBindTexture(GL_TEXTURE_2D, 0);
}

WaterPlate::~WaterPlate()
{}

void WaterPlate::Render(RENDER_MODE rmode)
{
    UseShader();

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    glEnable(GL_BLEND);
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ZERO, GL_ONE);

    // // Renderer::GetInstance()->MyBufferClearSystem(0,0,0);

    
    // // GLint TransformUniformPos;
    // // bool firstSprite = true;

    //Camera loading
    glm::mat4 viewProj = glm::mat4(1.0f);
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    viewProj = cam->GetViewProjectionMatrix();
    glUniformMatrix4fv(
        glGetUniformLocation(m_shader->Program, "view_proj_matrix"),
        1, GL_FALSE, glm::value_ptr(viewProj)
    );
    glUniform3fv(
        glGetUniformLocation(m_shader->Program, "camera_pos"), 
        1, glm::value_ptr(cam->GetPosition())
    );
    glUniform2fv(
        glGetUniformLocation(m_shader->Program, "z_near_far"), 
        1, glm::value_ptr(cam->GetZDists())
    );
    glUniformMatrix4fv(
        glGetUniformLocation(m_shader->Program, "inv_proj_matrix"),
        1, GL_FALSE, glm::value_ptr(cam->GetInvProjMatrix())
    );



    glUniform2fv(
        glGetUniformLocation(m_shader->Program, "window_size"),
        1, glm::value_ptr(glm::vec2 {SCREEN_WIDTH, SCREEN_HEIGHT})
    );


    timer += GameManager::GetInstance()->timeDelta;
    glUniform1f(glGetUniformLocation(m_shader->Program, "time_seconds"), timer);


    Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
            2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];
    TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
    SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();
    glUniform3fv(
        glGetUniformLocation(m_shader->Program, "sun_pos"), 
        1, glm::value_ptr(sunTransform->position)
    );
    glUniform3fv(
        glGetUniformLocation(m_shader->Program, "sun_color"), 
        1, glm::value_ptr(sunComp->color)
    );

    
    //     GLint ViewProjMatrixUniformPos, GlobalIlluminationPointUniformPos,
    //             GlobalIlluminationColorUniformPos, diffuseColorPos, isIsotropicColorPos,
    //             isDepthOnlyPos, GlobalIlluminationViewProjUniformPos, CameraPositionPos;
        
    //     ViewProjMatrixUniformPos = glGetUniformLocation(model->m_shader->Program, "view_proj_matrix");
    //     GlobalIlluminationPointUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_pos");
    //     GlobalIlluminationColorUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_color"); 
    //     GlobalIlluminationViewProjUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_view_proj_matrix");
    //     diffuseColorPos = glGetUniformLocation(model->m_shader->Program, "diffuse_color"); 
    //     isIsotropicColorPos = glGetUniformLocation(model->m_shader->Program, "is_isotropic_color"); 
    //     isDepthOnlyPos = glGetUniformLocation(model->m_shader->Program, "is_depth_only"); 
        // CameraPositionPos = glGetUniformLocation(model->m_shader->Program, "camera_pos"); 

    //     // std::vector<Entity*> sunLightEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
    //     //     1, 0, std::vector<unsigned>{SunLightComponent::GetCompId()}, std::vector<std::string>{}, true));
    //     // SunLightComponent* sunS = sunLightEntities[0]->GetComponent<SunLightComponent>();
    //     // TransformComponent* sunT = sunLightEntities[0]->GetComponent<TransformComponent>();
        

    //     model->LoadInShader_Scene3DObject();
    //     glUniformMatrix4fv(ViewProjMatrixUniformPos, 1, GL_FALSE, glm::value_ptr(viewProj));
    //     glUniformMatrix4fv(GlobalIlluminationViewProjUniformPos, 1, GL_FALSE, glm::value_ptr(sunViewProj));
    
    //     glUniform3f(CameraPositionPos, cam->position.x, cam->position.y, cam->position.z);
    //     glUniform1i(isIsotropicColorPos, !model->isUsingDiffuseTexture);
    //     glUniform1i(isDepthOnlyPos, false);
    //     glUniform4f(diffuseColorPos, model->diffuseColor.x, model->diffuseColor.y, model->diffuseColor.z, model->diffuseColor.z);
    //     glUniform3f(glGetUniformLocation(model->m_shader->Program, "sphere_center"),
    //         model->position.x, model->position.y, model->position.z);

    //TEXTURES
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, underWaterTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "under_water_texture"), 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, aboveWaterTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "above_water_texture"), 1);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, DUDVTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "dudv_texture"), 2);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, normalMapTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "normal_texture"), 3);
    glActiveTexture(GL_TEXTURE4);
    glBindTexture(GL_TEXTURE_2D, bloomTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "bloom_texture"), 4);
    glActiveTexture(GL_TEXTURE5);
    glBindTexture(GL_TEXTURE_2D, depthTexture);
    glUniform1i(glGetUniformLocation(m_shader->Program, "depth_under_water_texture"), 5);
        
    //Transform matrix loading
    LoadInShader_Scene3DObject();

    //DRAW CALL
    RenderableObjectDataHolder::Render(rmode);

    //RETURN OPTION
    glDisable(GL_BLEND);
}