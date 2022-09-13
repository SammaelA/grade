#include "cities_generator/Renderer.h"
#include "cities_generator/ECSClass.h"
#include <memory>

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

Renderer* Renderer::instance = nullptr;
std::vector<Renderer::DebugSphere*> Renderer::DebugSphere::spheres{};

Renderer* Renderer::GetInstance()
{
    if (instance == nullptr)
        instance = new Renderer();
    return instance;
}

void Renderer::MyBufferClearSystem(GLfloat r, GLfloat g, GLfloat b)
{
    glClearColor(r,g,b,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Renderer::RenderTexture(GLuint tex, GLfloat scaleFrom, GLfloat scaleTo, bool keepRatio)
{
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    
    const std::string tag = "TextureDisplayObject";
    Entity* displayer = Entity::AnyObjectWithTag(tag, true);
    if (displayer == nullptr)
    {
        displayer = new Entity();
        displayer->giveTag(tag);
        // displayer->components.push_back(new ImageVertexInfoComponent(displayer->id));
        displayer->AddComponent(new ImageVertexInfoComponent(displayer->id));
    }
    ImageVertexInfoComponent* displayComponent = displayer->GetComponent<ImageVertexInfoComponent>();
    
    // debug(tex, scaleFrom, scaleTo, keepRatio);

    float scaleRatio[] = {1.0, 1.0};
    if (keepRatio)
    {
        glBindTexture(GL_TEXTURE_2D, tex);
        int w, h;
        int miplevel = 0;
        glGetTexLevelParameteriv(GL_TEXTURE_2D, miplevel, GL_TEXTURE_WIDTH, &w);
        glGetTexLevelParameteriv(GL_TEXTURE_2D, miplevel, GL_TEXTURE_HEIGHT, &h);
        scaleRatio[0] = std::max((float)SCREEN_WIDTH / SCREEN_HEIGHT * h / w, 1.0f);
        scaleRatio[1] = std::max((float)SCREEN_HEIGHT / SCREEN_WIDTH * w / h, 1.0f);
    } 


    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    Renderer::MyBufferClearSystem(0,0,0);
    displayComponent->shader->Use();
    glBindVertexArray(displayComponent->VAO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, displayComponent->EBO);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glUniform1i(glGetUniformLocation(displayComponent->shader->Program, "Texture"), 0);
    glUniform1f(glGetUniformLocation(displayComponent->shader->Program, "scaleFrom"), scaleFrom);
    glUniform1f(glGetUniformLocation(displayComponent->shader->Program, "scaleTo"), scaleTo);
    glUniform2f(glGetUniformLocation(displayComponent->shader->Program, "scaleRatio"), 
        scaleRatio[0], scaleRatio[1]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void Renderer::CreateDebugSphere(float f, glm::vec3 pos, glm::vec3 col, std::string name)
{
    new DebugSphere(f, pos, col, name);
}

Renderer::DebugSphere::DebugSphere(float Rad, glm::vec3 Pos, glm::vec3 Col, std::string Name)
{
    radius = Rad;
    name = Name;
    col = Col;

    if (name.length() == 0)
    {
        int counter = spheres.size();
        name = std::to_string(counter);
        int i = 0;
        while (i < spheres.size())
        {
            if (spheres[i]->name == name)
            {
                name = std::to_string(++counter);
                i = 0;
            }
            i++;
        }
    }
    else
    {
        for (DebugSphere* sp: spheres)
        {
            if (sp->name == Name)
            {
                debug("SUCH DEBUG SPHERE ALREADY EXISTS!");
                throw std::exception{};
            }
        }
    }

    spheres.push_back(this);
    Entity* sphere = new Entity();
    sphere->LoadGltfModel
        (std::string(ASSETS_FOLDER)+"Common/DebugSphere/DebugSphere.gltf", "Sphere", "debugSphereShader");
    
    if (col.x >= 0 || col.y >= 0 || col.z >= 0)
    {
        sphere->GetComponent<ModelInfoComponent>()->SetColor(glm::vec4(col, 1.0));
    }
    sphere->GetComponent<ModelInfoComponent>()->position = glm::vec3(Pos.x, Pos.y, Pos.z);
    sphere->GetComponent<ModelInfoComponent>()->scale = glm::vec3(radius,radius,radius);
    entityPointer = sphere;
}

Renderer::DebugSphere* Renderer::DebugSphere::GetSphereByName(std::string _name, bool safe)
{
    for (DebugSphere* sphere: spheres)
    {
        if (sphere->name == _name)
            return sphere;
    }
    if (!safe)
    {
        debug("NO SPHERE WITH THIS NAME! " + _name);
        throw std::exception{};
    }
    return nullptr;
}

void Renderer::DebugSphereMove(glm::vec3 _pos, std::string _name)
{
    DebugSphere* sphere = DebugSphere::GetSphereByName(_name);
    sphere->entityPointer->GetComponent<ModelInfoComponent>()->position = _pos;
}

void Renderer::DebugSphereSetColor(glm::vec3 _col, std::string _name)
{
    DebugSphere* sphere = DebugSphere::GetSphereByName(_name);
    sphere->col = _col;
    sphere->entityPointer->GetComponent<ModelInfoComponent>()->SetColor(glm::vec4(_col, 1.0));
}

void Renderer::Destroy()
{
    DebugSphere::Clear();
}

void Renderer::DebugSphere::Clear()
{
    for (auto sphere: spheres)
    {
        delete sphere;
    }
}

void Renderer::RenderDebugSpheres()
{
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // Renderer::GetInstance()->MyBufferClearSystem(0,0,0);

    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // GLint TransformUniformPos;
    // bool firstSprite = true;
                                            
    glm::mat4 viewProj = glm::mat4(1.0f);
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    
    viewProj = cam->GetViewProjectionMatrix();
             
    for (size_t i = 0; i < DebugSphere::spheres.size(); i++)
    {
        ModelInfoComponent *model = DebugSphere::spheres[i]->
            entityPointer->GetComponent<ModelInfoComponent>();
        
        Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
                2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];

        TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
        SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();

        model->m_shader->Use();

        GLint ViewProjMatrixUniformPos=0, GlobalIlluminationPointUniformPos=0,
                GlobalIlluminationColorUniformPos=0, diffuseColorPos=0, isIsotropicColorPos=0,
                isDepthOnlyPos=0, GlobalIlluminationViewProjUniformPos=0, CameraPositionPos=0;
        
        ViewProjMatrixUniformPos = glGetUniformLocation(model->m_shader->Program, "view_proj_matrix");
        GlobalIlluminationPointUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_pos");
        GlobalIlluminationColorUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_color"); 
        GlobalIlluminationViewProjUniformPos = glGetUniformLocation(model->m_shader->Program, "sun_view_proj_matrix");
        isDepthOnlyPos = glGetUniformLocation(model->m_shader->Program, "is_depth_only"); 
        CameraPositionPos = glGetUniformLocation(model->m_shader->Program, "camera_position"); 

        // std::vector<Entity*> sunLightEntities = std::move(Entity::AllEntitiesWithComponentsAndTags(
        //     1, 0, std::vector<unsigned>{SunLightComponent::GetCompId()}, std::vector<std::string>{}, true));
        // SunLightComponent* sunS = sunLightEntities[0]->GetComponent<SunLightComponent>();
        // TransformComponent* sunT = sunLightEntities[0]->GetComponent<TransformComponent>();
        glm::mat4 sunProjM = sunComp->projectionMatrix;
        glm::mat4 sunVeiwM = sunTransform->GetViewMatrixFromThisPoint();
        glm::mat4 sunViewProj = sunProjM * sunVeiwM;

        model->LoadInShader_Scene3DObject();
        glUniformMatrix4fv(ViewProjMatrixUniformPos, 1, GL_FALSE, glm::value_ptr(viewProj));
        glUniformMatrix4fv(GlobalIlluminationViewProjUniformPos, 1, GL_FALSE, glm::value_ptr(sunViewProj));
        glUniform3f(GlobalIlluminationPointUniformPos, sunTransform->position.x, sunTransform->position.y, sunTransform->position.z);
        glUniform3f(GlobalIlluminationColorUniformPos, sunComp->color.x, sunComp->color.y, sunComp->color.z);
        glUniform3fv(CameraPositionPos, 1, glm::value_ptr(cam->GetPosition()));
        glUniform1i(isIsotropicColorPos, !model->isUsingDiffuseTexture);
        glUniform1i(isDepthOnlyPos, false);
        glUniform4f(diffuseColorPos, model->diffuseColor.x, model->diffuseColor.y, model->diffuseColor.z, model->diffuseColor.z);
        glUniform3f(glGetUniformLocation(model->m_shader->Program, "sphere_center"),
            model->position.x, model->position.y, model->position.z);

        //TEXTURES
        int TEX_COUNT = 0;
        model->LoadInShader_BasicTexturedObject(TEX_COUNT);
        glActiveTexture(GL_TEXTURE0 + TEX_COUNT);
        // glBindTexture(GL_TEXTURE_2D, sunComp->shadowTexture);
        // glUniform1i(glGetUniformLocation(model->m_shader->Program, "shadow_texture"), TEX_COUNT);
        
        //DRAW CALL
        model->Render(RENDER_MODE::DEFUALT);
    }
}

Renderer::DebugSphere::~DebugSphere()
{}

void Renderer::SetFullViewport()
{
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
}

bool Renderer::DebugSphereDoesExist(std::string name)
{
    return DebugSphere::GetSphereByName(name, true) != nullptr;
}