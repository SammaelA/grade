#include "BasicGLTFModel.h"
#include "../../ECSClass.h"
#include "../dependencies/JSON/json/json.h"
#include "../dependencies/JSON/json/json-forwards.h"
#include "../dependencies/JSON/jsoncpp.cpp"
#include "glm/gtx/quaternion.hpp"
#include "third_party/stb_image.h"

extern GLuint SCREEN_WIDTH, SCREEN_HEIGHT;

BasicGLTFModel* BasicGLTFModel::LoadGltfModel
    (std::string path, std::string obj_name, std::string shaderName)
{
    BasicGLTFModel* model = new BasicGLTFModel();

    if (!IS_DERIVED(model, RenderableObjectDataHolder))
        throw std::exception{};

    model->m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer(shaderName);
    model->Set_IBO_element_type(GL_UNSIGNED_SHORT);
    
    const std::string gltf_file_name = path;

    std::ifstream model_file(gltf_file_name);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(model_file, obj);
    std::string bin_name = obj["buffers"][0]["uri"].asString();
    std::string bin_file_name = gltf_file_name;
    bin_file_name.erase(bin_file_name.rfind("/") + 1, bin_file_name.size());
    bin_file_name.append(bin_name);
    FILE *f = fopen(bin_file_name.c_str(), "rb");

    const Json::Value& nodes = obj["nodes"];
    int i;
    for (i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]["name"] == obj_name)
        {
            break;
        }
    }
    if (i == nodes.size())
    {
        debug("WRONG GLTF NAME");
        throw std::exception{};
    }
    
    if (IS_DERIVED(model, Scene3DObjectDataHolder))
    {
        if (nodes[i]["translation"].isArray())
        {
            const Json::Value& translation = nodes[i]["translation"];
            for (int j = 0; j < 3; j++)
            {
                model->position[j] = translation[j].asFloat();
            }  
        }
        
        if (nodes[i]["rotation"].isArray())
        {
            const Json::Value& rotation = nodes[i]["rotation"];
            glm::quat Qu;
            Qu.w = rotation[3].asFloat();
            Qu.x = rotation[0].asFloat();
            Qu.y = rotation[1].asFloat();
            Qu.z = rotation[2].asFloat();
            model->rotation = glm::mat4_cast(Qu); 
        }

        if (nodes[i]["scale"].isArray())
        {
            const Json::Value& scale = nodes[i]["scale"];
            for (int j = 0; j < 3; j++)
            {
                model->scale[j] = scale[j].asFloat(); 
            }
        }
    }   

    const Json::Value& meshes = obj["meshes"];
    const Json::Value& mesh = meshes[nodes[i]["mesh"].asUInt()];
    const Json::Value& primitive = mesh["primitives"][0];
    const Json::Value& material = obj["materials"][primitive["material"].asUInt()];

    unsigned Vertex_accessor, Normal_accessor, Texcoord_accessor, Index_accessor;
    Vertex_accessor = primitive["attributes"]["POSITION"].asUInt();
    Normal_accessor = primitive["attributes"]["NORMAL"].asUInt();
    Texcoord_accessor = primitive["attributes"]["TEXCOORD_0"].asUInt();
    Index_accessor = primitive["indices"].asUInt();


    const Json::Value& accessors = obj["accessors"];
    const Json::Value& bufferViews = obj["bufferViews"];
    //5126 - FLOAT 5123 - UNSIGNE SHORT

    //VERTEXES
    const Json::Value* accessor = &accessors[Vertex_accessor];
    const Json::Value* bufferView = &bufferViews[(*accessor)["bufferView"].asUInt()];
    model->vertexesNum = (*accessor)["count"].asUInt();
    GLfloat* vertexes = new GLfloat[model->vertexesNum * 3];
    unsigned offset = (*bufferView)["byteOffset"].asUInt() + (*accessor)["byteOffset"].asUInt();//in bytes
    unsigned length = model->vertexesNum * 3; // in elements
    fseek(f, offset, SEEK_SET);
    fread(vertexes, sizeof(vertexes[0]), length, f);
    
    // std::cout << offset << " " << length << std::endl;

    //NORMALS
    accessor = &accessors[Normal_accessor];
    bufferView = &bufferViews[(*accessor)["bufferView"].asUInt()];
    GLfloat* normals = new GLfloat[model->vertexesNum * 3];
    offset = (*bufferView)["byteOffset"].asUInt() + (*accessor)["byteOffset"].asUInt(); //in bytes
    length = model->vertexesNum * 3; // in elements
    fseek(f, offset, SEEK_SET);
    fread(normals, sizeof(normals[0]), length, f);

    //Texture Coordinates
    accessor = &accessors[Texcoord_accessor];
    bufferView = &bufferViews[(*accessor)["bufferView"].asUInt()];
    GLfloat* TCoords = new GLfloat[model->vertexesNum * 2];
    offset = (*bufferView)["byteOffset"].asUInt() + (*accessor)["byteOffset"].asUInt(); //in bytes
    length = model->vertexesNum * 2; // in elements
    fseek(f, offset, SEEK_SET);
    fread(TCoords, sizeof(TCoords[0]), length, f);

    //INDICES
    accessor = &accessors[Index_accessor];
    bufferView = &bufferViews[(*accessor)["bufferView"].asUInt()];
    model->indicesNum = (*accessor)["count"].asUInt();
    GLushort* indices = new GLushort[model->indicesNum];
    offset = (*bufferView)["byteOffset"].asUInt() + (*accessor)["byteOffset"].asUInt(); //in bytes
    length = model->indicesNum; // in elements
    fseek(f, offset, SEEK_SET);
    fread(indices, sizeof(indices[0]), length, f);

    //INDICES AND VERTEXES
    GLfloat* generalBuffer = new GLfloat[model->vertexesNum * 8];
    for (int i = 0; i < model->vertexesNum; i++)
    {
        memcpy(&(generalBuffer[i*8]), &(vertexes[i*3]), sizeof(generalBuffer[0]) * 3);
        memcpy(&(generalBuffer[i*8 + 3]), &(normals[i*3]), sizeof(generalBuffer[0]) * 3);
        memcpy(&(generalBuffer[i*8 + 6]), &(TCoords[i*2]), sizeof(generalBuffer[0]) * 2);
    }
    
    //LOADING TEXTURES
    std::string path_to_models = path;
    path_to_models.erase(path_to_models.rfind("/") + 1, path_to_models.size());

    const Json::Value& images = obj["images"];
    if (material["pbrMetallicRoughness"].isMember("baseColorTexture"))
    {
        model->isUsingDiffuseTexture = true;
        unsigned diffuseId = obj["textures"][material["pbrMetallicRoughness"]["baseColorTexture"]["index"].asUInt()]["source"].asUInt();
        const Json::Value& image = images[diffuseId];

        int imageW, imageH;
        std::string textureName = path_to_models;
        textureName.append(image["uri"].asString());
        unsigned char* diffuseTex = stbi_load(textureName.c_str(), &imageW, &imageH, 0, 4);
        if (!diffuseTex)
        {
            fprintf(stderr,"Diffuse texture load failed");
        }
        
        glGenTextures(1, &model->diffuseTexture);
        glBindTexture(GL_TEXTURE_2D, model->diffuseTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageW, imageH, 0, GL_RGBA, GL_UNSIGNED_BYTE, diffuseTex);
        glGenerateMipmap(GL_TEXTURE_2D);
        stbi_image_free(diffuseTex);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    else if (material["pbrMetallicRoughness"].isMember("baseColorFactor"))
    {
        model->isUsingDiffuseTexture = false;
        const Json::Value& color = material["pbrMetallicRoughness"]["baseColorFactor"];
        for (int i = 0; i < 4; i++)
        {
            model->diffuseColor[i] = color[i].asFloat();
        }
    }

    if (material["pbrMetallicRoughness"].isMember("roughnessFactor"))
    {
        model->roughnessFactor = material["pbrMetallicRoughness"]["roughnessFactor"].asFloat();
    }


    //LOADING BUFFERS
    glBindVertexArray(model->VAO);
    glBindBuffer(GL_ARRAY_BUFFER, model->VBO);
    glBufferData(GL_ARRAY_BUFFER, model->vertexesNum * sizeof(generalBuffer[0]) * 8, generalBuffer, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model->IBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, model->indicesNum * sizeof(indices[0]), indices, GL_STATIC_DRAW);

    // glBindBuffer(GL_ARRAY_BUFFER, model->VBO);
    // glBufferData(GL_ARRAY_BUFFER, model->vertexesNum * sizeof(generalBuffer[0]) * 3, generalBuffer, GL_STATIC_DRAW);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, model->IBO);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, model->indicesNum * sizeof(indices[0]), indices, GL_STATIC_DRAW);
    
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    // glEnableVertexAttribArray(0);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*)(6 * sizeof(GLfloat)));
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    delete [] generalBuffer;
    delete [] vertexes;
    delete [] normals;
    delete [] indices;
    delete [] TCoords;

    return model;
}

std::vector<glm::vec3> BasicGLTFModel::GetGltfVertexes(std::string path, std::string obj_name)
{
    const std::string gltf_file_name = path;

    std::ifstream model_file(gltf_file_name);
    Json::Reader reader;
    Json::Value obj;
    reader.parse(model_file, obj);
    std::string bin_name = obj["buffers"][0]["uri"].asString();
    std::string bin_file_name = gltf_file_name;
    bin_file_name.erase(bin_file_name.rfind("/") + 1, bin_file_name.size());
    bin_file_name.append(bin_name);
    FILE *f = fopen(bin_file_name.c_str(), "rb");

    const Json::Value& nodes = obj["nodes"];
    int i;
    for (i = 0; i < nodes.size(); i++)
    {
        if (nodes[i]["name"] == obj_name)
        {
            break;
        }
    }
    if (i == nodes.size())
    {
        debug("WRONG GLTF NAME");
        throw std::exception{};
    } 

    const Json::Value& meshes = obj["meshes"];
    const Json::Value& mesh = meshes[nodes[i]["mesh"].asUInt()];
    const Json::Value& primitive = mesh["primitives"][0];

    unsigned Vertex_accessor;
    Vertex_accessor = primitive["attributes"]["POSITION"].asUInt();

    const Json::Value& accessors = obj["accessors"];
    const Json::Value& bufferViews = obj["bufferViews"];
    //5126 - FLOAT 5123 - UNSIGNE SHORT

    //VERTEXES
    const Json::Value* accessor = &accessors[Vertex_accessor];
    const Json::Value* bufferView = &bufferViews[(*accessor)["bufferView"].asUInt()];
    int vertexesNum = (*accessor)["count"].asUInt();
    GLfloat* vertexes = new GLfloat[vertexesNum * 3];
    unsigned offset = (*bufferView)["byteOffset"].asUInt() + (*accessor)["byteOffset"].asUInt();//in bytes
    unsigned length = vertexesNum * 3; // in elements
    fseek(f, offset, SEEK_SET);
    fread(vertexes, sizeof(vertexes[0]), length, f);

    std::vector<glm::vec3> result;
    for (int i = 0; i < vertexesNum; i++)
    {
        result.emplace_back(vertexes[3 * i], vertexes[3 * i + 1], vertexes[3 * i + 2]);
    }
    
    return result;
}

void BasicGLTFModel::Render(RENDER_MODE mode)
{
    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);


    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    glm::mat4 viewProj = cam->GetViewProjectionMatrix();

    Entity* sun = Entity::AllEntitiesWithComponentsAndTags(
            2,0,{SunLightComponent::GetCompId(),TransformComponent::GetCompId()},{})[0];

    TransformComponent* sunTransform = sun->GetComponent<TransformComponent>();
    SunLightComponent* sunComp = sun->GetComponent<SunLightComponent>();

    glm::mat4 sunProjM = sunComp->projectionMatrix;
    glm::mat4 sunVeiwM = sunTransform->GetViewMatrixFromThisPoint();
    glm::mat4 sunViewProj = sunProjM * sunVeiwM;

    UseShader();

    Scene3DObjectDataHolder::VerifyUnoformPos(viewProjMatrix_Uniform, "view_proj_matrix");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunPos_Uniform, "sun_pos");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunColor_Uniform, "sun_color");
    Scene3DObjectDataHolder::VerifyUnoformPos(sunViewProjMatrix_Uniform, "sun_view_proj_matrix");
    Scene3DObjectDataHolder::VerifyUnoformPos(isDepthOnly_Uniform, "is_depth_only");
    Scene3DObjectDataHolder::VerifyUnoformPos(cameraPosition_Uniform, "camera_position");
    Scene3DObjectDataHolder::VerifyUnoformPos(shadowTexture_Unifrom, "shadow_texture");

    LoadInShader_Scene3DObject();
    glUniformMatrix4fv(viewProjMatrix_Uniform, 1, GL_FALSE, glm::value_ptr(viewProj));
    glUniformMatrix4fv(sunViewProjMatrix_Uniform, 1, GL_FALSE, glm::value_ptr(sunViewProj));
    glUniform3fv(sunPos_Uniform, 1, glm::value_ptr(sunTransform->position));
    glUniform3fv(sunColor_Uniform, 1, glm::value_ptr(sunComp->color));
    glUniform3fv(isDepthOnly_Uniform, 1, glm::value_ptr(cam->GetPosition()));
    glUniform1i(cameraPosition_Uniform, false);
    
    //TEXTURES
    int TEX_COUNT = 0;
    LoadInShader_BasicTexturedObject(TEX_COUNT);
    glActiveTexture(GL_TEXTURE0 + TEX_COUNT);
    glBindTexture(GL_TEXTURE_2D, sunComp->shadowTexture);
    glUniform1i(shadowTexture_Unifrom, TEX_COUNT);
    
    //DRAW CALL
    RenderableObjectDataHolder::Render(mode);
}