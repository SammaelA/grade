#include "Skybox.h"
#include "cities_generator/ECSClass.h"
#include "cities_generator/Renderer.h"
#include "third_party/stb_image.h"

Skybox::Skybox(std::string path, glm::vec3 _sunPos, float _rotation)
{
    rotationRadians = _rotation;
    sunPos = _sunPos;
    checkForGlErrors("PRE-CUBEMAP CREATION");
    m_shader = ShaderLibrary::GetLibrary()->GetShaderPointer("skyboxShader");
    GLfloat vertexes[] = {
        -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        -1.0f,  1.0f, -1.0f,
        1.0f,  1.0f, -1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
        1.0f, -1.0f,  1.0f
    };

    sunPos = glm::rotateY(sunPos, rotationRadians);
    rotaionMatrix = glm::eulerAngleZY(0.f, -rotationRadians);

    vertexesNum = 36;
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertexes), vertexes, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    std::vector<std::string> faces =
    {
        "right.png",
        "left.png",
        "top.png",
        "bottom.png",
        "front.png",
        "back.png"
    };
    glGenTextures(1, &cubeMap);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubeMap);
    int width, height;
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        unsigned char *data = stbi_load((path + faces[i]).c_str(), &width, &height, 0, 4);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 
                         0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data
            );
            stbi_image_free(data);
        }
        else
        {
            std::cout << "Cubemap texture failed to load at path: " << faces[i] << std::endl;
            stbi_image_free(data);
        }
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    checkForGlErrors("CUBEMAP CREATION");
    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
}

void Skybox::Render(RENDER_MODE mode)
{
    UseShader();
    glEnable(GL_DEPTH_TEST);
    Renderer::GetInstance()->SetFullViewport();
    glDepthFunc(GL_LEQUAL);

    //Camera loading
    CameraComponent* cam = Entity::AnyObjectWithTag("mainCamera")->GetComponent<CameraComponent>();
    glm::mat4 viewProj = cam->GetViewMatrix();
    viewProj = glm::mat4(glm::mat3(viewProj));
    viewProj = cam->GetProjectionMatrix() * viewProj;
    glUniformMatrix4fv(
        glGetUniformLocation(m_shader->Program, "view_proj_matrix"),
        1, GL_FALSE, glm::value_ptr(viewProj)
    );
    glUniformMatrix4fv(
        glGetUniformLocation(m_shader->Program, "skybox_rotation_matrix"),
        1, GL_FALSE, glm::value_ptr(rotaionMatrix)
    );
    
    //TEXTURES
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(m_shader->Program, "skybox"), 0);
    glBindTexture(GL_TEXTURE_CUBE_MAP, cubeMap);

    //DRAW CALL
    RenderableObjectDataHolder::Render(mode);

    //RETURN OPENGL STATE
    glDepthFunc(GL_LESS);
}

glm::vec3 Skybox::GetSunPosition()
{
    return sunPos;
}