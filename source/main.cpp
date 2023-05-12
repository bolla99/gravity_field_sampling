#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

//#define SDL_MAIN_HANDLED

#include <iostream>
#include <sstream>
#include <cmath>
#include <SDL.h>

#ifdef WIN32
#include <GL/glew.h>
#endif

#include <imgui.h>
#include <imfilebrowser.h>
#include <imgui_impl_sdl2.h>
#include <imgui_impl_opengl3.h>
#include <SDL_opengl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <Shader.hpp>
#include <omp.h>
#include <Mesh.hpp>
#include <stb_image.h>
#include <random>
#include <string>

#ifdef METAL_READY
// METAL
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <simd/simd.h>
#endif

// LOAD FILE AS STRING FOR METAL SHADERS
std::string loadMetalShader(const std::string& path);

int main(int argv, char** args) {

    // SDL INIT
    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) < 0) {
        SDL_Log("INIT FAILED: %s", SDL_GetError());
        return -1;
    }

    // SET VERSION ATTRIBUTES
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG); // Always required on Mac
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);

    // SET GL ATTRIBUTES
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);

    // CREATE WINDOW AND CONTEXT
    auto window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);

    SDL_Window* window = SDL_CreateWindow(
            "",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED, 800, 600, window_flags
            );
    // CREATE GL CONTEXT
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);

#ifdef WIN32
    glewInit();
#endif

    std::cout << "gl version: " << glGetString(GL_VERSION);

    // SET GL CONTEXT
    SDL_GL_MakeCurrent(window, gl_context);

    // ENABLE VSYNC
    SDL_GL_SetSwapInterval(1);

    // IMGUI SETUP
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();

    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init("#version 410");

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // SHADER LOADING AND SETUP
    Shader shader = Shader("../shaders/vs.glsl", "../shaders/fs.glsl");
    Shader tetra_shader = Shader("../shaders/vs_tetra.glsl", "../shaders/fs_tetra.glsl");


    // IMGUI FIlE BROWSER
    ImGui::FileBrowser meshBrowser;
    meshBrowser.SetTitle("mesh .obj browser");
    meshBrowser.SetTypeFilters({".obj"});

    ImGui::FileBrowser textureBrowser;
    textureBrowser.SetTitle("texture browser");
    textureBrowser.SetTypeFilters({".png", ".jpg", ".jpeg"});

    // MESH
    Mesh mesh = Mesh();

    // TEXTURE -> discorso texture da approfondire
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // SET DEFAULT TEXTURE
    auto* data = new unsigned char[3 * sizeof(char)];
    data[0] = 127; data[1] = 127; data[2] = 127;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    delete[](data);

    float f[3] = {0, 0, 0};
    std::string isInside;

    glm::vec3 cam_position = glm::vec3{0.f, 0.f, 10.f};
    float h_cam_angle = 0.f;
    float v_cam_angle = 0.f;

    float clear_color[3] = {0.5, 0.5, 0.5};
    float default_texture_color[3] = {0.5, 0.5, 0.5};
    float ambient_color[3] = {0.5, 0.5, 0.5};
    float diffuse_color[3] = {0.5, 0.5, 0.5};
    float diffuse_position[3] = {10.0, 0.0, 0.0};
    float ambient_intensity = 0.5;
    float specularStrenght = 0.5;
    int shininess = 32;

    // model matrix
    float scale = 1.0f;
    float position[3] = {0.0f, 0.0f, 0.0f};
    float rotation[3] = {0.0f, 0.0f, 0.0f};

    float tetrahedron_vertex[3] = {0.f, 0.f, 0.f};
    float potential_point[3] = {0.f, 0.f, 0.f};
    glm::vec3 gravity = {0.f, 0.f, 0.f};
    float volume = 0.f;

    // tetrahedron data
    /*
    unsigned int tVBOs[1000];
    unsigned int tVAOs[1000];
    unsigned int tEBOs[1000];
     */

    // RAY
    float origin[3] = {0, 0, 0};
    float direction[3] = {0, 0, 0};

    int gravity_resolution = 0;
    std::vector<tube> tubes = {};
    std::vector<mass> masses = {};
    std::vector<glm::vec3> gravity_output;
    auto gravity_output_p = (float *)&gravity_output.front();
    int masses_size = 0;
    float *gravity_output_heap;

    unsigned int rayVBO, rayVAO;
    glGenBuffers(1, &rayVBO);
    glGenVertexArrays(1, &rayVAO);
    glBindVertexArray(rayVAO);
    tetra_shader.use();
    glBindBuffer(GL_ARRAY_BUFFER, rayVBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);

      // ******************************************************** //
     // ************************ MAIN LOOP ********************* //
    // ***************************************************O**** //
    bool quit = false;
    while(!quit) {


          // ******************************************************** //
         // ************************** UPDATE ********************** //
        // ******************************************************** //

        // PROJECTION MATRIX
        glm::mat4 projectionMatrix = glm::perspective(
                glm::radians(45.0f),
                io.DisplaySize.x / io.DisplaySize.y,
                1.f,
                100.0f
                );

        // VIEW MATRIX
        glm::mat4 camTranslation = glm::translate(glm::mat4(1.f), cam_position);
        glm::mat4 camRotationX = glm::rotate(glm::mat4(1.f), v_cam_angle, glm::vec3{-1.f, 0.f, 0.f});
        glm::mat4 camRotationY = glm::rotate(glm::mat4(1.f), h_cam_angle, glm::vec3{0.f, 1.f, 0.f});
        glm::mat4 viewMatrix = glm::inverse(camRotationX * camRotationY * camTranslation);

        glm::vec4 camPosition4 = camRotationX * camRotationY * camTranslation * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
        glm::vec3 cameraPosition3 = glm::vec3{camPosition4.x, camPosition4.y, camPosition4.z};
        shader.use();
        glUniform3fv(glGetUniformLocation(shader.programID, "viewerPosition"), 1, glm::value_ptr(cameraPosition3));

        // MODEL MATRIX
        glm::mat4 modelTranslation = glm::translate(glm::mat4(1.f), glm::vec3{position[0], position[1], position[2]});
        glm::mat4 modelRotationX = glm::rotate(glm::mat4(1.f), rotation[0], glm::vec3{1.f, 0.f, 0.f});
        glm::mat4 modelRotationY = glm::rotate(glm::mat4(1.f), rotation[1], glm::vec3{0.f, 1.f, 0.f});
        glm::mat4 modelRotationZ = glm::rotate(glm::mat4(1.f), rotation[2], glm::vec3{-1.f, 0.f, 0.f});
        glm::mat4 modelScale = glm::scale(glm::mat4(1.f), glm::vec3{scale});
        glm::mat4 modelMatrix = modelScale * modelTranslation * modelRotationZ * modelRotationY * modelRotationX;

        // GLM::LOOKAT METHOD
        /* GLM::LOOK AT METHOD
        glm::vec4 camTransform = glm::vec4(cam_position.x, cam_position.y, cam_position.z, 1.0f);
        glm::vec3 v_up = glm::vec3{0.f, 1.f, 0.f};
        camTransform = camRotationY * camRotationX * camTransform;
        auto viewMatrix = glm::lookAt(
                glm::vec3{camTransform.x, camTransform.y, camTransform.z},
                glm::vec3{0.f, 0.f, 0.f},
                v_up
                );
        */
        shader.use();
        glUniformMatrix4fv(glGetUniformLocation(shader.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(glGetUniformLocation(shader.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(viewMatrix));
        glUniformMatrix4fv(glGetUniformLocation(shader.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(modelMatrix));
        tetra_shader.use();
        glUniformMatrix4fv(glGetUniformLocation(tetra_shader.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(glGetUniformLocation(tetra_shader.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(viewMatrix));

          // ******************************************************* //
         // ************************** INPUT ********************** //
        // ******************************************************* //
        SDL_Event event;

        while (SDL_PollEvent(&event)) {
            if(io.WantCaptureMouse || io.WantCaptureKeyboard) {
                ImGui_ImplSDL2_ProcessEvent(&event);
                continue;
            }
            if((event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) || event.type == SDL_QUIT) quit = true;
            // MOUSE CAMERA MOVEMENT
            if(event.type == SDL_MOUSEMOTION) {
                if(event.motion.state & SDL_BUTTON_LMASK) {
                    h_cam_angle -= io.DeltaTime * (float)event.motion.xrel;
                    if(
                            (v_cam_angle > (-M_PI / 2 + 0.1) && (float)event.motion.yrel < 0.f)
                            || (v_cam_angle < (M_PI / 2 - 0.1) && (float)event.motion.yrel > 0.f)
                            )
                    {
                        v_cam_angle += io.DeltaTime * (float)event.motion.yrel;
                    }
                } else if(SDL_GetKeyboardState(nullptr)[SDL_SCANCODE_LGUI]) { //if(event.motion.state & SDL_BUTTON_RMASK) {
                    cam_position.z += io.DeltaTime * (float)event.motion.yrel * 10.f;
                }
            }
        }

        // PROCESS KEYBOARD INPUT FOR CAMERA MOVEMENT
        const Uint8* numkey = SDL_GetKeyboardState(nullptr);
        if(numkey[SDL_SCANCODE_DOWN]) cam_position.z += io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_UP]) cam_position.z -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_LEFT]) h_cam_angle += io.DeltaTime * 3.f;
        if(numkey[SDL_SCANCODE_RIGHT]) h_cam_angle -= io.DeltaTime * 3.f;
        if(numkey[SDL_SCANCODE_L] && v_cam_angle > (-M_PI / 2 + 0.1)) v_cam_angle -= io.DeltaTime * 3.f;
        if(numkey[SDL_SCANCODE_O] && v_cam_angle < (M_PI / 2 - 0.1)) v_cam_angle += io.DeltaTime * 3.f;

        if(numkey[SDL_SCANCODE_W]) cam_position.z -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_S]) cam_position.z += io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_A]) cam_position.x -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_D]) cam_position.x += io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_Q]) cam_position.y -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_E]) cam_position.y += io.DeltaTime * 10.f;

         // ******************************************************** //
         // ************************** IMGUI ********************** //
        // ******************************************************* //
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::Begin("General");
        if(ImGui::Button("WIREFRAME")) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        if(ImGui::Button("FILL")) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        ImGui::Text("Framerate: %f", io.Framerate);

        ImGui::Spacing();ImGui::Spacing();ImGui::Spacing();
        if(!mesh.getVertices().empty()) {
            ImGui::Text("Mesh loaded");
            if(ImGui::Button("Save Mesh")) {
                mesh.writeOnDiskAsObj("saved_mesh.obj");
            }
        }
        if(mesh.isLoaded()) {
            ImGui::SliderFloat3("tetrahedron vertex", tetrahedron_vertex,-1000, 1000);
            ImGui::Text("Volume %f", volume);

            ImGui::SliderInt("ray gravity resolution", &gravity_resolution, 0, 1000);
            ImGui::InputFloat3("potential point", potential_point);
            if(ImGui::Button("calculate gravity rt")) {
                gravity = mesh.getGravityRT(gravity_resolution, {potential_point[0], potential_point[1], potential_point[2]});
            }
            if(ImGui::Button("set up tubes")) { tubes = mesh.getTubes(gravity_resolution); }
            if(ImGui::Button("calculate gravity with tubes")) {
                gravity = mesh.getGravityFromTubes(gravity_resolution, tubes, {potential_point[0], potential_point[1], potential_point[2]});
            }
            if(ImGui::Button("set up masses")) { masses = mesh.getMasses(gravity_resolution); }
            if(ImGui::Button("calculate gravity with masses")) {
                gravity = Mesh::getGravityFromMasses(masses, 10, {potential_point[0], potential_point[1], potential_point[2]});
            }
#ifdef METAL_READY
            if(ImGui::Button("GPU COMPUTING")) {
                auto masses_as_float = (float *)&masses.front();
                std::vector<glm::vec3> space = mesh.getDiscreteSpace(gravity_resolution);
                for(int i = 0; i < 100; i++) {
                    std::cout << "{" << space[i].x << " " << space[i].y << " " << space[i].z << "}" << std::endl;
                }
                auto space_as_float = (float *)(glm::value_ptr(space.front()));
                gravity_output.reserve((int)pow(gravity_resolution, 3));
                gravity_output_heap = (float *)malloc(sizeof(glm::vec3) * (int)space.size());

                UInt64 start = SDL_GetTicks64();

                std::string metal_shader = loadMetalShader("../shaders/add.metal");
                const char* metal_shader_c = metal_shader.c_str();
                auto source = NS::String::string(metal_shader_c, NS::UTF8StringEncoding);

                NS::Error* error;
                MTL::Device *device = MTL::CreateSystemDefaultDevice();
                MTL::Library *library = device->newLibrary(source, nullptr, &error);

                auto str = NS::String::string("add_arrays", NS::ASCIIStringEncoding);
                MTL::Function *function = library->newFunction(str);
                if(function == nullptr) return 1;

                NS::Error* e;
                MTL::ComputePipelineState* pso = device->newComputePipelineState(function, &e);
                function->release();
                MTL::CommandQueue* commandQueue = device->newCommandQueue();

                masses_size = (int)(sizeof(mass) * masses.size());
                std::cout << "MASSES SIZE " << masses_size;

                MTL::Buffer* space_buffer = device->newBuffer(space_as_float, (int)space.size() * sizeof(glm::vec3), MTL::ResourceStorageModeShared);
                MTL::Buffer* masses_buffer = device->newBuffer(masses_as_float, masses.size() * sizeof(mass), MTL::ResourceStorageModeShared);
                MTL::Buffer* gravity_buffer = device->newBuffer(gravity_output_heap, (int)space.size() * sizeof(glm::vec3), MTL::ResourceStorageModeShared);
                MTL::Buffer* masses_size_buffer = device->newBuffer(&masses_size, sizeof(int), MTL::ResourceStorageModeShared);

                MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
                MTL::ComputeCommandEncoder* computeCommandEncoder = commandBuffer->computeCommandEncoder();

                computeCommandEncoder->setComputePipelineState(pso);
                computeCommandEncoder->setBuffer(space_buffer, 0, 0);
                computeCommandEncoder->setBuffer(masses_buffer, 0, 1);
                computeCommandEncoder->setBuffer(gravity_buffer, 0, 2);
                computeCommandEncoder->setBuffer(masses_size_buffer, 0, 3);

                auto size = (int)space.size();
                MTL::Size gridsize = MTL::Size::Make(size, 1, 1);
                NS::UInteger threadGroupSize = pso->maxTotalThreadsPerThreadgroup();
                if(threadGroupSize > size) {
                    threadGroupSize = size;
                }
                MTL::Size threadsize = MTL::Size::Make(threadGroupSize, 1, 1);
                computeCommandEncoder->dispatchThreads(gridsize, threadsize);

                computeCommandEncoder->endEncoding();
                commandBuffer->commit();
                commandBuffer->waitUntilCompleted();

                std::cout << "time elapsed: " << ((float)SDL_GetTicks64() - (float)start)/ 1000.f << std::endl;

                for(int i = 0; i < 1000; i+=3) {
                    std::cout << ((float *)gravity_buffer->contents())[i] << " "
                    << ((float *)gravity_buffer->contents())[i + 1] << " "
                    << ((float *)gravity_buffer->contents())[i + 2] << std::endl;
                }
            }
#endif

            //gravity = mesh.getGravity(gravity_resolution, {potential_point[0], potential_point[1], potential_point[2]});
            ImGui::Text("gravity: %f %f %f", gravity.x, gravity.y, gravity.z);
            ImGui::Text("gravity force: %f", glm::length(gravity));

            ImGui::InputFloat3("ray origin", origin);
            ImGui::InputFloat3("ray direction", direction);
            std::vector<glm::vec3> intersections = {};/*mesh.rayMeshIntersections({
                                              {origin[0], origin[1], origin[2]},
                                              {direction[0], direction[1], direction[2]}});*/
            ImGui::Text("number of intersection: %d", (int)intersections.size());
            for(auto & intersection : intersections) {
                ImGui::Text("intersection: ( %f, %f, %f )", intersection.x, intersection.y, intersection.z);
            }
        }
        /*gravity = mesh.gravity(
                glm::vec3{potential_point[0], potential_point[1], potential_point[2]},
                     glm::vec3{tetrahedron_vertex[0], tetrahedron_vertex[1], tetrahedron_vertex[2]});*/
        volume = mesh.volume({tetrahedron_vertex[0], tetrahedron_vertex[1], tetrahedron_vertex[2]});
        if(ImGui::Button("OPEN MESH")) meshBrowser.Open();
        if(ImGui::Button("OPEN TEXTURE")) textureBrowser.Open();
        if(ImGui::Button("SET COLORMODE VERTEX COLOR")) {
            shader.use();
            glUniform1i(glGetUniformLocation(shader.programID, "colorMode"), 1);
        }
        if(ImGui::Button("SET COLORMODE TEXTURE")) {
            shader.use();
            glUniform1i(glGetUniformLocation(shader.programID, "colorMode"), 0);
        }
        if(ImGui::Button("RESET TO DEFAULT TEXTURE")) {
            unsigned char* data = new unsigned char[3];
            data[0] = (char)(default_texture_color[0]*255.0f);
            data[1] = (char)(default_texture_color[1]*255.0f);
            data[2] = (char)(default_texture_color[2]*255.0f);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            delete[](data);
        }

        // WINDOW SIZE
        ImGui::Text("DisplaySize for imgui: %f/%f", io.DisplaySize.x, io.DisplaySize.y);

        ImGui::ColorEdit3("CLEAR COLOR", clear_color);
        ImGui::ColorEdit3("DEFAULT TEXTURE COLOR", default_texture_color);
        ImGui::End();

        ImGui::Begin("Model Transform");
        ImGui::SliderFloat("model scale", &scale, 0.1f, 5.0f);
        ImGui::SliderFloat3("model position", position, -20.0, 20.0);
        ImGui::SliderFloat3("model rotation", rotation, 0.0f, 2.0f * M_PI);
        ImGui::End();

        ImGui::SetWindowPos(ImVec2(io.DisplaySize.x - io.DisplaySize.x / 5.0, io.DisplaySize.y - io.DisplaySize.y / 5.0));
        ImGui::Begin("Light");
        ImGui::ColorEdit3("Ambient Color", ambient_color);
        ImGui::ColorEdit3("Diffuse Color", diffuse_color);
        ImGui::SliderFloat3("Diffuse position", diffuse_position, -20.0, 20.0);
        ImGui::SliderFloat("ambient intensity", &ambient_intensity, 0.0f, 5.0f);
        ImGui::SliderFloat("specular strenght", &specularStrenght, 0.0f, 1.0f);
        ImGui::SliderInt("shininess", &shininess, 0, 512);
        ImGui::End();

        shader.use();
        glUniform1f(glGetUniformLocation(shader.programID, "ambientLightIntensity"), ambient_intensity);
        glUniform3fv(glGetUniformLocation(shader.programID, "ambientLightColor"), 1, ambient_color);
        glUniform3fv(glGetUniformLocation(shader.programID, "diffuseLightColor"), 1, diffuse_color);
        glUniform3fv(glGetUniformLocation(shader.programID, "diffuseLightPosition"), 1, diffuse_position);
        glUniform1f(glGetUniformLocation(shader.programID, "specularStrenght"), specularStrenght);
        glUniform1i(glGetUniformLocation(shader.programID, "shininess"), shininess);

        meshBrowser.Display();
        textureBrowser.Display();

        if(meshBrowser.HasSelected())
        {
            /*
            glDeleteBuffers(1000, tVBOs);
            glDeleteBuffers(1000, tEBOs);
            glDeleteVertexArrays(1000, tVAOs);*/

            std::cout << "Selected filename" << meshBrowser.GetSelected().string() << std::endl;
            mesh.loadFromObj(meshBrowser.GetSelected().string());
            if(mesh.hasColor()) {
                shader.use();
                glUniform1i(glGetUniformLocation(shader.programID, "colorMode"), 1);
            } else {
                shader.use();
                glUniform1i(glGetUniformLocation(shader.programID, "colorMode"), 0);
            }
            meshBrowser.ClearSelected();
            /*
            glGenVertexArrays(mesh.getFaces().size(), tVAOs);
            glGenBuffers(mesh.getFaces().size(), tVBOs);
            glGenBuffers(mesh.getFaces().size(), tEBOs);*/
        }

        if(textureBrowser.HasSelected()) {
            int width, height, nrChannels;
            unsigned char *data = stbi_load(reinterpret_cast<const char *>(textureBrowser.GetSelected().c_str()), &width, &height, &nrChannels, 3);
            if(data) {
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
                //glGenerateMipmap(GL_TEXTURE_2D);
            } else { std::cout << "faild to load texture"; }
            textureBrowser.ClearSelected();
            stbi_image_free(data);
        }


          // ******************************************************** //
         // ************************** RENDER ********************** //
        // ******************************************************** //
        ImGui::Render();
        int w, h;
        SDL_GL_GetDrawableSize(window, &w, &h);
        glViewport(0, 0, w, h);

        glClearColor(clear_color[0], clear_color[1], clear_color[2], 1);

        glStencilMask(0xFF);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if(!mesh.getVertices().empty()) {
            shader.use();
            glBindVertexArray(mesh.getVAO());
            glDrawElements(GL_TRIANGLES, (int)mesh.getElements().size(), GL_UNSIGNED_INT, nullptr);
        }

        tetra_shader.use();
        glBindVertexArray(rayVAO);
        glBindBuffer(GL_ARRAY_BUFFER, rayVBO);
        float ray_data[6] = {origin[0], origin[1], origin[2],
                             origin[0] + 100 * direction[0],
                             origin[1] + 100 * direction[1],
                             origin[2] + 100 * direction[2]};
        glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(float), ray_data, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
        glEnableVertexAttribArray(0);
        glDrawArrays(GL_LINES, 0, 2);
        /*
        for(int i = 0; i < mesh.getFaces().size(); i++) {
            glm::vec3 t[4] = {mesh.getVertices()[mesh.getFaces()[i].x - 1],
                          mesh.getVertices()[mesh.getFaces()[i].y - 1],
                          mesh.getVertices()[mesh.getFaces()[i].z - 1],
                          {tetrahedron_vertex[0], tetrahedron_vertex[1], tetrahedron_vertex[2] }};
            tetrahedron tetrahedron = {t[0], t[1], t[2], t[3]};
            float volume = Mesh::tetrahedronVolume(tetrahedron);

            const float green[4] = {0.f, 1.f, 0.f, 1.f};
            const float red[4] = {1.f, 0.f, 0.f, 1.f};

            if(volume >= 0) glUniform4fv(glGetUniformLocation(tetra_shader.programID, "color"), 1, green);
            else glUniform4fv(glGetUniformLocation(tetra_shader.programID, "color"), 1, red);

            unsigned int f[12] = {0, 1, 2, 0, 1, 3, 1, 2, 3, 2, 0, 3};
            glBindVertexArray(tVAOs[i]);
            glBindBuffer(GL_ARRAY_BUFFER, tVBOs[i]);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tEBOs[i]);
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, 12 * sizeof(int), f, GL_DYNAMIC_DRAW);
            glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * 4, t, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
            glEnableVertexAttribArray(0);
            glDrawElements(GL_TRIANGLES, 12, GL_UNSIGNED_INT, nullptr);
        }
         */

        SDL_GL_SwapWindow(window);
    }
    std::cout << "quitting gracefully";
}


std::string loadMetalShader(const std::string& path) {
    std::ifstream file;
    file.open(path, std::ios::in);
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}
