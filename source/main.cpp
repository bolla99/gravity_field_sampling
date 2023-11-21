  // ******************************************************** //
 // ********************* PRE PROCESSOR ******************** //
// ******************************************************** //

/* AUTHOR: GIOVANNI BOLLATI
 * UNIMI THESIS PROJECT
 */

/* definition needed to activate gl functions
 * declared in sdl_opengl_glext.h,
 * which is included by sdl_opengl.h */
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

// user defined types
#include <mesh.hpp>
#include <util.hpp>
#include <shader.hpp>
#include <GPUComputing.hpp>
#include <gravity.hpp>
#include <timer.hpp>
#include <octree.hpp>

// stdlib
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>

// sdl
#include <SDL.h>

#ifdef WIN32
#include <GL/glew.h>
#endif

// imgui
#include <imgui.h>
#include <imfilebrowser.h>
#include <imgui_impl_sdl2.h>
#include <imgui_impl_opengl3.h>

// opengl
#include <SDL_opengl.h>

// glm
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL

// image loading
#include <stb_image.h>

  // ******************************************************** //
 // ******************** GLOBAL FUNCTIONS ****************** //
// ******************************************************** //

// LOAD FILE AS STRING FOR METAL SHADERS
std::string loadMetalShader(const std::string& path);


  // ******************************************************** //
 // ************************** MAIN ************************ //
// ******************************************************** //
int main(int argv, char** args) {

    // SDL INIT
    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) < 0) {
        SDL_Log("INIT FAILED: %s", SDL_GetError());
        return -1;
    }

    // SET VERSION ATTRIBUTES

    // Always required on Mac
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);

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

    // IF ON WINDOWS INIT GLEW
#ifdef WIN32
    glewInit();
#endif

    // LOG GL VERSION
    SDL_Log("GL VERSION: %s", glGetString(GL_VERSION));

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

    // ENABLE DEPTH TEST AND SET DEPTH FUNCTION
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // GL_LESS: minor z passes

    // SHADER LOADING
    shader shdr = shader("../shaders/vs.glsl", "../shaders/fs.glsl");
    shader tetra_shdr = shader("../shaders/vs_tetra.glsl", "../shaders/fs_tetra.glsl");


    // CREATE MESH FILE BROWSER
    ImGui::FileBrowser mesh_browser;
    mesh_browser.SetTitle("mesh .obj browser");
    mesh_browser.SetTypeFilters({".obj"});
    // CREATE TEXTURE FILE BROWSER
    ImGui::FileBrowser texture_browser;
    texture_browser.SetTitle("texture browser");
    texture_browser.SetTypeFilters({".png", ".jpg", ".jpeg"});

    // CREATE MESH OBJECTS
    mesh msh = mesh();
    mesh arrow = mesh();

    // CREATE TEXTURE AND BIND
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    // SET WRAPPING MODE
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // MIPMAP (CURRENTLY NOT WORKING)
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    // SET FILTER MODE
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // SET DEFAULT TEXTURE DATA
    auto* data = new unsigned char[3 * sizeof(char)];
    data[0] = 127; data[1] = 127; data[2] = 127; // GRAY COLOR
    // set data for binded texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    delete[](data);

    // CAMERA PARAMETERS
    glm::vec3 cam_position = glm::vec3{0.f, 0.f, 10.f};
    float h_cam_angle = 0.f;
    float v_cam_angle = 0.f;

    // COLOR PARAMETERS
    float clear_color[3] = {0.5, 0.5, 0.5};
    float default_texture_color[3] = {0.5, 0.5, 0.5};

    // LIGHT UNIFORM PARAMETERS
    float ambient_color[3] = {0.5, 0.5, 0.5};
    float ambient_intensity = 0.5;
    float diffuse_color[3] = {0.5, 0.5, 0.5};
    float diffuse_position[3] = {10.0, 0.0, 0.0};
    float specular_strength = 0.5;
    int shininess = 32;

    // MODEL MATRIX PARAMETERS
    float scale = 1.0f;
    float position[3] = {0.0f, 0.0f, 0.0f};
    float rotation[3] = {0.0f, 0.0f, 0.0f};

    // GRAVITY VALUES

    // tetrahedron vertex for gravity and volume calculation; chosen arbitrarily
    float tetrahedron_vertex[3] = {0.f, 0.f, 0.f};
    // scale value for gravity calculated with tetrahedrons method
    float tetrahedrons_gravity_scale = 60.f;

    // point where gravity is calculated
    float potential_point[3] = {0.f, 0.f, 0.f};

    // output holder for gravity calculated by alternative method (tubes, masses, RT, GPU)
    glm::vec3 gravity = {0.f, 0.f, 0.f};

    // output holders for gravity vector calculated by GPUComputing and space vector
    std::vector<glm::vec3> discrete_space{};
    std::vector<glm::vec3> gpu_output_gravity{};

    // output holders for gravity with tetrahedrons metho
    glm::vec3 gravity_with_tetrahedrons = {0.f, 0.f, 0.f};
    glm::vec3 gravity_with_tetrahedrons_corrected = {0.f, 0.f, 0.f};

    // parameters used where function parameter "resolution" is required.
    int gravity_resolution = 0;

    // tubes and masses containers
    std::vector<gravity::tube> tubes = {};
    std::vector<gravity::mass> masses = {};

    // mesh volume
    float volume = 0.f;

    // DEBUG RAY; used for ray - mesh intersection
    float origin[3] = {0, 0, 0};
    float direction[3] = {0, 0, 0};

    // SET UP RAY GL
    unsigned int rayVBO, rayVAO;
    glGenBuffers(1, &rayVBO);
    glGenVertexArrays(1, &rayVAO);
    glBindVertexArray(rayVAO);
    tetra_shdr.use();
    glBindBuffer(GL_ARRAY_BUFFER, rayVBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);

      // ******************************************************** //
     // ************************ MAIN LOOP ********************* //
    // ******************************************************** //
    bool quit = false;
    while(!quit) {

          // ******************************************************** //
         // *********************** UPDATE ************************* //
        // ******************************************************** //

        // PROJECTION MATRIX
        glm::mat4 projection_matrix = glm::perspective(
                glm::radians(45.0f),
                io.DisplaySize.x / io.DisplaySize.y,
                1.f,
                100.0f
                );

        // VIEW MATRIX
        glm::mat4 cam_translation = glm::translate(glm::mat4(1.f), cam_position);
        glm::quat quatY = {cos(h_cam_angle/2.f), 0.f, sin(h_cam_angle/2.f), 0.f};
        glm::quat quatX = glm::quat(cos(v_cam_angle/2.f), sin(v_cam_angle/2.f)*(glm::toMat4(quatY)*glm::vec4(-1.f, 0.f, 0.f, 0.f)));
        glm::quat quatXY = quatX * quatY;

        glm::mat4 view_matrix = glm::inverse(glm::toMat4(quatXY) * cam_translation);

        glm::vec4 cam_position4 = glm::toMat4(quatXY) * cam_translation * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
        glm::vec3 camera_position3 = glm::vec3{cam_position4.x, cam_position4.y, cam_position4.z};

        // UPDATE SHADER CAMERA POSITION VALUE FOR SPECULAR LIGHT
        shdr.use();
        glUniform3fv(glGetUniformLocation(shdr.programID, "viewerPosition"), 1, glm::value_ptr(camera_position3));

        // MODEL MATRIX
        glm::mat4 model_translation = glm::translate(glm::mat4(1.f), {position[0], position[1], position[2]});
        glm::mat4 model_rotation_x = glm::rotate(glm::mat4(1.f), rotation[0], {1.f, 0.f, 0.f});
        glm::mat4 model_rotation_y = glm::rotate(glm::mat4(1.f), rotation[1], {0.f, 1.f, 0.f});
        glm::mat4 model_rotation_z = glm::rotate(glm::mat4(1.f), rotation[2], {0.f, 0.f, 1.f});
        glm::mat4 model_scale = glm::scale(glm::mat4(1.f), glm::vec3{scale});
        glm::mat4 model_matrix = model_scale * model_translation * model_rotation_z * model_rotation_y * model_rotation_x;

        glm::mat4 arrow_model_matrix = glm::translate(glm::mat4(1.f), glm::make_vec3(potential_point))
                * glm::toMat4(util::rotation_between_vectors({0.f, 0.f, -1.f}, glm::normalize(gravity)))
                * glm::scale(glm::mat4(1.f), {glm::length(gravity) / 100, glm::length(gravity) / 100, glm::length(gravity) / 100});

        // UPDATE SHADERS

        // matrices parameters
        shdr.use();
        glUniformMatrix4fv(glGetUniformLocation(shdr.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projection_matrix));
        glUniformMatrix4fv(glGetUniformLocation(shdr.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(view_matrix));
        glUniformMatrix4fv(glGetUniformLocation(shdr.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(model_matrix));
        tetra_shdr.use();
        glUniformMatrix4fv(glGetUniformLocation(tetra_shdr.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projection_matrix));
        glUniformMatrix4fv(glGetUniformLocation(tetra_shdr.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(view_matrix));

        // light parameters
        shdr.use();
        glUniform1f(glGetUniformLocation(shdr.programID, "ambientLightIntensity"), ambient_intensity);
        glUniform3fv(glGetUniformLocation(shdr.programID, "ambientLightColor"), 1, ambient_color);
        glUniform3fv(glGetUniformLocation(shdr.programID, "diffuseLightColor"), 1, diffuse_color);
        glUniform3fv(glGetUniformLocation(shdr.programID, "diffuseLightPosition"), 1, diffuse_position);
        glUniform1f(glGetUniformLocation(shdr.programID, "specularStrenght"), specular_strength);
        glUniform1i(glGetUniformLocation(shdr.programID, "shininess"), shininess);


          // ******************************************************* //
         // ************** UPDATE.(GRAVITY PROCESSING) ************ //
        // ******************************************************* //

        /*
         * Gravity processing is done real time with tetrahedrons method;
         * volume processing is done real time
         * other gravity processing methods are executed by gui interactions
         */

        // UPDATE MESH VOLUME
        volume = gravity::volume(msh.get_vertices(), msh.get_faces(), {tetrahedron_vertex[0], tetrahedron_vertex[1], tetrahedron_vertex[2]});

        // UPDATE GRAVITY WITH TETRAHEDRONS (REAL TIME - done every frame)
        gravity_with_tetrahedrons = gravity::get_gravity_from_tetrahedrons(msh.get_vertices(), msh.get_faces(), glm::make_vec3(potential_point), glm::make_vec3(tetrahedron_vertex)) * tetrahedrons_gravity_scale;
        gravity_with_tetrahedrons_corrected = gravity::get_gravity_from_tetrahedrons_corrected(msh.get_vertices(), msh.get_faces(), glm::make_vec3(potential_point), glm::make_vec3(tetrahedron_vertex)) * tetrahedrons_gravity_scale;



          // ******************************************************* //
         // ************************** INPUT ********************** //
        // ******************************************************* //
        SDL_Event event;

        while (SDL_PollEvent(&event)) {
            // IF USING GUI SKIP INPUT PROCESSING
            if(io.WantCaptureMouse || io.WantCaptureKeyboard) {
                ImGui_ImplSDL2_ProcessEvent(&event);
                continue;
            }
            // QUIT EVENT
            if((event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) || event.type == SDL_QUIT) quit = true;

            // MOUSE CAMERA MOVEMENT
            if(event.type == SDL_MOUSEMOTION) {
                // if motion + left button pressed
                if(event.motion.state & SDL_BUTTON_LMASK) {
                    h_cam_angle -= io.DeltaTime * (float)event.motion.xrel;
                    v_cam_angle += io.DeltaTime * (float)event.motion.yrel;
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
        if(numkey[SDL_SCANCODE_L]) v_cam_angle -= io.DeltaTime * 3.f;
        if(numkey[SDL_SCANCODE_O]) v_cam_angle += io.DeltaTime * 3.f;

        if(numkey[SDL_SCANCODE_W]) cam_position.z -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_S]) cam_position.z += io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_A]) cam_position.x -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_D]) cam_position.x += io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_Q]) cam_position.y -= io.DeltaTime * 10.f;
        if(numkey[SDL_SCANCODE_E]) cam_position.y += io.DeltaTime * 10.f;

         // ******************************************************** //
         // ************************** IMGUI ********************** //
        // ******************************************************* //

        // NEW FRAME
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame();
        ImGui::NewFrame();

        // MAIN IMGUI TAB (MESH LOADING AND GRAVITY PROCESSING
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::Begin("IO");

        // WIREFRAME - FILL CHOICE
        if(ImGui::Button("WIREFRAME")) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        if(ImGui::Button("FILL")) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // MESH LOADING OPERATIONS
        ImGui::Spacing();ImGui::Spacing();ImGui::Spacing();
        if(ImGui::Button("OPEN MESH")) mesh_browser.Open();
        if(!msh.get_vertices().empty()) {
            ImGui::Text("Mesh loaded");
            if(ImGui::Button("Save Mesh")) {
                msh.write_on_disk_as_obj("saved_mesh.obj");
            }
        }

        // GRAVITY CALCULATION OPERATIONS
        if(msh.is_loaded()) {
            ImGui::Begin("GRAVITY PROCESSING");
            ImGui::InputFloat3("tetrahedron vertex", tetrahedron_vertex);
            ImGui::Text("Volume %f", volume);

            ImGui::InputInt("ray gravity resolution", &gravity_resolution);
            ImGui::InputFloat3("potential point", potential_point);
            if(ImGui::Button("calculate gravity rt")) {
                gravity = gravity::get_gravity_RT(msh.get_vertices(), msh.get_faces(), gravity_resolution, glm::make_vec3(potential_point));
            }
            if(ImGui::Button("set up tubes")) { tubes = gravity::get_tubes(msh.get_vertices(), msh.get_faces(), gravity_resolution); }
            if(ImGui::Button("calculate gravity with tubes")) {
                gravity = gravity::get_gravity_from_tubes(msh.get_vertices(), gravity_resolution, tubes, glm::make_vec3(potential_point));
            }
            if(ImGui::Button("set up masses")) { masses = gravity::get_masses(msh.get_vertices(), msh.get_faces(), gravity_resolution); }
            if(ImGui::Button("calculate gravity with masses")) {
                Timer t{};
                t.log();
                gravity = gravity::get_gravity_from_masses(masses, 10, glm::make_vec3(potential_point));
                t.log();
            }
            ImGui::InputFloat("tetrahedrons gravity scale", &tetrahedrons_gravity_scale);

            // GPU COMPUTING BUTTON -> calculate gravity with gpu computing
            if(ImGui::Button("GPU COMPUTING")) {
                Timer t{};
                t.log();
                auto masses_as_float = (float *) &masses.front();
                std::vector<glm::vec3> space = gravity::get_discrete_space(util::get_min(msh.get_vertices()), util::get_max(msh.get_vertices()), gravity_resolution);
                for (int i = 0; i < 100; i++) {
                    std::cout << "{" << space[i].x << " " << space[i].y << " " << space[i].z << "}" << std::endl;
                }
                auto space_as_float = (float *) (glm::value_ptr(space.front()));
                float *output_gravity = GPUComputing::get_gravity_from_point_masses_and_discrete_space(
                        masses_as_float,
                        masses.size() * sizeof(gravity::mass),
                        space_as_float,
                        space.size() * sizeof(glm::vec3)
                        );
                std::vector<glm::vec3> output_gravity_as_glm_vec{};
                int k = 0;
                for(int i = 0; i < space.size(); i++) {
                    output_gravity_as_glm_vec.emplace_back(output_gravity[k], output_gravity[k+1], output_gravity[k+2]);
                    k += 3;
                }
                discrete_space = space;
                gpu_output_gravity = output_gravity_as_glm_vec;
                t.log();
            }
            if(ImGui::Button("test discrete space as octree")) {
                Timer timer{};
                timer.log(); /*
                auto o = gravity::getDiscreteSpaceAsOctree(util::getMin(mesh.getVertices()), util::getMax(mesh.getVertices()), gravity_resolution);
                std::cout << "is leaf?: " << (o->isLeaf() ? "true" : "false") << std::endl;
                std::cout << "size of octree: " << o->bytesize() << std::endl;
                delete o;
                /*
                // test trilinear interpolation
                std::array<glm::vec3, 8> cube = {glm::vec3{-1.0, -1.0, -1.0},
                                                 glm::vec3{1.0, -1.0, -1.0},
                                                 glm::vec3{-1.0, 1.0, -1.0},
                                                 glm::vec3{1.0, 1.0, -1.0},
                                                 glm::vec3{-1.0, -1.0, 1.0},
                                                 glm::vec3{1.0, -1.0, 1.0},
                                                 glm::vec3{-1.0, 1.0, 1.0},
                                                 glm::vec3{1.0, 1.0, 1.0}};
                std::cout << "il punto è nel cubo? " << (util::isInsideCube(glm::make_vec3(potential_point), cube) ? "true\n" : "false\n");
                std::array<float, 8> trilinear_coords = util::trilinearCoordinates(glm::make_vec3(potential_point), cube);
                for(int i = 0; i < 8; i++) {
                    std::cout << "trilinear coords " << i << ": " << trilinear_coords[i] << std::endl;
                } */
                auto oct = gravity::get_gravity_octree_from_masses(util::get_min(msh.get_vertices()), util::get_max(msh.get_vertices()), gravity_resolution, masses);
                timer.log();
            }

            if(ImGui::Button("compute gravity from gpu vector")) {
                Timer t{};
                t.log();
                gravity = gravity::get_gravity_from_1D_precomputed_vector(
                    glm::make_vec3(potential_point),
                    gpu_output_gravity,
                    discrete_space,
                    util::get_min(msh.get_vertices()),
                    util::get_box(util::get_min(msh.get_vertices()), util::get_max(msh.get_vertices()))[3],
                    gravity_resolution
                    );
                t.log();
            }

            // GRAVITY OUTPUTS
            ImGui::Text("gravity: %f %f %f", gravity.x, gravity.y, gravity.z);
            ImGui::Text("gravity force: %f", glm::length(gravity));

            ImGui::Text("gravity with tetrahedrons: %f %f %f", gravity_with_tetrahedrons.x, gravity_with_tetrahedrons.y, gravity_with_tetrahedrons.z);
            ImGui::Text("gravity force with tetrahedrons: %f", glm::length(gravity_with_tetrahedrons));
            ImGui::Text("gravity with tetrahedrons corrected: %f %f %f", gravity_with_tetrahedrons_corrected.x, gravity_with_tetrahedrons_corrected.y, gravity_with_tetrahedrons_corrected.z);
            ImGui::Text("gravity force with tetrahedrons corrected: %f", glm::length(gravity_with_tetrahedrons_corrected));

            // RAY PARAMETERS INPUT
            ImGui::InputFloat3("ray origin", origin);
            ImGui::InputFloat3("ray direction", direction);
            std::vector<glm::vec3> intersections = {};

            // RAY MESH INTERSECTION OUTPUT
            ImGui::Text("number of intersection: %d", (int)intersections.size());
            for(auto & intersection : intersections) {
                ImGui::Text("intersection: ( %f, %f, %f )", intersection.x, intersection.y, intersection.z);
            }
            ImGui::End();
        }
        ImGui::End();

        // IMGUI TAB: TEXTURE OPERATIONS
        ImGui::Begin("Texture");
        if(ImGui::Button("OPEN TEXTURE")) texture_browser.Open();
        if(ImGui::Button("SET COLORMODE VERTEX COLOR")) {
            shdr.use();
            glUniform1i(glGetUniformLocation(shdr.programID, "colorMode"), 1);
        }
        if(ImGui::Button("SET COLORMODE TEXTURE")) {
            shdr.use();
            glUniform1i(glGetUniformLocation(shdr.programID, "colorMode"), 0);
        }
        if(ImGui::Button("RESET TO DEFAULT TEXTURE")) {
            unsigned char* data = new unsigned char[3];
            data[0] = (char)(default_texture_color[0]*255.0f);
            data[1] = (char)(default_texture_color[1]*255.0f);
            data[2] = (char)(default_texture_color[2]*255.0f);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            delete[](data);
        }

        // COLOR SETTINGS
        ImGui::ColorEdit3("CLEAR COLOR", clear_color);
        ImGui::ColorEdit3("DEFAULT TEXTURE COLOR", default_texture_color);
        ImGui::End();

        // IMGUI TAB FOR GENERIC INFO
        ImGui::Begin("GENERIC INFO");
        ImGui::Text("Framerate: %f", io.Framerate); // FRAMERATE
        ImGui::Text("DisplaySize for imgui: %f/%f", io.DisplaySize.x, io.DisplaySize.y); // DISPLAY SIZE
        ImGui::End();

        // IMGUI TAB FOR MODEL TRANSFORM
        ImGui::Begin("Model Transform");
        ImGui::SliderFloat("model scale", &scale, 0.1f, 5.0f);
        ImGui::SliderFloat3("model position", position, -20.0, 20.0);
        ImGui::SliderFloat3("model rotation", rotation, 0.0f, 2.0f * M_PI);
        ImGui::End();

        // IMGUI TAB FOR LIGHT SETTINGS
        ImGui::SetWindowPos(ImVec2(io.DisplaySize.x - io.DisplaySize.x / 5.0, io.DisplaySize.y - io.DisplaySize.y / 5.0));
        ImGui::Begin("Light");
        ImGui::ColorEdit3("Ambient Color", ambient_color);
        ImGui::ColorEdit3("Diffuse Color", diffuse_color);
        ImGui::SliderFloat3("Diffuse position", diffuse_position, -20.0, 20.0);
        ImGui::SliderFloat("ambient intensity", &ambient_intensity, 0.0f, 5.0f);
        ImGui::SliderFloat("specular strenght", &specular_strength, 0.0f, 1.0f);
        ImGui::SliderInt("shininess", &shininess, 0, 512);
        ImGui::End();

        // DISPLAY FILE BROWSERS
        mesh_browser.Display();
        texture_browser.Display();

        // IF A MESH IS SELECTED WITH THE FILE BROWSER
        if(mesh_browser.HasSelected())
        {
            std::cout << "Selected filename" << mesh_browser.GetSelected().string() << std::endl;

            // LOAD MESH
            msh.load_from_obj(mesh_browser.GetSelected().string());

            // UPDATE COLORMODE SHADER PARAMETER
            if(msh.has_color()) {
                shdr.use();
                glUniform1i(glGetUniformLocation(shdr.programID, "colorMode"), 1);
            } else {
                shdr.use();
                glUniform1i(glGetUniformLocation(shdr.programID, "colorMode"), 0);
            }
            mesh_browser.ClearSelected();

            // LOAD ARROW MESH FOR GRAVITY DEBUG
            arrow.load_from_obj("../resources/arrow.obj");
        }

        // IF A TEXTURE IS SELECTED
        if(texture_browser.HasSelected()) {
            int width, height, nrChannels;
            unsigned char *data = stbi_load(
                    reinterpret_cast<const char *>(texture_browser.GetSelected().c_str()),
                    &width,
                    &height,
                    &nrChannels,
                    3
                    );
            if(data) {
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            } else { std::cout << "faild to load texture"; }
            texture_browser.ClearSelected();
            stbi_image_free(data);
        }


          // ******************************************************** //
         // ************************ RENDERING ********************* //
        // ******************************************************** //
        ImGui::Render();
        int w, h;
        SDL_GL_GetDrawableSize(window, &w, &h);
        glViewport(0, 0, w, h);

        glClearColor(clear_color[0], clear_color[1], clear_color[2], 1);

        glStencilMask(0xFF);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        if(!msh.get_vertices().empty()) {
            // pick shader
            shdr.use();
            // bind VERTEX ARRAY OBJECT
            glBindVertexArray(msh.get_VAO());
            // DRAW
            glDrawElements(GL_TRIANGLES, (int)msh.get_elements().size(), GL_UNSIGNED_INT, nullptr);

            glBindVertexArray(arrow.get_VAO());
            // update model matrix
            glUniformMatrix4fv(glGetUniformLocation(shdr.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(arrow_model_matrix));
            glDrawElements(GL_TRIANGLES, (int)arrow.get_elements().size(), GL_UNSIGNED_INT, nullptr);
        }

        // RAY
        tetra_shdr.use();

        // make data (two points -> line)
        float ray_data[6] = {origin[0], origin[1], origin[2],
                             origin[0] + 100 * direction[0],
                             origin[1] + 100 * direction[1],
                             origin[2] + 100 * direction[2]};

        // update VBO
        glBindBuffer(GL_ARRAY_BUFFER, rayVBO);
        glBufferData(GL_ARRAY_BUFFER, 6 * sizeof(float), ray_data, GL_DYNAMIC_DRAW);

        // draw
        glBindVertexArray(rayVAO);
        glDrawArrays(GL_LINES, 0, 2);

        // AXIS DRAW
        float axis[18] = {
                0.f, 0.f, 0.f, 10.f, 0.f, 0.f,
                0.f, 0.f, 0.f, 0.f, 10.f, 0.f,
                0.f, 0.f, 0.f, 0.f, 0.f, 10.f
        };
        glBufferData(GL_ARRAY_BUFFER, 18 * sizeof(float), axis, GL_DYNAMIC_DRAW);
        glUniform4f(glGetUniformLocation(tetra_shdr.programID, "color"), 0.8, 0, 0, 1);
        glDrawArrays(GL_LINES, 0, 2);
        glUniform4f(glGetUniformLocation(tetra_shdr.programID, "color"), 0, 0.8, 0, 1);
        glDrawArrays(GL_LINES, 2, 2);
        glUniform4f(glGetUniformLocation(tetra_shdr.programID, "color"), 0, 0, 0.8, 1);
        glDrawArrays(GL_LINES, 4, 2);

        // SWAP WINDOWS
        SDL_GL_SwapWindow(window);
    }

    SDL_Log("quitting gracefully");
}


std::string loadMetalShader(const std::string& path) {
    std::ifstream file;
    file.open(path, std::ios::in);
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}
