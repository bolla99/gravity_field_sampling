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

// stdlib
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <unordered_map>
#include <queue>
#include <tuple>
#include <cstdio>
#include <algorithm>

// sdl
#include <SDL.h>

#ifdef WIN32
#include <GL/glew.h>
#endif

// bullet

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
#include <glm/gtx/hash.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL

// image loading
#include <stb_image.h>

#include <omp.h>

  // ******************************************************** //
 // ******************** GLOBAL FUNCTIONS ****************** //
// ******************************************************** //



struct DATA {
    struct GRAVITY {
        constexpr static glm::vec3 MIN{-5.f, -5.f, -5.f};
        constexpr static float EDGE = 10.f;
        constexpr static float TUBES_RESOLUTION = 128;
        constexpr static float PRECISION = 0.5;
    };
};
bool file_exists(const std::string& name);
void set_up_gl_attributes();
void set_up_gl();

  // ******************************************************** //
 // ************************** MAIN ************************ //
// ******************************************************** //
int main(int argv, char** args) {
    float ball_initial_velocity[3] = {0, 0, 0};

    // SDL INIT
    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) < 0) {
        SDL_Log("INIT FAILED: %s", SDL_GetError());
        return -1;
    }

    set_up_gl_attributes();

      // CREATE WINDOW
    auto window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_Window* window = SDL_CreateWindow(
            "GRAVITY PROCESSING",
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED, 800, 600, window_flags
            );
    // CREATE GL CONTEXT
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_MakeCurrent(window, gl_context);

    set_up_gl();
    // LOG GL VERSION
    SDL_Log("GL VERSION: %s", glGetString(GL_VERSION));

    // IF ON WINDOWS INIT GLEW
#ifdef WIN32
    glewInit();
#endif

    // IMGUI SETUP
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init("#version 410");

    // SHADER LOADING
    shader mesh_shader = shader("../shaders/vs_0.glsl", "../shaders/fs_0.glsl");
    shader debug_shader = shader("../shaders/vs_1.glsl", "../shaders/fs_1.glsl");
    shader octree_shader = shader("../shaders/octree_vs.glsl", "../shaders/octree_fs.glsl");

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

    // DEBUG BALL
    mesh debug_ball = mesh();
    glm::vec3 debug_ball_velocity{0.f};
    bool is_gravity_on = false;

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

    // DEBUG BALL MODEL MATRIX PARAMETERS
    float debug_ball_scale = 0.1f;
    float debug_ball_position[3] = {1.f, 0.f, 0.f};

    // GRAVITY VALUES

    // point where gravity is calculated
    float potential_point[3] = {1.f, 0.f, 0.f};

    // output holders for gravity
    glm::vec3 gravity{0.f, 0.f, 0.f};

    glm::vec3 gravity_from_octree_potential{0.f, 0.f, 0.f};

    glm::vec3 gravity_no_gpu{0.f, 0.f, 0.f};

    glm::vec3 gravity_from_octree{0.f, 0.f, 0.f};

    glm::vec3 gravity_from_potential_gradient{0.f, 0.f, 0.f};

    // output holders for gravity vector calculated by GPUComputing and space vector
    std::vector<glm::vec3> discrete_space{};
    std::vector<glm::vec3> gpu_output_gravity{};

    // output holder for potentials vector calculated by GPUComputing
    std::vector<float> gpu_output_potential{};

    float cylinder_R;

    // G constant
    float G = 10;

    // TUBES RESOLUTION
    int tube_resolution = 127;

    // benchmark parameters
    int n_test = 1000;

    int n_real_time = 1;
    int n_real_time_octree = 1;
    float real_time_performance = 0;
    float real_time_performance_octree = 0;

    // performance data
    int n_gravity = 1;
    int n_potential = 1;
    int space_test_resolution = 1;

    std::vector<gravity::tube> tubes{};

    // octree data
    // semantics: if octree[i] <= 0 then -octree[i] is the index of gravity value; the following
    // seven values have the same semantics; if octree[i] > 0 then octree[i] is the index (on
    // the same vector) of the first vertex of inner cube; the following values have the same semantics
    auto octree = std::vector<int>();
    auto potential_octree = std::vector<int>();

    // gravity values; this vector is accessed from indexes retrieved from octree vecto
    auto gravity_values = std::vector<glm::vec3>();

    auto octree_space_locations = std::vector<glm::vec3>();
    auto potential_values = std::vector<float>();

    auto potential_cached_values = std::unordered_map<glm::ivec3, float>();

    // if int value < 0 -> value is in tmp_gravity_value with index -i - 1
    // if int value >= 0 -> value is in gravity_value
    // when cached_values is filled with gpu_output_gravity values (precomputed grid)
    // integer coordinates depends on octree_depth and cached_gpu_depth (gpu_depth
    // when gpu grid was computed; that means that if octree_depth is changed,
    // cached_values must be refilled with same values but different coordinates;
    // workflow is: set gpu depth -> compute grid; gpu depth is cached; than fill cached_values;
    // than build octree; if gpu_depth is changed repeat whole process; is octree_depth is
    // changed repeat from fill cached values; if precision is changed repeat from build octree.
    // elements are integer coordinates -> index
    auto cached_values = std::unordered_map<glm::ivec3, int>();

    // this map contains values already used and added to gravity_values during building
    // this is cleared before octree building and it is filled during octree building, and
    // then cleared
    auto gravity_values_map = std::unordered_map<glm::ivec3, int>();

    // these are temporary gravity values: they are cached (already computed) gravity values,
    // they will be added to gravity values if they will be used, in fact gravity_values
    // size is always less than tmp_gravity_values size
    auto tmp_gravity_values = std::vector<glm::vec3>();

    // area of octree computation; changes to these values imply the repetition of octree
    // workflow (from grid computation)
    float min[3] {-5.f, -5.f, -5.f};
    float edge = 10.f;

    // octree parameters

    int should_divide_method = 2;

    // depth of grid computed by gpu
    int gpu_depth = 8;
    int cached_gpu_depth = gpu_depth;
    // max depth of octree
    int octree_depth = 20;
    int cached_octree_depth = octree_depth;
    float precision = 1.0f;
    float alpha = 1.0f;

    int min_depth = 0;

    // rendering octree data structure
    struct octree_rendering_unit {
        glm::vec3 p;
        float depth;
    };
    std::vector<octree_rendering_unit> octree_rendering{};
    std::vector<octree_rendering_unit> potential_octree_rendering{};

    int max_depth_reached = 0;
    int potential_max_depth_reached = 0;

    // DEBUG RAY; used for ray - mesh intersection
    float origin[3] = {0, 0, 0};
    float direction[3] = {0, 0, 0};

    // SET UP RAY GL
    unsigned int rayVBO, rayVAO;
    glGenBuffers(1, &rayVBO);
    glGenVertexArrays(1, &rayVAO);
    glBindVertexArray(rayVAO);
    debug_shader.use();
    glBindBuffer(GL_ARRAY_BUFFER, rayVBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);

    // setup octree representation shader
    unsigned int octreeVBO, octreeVAO;
    glGenBuffers(1, &octreeVBO);
    glGenVertexArrays(1, &octreeVAO);
    glBindVertexArray(octreeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, octreeVBO);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

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
        mesh_shader.use();
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "viewerPosition"), 1, glm::value_ptr(camera_position3));

        // MODEL MATRIX
        glm::mat4 model_translation = glm::translate(glm::mat4(1.f), {position[0], position[1], position[2]});
        glm::mat4 model_rotation_x = glm::rotate(glm::mat4(1.f), rotation[0], {1.f, 0.f, 0.f});
        glm::mat4 model_rotation_y = glm::rotate(glm::mat4(1.f), rotation[1], {0.f, 1.f, 0.f});
        glm::mat4 model_rotation_z = glm::rotate(glm::mat4(1.f), rotation[2], {0.f, 0.f, 1.f});
        glm::mat4 model_scale = glm::scale(glm::mat4(1.f), glm::vec3{scale});
        glm::mat4 model_matrix = model_scale * model_translation * model_rotation_z * model_rotation_y * model_rotation_x;

        /*
        // update ball position
        if(is_gravity_on) {
            debug_ball_position[0] += debug_ball_velocity.x * io.DeltaTime;
            debug_ball_position[1] += debug_ball_velocity.y * io.DeltaTime;
            debug_ball_position[2] += debug_ball_velocity.z * io.DeltaTime;
        }
        */

        glm::mat4 debug_ball_model_translation = glm::translate(
            glm::mat4(1.f), {debug_ball_position[0], debug_ball_position[1], debug_ball_position[2]}
            );
        glm::mat4 debug_ball_model_scale = glm::scale(glm::mat4(1.f), glm::vec3{debug_ball_scale});
        auto debug_ball_model_matrix = debug_ball_model_translation * debug_ball_model_scale;

        glm::mat4 arrow_model_matrix = glm::translate(glm::mat4(1.f), glm::make_vec3(debug_ball_position))
                * glm::toMat4(util::rotation_between_vectors({0.f, 0.f, -1.f}, glm::normalize(gravity)))
                * glm::scale(glm::mat4(1.f), {glm::length(gravity) / 100, glm::length(gravity) / 100, glm::length(gravity) / 100});

        // UPDATE SHADERS

        // matrices parameters
        mesh_shader.use();
        glUniformMatrix4fv(glGetUniformLocation(mesh_shader.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projection_matrix));
        glUniformMatrix4fv(glGetUniformLocation(mesh_shader.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(view_matrix));
        glUniformMatrix4fv(glGetUniformLocation(mesh_shader.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(model_matrix));
        debug_shader.use();
        glUniformMatrix4fv(glGetUniformLocation(debug_shader.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projection_matrix));
        glUniformMatrix4fv(glGetUniformLocation(debug_shader.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(view_matrix));
        octree_shader.use();
        glUniformMatrix4fv(glGetUniformLocation(octree_shader.programID, "projectionMatrix"), 1, GL_FALSE, glm::value_ptr(projection_matrix));
        glUniformMatrix4fv(glGetUniformLocation(octree_shader.programID, "viewMatrix"), 1, GL_FALSE, glm::value_ptr(view_matrix));


        // light parameters
        mesh_shader.use();
        glUniform1f(glGetUniformLocation(mesh_shader.programID, "ambientLightIntensity"), ambient_intensity);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "ambientLightColor"), 1, ambient_color);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "diffuseLightColor"), 1, diffuse_color);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "diffuseLightPosition"), 1, diffuse_position);
        glUniform1f(glGetUniformLocation(mesh_shader.programID, "specularStrenght"), specular_strength);
        glUniform1i(glGetUniformLocation(mesh_shader.programID, "shininess"), shininess);


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
        ImGui::Begin("SETTINGS");

        // WIREFRAME - FILL CHOICE
        if(ImGui::Button("WIREFRAME")) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        if(ImGui::Button("FILL")) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        bool showMesh;
        bool showTubes;
        bool showOctree;
        bool showPotentialOctree;
        ImGui::Checkbox("show mesh", &showMesh);
        ImGui::Checkbox("show tubes", &showTubes);
        ImGui::Checkbox("show gravity force octree", &showOctree);
        ImGui::Checkbox("show potential octree", &showPotentialOctree);
        ImGui::Checkbox("gravity on", &is_gravity_on);

        // MESH LOADING OPERATIONS
        ImGui::Spacing();ImGui::Spacing();ImGui::Spacing();
        if(ImGui::Button("OPEN MESH")) mesh_browser.Open();
        if(!msh.get_vertices().empty()) {
            ImGui::Text("mesh loaded");
        }

        // GRAVITY CALCULATION OPERATIONS
        if(msh.is_loaded()) {
            ImGui::Begin("GRAVITY PROCESSING");
            ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.45f);

            ImGui::InputFloat("G", &G);
            /*
            ImGui::SliderFloat("point x", &potential_point[0], -5, 5);
            ImGui::SliderFloat("point y", &potential_point[1], -5, 5);
            ImGui::SliderFloat("point z", &potential_point[2], -5, 5);
            ImGui::InputFloat3("point", potential_point);
            */
            ImGui::InputInt("tubes res", &tube_resolution);

            if(ImGui::Button("set up tubes")) {
                tubes.erase(tubes.begin(), tubes.end());
                tubes = gravity::get_tubes(msh.get_vertices(), msh.get_faces(), tube_resolution, &cylinder_R);
            }

            // REAL TIME GRAVITY UPDATE WITH GPU
            Timer t3{};
            if(!tubes.empty()) {
                // update tubes
                for(int i = 0; i < n_real_time; i++) {
                    auto transformed_point = glm::inverse(model_matrix) *
                                             glm::vec<4, float>{debug_ball_position[0], debug_ball_position[1],
                                                                debug_ball_position[2], 1};
                    auto output = GPUComputing::get_gravity_from_tubes_with_integral(glm::value_ptr(tubes.front().t1),
                                                                                     tubes.size(),
                                                                                     glm::value_ptr(transformed_point),
                                                                                     cylinder_R, G);


                    glm::vec3 output_gravity{0.f, 0.f, 0.f};
                    omp_set_num_threads(omp_get_max_threads());
                    glm::vec3 thread_gravity[omp_get_max_threads()];

                    for (int i = 0; i < omp_get_max_threads(); i++) {
                        thread_gravity[i] = {0, 0, 0};
                    }


#pragma omp parallel for default(none) shared(output, thread_gravity, tubes, std::cout)
                    for (int i = 0; i < tubes.size(); i++) {
                        auto o = glm::vec3(output[3 * i], output[3 * i + 1], output[3 * i + 2]);
                        if (glm::length(o) > 1000) std::cout << "    " << i << "   ";
                        thread_gravity[omp_get_thread_num()] += glm::vec3(output[3 * i], output[3 * i + 1],
                                                                          output[3 * i + 2]);
                    }

                    for (int i = 0; i < omp_get_max_threads(); i++) {
                        output_gravity += thread_gravity[i];
                    }
                    gravity = output_gravity;
                    delete output;
                }
                //gravity = gravity::get_gravity_from_tubes_with_integral_with_gpu(transformed_point, tubes, G, cylinder_R);
                gravity = model_rotation_z * model_rotation_y * model_rotation_x * glm::vec<4, float>{gravity, 1};

                if(is_gravity_on && octree.size() == 0) {
                    debug_ball_velocity[0] += gravity.x * io.DeltaTime;
                    debug_ball_velocity[1] += gravity.y * io.DeltaTime;
                    debug_ball_velocity[2] += gravity.z * io.DeltaTime;
                } else if(!is_gravity_on) {
                    debug_ball_velocity = glm::make_vec3(ball_initial_velocity);
                }
            }

            real_time_performance = t3.time();

            ImGui::Text("OCTREE PARAMETERS");
            ImGui::InputInt("octree depth", &octree_depth);
            ImGui::InputInt("octree min depth", &min_depth, 0, 10);
            ImGui::InputFloat("accuracy [%]", &precision);
            //ImGui::SliderFloat("alpha", &alpha, 0.f, 10.f);
            ImGui::InputFloat3("min", min);
            ImGui::InputFloat("edge", &edge);
            ImGui::Text("BENCHMARK PARAMETERS");
            ImGui::InputInt("#tests for random benchmark", &n_test);
            ImGui::InputInt("grav. comp. per frame", &n_gravity);
            ImGui::InputInt("pot. comp. per frame", &n_potential);
            //ImGui::SliderInt("space test resolution", &space_test_resolution, 0, 1000);
            ImGui::InputInt("direct comp. per cycle", &n_real_time);
            ImGui::InputInt("octree comp. per cycle", &n_real_time_octree);

            /*
            if(ImGui::Button("compute gpu grid")) {
                Timer t{};
                // cached_gpu_depth: gpu_depth for last gpu grid computed
                cached_gpu_depth = gpu_depth;
                gpu_output_gravity.clear();
                discrete_space.clear();
                // compute discrete space from min, edge and gpu_deph
                discrete_space = gravity::get_discrete_space(glm::make_vec3(min), edge, std::pow(2, gpu_depth));
                // compute gravity grid
                auto output = GPUComputing::get_gravities_from_tubes_with_integral(
                        glm::value_ptr(tubes.front().t1),
                        (int)tubes.size(),
                        glm::value_ptr(discrete_space.front()),
                        (int)discrete_space.size(),
                        cylinder_R, G
                        );
                // fill gpu_output_gravity
                int k = 0;
                for(int i = 0; i < discrete_space.size(); i++) {
                    gpu_output_gravity.emplace_back(output[k], output[k+1], output[k+2]);
                    k += 3;
                }

                delete output;
                std::cout << "gpu grid duration: "; t.log(); std::cout << std::endl;
            }

            if(ImGui::Button("compute potential grid")) {
                cached_gpu_depth = gpu_depth;
                gpu_output_potential.clear();
                discrete_space.clear();
                discrete_space = gravity::get_discrete_space(glm::make_vec3(min), edge, std::pow(2, gpu_depth));
                Timer t{};
                auto output = GPUComputing::get_potentials_from_tubes_with_integral(
                        glm::value_ptr(tubes.front().t1),
                        (int)tubes.size(),
                        glm::value_ptr(discrete_space.front()),
                        (int)discrete_space.size(),
                        cylinder_R, G
                );
                for(int i = 0; i < discrete_space.size(); i++) {
                    gpu_output_potential.emplace_back(output[i]);
                }
                delete output;
                std::cout << "potential gpu grid duration: "; t.log(); std::cout << std::endl;
            }

            // can call it gpu_output_gravity is not empty; it depends on octree_depth
            // that means that if octree_depth changes and gpu_depth doesn't only load grid values
            // must be called again to update integer coordinates -> gravity value mapping done by
            // cache data structures
            if(ImGui::Button("load grid values") && !gpu_output_gravity.empty()) {
                Timer t{};

                // clear maps
                // cached values and tmp_gravity_values contains gpu computed data
                cached_values.clear();
                tmp_gravity_values.clear();

                // cached_octree_depth
                cached_octree_depth = octree_depth;

                // integer distance between gpu computed values in the
                // grid with increased resolution, which consists in the difference between the
                // target depth (octree_depth) and the gpu_depth plus 2, which is for address
                // further computation during should divide procedure.
                // (when computing should divide, points examined belongs to current depth level + 2,
                // and in order to cache them they need to be addressed with integer indexes.
                auto m = (int) std::pow(2, cached_octree_depth + 2 - cached_gpu_depth);

                // here res means number of values per axis
                // for example, if pre_gpu_depth is 1, that means 1 division, so 3 values = 2^1 + 1
                // values per axis = segments per axis + 1
                auto res = (int) std::pow(2, cached_gpu_depth) + 1;

                // cached_values maps integer coordinates -> index for tmp_gravity_values
                for (int i = 0; i < discrete_space.size(); i++) {
                    cached_values.emplace(
                            glm::ivec3{m * (i / (res * res)), m * ((i % (res * res)) / res), m * (i % res)},
                            i
                    );
                    tmp_gravity_values.emplace_back(gpu_output_gravity[i]);
                }

                std::cout << "load grid values duration: ";
                t.log();
                std::cout << " ms" << std::endl;
            }

            if(ImGui::Button("load potential grid values")) {
                Timer t{};
                cached_octree_depth = octree_depth;
                potential_cached_values.clear();
                auto m = (int) std::pow(2, cached_octree_depth - cached_gpu_depth);

                // here res means number of values per axis
                // for example, if pre_gpu_depth is 1, that means 1 division, so 3 values = 2^1 + 1
                // values per axis = segments per axis + 1
                auto res = (int) std::pow(2, cached_gpu_depth) + 1;

                // cached_values maps integer coordinates -> index for tmp_gravity_values
                for (int i = 0; i < discrete_space.size(); i++) {
                    potential_cached_values.emplace(
                            glm::ivec3{m * (i / (res * res)), m * ((i % (res * res)) / res), m * (i % res)},
                            gpu_output_potential[i]
                    );
                }

                std::cout << "load potential grid values duration: ";
                t.log();
                std::cout << " ms" << std::endl;
            }
            */
            //ImGui::SliderInt("should divide method", &should_divide_method, 0, 3);
            if(ImGui::Button("build octree")) {
                Timer t{};
                octree.clear();

                // gravity_values and gravity_values_map contain values computed during octree building
                // while cached_values and tmp_gravity_values contain gpu computed values; these data structures
                // are not modified during octree building, while gravitu_values and gravity_values_map are
                // cleared and filled each time octree building is called
                gravity_values.clear();
                gravity_values_map.clear();

                // reset cached values only when change resolution
                if(cached_octree_depth != octree_depth) {
                    tmp_gravity_values.clear();
                    cached_values.clear();
                    cached_octree_depth = octree_depth;
                }

                gravity::build_octree(
                        should_divide_method,
                        precision,
                        octree,
                        0,
                        cached_octree_depth,
                        cached_octree_depth,
                        min_depth,
                        glm::make_vec3(min),
                        edge,
                        glm::ivec3{0, 0, 0},
                        (int)std::pow(2, cached_octree_depth) * 4,
                        gravity_values,
                        tmp_gravity_values,
                        cached_values,
                        gravity_values_map,
                        tubes, G, cylinder_R,
                        msh.get_vertices(), msh.get_faces()
                        );
                std::cout << "max size: " << (32.f/7.f)*std::pow(8, cached_octree_depth + 1) - 32.f/7.f + 12.f*std::pow((int)std::pow(2, cached_octree_depth) + 1, 3)<< " byte" << std::endl;
                std::cout << "true size: " << octree.size()*4 + gravity_values.size()*12 << " byte" << std::endl;
                //std::cout << "cached values size " << cached_values.size() << std::endl;
                std::cout << "octree size " << 4*octree.size() << std::endl;
                std::cout << "gravity values size " << 12*gravity_values.size() << std::endl;
                //std::cout << "tmp gravity values size " << tmp_gravity_values.size() << std::endl;
                std::cout << "tempo totale: "; t.log(); std::cout << std::endl;

                // build octree rendering data structure
                octree_rendering.clear();
                std::queue<std::tuple<int, glm::vec3, float, float>> q{};
                q.emplace(0, glm::make_vec3(min), edge, -1.f);
                while(!q.empty()) {
                    auto k = q.front(); q.pop();
                    if(octree[std::get<0>(k)] > 0) {
                        auto box = util::get_box(std::get<1>(k), std::get<2>(k) / 2.f);
                        for(int i = 0; i < 8; i++)
                            q.emplace(octree[std::get<0>(k) + i], box[i], std::get<2>(k) / 2.f, std::get<3>(k) + 1);
                    } else {
                        auto m = std::get<1>(k);
                        auto e = std::get<2>(k);
                        auto box = util::get_box(m, e);
                        float d = std::get<3>(k);
                        max_depth_reached = std::max(max_depth_reached, (int)d);
                        octree_rendering.push_back({box[0], d});
                        octree_rendering.push_back({box[1], d});
                        octree_rendering.push_back({box[0], d});
                        octree_rendering.push_back({box[2], d});
                        octree_rendering.push_back({box[1], d});
                        octree_rendering.push_back({box[3], d});
                        octree_rendering.push_back({box[2], d});
                        octree_rendering.push_back({box[3], d});
                        octree_rendering.push_back({box[4], d});
                        octree_rendering.push_back({box[5], d});
                        octree_rendering.push_back({box[4], d});
                        octree_rendering.push_back({box[6], d});
                        octree_rendering.push_back({box[5], d});
                        octree_rendering.push_back({box[7], d});
                        octree_rendering.push_back({box[6], d});
                        octree_rendering.push_back({box[7], d});
                        octree_rendering.push_back({box[0], d});
                        octree_rendering.push_back({box[4], d});
                        octree_rendering.push_back({box[1], d});
                        octree_rendering.push_back({box[5], d});
                        octree_rendering.push_back({box[2], d});
                        octree_rendering.push_back({box[6], d});
                        octree_rendering.push_back({box[3], d});
                        octree_rendering.push_back({box[7], d});
                    }
                }
            }

            if(ImGui::Button("build potential octree")) {
                std::cout << "alpha: " << alpha << std::endl;
                Timer t{};

                if(octree_depth != cached_octree_depth) {
                    potential_cached_values.clear();
                    cached_octree_depth = octree_depth;
                }

                potential_octree.clear();

                gravity::potential::build_octree(
                        should_divide_method,
                        alpha,
                        precision,
                        potential_octree,
                        0,
                        cached_octree_depth,
                        cached_octree_depth,
                        glm::make_vec3(min),
                        edge,
                        glm::ivec3{0, 0, 0},
                        (int)std::pow(2, cached_octree_depth),
                        potential_cached_values,
                        tubes, G, cylinder_R, msh.get_vertices(), msh.get_faces()
                );
                std::cout << "max size without potential: " << (32.f/7.f)*std::pow(8, cached_octree_depth + 1) - 32.f/7.f + 12.f*std::pow((int)std::pow(2, cached_octree_depth) + 1, 3)<< " byte" << std::endl;
                std::cout << "potential max size: " << (32.f/7.f)*std::pow(8, cached_octree_depth + 1) - 32.f/7.f << " byte" << std::endl;
                std::cout << "potential true size: " << potential_octree.size()*4 << " byte" << std::endl;
                std::cout << "tempo totale: "; t.log(); std::cout << std::endl;

                // build octree rendering data structure
                potential_octree_rendering.clear();
                std::queue<std::tuple<int, glm::vec3, float, float>> q{};
                q.emplace(0, glm::make_vec3(min), edge, -1.f);
                while(!q.empty()) {
                    auto k = q.front(); q.pop();
                    if(potential_octree[std::get<0>(k)] > 0) {
                        auto box = util::get_box(std::get<1>(k), std::get<2>(k) / 2.f);
                        for(int i = 0; i < 8; i++)
                            q.emplace(potential_octree[std::get<0>(k) + i], box[i], std::get<2>(k) / 2.f, std::get<3>(k) + 1);
                    } else {
                        auto m = std::get<1>(k);
                        auto e = std::get<2>(k);
                        auto box = util::get_box(m, e);
                        float d = std::get<3>(k);
                        potential_max_depth_reached = std::max(potential_max_depth_reached, (int)d);
                        potential_octree_rendering.push_back({box[0], d});
                        potential_octree_rendering.push_back({box[1], d});
                        potential_octree_rendering.push_back({box[0], d});
                        potential_octree_rendering.push_back({box[2], d});
                        potential_octree_rendering.push_back({box[1], d});
                        potential_octree_rendering.push_back({box[3], d});
                        potential_octree_rendering.push_back({box[2], d});
                        potential_octree_rendering.push_back({box[3], d});
                        potential_octree_rendering.push_back({box[4], d});
                        potential_octree_rendering.push_back({box[5], d});
                        potential_octree_rendering.push_back({box[4], d});
                        potential_octree_rendering.push_back({box[6], d});
                        potential_octree_rendering.push_back({box[5], d});
                        potential_octree_rendering.push_back({box[7], d});
                        potential_octree_rendering.push_back({box[6], d});
                        potential_octree_rendering.push_back({box[7], d});
                        potential_octree_rendering.push_back({box[0], d});
                        potential_octree_rendering.push_back({box[4], d});
                        potential_octree_rendering.push_back({box[1], d});
                        potential_octree_rendering.push_back({box[5], d});
                        potential_octree_rendering.push_back({box[2], d});
                        potential_octree_rendering.push_back({box[6], d});
                        potential_octree_rendering.push_back({box[3], d});
                        potential_octree_rendering.push_back({box[7], d});
                    }
                }
            }
            /*
            if(!tubes.empty()) {
                float volume = 0;
                for(int i = 0; i < tubes.size(); i++) {
                    volume += std::pow(cylinder_R, 2)*M_PI*(glm::distance(tubes[i].t1, tubes[i].t2));
                }
                std::cout << "volume sfera: " << volume << std::endl;
            }*/

            if(!octree.empty() && ImGui::Button("gravity random benchmark")) {
                Timer t{};
                auto e = gravity::benchmark(
                        util::random_locations(glm::make_vec3(min), edge, n_test),
                        tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                        );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "gravity benchmark: " << e * 100 << std::endl;
            }

            if(!potential_octree.empty() && ImGui::Button("potential random benchmark")) {
                Timer t{};
                auto e = gravity::potential::benchmark(
                    util::random_locations(glm::make_vec3(min), edge, n_test),
                    tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "potential benchmark: " << e * 100 << std::endl;
            }

            if(!octree.empty() && ImGui::Button("near mesh gravity benchmark")) {
                Timer t{};
                auto e = gravity::benchmark(
                        util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                        tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "gravity benchmark: " << e * 100 << std::endl;
            }

            if(!potential_octree.empty() && ImGui::Button("near mesh potential benchmark")) {
                Timer t{};
                auto e = gravity::potential::benchmark(
                        util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                        tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "potential benchmark: " << e * 100 << std::endl;
            }

            if(!octree.empty() && ImGui::Button("outside mesh gravity benchmark")) {
                Timer t{};
                auto e = gravity::benchmark(
                        util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min), edge),
                        tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "gravity benchmark: " << e * 100 << std::endl;
            }

            if(!potential_octree.empty() && ImGui::Button("outside mesh potential benchmark")) {
                Timer t{};
                auto e = gravity::potential::benchmark(
                        util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min), edge),
                        tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                );
                std::cout << "tempo impiegato per il benchmark con " << n_test << " campioni: " << t.time() << std::endl;
                std::cout << "potential benchmark: " << e * 100 << std::endl;
            }

            if(ImGui::Button("gravity octree log data")) {
                int i = 0;
                std::string file_name = "./gravity_octree_log_" + std::to_string(i) + ".txt";
                while(file_exists(file_name)) {
                    i++; file_name = "./gravity_octree_log_" + std::to_string(i) + ".txt";
                }

                std::vector<float> precisions = {
                        10, 9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3,
                        2.5, 2, 1.8, 1.5, 1.3, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
                        0.1, 0.09, 0.08, 0.07, 0.06, 0.05
                };
                auto file = std::ofstream(file_name, std::ios::app);
                file << "method: " << should_divide_method << std::endl;
                cached_values.clear();
                tmp_gravity_values.clear();
                for(auto p : precisions) {
                    octree.clear();
                    gravity_values.clear();
                    gravity_values_map.clear();
                    gravity::build_octree(
                            should_divide_method,
                            p,
                            octree,
                            0,
                            20,
                            20,
                            0,
                            glm::make_vec3(min),
                            edge,
                            glm::ivec3{0, 0, 0},
                            (int) std::pow(2, 20) * 4,
                            gravity_values,
                            tmp_gravity_values,
                            cached_values,
                            gravity_values_map,
                            tubes, G, cylinder_R,
                            msh.get_vertices(), msh.get_faces()
                    );

                    file << (octree.size() * 4 + gravity_values.size() * 12)/1000.f << ";";
                    file << gravity::benchmark(
                            util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                            tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                    ) * 100 << ";";
                    file << gravity::benchmark(
                            util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min), edge),
                            tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                    ) * 100 << ";";
                    file << gravity::benchmark(
                            util::random_locations(glm::make_vec3(min), edge, 10000),
                            tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                    ) * 100 << std::endl;
                    /*
                    file << "precision: " << p << std::endl;
                    file << "storage: " << octree.size() * 4 + gravity_values.size() * 12 << " byte" << std::endl;
                    file << "benchmark near mesh: " << gravity::benchmark(
                            util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                            tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                    ) * 100 << "%" << std::endl;
                    file << "benchmark outside mesh: " << gravity::benchmark(
                            util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min),
                                                         edge),
                            tubes, cylinder_R, G, octree, gravity_values, glm::make_vec3(min), edge
                    ) * 100 << "%" << std::endl;
                    file << std::endl << std::endl << std::endl;
                    */
                }
                file.close();
            }

            if(ImGui::Button("potential octree log data")) {
                int i = 0;
                std::string file_name = "./potential_octree_log_" + std::to_string(i) + ".txt";
                while(file_exists(file_name)) {
                    i++; file_name = "./potential_octree_log_" + std::to_string(i) + ".txt";
                }

                std::vector<float> precisions = {
                        10, 9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3,
                        2.5, 2, 1.8, 1.5, 1.3, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2,
                        0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01
                };

                auto file = std::ofstream(file_name, std::ios::app);
                file << "method: " << should_divide_method << std::endl;
                potential_cached_values.clear();
                for(auto p : precisions) {
                    potential_octree.clear();
                    gravity::potential::build_octree(
                            should_divide_method,
                            p,
                            p,
                            potential_octree,
                            0,
                            20,
                            20,
                            glm::make_vec3(min),
                            edge,
                            glm::ivec3{0, 0, 0},
                            (int)std::pow(2, 20),
                            potential_cached_values,
                            tubes, G, cylinder_R, msh.get_vertices(), msh.get_faces()
                    );
                    file << (float(potential_octree.size()) * 4.f)/1000.f << ";";
                    file << gravity::potential::benchmark(
                            util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                            tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                    ) * 100 << ";";
                    file << gravity::potential::benchmark(
                            util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min),
                                                         edge),
                            tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                    ) * 100 << ";";
                    file << gravity::potential::benchmark(
                            util::random_locations(glm::make_vec3(min), edge, 10000),
                            tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                    ) * 100 << std::endl;
                    /*
                    file << "p: " << p << std::endl;
                    file << "storage: " << potential_octree.size() * 4 << " byte" << std::endl;
                    file << "benchmark near mesh: " << gravity::potential::benchmark(
                            util::near_mesh_locations(msh.get_vertices(), msh.get_faces()),
                            tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                    ) * 100 << "%" << std::endl;
                    file << "benchmark outside mesh: " << gravity::potential::benchmark(
                            util::outside_mesh_locations(msh.get_vertices(), msh.get_faces(), glm::make_vec3(min),
                                                         edge),
                            tubes, cylinder_R, G, potential_octree, glm::make_vec3(min), edge
                    ) * 100 << "%" << std::endl;
                    file << std::endl << std::endl << std::endl;
                    */

                }
                file.close();
            }

            /*

            // get locations from octree
            // locations are placed in the same order of sampled values vector
            if(ImGui::Button("set up locations vector for potential computation")) {
                octree_space_locations.clear();
                octree_space_locations.resize(gravity_values.size());
                std::queue<std::tuple<int, glm::vec3, float>> q{};
                q.emplace(0, glm::make_vec3(min), edge);
                while(!q.empty()) {
                    auto k = q.front(); q.pop();
                    if(octree[std::get<0>(k)] > 0) {
                        auto box = util::get_box(std::get<1>(k), std::get<2>(k) / 2.f);
                        for(int i = 0; i < 8; i++)
                            q.emplace(octree[std::get<0>(k) + i], box[i], std::get<2>(k) / 2.f);
                    } else {
                        auto box = util::get_box(std::get<1>(k), std::get<2>(k));
                        for(int i = 0; i < 8; i++) {
                            octree_space_locations[-octree[std::get<0>(k) + i]] = box[i];
                        }
                    }
                }
            }

            // compute a potential values vector from locations vector
            if(ImGui::Button("fill potential values")) {
                potential_values.clear();
                auto output = GPUComputing::get_potentials_from_tubes_with_integral(
                        glm::value_ptr(tubes.front().t1), tubes.size(),
                        glm::value_ptr(octree_space_locations.front()),
                        octree_space_locations.size(), cylinder_R, G
                        );
                for(int i = 0; i < octree_space_locations.size(); i++) {
                    potential_values.emplace_back(output[i]);
                }
            }

            */

            // WRITE FILE
            // file structure
            // min (3 * float) edge (float) octree_size (32bit int) octree values
            // gravity values size (32bit int) gravity values
            // endl endl endl "file info: " ...
            if(!octree.empty()) {
                char file_name[20];
                ImGui::InputText("file name", file_name, 20);
                if(ImGui::Button("write octree")) {
                    std::ofstream ofs(file_name, std::ofstream::binary);
                    int octree_type = 0;
                    ofs.write(reinterpret_cast<char*>(&octree_type), sizeof(int));
                    ofs.write(reinterpret_cast<char*>(min), 3*sizeof(float));
                    ofs.write(reinterpret_cast<char*>(&edge), sizeof(float));
                    int octree_size = octree.size();
                    ofs.write(reinterpret_cast<char*>(&octree_size), sizeof(int));
                    for(int i = 0; i < octree_size; i++) {
                        ofs.write(reinterpret_cast<char*>(&octree[i]), sizeof(int));
                    }
                    int gravity_values_size = gravity_values.size();
                    ofs.write(reinterpret_cast<char*>(&gravity_values_size), sizeof(int));
                    for(int i = 0; i < gravity_values_size; i++) {
                        ofs.write(reinterpret_cast<char*>(glm::value_ptr(gravity_values[i])), 3*sizeof(float));
                    }
                    ofs << std::endl << std::endl << std::endl;
                    ofs << "file info:" << std::endl;
                    ofs << "mesh: " << mesh_browser.GetSelected().string() << std::endl;
                    ofs << "min: " << min[0] << " " << min[1] << " " << min[2] << std::endl;
                    ofs << "edge: " << edge << std::endl;
                    ofs << "this is a gravity values file" << std::endl;
                    ofs << "last should divide iteration method: " << std::endl;
                    ofs << "should divide method: " << std::endl;
                    ofs << "octree depth: " << cached_octree_depth;
                }
            }

            if(!potential_octree.empty()) {
                char file_name[20];
                ImGui::InputText("file name", file_name, 20);
                if(ImGui::Button("write potential octree")) {
                    std::ofstream ofs(file_name, std::ofstream::binary);
                    int octree_type = 1;
                    ofs.write(reinterpret_cast<char*>(&octree_type), sizeof(int));
                    ofs.write(reinterpret_cast<char*>(min), 3*sizeof(float));
                    ofs.write(reinterpret_cast<char*>(&edge), sizeof(float));
                    int octree_size = potential_octree.size();
                    ofs.write(reinterpret_cast<char*>(&octree_size), sizeof(int));
                    for(int i = 0; i < octree_size; i++) {
                        ofs.write(reinterpret_cast<char*>(&potential_octree[i]), sizeof(int));
                    }

                    ofs << std::endl << std::endl << std::endl;
                    ofs << "file info:" << std::endl;
                    ofs << "mesh: " << mesh_browser.GetSelected().string() << std::endl;
                    ofs << "min: " << min[0] << " " << min[1] << " " << min[2] << std::endl;
                    ofs << "edge: " << edge << std::endl;
                    ofs << "this is a potential file" << std::endl;
                    ofs << "last should divide iteration method: " << std::endl;
                    ofs << "should divide method: " << std::endl;
                    ofs << "octree depth: " << cached_octree_depth;
                }
            }

            auto space_test = gravity::get_discrete_space(glm::make_vec3(min), edge, space_test_resolution);

            // GRAVITY COMPUTATION FROM OCTREE IN REAL TIME
            Timer t{};
            int d = 0;
            if(!octree.empty()) {
                auto p = debug_ball_position;
                if(util::is_inside_box(glm::make_vec3(p), util::get_box(glm::make_vec3(min), edge))) {
                    gravity_from_octree = gravity::get_gravity_from_octree(glm::make_vec3(p), octree, glm::make_vec3(min), edge, gravity_values, &d);
                }

                // repeat for performance mesurement
                for(int i = 1; i < n_real_time_octree; i++) {
                    gravity::get_gravity_from_octree(glm::make_vec3(p), octree, glm::make_vec3(min), edge, gravity_values, &d);
                }
            }
            ImGui::Text("gravity computation time: %d", t.time());
            real_time_performance_octree = t.time();

            // GRAVITY COMPUTATION FROM POTENTIAL OCTREE IN REAL TIME
            Timer t2{};
            int potential_d = 0;
            if(!potential_octree.empty()) {
                auto p = glm::make_vec3(debug_ball_position);
                if(util::is_inside_box(p, util::get_box(glm::make_vec3(min), edge))) {
                    gravity_from_octree_potential = gravity::potential::get_gravity_from_octree(p, potential_octree, glm::make_vec3(min), edge, &potential_d);
                }
            }
            ImGui::Text("potential computation time: %d", t2.time());

            ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

            // GRAVITY OUTPUTS
            ImGui::Text("gravity real time:");
            ImGui::Text("%.2f %.2f %.2f  length: %.2f", gravity.x, gravity.y, gravity.z, glm::length(gravity));
            ImGui::Text("gravity from octree:");
            ImGui::Text("%.2f %.2f %.2f  length: %.2f at depth: %d",
                    gravity_from_octree.x,
                    gravity_from_octree.y,
                    gravity_from_octree.z,
                    glm::length(gravity_from_octree),
                    d
                    );
            ImGui::Text("gravity from potential octree:");
            ImGui::Text("%.2f %.2f %.2f  length: %.2f at depth: %d",
                    gravity_from_octree_potential.x,
                    gravity_from_octree_potential.y,
                    gravity_from_octree_potential.z,
                    glm::length(gravity_from_octree_potential),
                    potential_d
            );

            ImGui::End();
        }
        ImGui::End();

        if(is_gravity_on && octree.size() > 0) {
            debug_ball_velocity[0] += gravity_from_octree.x * io.DeltaTime;
            debug_ball_velocity[1] += gravity_from_octree.y * io.DeltaTime;
            debug_ball_velocity[2] += gravity_from_octree.z * io.DeltaTime;
        }

        if(is_gravity_on) {
            debug_ball_position[0] += debug_ball_velocity.x * io.DeltaTime;
            debug_ball_position[1] += debug_ball_velocity.y * io.DeltaTime;
            debug_ball_position[2] += debug_ball_velocity.z * io.DeltaTime;
        }

        // IMGUI TAB: TEXTURE OPERATIONS
        ImGui::Begin("Texture");
        if(ImGui::Button("OPEN TEXTURE")) texture_browser.Open();
        if(ImGui::Button("SET COLORMODE VERTEX COLOR")) {
            mesh_shader.use();
            glUniform1i(glGetUniformLocation(mesh_shader.programID, "colorMode"), 1);
        }
        if(ImGui::Button("SET COLORMODE TEXTURE")) {
            mesh_shader.use();
            glUniform1i(glGetUniformLocation(mesh_shader.programID, "colorMode"), 0);
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
        ImGui::Text("REAL TIME PERFORMANCE: %f", real_time_performance);
        ImGui::Text("REAL TIME PERFORMANCE OCTREE: %f", real_time_performance_octree);
        ImGui::End();

        // IMGUI TAB FOR MODEL TRANSFORM
        ImGui::Begin("Model Transform");
        ImGui::SliderFloat("model scale", &scale, 0.1f, 5.0f);
        ImGui::SliderFloat3("model position", position, -20.0, 20.0);
        ImGui::SliderFloat3("model rotation", rotation, 0.0f, 2.0f * M_PI);

        ImGui::SliderFloat("debug ball scale", &debug_ball_scale, 0.1f, 1.f);
        ImGui::InputFloat3("debug ball position", debug_ball_position);
        ImGui::InputFloat3("debug_ball_position", debug_ball_position);
        ImGui::InputFloat3("initial velocity", ball_initial_velocity);
        if(ImGui::Button("reset ball velocity")) debug_ball_velocity = {0.f, 0.f, 0.f};
        ImGui::Text("ball velocity: %f %f %f",
                    debug_ball_velocity.x,
                    debug_ball_velocity.y,
                    debug_ball_velocity.z);
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
                mesh_shader.use();
                glUniform1i(glGetUniformLocation(mesh_shader.programID, "colorMode"), 1);
            } else {
                mesh_shader.use();
                glUniform1i(glGetUniformLocation(mesh_shader.programID, "colorMode"), 0);
            }
            mesh_browser.ClearSelected();

            // LOAD ARROW MESH FOR GRAVITY DEBUG
            arrow.load_from_obj("../resources/arrow.obj");
            debug_ball.load_from_obj("../resources/sphere.obj");
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
            mesh_shader.use();
            // bind VERTEX ARRAY OBJECT
            glBindVertexArray(msh.get_VAO());
            // DRAW
            if(showMesh)
                glDrawElements(GL_TRIANGLES, (int)msh.get_elements().size(), GL_UNSIGNED_INT, nullptr);

            glBindVertexArray(debug_ball.get_VAO());
            glUniformMatrix4fv(glGetUniformLocation(mesh_shader.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(debug_ball_model_matrix));
            glDrawElements(GL_TRIANGLES, (int)debug_ball.get_elements().size(), GL_UNSIGNED_INT, nullptr);

            glBindVertexArray(arrow.get_VAO());
            glUniformMatrix4fv(glGetUniformLocation(mesh_shader.programID, "modelMatrix"), 1, GL_FALSE, glm::value_ptr(arrow_model_matrix));
            glDrawElements(GL_TRIANGLES, (int)arrow.get_elements().size(), GL_UNSIGNED_INT, nullptr);
        }

        // RAY
        debug_shader.use();

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
        glUniform4f(glGetUniformLocation(debug_shader.programID, "color"), 0.8, 0, 0, 1);
        glDrawArrays(GL_LINES, 0, 2);
        glUniform4f(glGetUniformLocation(debug_shader.programID, "color"), 0, 0.8, 0, 1);
        glDrawArrays(GL_LINES, 2, 2);
        glUniform4f(glGetUniformLocation(debug_shader.programID, "color"), 0, 0, 0.8, 1);
        glDrawArrays(GL_LINES, 4, 2);

        // TUBES DRAW
        if(!tubes.empty() && showTubes) {
            glBufferData(GL_ARRAY_BUFFER, tubes.size() * 6 * sizeof(float), glm::value_ptr(tubes.front().t1), GL_DYNAMIC_DRAW);
            glDrawArrays(GL_LINES, 0, tubes.size() * 2);
        }

        // octree draw
        if(!octree.empty() && showOctree) {
            octree_shader.use();
            glUniform4f(glGetUniformLocation(octree_shader.programID, "color"), 0, 0, 0.8, 1);
            glUniform1i(glGetUniformLocation(octree_shader.programID, "max_depth"), max_depth_reached);
            glBindBuffer(GL_ARRAY_BUFFER, octreeVBO);
            glBufferData(GL_ARRAY_BUFFER, octree_rendering.size() * 4 * sizeof(float), glm::value_ptr(octree_rendering.front().p), GL_DYNAMIC_DRAW);

            glBindVertexArray(octreeVAO);
            glDrawArrays(GL_LINES, 0, octree_rendering.size());
        }

        // potential octree draw
        if(!potential_octree.empty() && showPotentialOctree) {
            octree_shader.use();
            glUniform4f(glGetUniformLocation(octree_shader.programID, "color"), 0, 0, 0.8, 1);
            glUniform1i(glGetUniformLocation(octree_shader.programID, "max_depth"), potential_max_depth_reached);
            glBindBuffer(GL_ARRAY_BUFFER, octreeVBO);
            glBufferData(GL_ARRAY_BUFFER, potential_octree_rendering.size() * 4 * sizeof(float), glm::value_ptr(potential_octree_rendering.front().p), GL_DYNAMIC_DRAW);

            glBindVertexArray(octreeVAO);
            glDrawArrays(GL_LINES, 0, potential_octree_rendering.size());
        }

        // SWAP WINDOWS
        SDL_GL_SwapWindow(window);
    }

    SDL_Log("quitting gracefully");
}

bool file_exists(const std::string& name) {
      std::ifstream f(name.c_str());
      return f.good();
}

void set_up_gl_attributes() {
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
}

void set_up_gl() {
    // ENABLE VSYNC
    SDL_GL_SetSwapInterval(1);

    // ENABLE DEPTH TEST AND SET DEPTH FUNCTION
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // GL_LESS: minor z passes
}