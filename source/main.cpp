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


  // ******************************************************** //
 // ************************** MAIN ************************ //
// ******************************************************** //
int main(int argv, char** args) {
    std::cout << "NUMERIC LIMIT" <<  std::numeric_limits<float>::epsilon();

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
    // this shader is used for meshes; only position is needed; normals are required for lightening
    // vs.glsl vertex shader support position, vertex color, uv and normals and out attributes
    // needed by fragment shader for lightening, texture and vertex color
    // fs.glsl fragment shader support vertex color, texture and lightning
    // for info about layout attributes and uniforms see glsl source code
    shader mesh_shader = shader("../shaders/vs_0.glsl", "../shaders/fs_0.glsl");

    // this shader is used for debug representations
    // vs_tetra.glsl shader support position without model matrix and doesn't pass any data to fragment
    // shader
    // fs_tetra.glsl support only color with uniform, no lightening nor texture
    // used for stuff like lines (axis, tubes), particles,
    shader debug_shader = shader("../shaders/vs_1.glsl", "../shaders/fs_1.glsl");


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

    // output holder for gravity
    glm::vec3 gravity{0.f, 0.f, 0.f};

    glm::vec3 gravity_from_octree{0.f, 0.f, 0.f};

    glm::vec3 gravity_from_potential_gradient{0.f, 0.f, 0.f};

    // output holders for gravity vector calculated by GPUComputing and space vector
    std::vector<glm::vec3> discrete_space{};
    std::vector<glm::vec3> gpu_output_gravity{};

    // parameters used where function parameter "resolution" is required.
    int gravity_resolution = 128;

    // mass radius
    float sphere_R;
    float cylinder_R;

    // G constant
    float G = 10;

    // masses container
    std::vector<gravity::mass> masses{};
    std::vector<gravity::tube> tubes{};

    // octree data
    auto octree = std::vector<int>();
    // gravity values; they are indexed with negatives valus; -octree.value -1
    auto gravity_values = std::vector<glm::vec3>();

    // area of octree computation
    float min[3] {-3.f, -3.f, -3.f};
    float edge = 6.f;

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
    debug_shader.use();
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
        mesh_shader.use();
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "viewerPosition"), 1, glm::value_ptr(camera_position3));

        // MODEL MATRIX
        glm::mat4 model_translation = glm::translate(glm::mat4(1.f), {position[0], position[1], position[2]});
        glm::mat4 model_rotation_x = glm::rotate(glm::mat4(1.f), rotation[0], {1.f, 0.f, 0.f});
        glm::mat4 model_rotation_y = glm::rotate(glm::mat4(1.f), rotation[1], {0.f, 1.f, 0.f});
        glm::mat4 model_rotation_z = glm::rotate(glm::mat4(1.f), rotation[2], {0.f, 0.f, 1.f});
        glm::mat4 model_scale = glm::scale(glm::mat4(1.f), glm::vec3{scale});
        glm::mat4 model_matrix = model_scale * model_translation * model_rotation_z * model_rotation_y * model_rotation_x;

        // update ball position
        if(is_gravity_on) {
            debug_ball_position[0] += debug_ball_velocity.x * io.DeltaTime;
            debug_ball_position[1] += debug_ball_velocity.y * io.DeltaTime;
            debug_ball_position[2] += debug_ball_velocity.z * io.DeltaTime;
        }

        glm::mat4 debug_ball_model_translation = glm::translate(
            glm::mat4(1.f), {debug_ball_position[0], debug_ball_position[1], debug_ball_position[2]}
            );
        glm::mat4 debug_ball_model_scale = glm::scale(glm::mat4(1.f), glm::vec3{debug_ball_scale});
        auto debug_ball_model_matrix = debug_ball_model_translation * debug_ball_model_scale;

        glm::mat4 arrow_model_matrix = glm::translate(glm::mat4(1.f), glm::make_vec3(potential_point))
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

        // light parameters
        mesh_shader.use();
        glUniform1f(glGetUniformLocation(mesh_shader.programID, "ambientLightIntensity"), ambient_intensity);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "ambientLightColor"), 1, ambient_color);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "diffuseLightColor"), 1, diffuse_color);
        glUniform3fv(glGetUniformLocation(mesh_shader.programID, "diffuseLightPosition"), 1, diffuse_position);
        glUniform1f(glGetUniformLocation(mesh_shader.programID, "specularStrenght"), specular_strength);
        glUniform1i(glGetUniformLocation(mesh_shader.programID, "shininess"), shininess);


          // ******************************************************* //
         // ************** UPDATE.(GRAVITY PROCESSING) ************ //
        // ******************************************************* //

        /*
         * Gravity processing is done real time with tetrahedrons method;
         * volume processing is done real time
         * other gravity processing methods are executed by gui interactions
         */

        // UPDATE MESH VOLUME
        // volume = gravity::volume(msh.get_vertices(), msh.get_faces(), {2.f, 2.f, 2.f});


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
        bool showMesh;
        bool showTubes;
        ImGui::Checkbox("show mesh", &showMesh);
        ImGui::Checkbox("show tubes", &showTubes);
        ImGui::Checkbox("gravity on", &is_gravity_on);
        if(ImGui::Button("reset ball velocity")) debug_ball_velocity = {0.f, 0.f, 0.f};

        // MESH LOADING OPERATIONS
        ImGui::Spacing();ImGui::Spacing();ImGui::Spacing();
        if(ImGui::Button("OPEN MESH")) mesh_browser.Open();
        if(!msh.get_vertices().empty()) {
            ImGui::Text("Mesh loaded");
            /*if(ImGui::Button("Save Mesh")) {
                msh.write_on_disk_as_obj("saved_mesh.obj");
            }*/
        }

        // GRAVITY CALCULATION OPERATIONS
        if(msh.is_loaded()) {
            ImGui::Begin("GRAVITY PROCESSING");
            ImGui::Text("Volume %f", volume);

            ImGui::SliderInt("ray gravity resolution", &gravity_resolution, 0, 2001);
            ImGui::SliderFloat("G", &G, 0, 100);
            ImGui::SliderFloat("point x", &potential_point[0], -5, 5);
            ImGui::SliderFloat("point y", &potential_point[1], -5, 5);
            ImGui::SliderFloat("point z", &potential_point[2], -5, 5);
            ImGui::InputFloat3("potential point", potential_point);

            if(ImGui::Button("set up masses")) {
                masses.erase(masses.begin(), masses.end());
                masses = gravity::get_masses(msh.get_vertices(), msh.get_faces(), gravity_resolution, &sphere_R);
            }
            if(ImGui::Button("set up tubes")) {
                tubes.erase(tubes.begin(), tubes.end());
                tubes = gravity::get_tubes(msh.get_vertices(), msh.get_faces(), gravity_resolution, &cylinder_R);
            }
            if(ImGui::Button("calculate gravity with masses")) {
                gravity = gravity::get_gravity_from_masses(masses, 10, sphere_R, glm::make_vec3(potential_point));
            }
            if(ImGui::Button("calculate gravity with tubes")) {
                gravity = gravity::get_gravity_from_tubes(msh.get_vertices(), gravity_resolution, tubes, glm::make_vec3(potential_point));
            }
            if(ImGui::Button("calculate gravity with tubes integral")) {
                gravity = gravity::get_gravity_from_tubes_with_integral(glm::make_vec3(potential_point), tubes, G, cylinder_R);
            }
            float offset;
            ImGui::InputFloat("gradient offset", &offset, 0, 0, "%.7f");
            if(ImGui::Button("calculate gravity from potential gradient")) {
                gravity_from_potential_gradient = util::get_gradient(
                    {
                        gravity::get_potential_from_tubes_with_integral({potential_point[0] - offset, potential_point[1], potential_point[2]}, tubes, G, cylinder_R),
                        gravity::get_potential_from_tubes_with_integral({potential_point[0] + offset, potential_point[1], potential_point[2]}, tubes, G, cylinder_R),
                        gravity::get_potential_from_tubes_with_integral({potential_point[0], potential_point[1] - offset, potential_point[2]}, tubes, G, cylinder_R),
                        gravity::get_potential_from_tubes_with_integral({potential_point[0], potential_point[1] + offset, potential_point[2]}, tubes, G, cylinder_R),
                        gravity::get_potential_from_tubes_with_integral({potential_point[0], potential_point[1], potential_point[2] - offset}, tubes, G, cylinder_R),
                        gravity::get_potential_from_tubes_with_integral({potential_point[0], potential_point[1], potential_point[2] + offset}, tubes, G, cylinder_R),
                    }, offset*2.0f
                    );
            }

            int pre_gpu_depth;
            int octree_depth;

            //if(!gpu_output_gravity.empty()) {
            if(ImGui::Button("compute gravity from gpu vector")) {
                gravity_from_octree = gravity::get_gravity_from_1D_precomputed_vector(
                    glm::make_vec3(debug_ball_position),
                    gpu_output_gravity,
                    discrete_space,
                    (int)std::pow(2, pre_gpu_depth)
                    );
            }


            // REAL TIME
            //if(ImGui::Button("compute gravity with tubes integral and GPU")) {
                //Timer t1{};
                //Timer t2{};
                //t1.log();
                //t2.log();
                if(!tubes.empty()) {
                    // update tubes

                    auto transformed_point = glm::inverse(model_matrix) * glm::vec<4, float>{debug_ball_position[0], debug_ball_position[1], debug_ball_position[2], 1};

                    auto output = GPUComputing::get_gravity_from_tubes_with_integral(glm::value_ptr(tubes.front().t1), tubes.size(), glm::value_ptr(transformed_point), cylinder_R, G);
                    //t1.log();

                    glm::vec3 output_gravity{0.f, 0.f, 0.f};
                    omp_set_num_threads(omp_get_max_threads());
                    glm::vec3 thread_gravity[omp_get_max_threads()];

                    for(int i = 0; i < omp_get_max_threads(); i++) {
                        thread_gravity[i] = {0, 0, 0};
                    }


#pragma omp parallel for default(none) shared(output, thread_gravity, tubes)
                    for(int i = 0; i < tubes.size(); i++) {
                        thread_gravity[omp_get_thread_num()] += glm::vec3(output[3*i], output[3*i + 1], output[3*i + 2]);
                    }

                    for(int i = 0; i < omp_get_max_threads(); i++) {
                        output_gravity += thread_gravity[i];
                    }
                    gravity = output_gravity;
                    gravity = model_rotation_z * model_rotation_y * model_rotation_x * glm::vec<4, float>{gravity, 1};

                    if(is_gravity_on) {
                        debug_ball_velocity[0] += gravity.x * io.DeltaTime;
                        debug_ball_velocity[1] += gravity.y * io.DeltaTime;
                        debug_ball_velocity[2] += gravity.z * io.DeltaTime;
                    }

                    delete output;

                    //t2.log();
                    //}
                }

            ImGui::SliderInt("octree max depth", &octree_depth, 0, 12);
            ImGui::SliderInt("pre gpu depth", &pre_gpu_depth, 0, 12);
            float precision;
            ImGui::SliderFloat("octree precision", &precision, 0.f, 10.f);
            ImGui::InputFloat3("min vertex of gravity sample cube", min);
            ImGui::InputFloat("edge of gravity sample cube", &edge);

            if(ImGui::Button("test new octree")) {
                gravity::node n = gravity::build_node_with_integral(
                    glm::make_vec3(min),
                    edge,
                    tubes, 10, sphere_R
                    );
                auto octree = std::vector<gravity::node>();
                octree.push_back(n);
                gravity::build_octree_with_integral(precision, octree, 0, octree_depth,
                     glm::make_vec3(min),
                    edge, tubes, 10, sphere_R);

                std::cout << "octree size: " << octree.size() << std::endl;
            }

            if(ImGui::Button("compute gravity grid with gpu given octree depth, min and edge")) {
                Timer t{};
                auto space = gravity::get_discrete_space(glm::make_vec3(min), edge, std::pow(2, pre_gpu_depth));
                auto output = GPUComputing::get_gravities_from_tubes_with_integral(
                        glm::value_ptr(tubes.front().t1),
                        (int)tubes.size(),
                        glm::value_ptr(space.front()),
                        (int)space.size(),
                        cylinder_R, G
                        );

                gpu_output_gravity.clear();
                gpu_output_gravity.resize(space.size());
                int k = 0;
                for(int i = 0; i < space.size(); i++) {
                    gpu_output_gravity[i] = glm::vec3{output[k], output[k+1], output[k+2]};
                    k += 3;
                }
                discrete_space.clear();
                discrete_space = space;

                delete output;
                t.log();
            }

            if(ImGui::Button("test new octree optimized")) {
                Timer t{};

                // here resolution stands for segments number for axis
                auto space = gravity::get_discrete_space(glm::make_vec3(min), edge, std::pow(2, pre_gpu_depth));
                auto output = GPUComputing::get_gravities_from_tubes_with_integral(
                        glm::value_ptr(tubes.front().t1),
                        (int)tubes.size(),
                        glm::value_ptr(space.front()),
                        (int)space.size(),
                        cylinder_R, G
                );
                std::cout << "init gpu";
                t.log();


                octree.clear();
                gravity_values.clear();


                // if int value < 0 -> value is in tmp_gravity_value with index -i - 1
                // if int value >= 0 -> value is in gravity_value
                auto cached_values = std::unordered_map<glm::ivec3, int>();

                // these are temporary gravity values: they are cached (already computed) gravity values,
                // they will be added to gravity values if they will be used, in fact gravity_values
                // size is always less than tmp_gravity_values size
                auto tmp_gravity_values = std::vector<glm::vec3>();
                tmp_gravity_values.reserve(space.size());
                std::cout << "space size" << space.size() << std::endl;

                // integer distance between gpu computed values in the
                // grid with increased resolution, which consists in the difference between the
                // target depth (octree_depth) and the pre computed depth (pre_gpu_depth) plus
                // 2, which is for address further computation during should divide procedure.
                // (when computing should divide, points examined belongs to current depth level + 2,
                // and in order to cache them they need to be addressed with integer indexes.
                auto m = (int)std::pow(2, octree_depth + 2 - pre_gpu_depth);

                // here res means number of values per axis
                // for example, if pre_gpu_depth is 1, that means 1 division, so 3 values = 2^1 + 1
                // values per axis = segments per axis + 1
                auto res = (int)std::pow(2, pre_gpu_depth) + 1;
                for(int i = 0; i < space.size(); i++) {
                    cached_values.emplace(
                            glm::ivec3{m*(i / (res*res)), m*((i % (res*res))/ res), m*(i % res)},
                            -i - 1
                            );
                    // WARNING -> cached values are pushed from index 0; cached index values are -index - 1
                    // to retrieve actual index, true_index = -cached_value - 1
                    tmp_gravity_values.push_back({output[3*i], output[3*i + 1], output[3*i + 2]});
                }

                std::cout << "preload values in tmp gravity values";
                t.log();

                gravity::build_octree_with_integral_optimized(
                        precision,
                        octree,
                        0,
                        octree_depth,
                        glm::make_vec3(min),
                        edge,
                        glm::ivec3{0, 0, 0},
                        (int)std::pow(2, octree_depth) * 4,
                        gravity_values,
                        tmp_gravity_values,
                        cached_values,
                        tubes, G, cylinder_R
                        );
                std::cout << "max size: " << (32.f/7.f)*std::pow(8, octree_depth + 1) - 32.f/7.f + 12.f*std::pow((int)std::pow(2, octree_depth) + 1, 3)<< " byte" << std::endl;
                std::cout << "true size: " << octree.size()*4 + gravity_values.size()*12 << " byte" << std::endl;
                std::cout << "cached values size " << cached_values.size() << std::endl;
                std::cout << "gravity values size " << gravity_values.size() << std::endl;
                std::cout << "tmp gravity values size " << tmp_gravity_values.size() << std::endl;
                std::cout << "tempo totale";
                t.log();

                delete output;
            }

            //if(ImGui::Button("compute gravity from octree")) {
            if(!octree.empty()) {
                auto p = debug_ball_position;
                if(!util::is_inside_box(glm::make_vec3(p), util::get_box(glm::make_vec3(min), edge))) {
                    std::cout << "point is outside octree" << std::endl;
                } else {
                    auto current_min = min;
                    auto current_edge = edge;

                    int i = 0;
                    for(int depth = 0;; depth++) {
                        if(octree[i] <= 0) {
                            std::array<glm::vec3, 8> values{};
                            for(int j = 0; j < 8; j++) values[j] = gravity_values[-octree[i + j]];
                            //std::cout << "is p inside testing box: " << util::is_inside_box(glm::make_vec3(p), util::get_box(glm::make_vec3(current_min), current_edge)) << std::endl;
                            //std::cout << "p tested at " << depth << " depth" << std::endl;
                            gravity_from_octree = util::interpolate(
                                    glm::make_vec3(p),
                                    util::get_box(glm::make_vec3(current_min), current_edge),
                                    values
                                    );
                            break;
                            // leaf: from i to i + 8 retreive gravity values and interpolate
                        } else {
                            int k = 0;
                            if(p[0] > current_min[0] + current_edge/2.f) k += 1;
                            if(p[1] > current_min[1] + current_edge/2.f) k += 2;
                            if(p[2] > current_min[2] + current_edge/2.f) k += 4;
                            i = octree[i + k];
                            current_edge /= 2.f;
                            current_min = glm::value_ptr(util::get_box(glm::make_vec3(current_min), current_edge)[k]);
                        }
                    }
                }
            }

            // GRAVITY OUTPUTS
            ImGui::Text("gravity: %f %f %f", gravity.x, gravity.y, gravity.z);
            ImGui::Text("gravity from octree: %f %f %f", gravity_from_octree.x, gravity_from_octree.y, gravity_from_octree.z);
            ImGui::Text("gravity as potential gradient: %f %f %f",
                gravity_from_potential_gradient.x,
                gravity_from_potential_gradient.y,
                gravity_from_potential_gradient.z
                );
            ImGui::Text("gravity force: %f", glm::length(gravity));

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
        ImGui::End();

        // IMGUI TAB FOR MODEL TRANSFORM
        ImGui::Begin("Model Transform");
        ImGui::SliderFloat("model scale", &scale, 0.1f, 5.0f);
        ImGui::SliderFloat3("model position", position, -20.0, 20.0);
        ImGui::SliderFloat3("model rotation", rotation, 0.0f, 2.0f * M_PI);

        ImGui::SliderFloat("debug ball scale", &debug_ball_scale, 0.1f, 1.f);
        ImGui::SliderFloat3("debug ball position", debug_ball_position, -20.0, 20.0);
        ImGui::InputFloat3("debug_ball_position", debug_ball_position);
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
            debug_ball.load_from_obj("../resources/ico.obj");
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

        // SWAP WINDOWS
        SDL_GL_SwapWindow(window);
    }

    SDL_Log("quitting gracefully");
}