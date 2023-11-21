//
//  Mesh.hpp
//  mesh_software
//
//  Created by Giovanni Bollati on 24/03/23.
//

// ********************************************************** //
// ********************************************************** //
// ********************************************************** //
// ********* WARNING WARNING WARNING WARNING WARNING ******** //
// ***** APP CRASHED WHILE LOADING MESH WITHOUT UV MAPPING ** //
// ********************************************************** //
// ********************************************************** //
// ********************************************************** //


#ifndef mesh_hpp
#define mesh_hpp

/* definition needed to activate gl functions
 * declared in sdl_opengl_glext.h,
 * which is included by sdl_opengl.h */
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include <gravity.hpp>

// stdlib
#include <cstdio>
#include <vector>
#include <string>

#ifdef WIN32
#include <GL/glew.h>
#endif

// SDL / GLM
#include <SDL_opengl.h>
#include <glm/glm.hpp>

// HELPER STRUCTS
using Face = glm::vec<3, unsigned int>;
using UV = glm::vec2;
struct Color {
    float r, g, b;
};

class mesh {

public:
    mesh() : vertices(), colors(), UVs(), normals(), faces(), faces_uv(), faces_normals(), attributes(), elements() {}

    // delete VAO, VBO and EBO if they have been created after object creation
    // no more cleanup needed
    ~mesh() {
        if(VAO != 0) glDeleteVertexArrays(1, &VAO);
        if(VBO != 0) glDeleteBuffers(1, &VBO);
        if(EBO != 0) glDeleteBuffers(1, &EBO);
    }

    // IO
    // fill every data structures except for rendering wise vectors
    void load_from_obj(const std::string& path);

    // to be modified
    void write_on_disk_as_obj(const std::string& path);

    // requires a shader with location 0 for position, 1 for color and 2 for uvs;
    unsigned int get_VAO();


    // CHECKERS
    [[nodiscard]] bool is_loaded() const;
    [[nodiscard]] bool has_color() const;
    [[nodiscard]] bool has_UVs() const;
    [[nodiscard]] bool has_normals() const;

    [[nodiscard]] bool is_inside(const glm::vec3& v) const;

    // GETTERS
    [[nodiscard]] const std::vector<glm::vec3>& get_vertices() const;
    [[nodiscard]] const std::vector<Color>& get_colors() const;
    [[nodiscard]] const std::vector<UV>& get_UVs() const;
    [[nodiscard]] const std::vector<glm::vec3>& get_normals() const;

    [[nodiscard]] const std::vector<Face>& get_faces() const;
    [[nodiscard]] const std::vector<Face>& get_faces_UV() const;
    [[nodiscard]] const std::vector<Face>& get_faces_normals() const;

    [[nodiscard]] const std::vector<float>& get_attributes() const;
    [[nodiscard]] const std::vector<int>& get_elements() const;

    // ABSTRACTION FUNCTION
    std::string to_string();
    
private:
    // REPRESENTATION
    // ATTRIBUTES
    std::vector<glm::vec3> vertices;
    std::vector<Color> colors;
    std::vector<UV> UVs; // uv values
    std::vector<glm::vec3> normals;

    // INDEXES
    std::vector<Face> faces; // position
    std::vector<Face> faces_uv; // uvs
    std::vector<Face> faces_normals; // normals

    // elements and attributes for rendering; MUST BE INITIALIZED WITH setUpWithRendering, which is called by
    // public api load function.
    std::vector<float> attributes; // position, color, UVs
    std::vector<int> elements;
    // buffers for opengl rendering (they are deleted by decostructor)
    unsigned int VAO = 0;
    unsigned int VBO = 0;
    unsigned int EBO = 0;

    // ********* HELPER FUNCTIONS CALLED BY LOAD FUNCTION ********
    // receives and parses a vector of string which represent a face in obj format, which is x/y/z;
    // while parsing it makes triangular any non-triangular face by fan method; it fills faces vectors
    // addIndices can be called only by load function in order to fill faces vectors with a single face line parsed.
    // CALLED BY LOAD FUNCTION; WHICH IS PUBLIC API
    void add_indices(std::vector<std::string>& indices);

    // fill attributes and elements data structure for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void set_up_for_rendering();

    // CLEAR: clears every and vector and delete buffers for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void clear();
    // *************************+*************************************
};
#endif /* Mesh_hpp */
