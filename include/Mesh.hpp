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


#ifndef Mesh_hpp
#define Mesh_hpp

/* definition needed to activate gl functions
 * declared in sdl_opengl_glext.h,
 * which is included by sdl_opengl.h */
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

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

struct Color {
    float r, g, b;
};
struct UV {
    float x, y;
};
struct tetrahedron {
    glm::vec3 b1, b2, b3, v;
};
struct ray {
    glm::vec3 origin;
    glm::vec3 dir; // to be normalized if needed
};
// tube -> line segment
struct tube {
    glm::vec3 t1, t2;
};
// point mass
struct mass {
    glm::vec3 p;
    float m;
};

class Mesh {

public:
    Mesh() : vertices(), colors(), UVs(), normals(), faces(), faces_uv(), faces_normals(), attributes(), elements() {}

    // delete VAO, VBO and EBO if they have been created after object creation
    // no more cleanup needed
    ~Mesh() {
        if(VAO != 0) glDeleteVertexArrays(1, &VAO);
        if(VBO != 0) glDeleteBuffers(1, &VBO);
        if(EBO != 0) glDeleteBuffers(1, &EBO);
    }

    // IO
    // fill every data structures except for rendering wise vectors
    void loadFromObj(const std::string& path);

    // to be modified
    void writeOnDiskAsObj(const std::string& path);

    // requires a shader with location 0 for position, 1 for color and 2 for uvs;
    unsigned int getVAO();


    // CHECKERS
    [[nodiscard]] bool isLoaded() const;
    [[nodiscard]] bool hasColor() const;
    [[nodiscard]] bool hasUVs() const;
    [[nodiscard]] bool hasNormals() const;

    [[nodiscard]] bool isInside(const glm::vec3& v) const;

    // GETTERS
    [[nodiscard]] const std::vector<glm::vec3>& getVertices() const;
    [[nodiscard]] const std::vector<Color>& getColors() const;
    [[nodiscard]] const std::vector<UV>& getUVs() const;
    [[nodiscard]] const std::vector<glm::vec3>& getNormals() const;

    [[nodiscard]] const std::vector<Face>& getFaces() const;
    [[nodiscard]] const std::vector<Face>& getFacesUV() const;
    [[nodiscard]] const std::vector<Face>& getFacesNormals() const;

    [[nodiscard]] const std::vector<float>& getAttributes() const;
    [[nodiscard]] const std::vector<int>& getElements() const;

    // for each face -> compute tetrahedron signed volume; find barycentre and add force to total
    [[deprecated]][[nodiscard]] glm::vec3 getGravityFromTetrahedrons(const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) const;
    [[deprecated]][[nodiscard]] glm::vec3 getGravityFromTetrahedronsCorrected(const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) const;

    [[nodiscard]] float volume(const glm::vec3& t) const;

    // ABSTRACTION FUNCTION
    std::string toString();

    // these are the function to read and modify
    [[nodiscard]] glm::vec3 getGravityRT(int resolution, glm::vec3 point) const;
    [[nodiscard]] glm::vec3 getGravityFromTubes(int resolution, const std::vector<tube>& tubes, glm::vec3 point) const;
    static glm::vec3 getGravityFromMasses(const std::vector<mass>& masses, float G, glm::vec3 point);

    [[nodiscard]] std::vector<tube> getTubes(int resolution) const;
    [[nodiscard]] std::vector<mass> getMasses(int resolution) const;

    // da capire
    std::vector<glm::vec3> getDiscreteSpace(int resolution);
    
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
    void addIndices(std::vector<std::string>& indices);

    // fill attributes and elements data structure for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void setUpForRendering();

    // CLEAR: clears every and vector and delete buffers for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void clear();
    // *************************+*************************************

    // returns a point which has the minimum coordinate along all the three axis
    [[nodiscard]] glm::vec3 getMin() const;
    // returns a point which has the maximum coordinate along the three axis
    [[nodiscard]] glm::vec3 getMax() const;
};
#endif /* Mesh_hpp */
