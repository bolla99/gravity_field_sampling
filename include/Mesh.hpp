//
//  Mesh.hpp
//  mesh_software
//
//  Created by Giovanni Bollati on 24/03/23.
//

#ifndef Mesh_hpp
#define Mesh_hpp


#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include <cstdio>
#include <vector>
#include <string>

#ifdef WIN32
#include <GL/glew.h>
#endif

#include <SDL_opengl.h>
#include <glm/glm.hpp>

struct Face {
    unsigned int x, y, z;
};
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
    glm::vec3 dir;
};
struct tube {
    glm::vec3 t1, t2;
};
struct mass {
    glm::vec3 p;
    float m;
};

class Mesh {

public:
    Mesh() {
        vertices = std::vector<glm::vec3>();
        colors = std::vector<Color>();
        UVs = std::vector<UV>();
        normals = std::vector<glm::vec3>();

        faces = std::vector<Face>();
        faces_uv = std::vector<Face>();
        faces_normals = std::vector<Face>();

        attributes = std::vector<float>();
        elements = std::vector<int>();
    }

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
    [[nodiscard]] glm::vec3 getGravityFromTetrahedrons(const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) const;
    [[nodiscard]] float volume(const glm::vec3& t) const;

    // ABSTRACTION FUNCTION
    std::string toString();

    // GEOMETRY
    static float pointEdgeDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2);
    static float pointTriangleDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2, const glm::vec3& e3);
    static float tetrahedronVolume(tetrahedron t);
    static glm::vec3 tetrahedronBarycentre(tetrahedron t);

    glm::vec3 getGravityRT(int resolution, glm::vec3 point);
    [[nodiscard]] glm::vec3 getGravityFromTubes(int resolution, const std::vector<tube>& tubes, glm::vec3 point) const;
    static glm::vec3 getGravityFromMasses(const std::vector<mass>& masses, float G, glm::vec3 point);
    std::vector<tube> getTubes(int resolution);
    std::vector<mass> getMasses(int resolution);

    [[maybe_unused]] std::vector<glm::vec3> rayMeshIntersections(ray r);

    std::vector<glm::vec3> rayMeshIntersectionsOptimized(ray r);
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
    std::vector<Face> faces_normals;

    // elements and attributes for rendering; MUST BE INITIALIZED WITH setUpWithRendering
    std::vector<float> attributes; // position, color, UVs
    std::vector<int> elements;
    unsigned int VAO = 0;
    unsigned int VBO = 0;
    unsigned int EBO = 0;

    // HELPER FUNCTIONS FOR LOADING FUNCTION

    // receives and parses a vector of string which represent a face in obj format, which is x/y/z;
    // while parsing it makes triangular any non-triangular face by fan method;
    // CALLED BY LOAD FUNCTION; WHICH IS PUBLIC API
    void addIndices(std::vector<std::string>& indices);

    // fill attributes and elements data structure for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void setUpForRendering();

    // CLEAR: clears every and vector and delete buffers for rendering; CALLED BY LOAD FUNCTION, WHICH IS PUBLIC API
    void clear();

    [[nodiscard]] glm::vec3 getMin() const;
    [[nodiscard]] glm::vec3 getMax() const;

    // RAYCAST GRAVITATIONAL FIELD CALCULATION
    static bool rayTriangleIntersection(ray r, glm::vec3 t1, glm::vec3 t2, glm::vec3 t3, float* parameter);
};
#endif /* Mesh_hpp */
