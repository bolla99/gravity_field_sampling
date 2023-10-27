//
//  Mesh.cpp
//  mesh_software
//
//  Created by Giovanni Bollati on 24/03/23.
//


/* definition needed to activate gl functions
 * declared in sdl_opengl_glext.h,
 * which is included by sdl_opengl.h */
#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include "Mesh.hpp"
#include "util.hpp"

// stdlib
#include <iostream>
#include <sstream>
#include <string>
#include <stack>
#include <fstream>

// parallel processing
#include <omp.h>

// SDL
#import <SDL.h>
#include <SDL_opengl.h>

// GLM
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

void Mesh::loadFromObj(const std::string& path) {
    // reset mesh -> can be done on existing mesh
    clear();
    
    // open obj file as input file stream
    std::ifstream obj;
    obj.open(path, std::ios::in);
    
    // current line
    std::string line;
    // line as string stream
    std::stringstream stream;
    // element of the line
    std::string element;
    
    // parsing
    while(obj.good()) {
        // reset line
        line = "";
        // reset element
        element = "";
        // get next line
        std::getline(obj, line);
        stream = std::stringstream(line);
        // get first element of line
        stream >> element;
        // discard comments
        if(element == "#") {
            continue;
        }
        // check for vertex
        if(element == "v") {
            glm::vec3 v = {0.f, 0.f, 0.f};
            stream >> v.x >> v.y >> v.z;
            vertices.push_back(v);
            if(stream.good()) {
                Color c = {0.f, 0.f, 0.f};
                stream >> c.r >> c.g >> c.b;
                colors.push_back(c);
            }
        }
        if(element == "vt") {
            UV uv = {0.f, 0.f};
            stream >> uv.x >> uv.y;
            UVs.push_back(uv);
            // check for face
        }
        if(element == "vn") {
            glm::vec3 normal = {0.f, 0.f, 0.f};
            stream >> normal.x >> normal.y >> normal.z;
            normals.push_back(normal);
            // check for face
        }
        if(element == "f") {
            std::vector<std::string> indices;
            while(stream.good()) {
                std::string s;
                if(stream >> s) indices.push_back(s);
            }
            addIndices(indices);
        }
    }
    // close obj file
    obj.close();
    std::cout << "mesh at path " << path << " successfully completed" << std::endl;
    setUpForRendering();
}

// TO BE ADAPTER TO NEW LOADING FEATURES ( UV )
void Mesh::writeOnDiskAsObj(const std::string& path) {
    std::ofstream file(path, std::ios::trunc);
    file << toString();
    file.close();
}

// TO BE ADAPTED
// TO STRING
[[maybe_unused]] std::string Mesh::toString() {
    std::stringstream stream;
    for(int i = 0; i < vertices.size(); i++) {
        stream << "v " + std::to_string(vertices[i].x) + " "
                + std::to_string(vertices[i].y) + " "
                + std::to_string(vertices[i].z);
        if(hasColor())
            stream << " " + std::to_string(colors[i].r)
                    + " " + std::to_string(colors[i].g)
                    + " " + std::to_string(colors[i].b);
        stream << std::endl;
    }
    if(hasUVs()) {
        for(auto & UV : UVs) {
            stream << "vt " + std::to_string(UV.x) + " " + std::to_string(UV.y) << std::endl;
        }
    }
    if(hasNormals()) {
        for(auto & normal : normals) {
            stream << "vn " + std::to_string(normal.x)
                    + " " + std::to_string(normal.y)
                    + " " + std::to_string(normal.z) << std::endl;
        }
    }
    for(int i = 0; i < faces.size(); i++) {
        stream << "f ";
        std::string s;
        s = "";
        s += std::to_string(faces[i].x);
        if(hasUVs()) s += "/" + std::to_string(faces_uv[i].x);
        if(hasNormals()) {
            s += "/"; if(!hasUVs()) s += "/";
            s += std::to_string(faces_normals[i].x);
        }
        s += " ";
        stream << s;

        s = "";
        s += std::to_string(faces[i].y);
        if(hasUVs()) s += "/" + std::to_string(faces_uv[i].y);
        if(hasNormals()) {
            s += "/"; if(!hasUVs()) s += "/";
            s += std::to_string(faces_normals[i].y);
        }
        s += " ";
        stream << s;

        s = "";
        s += std::to_string(faces[i].z);
        if(hasUVs()) s += "/" + std::to_string(faces_uv[i].z);
        if(hasNormals()) {
            s += "/"; if(!hasUVs()) s += "/";
            s += std::to_string(faces_normals[i].z);
        }
        s += " ";
        stream << s;
        stream << std::endl;
    }
    return stream.str();
}

void Mesh::setUpForRendering() {

    int j = 0;
    for(int i = 0; i < faces.size(); i++) {
        elements.push_back(j);
        attributes.push_back(vertices[faces[i].x - 1].x);
        attributes.push_back(vertices[faces[i].x - 1].y);
        attributes.push_back(vertices[faces[i].x - 1].z);
        if(hasColor()) {
            attributes.push_back(colors[faces[i].x - 1].r);
            attributes.push_back(colors[faces[i].x - 1].g);
            attributes.push_back(colors[faces[i].x - 1].b);
        }
        if(hasUVs()) {
            attributes.push_back(UVs[faces_uv[i].x - 1].x);
            attributes.push_back(UVs[faces_uv[i].x - 1].y);
        }
        if(hasNormals()) {
            attributes.push_back(normals[faces_normals[i].x - 1].x);
            attributes.push_back(normals[faces_normals[i].x - 1].y);
            attributes.push_back(normals[faces_normals[i].x - 1].z);
        }

        elements.push_back(j+1);
        attributes.push_back(vertices[faces[i].y - 1].x);
        attributes.push_back(vertices[faces[i].y - 1].y);
        attributes.push_back(vertices[faces[i].y - 1].z);
        if(hasColor()) {
            attributes.push_back(colors[faces[i].y - 1].r);
            attributes.push_back(colors[faces[i].y - 1].g);
            attributes.push_back(colors[faces[i].y - 1].b);
        }
        if(hasUVs()) {
            attributes.push_back(UVs[faces_uv[i].y - 1].x);
            attributes.push_back(UVs[faces_uv[i].y - 1].y);
        }
        if(hasNormals()) {
            attributes.push_back(normals[faces_normals[i].y - 1].x);
            attributes.push_back(normals[faces_normals[i].y - 1].y);
            attributes.push_back(normals[faces_normals[i].y - 1].z);
        }

        elements.push_back(j+2);
        attributes.push_back(vertices[faces[i].z - 1].x);
        attributes.push_back(vertices[faces[i].z - 1].y);
        attributes.push_back(vertices[faces[i].z - 1].z);
        if(hasColor()) {
            attributes.push_back(colors[faces[i].z - 1].r);
            attributes.push_back(colors[faces[i].z - 1].g);
            attributes.push_back(colors[faces[i].z - 1].b);
        }
        if(hasUVs()) {
            attributes.push_back(UVs[faces_uv[i].z - 1].x);
            attributes.push_back(UVs[faces_uv[i].z - 1].y);
        }
        if(hasNormals()) {
            attributes.push_back(normals[faces_normals[i].z - 1].x);
            attributes.push_back(normals[faces_normals[i].z - 1].y);
            attributes.push_back(normals[faces_normals[i].z - 1].z);
        }
        j += 3;
    }
    std::cout << "mesh rendering setup successfully complete" << std::endl;
}

unsigned int Mesh::getVAO() {
    if(VAO == 0) {
        // DECLARE BUFFERS
        glGenVertexArrays(1, &VAO);
        // GEN BUFFERS
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
        glBindVertexArray(VAO);
        // BIND BUFFERS
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        // BIND DATA
        glBufferData(GL_ARRAY_BUFFER, attributes.size() * sizeof(float), &(attributes.front()), GL_DYNAMIC_DRAW);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, elements.size() * sizeof(int), &(elements.front()), GL_DYNAMIC_DRAW);

        int stride = 3 * sizeof(float);
        int uvs_pointer = 3 * sizeof(float);
        int normals_pointer = 3 * sizeof(float);
        if(!colors.empty()) {
            stride += 3 * sizeof(float);
            uvs_pointer += 3 * sizeof(float);
            normals_pointer += 3 * sizeof(float);
        }
        if(!UVs.empty()) {
            stride += 2 * sizeof(float);
            normals_pointer += 2 * sizeof(float);
        }
        if(!normals.empty()) stride += 3 * sizeof(float);
        // SET POSITION ATTRIBUTE
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, nullptr);
        glEnableVertexAttribArray(0);
        // SET COLOR ATTRIBUTE ON LOCATION 1
        if(!colors.empty()) {
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(1);
        }
        // SET COLOR ATTRIBUTE ON LOCATION 2
        if(!UVs.empty()) {
            glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride, (void *)(size_t)uvs_pointer);
            glEnableVertexAttribArray(2);
        }
        if(!normals.empty()) {
            glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, stride, (void *)(size_t)normals_pointer);
            glEnableVertexAttribArray(3);
        }
    }
    return VAO;
}

void Mesh::addIndices(std::vector<std::string>& indices) {
    for(int i = 0; i < indices.size() - 2; i++) {
        // triangle
        Face face{};
        Face face_uv{};
        Face face_normals{};
        size_t slash;
        size_t slash2;

        slash = indices[0].find('/');
        face.x = std::stoi(indices[0].substr(0, slash), nullptr);
        slash2 = indices[0].find('/', slash + 1);
        face_uv.x = std::stoi(indices[0].substr(slash + 1, slash2 - slash - 1), nullptr);
        face_normals.x = std::stoi(indices[0].substr( slash2 + 1), nullptr);

        slash = indices[i + 1].find('/');
        face.y = std::stoi(indices[i + 1].substr(0, slash), nullptr);
        slash2 = indices[i + 1].find('/', slash + 1);
        face_uv.y = std::stoi(indices[i + 1].substr(slash + 1, slash2 - slash - 1), nullptr);
        face_normals.y = std::stoi(indices[i + 1].substr( slash2 + 1), nullptr);

        slash = indices[i + 2].find('/');
        face.z = std::stoi(indices[i + 2].substr(0, slash), nullptr);
        slash2 = indices[i + 2].find('/', slash + 1);
        face_uv.z = std::stoi(indices[i + 2].substr(slash + 1, slash2 - slash - 1), nullptr);
        face_normals.z = std::stoi(indices[i + 2].substr( slash2 + 1), nullptr);

        // uvs and normals checkers can be called since faces parsing comes after vertex parsing,
        // that's why uvs and normals vector are already filled
        faces.push_back(face);
        if(hasUVs()) faces_uv.push_back(face_uv);
        if(hasNormals()) faces_normals.push_back(face_normals);
    }
}

// GETTERS
const std::vector<glm::vec3>& Mesh::getVertices() const {
    return vertices;
}
const std::vector<Color>& Mesh::getColors() const {
    return colors;
}
const std::vector<UV>& Mesh::getUVs() const {
    return UVs;
}
const std::vector<glm::vec3>& Mesh::getNormals() const {
    return normals;
}
const std::vector<Face>& Mesh::getFaces() const {
    return faces;
}
const std::vector<Face>& Mesh::getFacesUV() const {
    return faces_uv;
}
const std::vector<Face>& Mesh::getFacesNormals() const {
    return faces_normals;
}
const std::vector<float>& Mesh::getAttributes() const {
    return attributes;
}
const std::vector<int>& Mesh::getElements() const {
    return elements;
}

void Mesh::clear() {
    vertices.clear();
    colors.clear();
    UVs.clear();
    normals.clear();

    faces.clear();
    faces_uv.clear();
    faces_normals.clear();

    attributes.clear();
    elements.clear();

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    VAO = 0;
    VBO = 0;
    EBO = 0;
}

bool Mesh::isLoaded() const { return !vertices.empty(); }
bool Mesh::hasColor() const { return !colors.empty(); }
bool Mesh::hasUVs() const { return !UVs.empty(); }
bool Mesh::hasNormals() const { return !normals.empty(); }

// IS INSIDE CHECKER + HELPERS
bool Mesh::isInside(const glm::vec3& v) const {
    double start = omp_get_wtime();
    float minDistance = util::pointTriangleDistance(
            v,
            vertices[faces[0].x],
            vertices[faces[0].y],
            vertices[faces[0].z]
            );
#pragma omp parallel for default(none) shared(minDistance, v)
    for(int i = 1; i < (int)faces.size(); i++) {
        float distance = util::pointTriangleDistance(
                v,
                vertices[faces[i].x],
                vertices[faces[i].y],
                vertices[faces[i].z]
        );
#pragma omp critical
        {
            if (abs(distance) < abs(minDistance)) {
                minDistance = distance;
            }
        }
    }
    std::cout << "time elapsed: " << (omp_get_wtime() - start) << "\n";
    return minDistance < 0;
}

glm::vec3 Mesh::getMin() const {
    glm::vec3 min = vertices[0];
    for(int i = 1; i < vertices.size(); i++) {
        if(vertices[i].x < min.x) min.x = vertices[i].x;
        if(vertices[i].y < min.y) min.y = vertices[i].y;
        if(vertices[i].z < min.z) min.z = vertices[i].z;
    }
    return min;
}
glm::vec3 Mesh::getMax() const {
    glm::vec3 max = vertices[0];
    for(int i = 1; i < vertices.size(); i++) {
        if(vertices[i].x > max.x) max.x = vertices[i].x;
        if(vertices[i].y > max.y) max.y = vertices[i].y;
        if(vertices[i].z > max.z) max.z = vertices[i].z;
    }
    return max;
}

glm::vec3 Mesh::getGravityFromTetrahedrons(const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) const {
    glm::vec3 gravity = {0.f, 0.f, 0.f};
    for(auto & face : faces) {
        tetrahedron t = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                tetrahedrons_vertex};
        float mass = util::tetrahedronVolume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedronBarycentre(t.b1, t.b2, t.b3, t.v);
        float distance = glm::length(p - barycentre);
        gravity = gravity - ((p - barycentre)*mass) / (float)pow(distance, 3);
    }
    return gravity;
}

glm::vec3 Mesh::getGravityFromTetrahedronsCorrected(const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) const {
    glm::vec3 gravity = {0.f, 0.f, 0.f};
    for(auto & face : faces) {
        tetrahedron t = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                tetrahedrons_vertex};
        float mass = util::tetrahedronVolume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedronBarycentre(t.b1, t.b2, t.b3, t.v);
        gravity = gravity + ((barycentre - p)*mass) / (float)pow(glm::length(barycentre - p), 3);
    }
    return gravity;
}

float Mesh::volume(const glm::vec3& t) const {
    float volume = 0.f;
    for(auto & face : faces) {
        tetrahedron tetrahedron = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                t};
        volume += util::tetrahedronVolume(tetrahedron.b1, tetrahedron.b2, tetrahedron.b3, tetrahedron.v) * 10;
    }
    return volume;
}

glm::vec3 Mesh::getGravityRT(int resolution, glm::vec3 point) const {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }

    Uint64 start = SDL_GetTicks64();
    // FIND XY PLANE
    glm::vec3 min = getMin();
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = getMax();

    max.x -= 0.0001; max.y -= 0.0001;
    min.z = min.z - 10.f;
    max.z = min.z;

    glm::vec3 center = (min + max) / 2.0f;

    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float max_extent = x_width;
    if(y_width > max_extent) max_extent = y_width;

    min = {center.x - (max_extent / 2.0f), center.y - (max_extent / 2.0f), center.z};

    // CUBE EDGE LENGTH
    float cube_edge = max_extent / (float)resolution;
    auto cube_volume = (float)pow(cube_edge, 3);
    float G = 10.f;

    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(resolution, cube_edge, ray_dir, min, cube_volume, point, thread_gravity, G)
    for(int i = 0; i < resolution + 1; i++) {
        // private vars for thread
        float units_per_tubef, tube_unit;
        float new_cube_volume, r3;
        int units_per_tube;
        tube t{};

        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = rayMeshIntersectionsOptimized(r);

            // FOR EACH INTERSECTIONS COUPLE
            for(int k = 0; k < intersections.size(); k += 2) {
                // TUBE
                t = {intersections[k], intersections[k + 1]};

                // TUBE LENGTH / UNIT
                units_per_tubef = glm::length(t.t2 - t.t1) / cube_edge;
                // ROUNDED
                units_per_tube = (int)std::round((double)units_per_tubef);
                // NEW CUBE MASS
                new_cube_volume = (cube_volume * units_per_tubef) / (float)units_per_tube;

                if(units_per_tube == 0) continue;
                // TUBE VERTEX OFFSET
                tube_unit = glm::length(t.t2 - t.t1) / (float)units_per_tube;
                for(int n = 0; n < units_per_tube + 1; n++) {
                    glm::vec3 mass = t.t1 + glm::normalize(t.t2 - t.t1) * (float)n * tube_unit;
                    glm::vec3 dir = mass - point;
                    r3 = (float)pow(glm::length(dir), 3);
                    if(r3 > -std::numeric_limits<float>::epsilon() && r3 < std::numeric_limits<float>::epsilon()) {
                        continue;
                    }
                    thread_gravity[omp_get_thread_num()] = thread_gravity[omp_get_thread_num()] + (dir * new_cube_volume * G) / r3;
                }
            }
        }
    }
    std::cout << "time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    glm::vec3 gravity = {0, 0, 0};
    for(int i = 0; i < omp_get_max_threads(); i++) {
        gravity = gravity + thread_gravity[i];
    }
    return gravity;
}
glm::vec3 Mesh::getGravityFromTubes(int resolution, const std::vector<tube>& tubes, glm::vec3 point) const {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }

    Uint64 start = SDL_GetTicks64();
    // FIND XY PLANE
    glm::vec3 min = getMin();
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = getMax();

    max.x -= 0.0001; max.y -= 0.0001;
    min.z = min.z - 10;
    max.z = min.z;

    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float max_extent = x_width;
    if(y_width > max_extent) max_extent = y_width;

    // CUBE EDGE LENGTH
    float cube_edge = max_extent / (float)resolution;
    auto cube_volume = (float)pow(max_extent / (float)resolution, 3);
    float G = 10.0;

    glm::vec3 ray_dir = {0.f, 0.f, 1.f};
#pragma omp parallel default(none), shared(tubes, cube_edge, cube_volume, point, thread_gravity, G)
    {
        float units_per_tubef, tube_unit;
        float new_cube_volume, r3;
        int units_per_tube;
#pragma omp for
        for (auto &t: tubes) {
            // TUBE LENGTH / UNIT
            units_per_tubef = glm::length(t.t2 - t.t1) / cube_edge;
            // ROUNDED
            units_per_tube = (int) std::round((double) units_per_tubef);
            // NEW CUBE MASS
            new_cube_volume = (cube_volume * units_per_tubef) / (float)units_per_tube;

            if (units_per_tube == 0) continue;
            // TUBE VERTEX OFFSET
            tube_unit = glm::length(t.t2 - t.t1) / (float) units_per_tube;
            for (int n = 0; n < units_per_tube + 1; n++) {
                glm::vec3 mass = t.t1 + glm::normalize(t.t2 - t.t1) * (float) n * tube_unit;
                glm::vec3 dir = mass - point;
                r3 = (float)pow(glm::length(dir), 3);
                if (r3 > -std::numeric_limits<float>::epsilon() && r3 < std::numeric_limits<float>::epsilon()) {
                    continue;
                }
                thread_gravity[omp_get_thread_num()] =
                        thread_gravity[omp_get_thread_num()] + (dir * new_cube_volume * G) / r3;
            }
        }
    }
    std::cout << "gravity tubes time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    glm::vec3 gravity = {0, 0, 0};
    for(int i = 0; i < omp_get_max_threads(); i++) {
        gravity = gravity + thread_gravity[i];
    }
    return gravity;
}
glm::vec3 Mesh::getGravityFromMasses(const std::vector<mass>& masses, float G, glm::vec3 point) {
    Uint64 start = SDL_GetTicks64();
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }
#pragma omp parallel for default(none) shared(masses, point, G, thread_gravity)
    for(auto & mass : masses) {
        glm::vec3 dir = mass.p - point;
        auto r3 = (float)pow(glm::length(dir), 3);
        if (r3 > -std::numeric_limits<float>::epsilon() && r3 < std::numeric_limits<float>::epsilon()) {
            continue;
        }
        thread_gravity[omp_get_thread_num()] = thread_gravity[omp_get_thread_num()] + (dir * mass.m * G) / r3;
    }
    glm::vec3 gravity{0, 0, 0};
    for(int i = 0; i < omp_get_max_threads(); i++) {
        gravity = gravity + thread_gravity[i];
    }
    std::cout << "gravity masses time elapsed: " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    return gravity;
}



std::vector<tube> Mesh::getTubes(int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = getMin();
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = getMax();

    max.x -= 0.0001; max.y -= 0.0001;
    min.z = min.z - 10;
    max.z = min.z;

    glm::vec3 center = (min + max) / 2.0f;

    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float max_extent = x_width;
    if(y_width > max_extent) max_extent = y_width;

    min = {center.x - (max_extent / 2.0f), center.y - (max_extent / 2.0f), center.z};

    // CUBE EDGE LENGTH
    float cube_edge = max_extent / (float)resolution;
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(cube_edge, resolution, ray_dir, tubes, min)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = rayMeshIntersectionsOptimized(r);

            // FOR EACH INTERSECTIONS COUPLE
            for(int k = 0; k < intersections.size(); k += 2) {
                // TUBE
#pragma omp critical
                tubes.push_back({intersections[k], intersections[k + 1]});
            }
        }
    }
    std::cout << "get tubes time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    return tubes;
}


std::vector<mass> Mesh::getMasses(int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};
    std::vector<mass> volumes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = getMin();
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = getMax();

    max.x -= 0.0001; max.y -= 0.0001;
    min.z = min.z - 10;
    max.z = min.z;

    glm::vec3 center = (min + max) / 2.0f;

    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float max_extent = x_width;
    if(y_width > max_extent) max_extent = y_width;

    min = {center.x - (max_extent / 2.0f), center.y - (max_extent / 2.0f), center.z};

    // CUBE EDGE LENGTH
    float cube_edge = max_extent / (float)resolution;
    auto cube_volume = (float)std::pow(cube_edge, 3);
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(resolution, cube_edge, cube_volume, ray_dir, tubes, min, volumes)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = rayMeshIntersectionsOptimized(r);

            // FOR EACH INTERSECTIONS COUPLE
            for(int k = 0; k < intersections.size(); k += 2) {
                // TUBE
                tube t = {intersections[k], intersections[k + 1]};
                float units_per_tubef, new_cube_volume;
                int units_per_tube;
                // TUBE LENGTH / UNIT
                units_per_tubef = glm::length(t.t2 - t.t1) / cube_edge;
                // ROUNDED
                units_per_tube = (int) std::round((double) units_per_tubef);
                // NEW CUBE MASS
                new_cube_volume = (cube_volume * units_per_tubef) / (float) units_per_tube;

                if (units_per_tube == 0) continue;
                // TUBE VERTEX OFFSET
                float tube_unit = glm::length(t.t2 - t.t1) / (float) units_per_tube;
                for (int n = 0; n < units_per_tube + 1; n++) {
#pragma omp critical
                    volumes.push_back({t.t1 + glm::normalize(t.t2 - t.t1) * (float) n * tube_unit, new_cube_volume});
                }
            }
        }
    }
    std::cout << "get masses time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    return volumes;
}
std::vector<glm::vec3> Mesh::rayMeshIntersections(ray r) const {
    // parameter: intersection = origin + direction * parameter
    std::vector<glm::vec3> intersections = {};
    std::vector<float> parameters = {};
    float parameter;
    for(auto & face : faces) {
        if(util::rayTriangleIntersection(r.origin, r.dir, vertices[face.x - 1], vertices[face.y - 1], vertices[face.z - 1], &parameter)) {
            parameters.push_back(parameter);
        }
    }
    if(parameters.empty()) return intersections;
    // sort to return intersection in order of distance from ray origin
    std::sort(parameters.begin(), parameters.end());

    // discard unique values in order to avoid double intersection through edges
    auto new_last = std::unique(parameters.begin(), parameters.end(),
                [](float f1, float f2) { return std::abs(f1 - f2) < std::numeric_limits<float>::epsilon(); });
    parameters.resize(std::distance(parameters.begin(), new_last));
    intersections.reserve(parameters.size());
    for(float p : parameters) {
        intersections.push_back(r.origin + r.dir * p);
    }
    return intersections;
}
std::vector<glm::vec3> Mesh::rayMeshIntersectionsOptimized(ray r) const {
    // parameter: intersection = origin + direction * parameter
    std::vector<glm::vec3> intersections = {};
    std::vector<float> parameters = {};
    float parameter;
    for(auto & face : faces) {
        glm::vec3 t0 = vertices[face.x - 1];
        glm::vec3 t1 = vertices[face.y - 1];
        glm::vec3 t2 = vertices[face.z - 1];
        glm::vec3 o = r.origin;
        // CUSTOM RAY MESH INTERSECTION: since ray is parallel to z axis, its possible
        // to do a prececk and exclude triangles that are out of ray vision;

        if(t0.x > o.x && t1.x > o.x && t2.x > o.x) continue;
        if(t0.x < o.x && t1.x < o.x && t2.x < o.x) continue;
        if(t0.y > o.y && t1.y > o.y && t2.y > o.y) continue;
        if(t0.y < o.y && t1.y < o.y && t2.y < o.y) continue;
        float Do = -o.x + o.y;
        if(t0.y < t0.x + Do && t1.y < t1.x + Do && t2.y < t2.x + Do) continue;
        if(t0.y > t0.x + Do && t1.y > t1.x + Do && t2.y > t2.x + Do) continue;

        if(util::rayTriangleIntersection(r.origin, r.dir, t0, t1, t2, &parameter)) {
            parameters.push_back(parameter);
        }
    }
    if(parameters.empty()) return intersections;
    // sort to return intersection in order of distance from ray origin
    std::sort(parameters.begin(), parameters.end());

    // discard unique values in order to avoid double intersection through edges
    auto new_last = std::unique(parameters.begin(), parameters.end(),
                                [](float f1, float f2) { return std::abs(f1 - f2) < std::numeric_limits<float>::epsilon(); });
    parameters.resize(std::distance(parameters.begin(), new_last));
    intersections.reserve(parameters.size());
    for(float p : parameters) {
        intersections.push_back(r.origin + r.dir * p);
    }
    return intersections;
}

std::vector<glm::vec3> Mesh::getDiscreteSpace(int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 min = getMin();
    glm::vec3 max = getMax();
    glm::vec3 center = (max + min) * 0.5f;
    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float z_width = max.z - min.z;
    float max_extent = x_width;
    if (y_width > max_extent) max_extent = y_width;
    if (z_width > max_extent) max_extent = z_width;

    min = center - (max_extent / 2.0f);
    max = center + (max_extent / 2.0f);

    float unit = max_extent / (float)resolution;

    std::vector<glm::vec3> v;
    for (int i = 0; i < resolution + 1; i++) {
        for (int j = 0; j < resolution + 1; j++) {
            for (int k = 0; k < resolution + 1; k++) {
                v.emplace_back(min.x + (float)i * unit, min.y + (float)j * unit, min.z + (float)k * unit);
            }
        }
    }
    return v;
}
