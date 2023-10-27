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