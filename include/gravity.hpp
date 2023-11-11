//
// Created by Giovanni Bollati on 27/10/23.
//

#ifndef GL_TEST_PROJECT_GRAVITY_HPP
#define GL_TEST_PROJECT_GRAVITY_HPP

#include <vector>
#include <glm/glm.hpp>
#include <octree.hpp>

namespace gravity {
    // tube -> line segment
    struct tube {
        glm::vec3 t1, t2;
    };
    // point mass
    struct mass {
        glm::vec3 p;
        float m;
    };
    struct ray {
        glm::vec3 origin;
        glm::vec3 dir; // to be normalized if needed
    };
    struct tetrahedron {
        glm::vec3 b1, b2, b3, v;
    };
    struct cube {
        glm::vec3 center;
        float extent;
    };
    struct gravity_cube {
        gravity::cube c;
        std::array<glm::vec3, 8> g;
    };

    std::vector<tube> getTubes(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution);
    std::vector<mass> getMasses(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution);

    glm::vec3 getGravityFromTubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point);
    glm::vec3 getGravityFromMasses(const std::vector<gravity::mass>& masses, float G, glm::vec3 point);

    octree<gravity::gravity_cube>* getGravityOctreeFromMasses(
            glm::vec3 min, glm::vec3 max, int resolution, const std::vector<gravity::mass>& masses);

    // vettore monodimensionale for(x) {for(y) {for(z)}}}
    std::vector<glm::vec3> getDiscreteSpace(glm::vec3 min, glm::vec3 max, int resolution);

    octree<gravity::cube>* getDiscreteSpaceAsOctree(glm::vec3 min, glm::vec3 max, int resolution);

    glm::vec3 getGravityRT(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces, int resolution, glm::vec3 point
    );

    // for each face -> compute tetrahedron signed volume; find barycentre and add force to total
    glm::vec3 getGravityFromTetrahedrons(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& p, const glm::vec3& tetrahedrons_vertex);
    glm::vec3 getGravityFromTetrahedronsCorrected(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& p, const glm::vec3& tetrahedrons_vertex);

    float volume(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& t);
}

#endif //GL_TEST_PROJECT_GRAVITY_HPP
