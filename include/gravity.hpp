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

    std::vector<tube> get_tubes(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution);
    std::vector<mass> get_masses(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution);

    glm::vec3 get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point);
    glm::vec3 get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, glm::vec3 point);

    // get gravity given 3d space and gravity vector (with related min vector, range and resolution) and point
    // space and gravity vector are meant to be precomputed during a non real-time phase, while
    // this function is meant to be called in real-time by interpolating precomputed values and compute
    // gravity for a given point;
    // min is the min space vertex of the box of which gravity has been previoudly computed
    // min, range and resolution are needed to find indices of the bounding box, then
    // the spatial bounding box coordinates are retrived and used to interpolate the gravity values of each
    // bounding box vertex;
    glm::vec3 get_gravity_from_1D_precomputed_vector(glm::vec3 point, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, glm::vec3 min, float range, int resolution);

    octree<gravity::gravity_cube>* get_gravity_octree_from_masses(
            glm::vec3 min, glm::vec3 max, int resolution, const std::vector<gravity::mass>& masses);

    // vettore monodimensionale for(x) {for(y) {for(z)}}}
    std::vector<glm::vec3> get_discrete_space(glm::vec3 min, glm::vec3 max, int resolution);

    octree<gravity::cube>* get_discrete_space_as_octree(glm::vec3 min, glm::vec3 max, int resolution);

    glm::vec3 get_gravity_RT(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces, int resolution, glm::vec3 point
    );

    // for each face -> compute tetrahedron signed volume; find barycentre and add force to total
    glm::vec3 get_gravity_from_tetrahedrons(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& p, const glm::vec3& tetrahedrons_vertex);

    glm::vec3 get_gravity_from_tetrahedrons_corrected(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& p, const glm::vec3& tetrahedrons_vertex);

    float volume(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& t);
}

#endif //GL_TEST_PROJECT_GRAVITY_HPP
