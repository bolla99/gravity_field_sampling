//
// Created by Giovanni Bollati on 27/10/23.
//

#ifndef GL_TEST_PROJECT_GRAVITY_HPP
#define GL_TEST_PROJECT_GRAVITY_HPP

#include <vector>
#include <glm/glm.hpp>
#include <octree.hpp>
#include <util.hpp>

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
            int resolution,
            float* mass_R
            );

    glm::vec3 get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point);
    glm::vec3 get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, float R, glm::vec3 point);

    // get gravity given 3d space and gravity vector (with related min vector, range and resolution) and point
    // space and gravity vector are meant to be precomputed during a non real-time phase, while
    // this function is meant to be called in real-time by interpolating precomputed values and compute
    // gravity for a given point;
    // min is the min space vertex of the box of which gravity has been previoudly computed
    // min, range and resolution are needed to find indices of the bounding box, then
    // the spatial bounding box coordinates are retrived and used to interpolate the gravity values of each
    // bounding box vertex;
    glm::vec3 get_gravity_from_1D_precomputed_vector(glm::vec3 point, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution);

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



    // OCTREE CONSTRUCTION STUFF
    struct node {
        int next_id;
        std::array<glm::vec3, 8> octant;
    };
    struct gravity_field {
        glm::vec3 min;
        float edge;
        std::vector<node> gravity_octree;
    };

    inline node build_node(glm::vec3 min, float edge, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution) {
        return node {
            0,
            {
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x, min.y, min.z}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x + edge, min.y, min.z}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x, min.y + edge, min.z}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x + edge, min.y + edge, min.z}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x, min.y, min.z + edge}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x + edge, min.y, min.z + edge}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x, min.y + edge, min.z + edge}, gravity, space, resolution
                    ),
                gravity::get_gravity_from_1D_precomputed_vector(
                    glm::vec3{min.x + edge, min.y + edge, min.z + edge}, gravity, space, resolution
                    ),
            }
        };
    }

    // always true -> no optimization
    inline bool should_divide(float precision, node n, glm::vec3 local_min, float edge, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution) {
        std::array<glm::vec3, 8> cube = util::get_box(local_min, edge);
        glm::vec3 min = glm::vec3{local_min.x + edge/4.0, local_min.y + edge/4.0, local_min.z + edge/4.0};
        return                   !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y + edge/2.0, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y + edge/2.0, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z + edge/2.0}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z + edge/2.0}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z + edge/2.0}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y + edge/2.0, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, cube, n.octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                    );

    }

    int build_octree(
        float precision,
        std::vector<node>& octree,
        int id,
        int next_id,
        int max_res,
        glm::vec3 min,
        float edge,
        const std::vector<glm::vec3>& gravity,
        const std::vector<glm::vec3>& space, int resolution,
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces
        );
}


#endif //GL_TEST_PROJECT_GRAVITY_HPP
