//
// Created by Giovanni Bollati on 27/10/23.
//

#ifndef GL_TEST_PROJECT_GRAVITY_HPP
#define GL_TEST_PROJECT_GRAVITY_HPP

#include <vector>
#include <glm/glm.hpp>
#include <util.hpp>
#include <iostream>

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
            int resolution,
            float* mass_R
            );
    std::vector<mass> get_masses_from_tube(tube t, int resolution, const std::vector<glm::vec3>& vertices);

    std::vector<mass> get_masses(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution,
            float* mass_R
            );

    // for tubes parallel to z axis
    inline glm::vec3 get_gravity_from_tube(glm::vec3 p, tube t, float G, float tube_radius) {
        // t1 -> p
        auto t1_p = p - t.t1;

        // new coordinate system x_axis
        auto x = glm::normalize(t.t2 - t.t1);
        if(glm::distance2(t.t2, p) < glm::distance2(t.t1, p))
            x = glm::normalize(t.t1 - t.t2);
        // new coordinate system origin
        auto o = t.t1 + x * glm::dot(x, t1_p);
        //std::cout << "ORIGIN " << o.x << " " << o.y << " " << o.z << std::endl;
        // new coordinate system y and z axis
        auto y = glm::normalize(p - o);
        auto z = glm::normalize(glm::cross(x, y));

        auto a = glm::distance(p, o);

        // potenzialmente da togliere -> ridondante
        //if(a < std::numeric_limits<float>::epsilon() && a > -std::numeric_limits<float>::epsilon()) {
        //    return {0.0, 0.0, 0.0};
        //};
        // punto dentro al cilindro -> da gestire

        if(a < tube_radius) return {0.0, 0.0, 0.0};
        auto b = glm::distance(t.t1, o);
        auto c = glm::distance(t.t2, o);

        if(b > c) {
            auto tmp = b; b = c; c = tmp;
        }

        // check potenzialmente da togliere -> ridondanti
        //auto denom_1 = std::sqrt(std::pow(a, 2) + std::pow(c, 2));
        //auto denom_2 = std::sqrt(std::pow(a, 2) + std::pow(b, 2));

        //if(denom_1 < std::numeric_limits<float>::epsilon() && denom_1 > -std::numeric_limits<float>::epsilon()) {
        //    return {0.0, 0.0, 0.0};
        //};

        //if(denom_2 < std::numeric_limits<float>::epsilon() && denom_2 > -std::numeric_limits<float>::epsilon()) {
        //    return {0.0, 0.0, 0.0};
        //};

        float fy = ((G*M_PI*std::pow(tube_radius,2)) *
                    ((c / std::sqrt(std::pow(a, 2) + std::pow(c, 2))) - (b / std::sqrt(std::pow(a, 2) + std::pow(b, 2))))
                    ) / a;
        float fx = (G*M_PI*std::pow(tube_radius,2)) *
                    ((1 / std::sqrt(std::pow(a, 2) + std::pow(b, 2))) - (1 / std::sqrt(std::pow(a, 2) + std::pow(c, 2))));

        fy = -fy; //change fy sign

        // p in mezzo, incrementare fy
        if(p.z > t.t1.z && p.z < t.t2.z) {
            fy -= 2.f * ((G*M_PI*std::pow(tube_radius,2)) *
                    ((b / std::sqrt(std::pow(a, 2) + std::pow(b, 2))))
                    ) / a;
        }

        // potenzialmente da togliere per ridondanza, comunque rende l'algoritmo che usa questa funzione
        // più robusto
        if(isnan(fx) || isnan(fy)) return {0.0, 0.0, 0.0};

        //std::cout << "RELATIVE FORCE " << force.x << " " << force.y << std::endl;
        /*
        std::cout << "x axis: " << x.x << " " << x.y << " " << x.z << std::endl;
        std::cout << "y axis: " << y.x << " " << y.y << " " << y.z << std::endl;
        std::cout << "z axis: " << z.x << " " << z.y << " " << z.z << std::endl;
        */

        return glm::mat3{x, y, z} * glm::vec3{fx, fy, 0};
    }

    glm::vec3 get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point);
    glm::vec3 get_gravity_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float R);

    glm::vec3 get_gravity_from_mass(gravity::mass m, float G, float R, glm::vec3 point);
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


    // vettore monodimensionale for(x) {for(y) {for(z)}}}
    std::vector<glm::vec3> get_discrete_space(glm::vec3 min, glm::vec3 max, int resolution);

    glm::vec3 get_gravity_RT(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces, int resolution, glm::vec3 point
    );

    float volume(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& t);

    // OCTREE CONSTRUCTION STUFF
    struct node {
        // first_child_id = -1 -> leaf node
        int first_child_id;
        std::array<glm::vec3, 8>* gravity_octant;
    };

    struct gravity_field {
        glm::vec3 min;
        float edge;
        std::vector<node> gravity_octree;
    };

    // create leaf node ( leaf: next_id = 0 ) next_id should be renamed first_child_id
    inline node build_node(glm::vec3 min, float edge, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution) {
        return node {
            -1,
            new std::array<glm::vec3, 8>{
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

    inline node build_node_with_integral(glm::vec3 min, float edge, const std::vector<tube>& tubes, float G, float R) {
        return node {
            -1,
            new std::array<glm::vec3, 8>{
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x, min.y, min.z}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x + edge, min.y, min.z}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x, min.y + edge, min.z}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x + edge, min.y + edge, min.z}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x, min.y, min.z + edge}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x + edge, min.y, min.z + edge}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x, min.y + edge, min.z + edge}, tubes, G, R
                    ),
                gravity::get_gravity_from_tubes_with_integral(
                    glm::vec3{min.x + edge, min.y + edge, min.z + edge}, tubes, G, R
                    ),
            }
        };
    }

    // check if a leaf node should be divided
    inline bool should_divide(float precision, node n, glm::vec3 local_min, float edge, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution) {
        std::array<glm::vec3, 8> cube = util::get_box(local_min, edge);
        glm::vec3 min = glm::vec3{local_min.x + edge/4.0, local_min.y + edge/4.0, local_min.z + edge/4.0};
        return                   !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y + edge/2.0, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y + edge/2.0, min.z}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x, min.y + edge/2.0, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_1D_precomputed_vector({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, gravity, space, resolution),
                            precision
                    );
    }

        inline bool should_divide_with_integral(float precision, node n, glm::vec3 local_min, float edge, const std::vector<tube>& tubes, float G, float R) {
        std::array<glm::vec3, 8> cube = util::get_box(local_min, edge);
        glm::vec3 min = glm::vec3{local_min.x + edge/4.0, local_min.y + edge/4.0, local_min.z + edge/4.0};
        return                   !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x, min.y, min.z}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x + edge/2.0, min.y, min.z}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x, min.y + edge/2.0, min.z}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x + edge/2.0, min.y + edge/2.0, min.z}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x, min.y, min.z + edge/2.0}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x + edge/2.0, min.y, min.z + edge/2.0}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x, min.y + edge/2.0, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x, min.y + edge/2.0, min.z + edge/2.0}, tubes, G, R),
                            precision
                            ) || !util::vectors_are_p_equal(
                            util::interpolate({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, cube, *n.gravity_octant),
                            get_gravity_from_tubes_with_integral({min.x + edge/2.0, min.y + edge/2.0, min.z + edge/2.0}, tubes, G, R),
                            precision
                    );
    }

    void build_octree(
        float precision,
        std::vector<node>& octree,
        int id,
        int max_res,
        glm::vec3 min,
        float edge,
        const std::vector<glm::vec3>& gravity,
        const std::vector<glm::vec3>& space, int resolution,
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        const std::vector<tube>& tubes,
        float G, float R
        );

    void build_octree_with_integral(
        float precision,
        std::vector<node>& octree,
        int id,
        int max_res,
        glm::vec3 min,
        float edge,
        const std::vector<tube>& tubes,
        float G, float R
        );
}


#endif //GL_TEST_PROJECT_GRAVITY_HPP
