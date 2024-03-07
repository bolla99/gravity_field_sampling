//
// Created by Giovanni Bollati on 27/10/23.
//

#ifndef GL_TEST_PROJECT_GRAVITY_HPP
#define GL_TEST_PROJECT_GRAVITY_HPP

#include <vector>
#include <glm/glm.hpp>
#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif
#include <glm/gtx/hash.hpp>
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
    struct box {
        glm::vec3 center;
        float edge;
    };
    struct gravity_cube {
        box b;
        std::array<glm::vec3, 8> g;
    };

    std::vector<tube> get_tubes(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution,
            float* cylinder_R
            );
    std::vector<mass> get_masses_from_tube(tube t, int resolution, const std::vector<glm::vec3>& vertices);

    std::vector<mass> get_masses(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution,
            float* sphere_R
            );

    // for tubes parallel to z axis
    inline glm::vec3 get_gravity_from_tube(glm::vec3 p, tube t, float G, float cylinder_R) {
        // t1 -> p
        auto t1_p = p - t.t1;

        // new coordinate system x_axis
        auto x = glm::normalize(t.t2 - t.t1);
        if(glm::distance2(t.t2, p) < glm::distance2(t.t1, p))
            x = glm::normalize(t.t1 - t.t2);

        // new coordinate system origin
        auto o = t.t1 + x * glm::dot(x, t1_p);

        // new coordinate system y and z axis
        auto y = glm::normalize(p - o);
        auto z = glm::normalize(glm::cross(x, y));

        auto a = glm::distance(p, o);

        // potenzialmente da togliere -> ridondante
        //if(a < std::numeric_limits<float>::epsilon() && a > -std::numeric_limits<float>::epsilon()) {
        //    return {0.0, 0.0, 0.0};
        //};
        // punto dentro al cilindro -> da gestire

        if(a < cylinder_R) return {0.0, 0.0, 0.0};
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

        float fy = ((G*M_PI*std::pow(cylinder_R,2)) *
                    ((c / std::sqrt(std::pow(a, 2) + std::pow(c, 2))) - (b / std::sqrt(std::pow(a, 2) + std::pow(b, 2))))
                    ) / a;
        float fx = (G*M_PI*std::pow(cylinder_R,2)) *
                    ((1 / std::sqrt(std::pow(a, 2) + std::pow(b, 2))) - (1 / std::sqrt(std::pow(a, 2) + std::pow(c, 2))));

        fy = -fy; //change fy sign

        // p in mezzo, incrementare fy
        if(p.z > t.t1.z && p.z < t.t2.z) {
            fy -= 2.f * ((G*M_PI*std::pow(cylinder_R,2)) *
                    ((b / std::sqrt(std::pow(a, 2) + std::pow(b, 2))))
                    ) / a;
        }

        // potenzialmente da togliere per ridondanza, comunque rende l'algoritmo che usa questa funzione
        // più robusto
        if(isnan(fx) || isnan(fy)) return {0.0, 0.0, 0.0};

        return glm::mat3{x, y, z} * glm::vec3{fx, fy, 0};
    }

    inline float get_potential_from_tube(glm::vec3 p, tube t, float G, float cylinder_R) {
        auto t1_p = p - t.t1;

        int left = 1;

        // new coordinate system x_axis
        auto x = glm::normalize(t.t2 - t.t1);
        if(glm::distance2(t.t2, p) < glm::distance2(t.t1, p)) {
            x = glm::normalize(t.t1 - t.t2);
            left = -1;
        }

        // new coordinate system origin
        auto o = t.t1 + x * glm::dot(x, t1_p);

        auto a = std::max(glm::distance(p, o), cylinder_R);
        auto b = glm::distance(t.t1, o);
        auto c = glm::distance(t.t2, o);

        auto integral_c = 0.5f*(float)std::log(2*c*(sqrt(pow(a, 2) + pow(c, 2)) + c)/pow(a, 2) + 1);
        auto integral_b = 0.5f*(float)std::log(2*b*(sqrt(pow(a, 2) + pow(b, 2)) + b)/pow(a, 2) + 1);

        // se p + in mezzo al tubo, devo sommare integrale da b a c a 2*integrale da 0 a b
        if(p.z >= t.t1.z && p.z <= t.t2.z) {
            return G*(float)std::pow(cylinder_R, 2)*(float)M_PI*(integral_c + integral_b);
        }

        return G*(float)std::pow(cylinder_R, 2)*(float)M_PI*((float)left*integral_c - (float)left*integral_b);
    }

    glm::vec3 get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<tube>& tubes, glm::vec3 point);
    glm::vec3 get_gravity_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);
    glm::vec3 get_gravity_from_tubes_with_integral_with_gpu(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

    float get_potential_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

    glm::vec3 get_gravity_from_mass(gravity::mass m, float G, float sphere_R, glm::vec3 point);
    glm::vec3 get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, float sphere_R, glm::vec3 point);

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
    // divite lo spazio su un asse in resolution segmenti (resolution + 1 campioni)
    std::vector<glm::vec3> get_discrete_space(glm::vec3 min, float edge, int resolution);

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

    // if child > 0 then child represent the index of the first of the eight children, so the node
    // is an internal node; if child <= 0, then abs(child) represents the index of an array from which you
    // can retreive the gravity value for that node -> the node is a leaf
    struct octree {
        int child;
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

    inline bool should_divide_with_integral_optimized(
        float precision,
        int id,
        const std::vector<int>& octree,
        int max_depth,
        glm::vec3 box_min_position,
        float edge,
        glm::ivec3 int_box_min_position,
        int int_edge,
        const std::vector<glm::vec3>& gravity_values,
        std::vector<glm::vec3>& tmp_gravity_values,
        std::unordered_map<glm::vec<3, int>, int>& cached_values,
        std::unordered_map<glm::vec<3, int>, int>& gravity_values_map,
        const std::vector<tube>& tubes,
        float G,
        float R) {

        if(max_depth == 1) {
            return true;
        }

        // box that is being tested (if it should be divided in eight boxes)
        std::array<glm::vec3, 8> box = util::get_box(box_min_position, edge);

        std::array<glm::vec3, 8> values{};
        for(int i = 0; i < 8; i++) {
            values[i] = gravity_values[-octree[id + i]];
        }

        // min of the eight test location where interpolated and real gravity are compared
        glm::vec3 min = glm::vec3{box_min_position.x + edge/4.0, box_min_position.y + edge/4.0, box_min_position.z + edge/4.0};
        glm::ivec3 int_min = glm::ivec3{int_box_min_position.x + int_edge/4, int_box_min_position.y + int_edge/4, int_box_min_position.z + int_edge/4};

        auto locations = util::get_box(min, edge/2.f);
        auto locations_vec = std::vector<glm::vec3>{};
        for(int i = 0; i < 8; i++) locations_vec.emplace_back(locations[i]);

        auto ilocations = util::get_int_box(int_min, int_edge/2);
        auto ilocations_vec = std::vector<glm::ivec3>{};
        for(int i = 0; i < 8; i++) ilocations_vec.emplace_back(ilocations[i]);

        glm::vec3 test_value;

        // add center and half edges to locations / i locations
        locations_vec.emplace_back(min.x + edge/2.f, min.y + edge/2.f, min.z + edge/2.f);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y + int_edge/2, int_min.z + int_edge/2);

        locations_vec.emplace_back(min.x + edge/2.f, min.y, min.z);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y, int_min.z);
        locations_vec.emplace_back(min.x, min.y + edge/2.f, min.z);
        ilocations_vec.emplace_back(int_min.x, int_min.y  + int_edge/2, int_min.z);
        locations_vec.emplace_back(min.x + edge/2.f, min.y + edge, min.z);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y + int_edge, int_min.z);
        locations_vec.emplace_back(min.x + edge, min.y + edge/2.f, min.z);
        ilocations_vec.emplace_back(int_min.x + int_edge, int_min.y + int_edge/2, int_min.z);

        locations_vec.emplace_back(min.x + edge/2.f, min.y, min.z + edge/2.f);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y, int_min.z + int_edge/2);
        locations_vec.emplace_back(min.x, min.y + edge/2.f, min.z + edge/2.f);
        ilocations_vec.emplace_back(int_min.x, int_min.y  + int_edge/2, int_min.z + int_edge/2);
        locations_vec.emplace_back(min.x + edge/2.f, min.y + edge, min.z + edge/2.f);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y + int_edge, int_min.z + int_edge/2);
        locations_vec.emplace_back(min.x + edge, min.y + edge/2.f, min.z + edge/2.f);
        ilocations_vec.emplace_back(int_min.x + int_edge, int_min.y + int_edge/2, int_min.z + int_edge/2);

        locations_vec.emplace_back(min.x + edge/2.f, min.y, min.z + edge);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y, int_min.z + int_edge);
        locations_vec.emplace_back(min.x, min.y + edge/2.f, min.z + edge);
        ilocations_vec.emplace_back(int_min.x, int_min.y  + int_edge/2, int_min.z + int_edge);
        locations_vec.emplace_back(min.x + edge/2.f, min.y + edge, min.z + edge);
        ilocations_vec.emplace_back(int_min.x + int_edge/2, int_min.y + int_edge, int_min.z + int_edge);
        locations_vec.emplace_back(min.x + edge, min.y + edge/2.f, min.z + edge);
        ilocations_vec.emplace_back(int_min.x + int_edge, int_min.y + int_edge/2, int_min.z + int_edge);

        for(int i = 0; i < locations_vec.size(); i++) {
            if(auto j = cached_values.find(ilocations[i]); j != cached_values.end()) {
                test_value = tmp_gravity_values[j->second];
            } else {
                test_value = get_gravity_from_tubes_with_integral_with_gpu(locations[i], tubes, G, R);
                if(max_depth > 1) {
                    cached_values.emplace(ilocations[i], tmp_gravity_values.size());
                    tmp_gravity_values.push_back(test_value);
                }
            }
            if(!util::vectors_are_p_equal(util::interpolate(locations[i], box, values), test_value, precision)) return true;
        }
        return false;
    }

    void build_octree(
        float precision,
        std::vector<node>& octree,
        int id,
        int max_res,
        glm::vec3 min,
        float edge,
        const std::vector<glm::vec3>& gravity,
        const std::vector<glm::vec3>& space, int resolution
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

    void build_octree_with_integral_optimized(
        float precision,
        std::vector<int>& octree,
        int id,
        int max_res,
        glm::vec3 min,
        float edge,
        glm::vec<3, int> int_min,
        int int_edge,
        std::vector<glm::vec3>& gravity_values,
        std::vector<glm::vec3>& tmp_gravity_values,
        std::unordered_map<glm::vec<3, int>, int>& cached_values,
        std::unordered_map<glm::vec<3, int>, int>& gravity_values_map,
        const std::vector<tube>& tubes,
        float G, float R
        );
}


#endif //GL_TEST_PROJECT_GRAVITY_HPP
