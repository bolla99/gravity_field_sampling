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

// this namespace depends on util and gpucomputing namespaces
// functions that depends on gpucomputing end with _with_gpu
// almost every function depends on util functions
// subnamespace potential is defined
namespace gravity {
    // tube -> line segment
    struct tube {
        glm::vec3 t1, t2;
    };

    // returns vector of tubes that represents {vertices, faces} mesh;
    // resolution +1 -> number of ray per axis
    // after execution cylinder_R contains tubes radius
    std::vector<tube> get_tubes(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            int resolution,
            float* cylinder_R
            );

    // returns gravity inducted by a single tube t using integral method
    inline glm::vec3 get_gravity_from_tube_with_integral(glm::vec3 p, tube t, float G, float cylinder_R) {
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

    // returns gravity with integral method
    // DEPENDS on get_gravity_from_tube_with_integral
    glm::vec3 get_gravity_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

    // return gravity with integral method implemented with gpu computing
    glm::vec3 get_gravity_from_tubes_with_integral_with_gpu(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

    // monodimensional vector for(x) {for(y) {for(z)}}}
    // resolution stands for number of segments, which means resolution + 1 gravity samples
    std::vector<glm::vec3> get_discrete_space(glm::vec3 min, float edge, int resolution);

    // used by build_octree (gravity version)
    // sd_method = 0 -> center is used
    // sd_method = 1 -> eight values (centers of the eight inner cubes)
    // sd_method = 2 -> both
    inline bool should_divide(
        int sd_method,
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

        /*if(max_depth == 1) {
            return true;
        }*/

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

        auto ilocations = util::get_int_box(int_min, int_edge/2);
        auto ilocations_vec = std::vector<glm::ivec3>{};

        glm::vec3 test_value{0.f, 0.f, 0.f};

        min = box_min_position;
        int_min = int_box_min_position;

        // add center and half edges to locations / i locations
        if(sd_method == 0) {
            locations_vec.emplace_back(min.x + edge / 2.f, min.y + edge / 2.f, min.z + edge / 2.f);
            ilocations_vec.emplace_back(int_min.x + int_edge / 2, int_min.y + int_edge / 2, int_min.z + int_edge / 2);
        } else if(sd_method == 1) {
            for(int i = 0; i < 8; i++) locations_vec.emplace_back(locations[i]);
            for(int i = 0; i < 8; i++) ilocations_vec.emplace_back(ilocations[i]);
        }  else {
            locations_vec.emplace_back(min.x + edge / 2.f, min.y + edge / 2.f, min.z + edge / 2.f);
            ilocations_vec.emplace_back(int_min.x + int_edge / 2, int_min.y + int_edge / 2, int_min.z + int_edge / 2);
            for(int i = 0; i < 8; i++) locations_vec.emplace_back(locations[i]);
            for(int i = 0; i < 8; i++) ilocations_vec.emplace_back(ilocations[i]);
        }

        static int n = 0;
        for(int i = 0; i < locations_vec.size(); i++) {
            if(auto j = cached_values.find(ilocations_vec[i]); j != cached_values.end()) {
                test_value = tmp_gravity_values[j->second];
            } else {
                test_value = get_gravity_from_tubes_with_integral_with_gpu(locations_vec[i], tubes, G, R);
                if(max_depth > 1) {
                    cached_values.emplace(ilocations_vec[i], tmp_gravity_values.size());
                    tmp_gravity_values.push_back(test_value);
                }
            }
            if(!util::vectors_are_p_equal(util::interpolate(locations_vec[i], box, values), test_value, precision)) return true;
        }
        return false;
    }

    // build octree
    // data is stored in octree and gravity_values parameters
    // DEPENDS on should_divide
    // sd_method specification in should_divide comment
    void build_octree(
        int sd_method,
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

    // returns gravity in p from given gravity octree, min and edge and gravity_values
    // depth stores depth at which gravity is computed (depth of box leaf that contains locaiton p)
    glm::vec3 get_gravity_from_octree(glm::vec3 p, const std::vector<int>& octree, glm::vec3 min, float edge, const std::vector<glm::vec3>& gravity_values, int* depth);


    // namespace potential contains functions related to potential computation
    namespace potential {

        // returns potential inducted by a single tube with integral method
        inline float get_potential_from_tube_with_integral(glm::vec3 p, tube t, float G, float cylinder_R) {
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

        // returns potential from tubes with integral method
        // DEPENDS on get_potential_from_tube_with_integral
        float get_potential_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

        // return potential computed with integral method with gpu
        float get_potential_with_gpu(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R);

        // build octree of potential
        // octree is stored in octree parameter
        // DEPENDS on potential::should_divide
        // for sd_method specification see should_divide comment
        void build_octree(
                int sd_method,
                float alpha,
                float precision,
                std::vector<int>& octree,
                int id,
                int max_res,
                glm::vec3 min,
                float edge,
                glm::vec3 int_min,
                int int_edge,
                std::unordered_map<glm::ivec3, float>& cached_values,
                const std::vector<tube>& tubes,
                float G, float R
                );

        // should divide method for potential::build_octree
        // sd_method == 0 -> check gravity in center
        // sd_method == 1 -> check gravity in the centers of the eight innner cubes
        // sd_method == 2 -> both
        // sd_method == 3 -> ad hoc method for potential octree build
        inline bool should_divide(int sd_method, float alpha, float precision, std::array<float, 8> values, glm::vec3 min, float edge, const std::vector<tube>& tubes, float G, float R) {
            if(sd_method == 3) {
                auto x_values = util::get_x_derivative_from_cube(values);
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        if (abs(x_values[i] - x_values[j]) > alpha * abs(std::min(x_values[i], x_values[j])))
                            return true;
                    }
                }
                auto y_values = util::get_y_derivative_from_cube(values);
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        if (abs(y_values[i] - y_values[j]) > alpha * abs(std::min(y_values[i], y_values[j])))
                            return true;
                    }
                }
                auto z_values = util::get_z_derivative_from_cube(values);
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        if (abs(z_values[i] - z_values[j]) > alpha * abs(std::min(z_values[i], z_values[j])))
                            return true;
                    }
                }
                return false;
            } else if(sd_method < 3) {
                std::vector<glm::vec3> locations{};
                auto locations_8 = util::get_box(glm::vec3{min.x + edge/4.0, min.y + edge/4.0, min.z + edge/4.0}, edge/2.f);
                auto center = min + edge/2.f;
                if(sd_method == 0) {
                    locations.push_back(center);
                } else if(sd_method == 1) {
                    for(int i = 0; i < 8; i++) { locations.push_back(locations_8[i]); }
                } else if(sd_method == 2) {
                    locations.push_back(center);
                    for(int i = 0; i < 8; i++) { locations.push_back(locations_8[i]); }
                }
                for(auto location : locations) {
                    if(!util::vectors_are_p_equal(
                            get_gravity_from_tubes_with_integral_with_gpu(location, tubes, G, R),
                            util::get_gradient_from_box(location, util::get_box(min, edge), values), precision)
                            ) return true;
                }
                return false;
            }
            return false;
        }

        // returns gravity in p from given potential octree, min and edge
        // depth stores depth at which gravity is computed (depth of box leaf that contains locaiton p)
        glm::vec3 get_gravity_from_octree(glm::vec3 p, const std::vector<int>& octree, glm::vec3 min, float edge, int* depth);
    }
}


#endif //GL_TEST_PROJECT_GRAVITY_HPP
