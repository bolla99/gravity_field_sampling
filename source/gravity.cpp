//
// Created by Giovanni Bollati on 27/10/23.
//

#include "gravity.hpp"
#include "GPUComputing.hpp"
#include "util.hpp"
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <SDL.h>

#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif

#include "glm/gtx/hash.hpp"
#include "glm/glm.hpp"
#include <glm/gtc/type_ptr.hpp>

std::vector<gravity::tube> gravity::get_tubes(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution,
        float* cylinder_R
        ) {
    if(resolution % 2 == 0) resolution++;

    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = util::get_min(vertices);

    //min.x -= 0.0001; min.y -= 0.0001;
    glm::vec3 max = util::get_max(vertices);

    //max.x += 0.0001; max.y += 0.0001;
    min.z = min.z - 10;
    max.z = min.z;

    glm::vec3 center = (min + max) / 2.0f;

    /*
    std::cout << "get tube; min_1: " << min.x << " " << min.y << " " << min.z << std::endl;
    std::cout << "get tube; max_1: " << max.x << " " << max.y << " " << max.z << std::endl;
    std::cout << "get tube; center: " << center.x << " " << center.y << " " << center.z << std::endl;
    */

    double max_extent = std::max(max.x - min.x, max.y - min.y);

    //std::cout << "max_extent: " << max_extent << std::endl;

    min = {center.x - (max_extent / 2.0), center.y - (max_extent / 2.0), center.z};

    /*
    std::cout << "get tube; min: " << min.x << " " << min.y << " " << min.z << std::endl;
    std::cout << "get tube; max: " << max.x << " " << max.y << " " << max.z << std::endl;
    */

    // CUBE EDGE LENGTH
    double cube_edge = (double)max_extent / (double)resolution;

    *cylinder_R = (float)std::sqrt(std::pow(cube_edge, 2) / M_PI);

    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

    auto real_max = glm::vec3{
        min.x + (float)resolution * (float)cube_edge,
        min.y + (float)resolution * (float)cube_edge,
        min.z + (float)resolution * (float)cube_edge
    };
    std::cout << "get tube; max reale: " << real_max.x << " " << real_max.y << " " << real_max.z << std::endl;
    //std::cout << "SCARTO ERRORE TRA MAX REALE E MAX CALCOLATO: " << (max.x) - (min.x + (float)((double)(resolution - 25) * cube_edge)) << std::endl;

    int t_left = 0;
    int t_right = 0;

    int t_up = 0;
    float d_up = 0;
    int t_down = 0;
    float d_down = 0;

#pragma omp parallel for default(none), shared(cube_edge, resolution, ray_dir, tubes, min, vertices, faces, std::cout, t_left, t_right, t_up, t_down, d_up, d_down)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            glm::vec3 ray_origin{min.x + (float)((double)i * cube_edge), min.y + (float)((double)j * cube_edge), min.z};

            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, ray_origin, ray_dir);

            // FOR EACH INTERSECTIONS COUPLE
            for(int k = 0; k < intersections.size(); k += 2) {
                /*
                if(glm::distance(intersections[k], intersections[k + 1]) > 10) {
                    std::cout << "impossible tube second check" << std::endl;
                }
                */
#pragma omp critical
                tubes.push_back({intersections[k], intersections[k + 1]});
                if(i < (resolution + 1) / 2) t_left++; else t_right++;
                if(j < (resolution + 1) / 2) {
                    t_down++;
                    d_down += glm::distance(intersections[k], intersections[k + 1]);
                } else {
                    t_up++;
                    d_up += glm::distance(intersections[k], intersections[k + 1]);
                }
            }
        }
    }

    /*
    std::cout << "t up: " << t_up << std::endl;
    std::cout << "t down: " << t_down << std::endl;
    std::cout << "t left: " << t_left << std::endl;
    std::cout << "t right: " << t_right << std::endl;
    std::cout << "d up: " << d_up << std::endl;
    std::cout << "d down: " << d_down << std::endl;
    */

    std::cout << "get tubes time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    return tubes;
}

glm::vec3 gravity::get_gravity_from_tubes_with_integral(glm::vec3 point, const std::vector<tube>& tubes, float G, float cylinder_R) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }
    glm::vec3 force_from_integral{};

#pragma omp parallel for default(none) shared(tubes, point, G, thread_gravity, cylinder_R)
    for(auto t : tubes) {
        thread_gravity[omp_get_thread_num()] += gravity::get_gravity_from_tube_with_integral(point, t, G, cylinder_R);
    }
    for(int i = 0; i < omp_get_max_threads(); i++) force_from_integral += thread_gravity[i];
    return force_from_integral;
}

glm::vec3 gravity::get_gravity_from_tubes_with_integral_with_gpu(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R) {
    auto output = GPUComputing::get_gravity_from_tubes_with_integral(glm::value_ptr(tubes.front().t1), tubes.size(), glm::value_ptr(point), cylinder_R, G);

    glm::vec3 output_gravity{0.f, 0.f, 0.f};

    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];

    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0.f, 0.f, 0.f};
    }

#pragma omp parallel for default(none) shared(output, thread_gravity, tubes, std::cout)
    for(int i = 0; i < tubes.size(); i++) {
        thread_gravity[omp_get_thread_num()] += glm::vec3(output[3*i], output[3*i + 1], output[3*i + 2]);
    }

    for(int i = 0; i < omp_get_max_threads(); i++) {
        output_gravity += thread_gravity[i];
    }
    delete output;
    return output_gravity;
}

std::vector<glm::vec3> gravity::get_discrete_space(glm::vec3 min, float edge, int resolution) {
    omp_set_num_threads(omp_get_max_threads());

    float unit = edge / (float)resolution;

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

void gravity::build_octree(
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
        ) {

    // first eight values -> first eight values
    auto location = util::get_box(min, edge);
    auto int_location = util::get_int_box(int_min, int_edge);
    for(int i = 0; i < 8; i++) {
        // check if values are already cached; values index in negative integer
        if(auto k1 = gravity_values_map.find(int_location[i]); k1 != gravity_values_map.end()) {
            octree.push_back(-(k1->second));
        } else if(auto k2 = cached_values.find(int_location[i]); k2 != cached_values.end()) {
            int j = gravity_values.size();
            octree.push_back(-j);
            gravity_values.push_back(tmp_gravity_values[k2->second]);
            gravity_values_map.emplace(int_location[i], j);
        } else {
            int new_value_index = (int)gravity_values.size();
            gravity_values.push_back(get_gravity_from_tubes_with_integral_with_gpu(location[i], tubes, G, R));
            gravity_values_map.emplace(int_location[i], new_value_index);
            //tmp_gravity_values.push_back(gravity_values[new_value_index]);
            //cached_values.emplace(int_location[i], new_value_index);
            octree.push_back(-new_value_index);
        }
    }

    if(max_res > 0 && should_divide(sd_method, precision, id, octree, max_res,min, edge, int_min, int_edge, gravity_values, tmp_gravity_values, cached_values, gravity_values_map, tubes, G, R)) {
        // inspecting node became internal node -> it gets children, and first_child_id must be updated with first
        // child id; the children are all adjacent; then delete gravity_octant, since it is no longer needed

        auto mins = util::get_box(min, edge/2.f);
        auto int_mins = util::get_int_box(int_min, int_edge/2);

        for(int i = 0; i < 8; i++) {
            octree[id + i] = (int)octree.size();
            build_octree(
                    sd_method, precision, octree, octree[id + i], max_res - 1, mins[i],
                    edge/2.f, int_mins[i],
                    int_edge / 2, gravity_values, tmp_gravity_values, cached_values, gravity_values_map, tubes, G, R
                    );
        }
    }
}

glm::vec3 gravity::get_gravity_from_octree(glm::vec3 p, const std::vector<int>& octree, glm::vec3 min, float edge, const std::vector<glm::vec3>& gravity_values, int* depth) {
    auto current_min = min;
    auto current_edge = edge;

    std::array<glm::vec3, 8> values{};
    std::array<glm::vec3, 8> box{};

    int i = 0;
    for(int d = 0;; d++) {
        if(octree[i] <= 0) {
            // leaf: from i to i + 8 retreive gravity values and interpolate
            for(int j = 0; j < 8; j++) values[j] = gravity_values[-octree[i + j]];
            box = util::get_box(glm::make_vec3(current_min), current_edge);
            *depth = d;
            return util::interpolate(
                    glm::make_vec3(p),
                    box,
                    values);
        } else {
            int k = 0;
            if(p[0] > current_min[0] + current_edge/2.f) k += 1;
            if(p[1] > current_min[1] + current_edge/2.f) k += 2;
            if(p[2] > current_min[2] + current_edge/2.f) k += 4;
            i = octree[i + k];
            current_edge /= 2.f;
            current_min = util::get_box(glm::make_vec3(current_min), current_edge)[k];
        }
    }
}

float gravity::potential::get_potential_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R) {
    omp_set_num_threads(omp_get_max_threads());
    float thread_potential[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_potential[i] = 0.f;
    }
    float potential = 0.f;

#pragma omp parallel for default(none) shared(tubes, point, G, thread_potential, cylinder_R)
    for(auto t : tubes) {
        thread_potential[omp_get_thread_num()] += gravity::potential::get_potential_from_tube_with_integral(point, t, G, cylinder_R);
    }
    for(int i = 0; i < omp_get_max_threads(); i++) potential += thread_potential[i];
    return potential;
}

float gravity::potential::get_potential_with_gpu(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R) {
    auto output = GPUComputing::get_potential_from_tubes_with_integral(glm::value_ptr(tubes.front().t1), tubes.size(), glm::value_ptr(point), cylinder_R, G);

    float output_potential = 0.f;

    omp_set_num_threads(omp_get_max_threads());
    float thread_potential[omp_get_max_threads()];

    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_potential[i] = 0.f;
    }

#pragma omp parallel for default(none) shared(output, thread_potential, tubes, std::cout)
    for(int i = 0; i < tubes.size(); i++) {
        thread_potential[omp_get_thread_num()] += output[i];
    }

    for(int i = 0; i < omp_get_max_threads(); i++) {
        output_potential += thread_potential[i];
    }
    delete output;
    return output_potential;
}

void gravity::potential::build_octree(
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
) {
    auto locations = util::get_box(min ,edge);
    auto int_locations = util::get_int_box(int_min, int_edge);
    std::array<float, 8> potential_values{};

    for(int i = 0; i < 8; i++) {
        if (auto k = cached_values.find(int_locations[i]); k != cached_values.end()) {
            float f_value = k->second;
            potential_values[i] = f_value;
            int *value = reinterpret_cast<int *>(&f_value);
            octree.push_back(*value);
        } else {
            float f_value = gravity::potential::get_potential_with_gpu(locations[i], tubes, G, R);
            potential_values[i] = f_value;
            cached_values.emplace(int_locations[i], f_value);
            int *value = reinterpret_cast<int *>(&f_value);
            octree.push_back(*value);
        }
    }

    if(gravity::potential::should_divide(sd_method, alpha, precision, potential_values, min, edge, tubes, G, R) && max_res > 0) {
        auto mins = util::get_box(min, edge/2.f);
        auto int_mins = util::get_int_box(int_min, int_edge/2);

        for(int i = 0; i < 8; i++) {
            octree[id + i] = (int)octree.size();
            build_octree(
                        sd_method, alpha, precision, octree, octree[id + i], max_res - 1, mins[i],
                        edge/2.f, int_mins[i],
                        int_edge / 2, cached_values, tubes, G, R);
        }
    }
}

glm::vec3 gravity::potential::get_gravity_from_octree(glm::vec3 p, const std::vector<int>& octree, glm::vec3 min, float edge, int* depth) {
    auto current_min = min;
    auto current_edge = edge;

    int i = 0;
    for(int d = 0;; d++) {
        if(octree[i] <= 0) {
            // leaf: from i to i + 8 retreive gravity values and interpolate
            auto box = util::get_box(glm::make_vec3(current_min), current_edge);
            auto values = std::array<float, 8>{};
            for(int j = 0; j < 8; j++) {
                auto pot = reinterpret_cast<const float*>(&octree[i+j]);
                values[j] = *pot;
            }
            *depth = d;
            return util::get_gradient_from_box(glm::make_vec3(p), box, values);
        } else {
            int k = 0;
            if(p[0] > current_min[0] + current_edge/2.f) k += 1;
            if(p[1] > current_min[1] + current_edge/2.f) k += 2;
            if(p[2] > current_min[2] + current_edge/2.f) k += 4;
            i = octree[i + k];
            current_edge /= 2.f;
            current_min = util::get_box(glm::make_vec3(current_min), current_edge)[k];
        }
    }
}


