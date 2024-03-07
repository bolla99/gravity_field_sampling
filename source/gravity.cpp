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
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = util::get_min(vertices);
    //min.x += 0.0001; min.y += 0.0001;
    min.x -= 0.0001; min.y -= 0.0001;
    glm::vec3 max = util::get_max(vertices);

    //max.x -= 0.0001; max.y -= 0.0001;
    max.x += 0.0001; max.y += 0.0001;
    min.z = min.z - 10;
    max.z = min.z;

    glm::vec3 center = (min + max) / 2.0f;

    double max_extent = std::max(max.x - min.x, max.y - min.y);

    min = {center.x - (max_extent / 2.0), center.y - (max_extent / 2.0f), center.z};

    // CUBE EDGE LENGTH
    double cube_edge = (double)max_extent / (double)resolution;

    *cylinder_R = (float)std::sqrt(std::pow(cube_edge, 2) / M_PI);

    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

    std::cout << "SCARTO ERRORE TRA MAX REALE E MAX CALCOLATO: " << (max.x) - (min.x + (float)((double)(resolution - 25) * cube_edge)) << std::endl;

#pragma omp parallel for default(none), shared(cube_edge, resolution, ray_dir, tubes, min, vertices, faces, std::cout)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)((double)i * cube_edge), min.y + (float)((double)j * cube_edge), min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, r.origin, r.dir);

            // FOR EACH INTERSECTIONS COUPLE
            for(int k = 0; k < intersections.size(); k += 2) {
                /*
                if(glm::distance(intersections[k], intersections[k + 1]) > 10) {
                    std::cout << "impossible tube second check" << std::endl;
                }*/
                // TUBE
#pragma omp critical
                tubes.push_back({intersections[k], intersections[k + 1]});
            }
        }
    }
    std::cout << "get tubes time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    return tubes;
}

std::vector<gravity::mass> gravity::get_masses_from_tube(tube t, int resolution, const std::vector<glm::vec3>& vertices) {
    auto cube_edge = util::xy_mesh_max_extent(vertices) / (float)resolution;
    auto cube_volume = (float)std::pow(cube_edge, 3);

    float units_per_tubef, new_cube_volume;
    int units_per_tube;

    // TUBE LENGTH / UNIT
    units_per_tubef = glm::length(t.t2 - t.t1) / cube_edge;
    // ROUNDED
    units_per_tube = (int) std::round((double) units_per_tubef);
    // NEW CUBE MASS
    new_cube_volume = (cube_volume * units_per_tubef) / (float) units_per_tube;

    // TUBE VERTEX OFFSET
    float tube_unit = glm::length(t.t2 - t.t1) / (float) units_per_tube;

    std::vector<gravity::mass> masses(units_per_tube);

    for (int n = 0; n < units_per_tube + 1; n++) {
        masses.push_back({t.t1 + glm::normalize(t.t2 - t.t1) * (float) n * tube_unit, new_cube_volume});
    }
    return masses;
}

std::vector<gravity::mass> gravity::get_masses(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution,
        float* sphere_R
        ) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};
    std::vector<mass> volumes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = util::get_min(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::get_max(vertices);

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

    // SET MASS SPHERE RADIUS
    *sphere_R = (float)std::cbrt(std::pow(cube_edge, 3) * (3.f / 4.f) / M_PI);

    auto cube_volume = (float)std::pow(cube_edge, 3);
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(resolution, cube_edge, cube_volume, ray_dir, tubes, min, volumes, vertices, faces, sphere_R)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, r.origin, r.dir);
            if(intersections.size() == 3) intersections.erase(intersections.end() - 1);

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
#pragma omp critical
                *sphere_R = std::min(*sphere_R, (float)std::cbrt(std::pow(tube_unit, 3) * (3.f / 4.f) / M_PI));
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

glm::vec3 gravity::get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<tube>& tubes, glm::vec3 point) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }

    Uint64 start = SDL_GetTicks64();
    // FIND XY PLANE
    glm::vec3 min = util::get_min(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::get_max(vertices);

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
    glm::vec3 gravity = {0, 0, 0};
    for(int i = 0; i < omp_get_max_threads(); i++) {
        gravity = gravity + thread_gravity[i];
    }
    return gravity;
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
        thread_gravity[omp_get_thread_num()] += gravity::get_gravity_from_tube(point, t, G, cylinder_R);
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
        thread_gravity[i] = {0, 0, 0};
    }

#pragma omp parallel for default(none) shared(output, thread_gravity, tubes)
    for(int i = 0; i < tubes.size(); i++) {
        thread_gravity[omp_get_thread_num()] += glm::vec3(output[3*i], output[3*i + 1], output[3*i + 2]);
    }

    for(int i = 0; i < omp_get_max_threads(); i++) {
        output_gravity += thread_gravity[i];
    }
    delete output;
    return output_gravity;
}

float gravity::get_potential_from_tubes_with_integral(glm::vec3 point, const std::vector<gravity::tube>& tubes, float G, float cylinder_R) {
    omp_set_num_threads(omp_get_max_threads());
    float thread_potential[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_potential[i] = 0.f;
    }
    float potential = 0.f;

#pragma omp parallel for default(none) shared(tubes, point, G, thread_potential, cylinder_R)
    for(auto t : tubes) {
        thread_potential[omp_get_thread_num()] += gravity::get_potential_from_tube(point, t, G, cylinder_R);
    }
    for(int i = 0; i < omp_get_max_threads(); i++) potential += thread_potential[i];
    return potential;
}

glm::vec3 gravity::get_gravity_from_mass(gravity::mass m, float G, float sphere_R, glm::vec3 point) {
    glm::vec3 dir = m.p - point;
    auto r3 = (float)pow(glm::length(dir), 3);
    if(glm::length(dir) < sphere_R) {
        r3 = (float)std::pow(sphere_R, 3);
    }
    if (r3 > -std::numeric_limits<float>::epsilon() && r3 < std::numeric_limits<float>::epsilon()) {
        return {0.0, 0.0, 0.0};
    }
    return (dir * m.m * G) / r3;
}

glm::vec3 gravity::get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, float sphere_R, glm::vec3 point) {
    Uint64 start = SDL_GetTicks64();
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }
#pragma omp parallel for default(none) shared(masses, point, G, thread_gravity, sphere_R)
    for(auto & mass : masses) {
        glm::vec3 dir = mass.p - point;
        auto r3 = (float)pow(glm::length(dir), 3);
        if(glm::length(dir) < sphere_R) {
#pragma omp critical
            r3 = (float)std::pow(sphere_R, 3);
        }
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

glm::vec3 gravity::get_gravity_from_1D_precomputed_vector(glm::vec3 point, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, int resolution) {
    // get indices of bounding box
    auto indices = util::get_box_indices(space.front(), abs(space.back().x - space.front().x), resolution, point);

    // translate 3d indices to 1d indices
    std::array<int, 8> indices_1d{};
    for(int i = 0; i < 8; i++) {
        indices_1d[i] = util::from_3d_indices_to_1d(indices[i], resolution);
    }

    // get space coordinates of bounding box
    std::array<glm::vec3, 8> space_cube{};
    for(int i = 0; i < 8; i++) {
        space_cube[i] = space[indices_1d[i]];
    }

    // obtain trilinear coordinates for interpolation
    auto trilinear_coordinates = util::trilinear_coordinates(point, space_cube);
    glm::vec3 output_gravity{};
    // interpolation
    for(int i = 0; i < 8; i++) {
        output_gravity +=  trilinear_coordinates[i] * gravity[indices_1d[i]];
    }
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

float gravity::volume(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        const glm::vec3& t) {
    float volume = 0.f;
    for(auto & face : faces) {
        tetrahedron tetrahedron = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                t};
        volume += util::tetrahedron_volume(tetrahedron.b1, tetrahedron.b2, tetrahedron.b3, tetrahedron.v) * 10;
    }
    return volume;
}

void gravity::build_octree(
    float precision,
    std::vector<node>& octree,
    int id,
    int max_res,
    glm::vec3 min,
    float edge,
    const std::vector<glm::vec3>& gravity,
    const std::vector<glm::vec3>& space, int resolution
    ) {
    if(should_divide(precision, octree[id], min, edge, gravity, space, resolution) && max_res > 0) {
        // inspecting node became internal node -> it gets children, and first_child_id must be updated with first
        // child id; the children are all adjacent; then delete gravity_octant, since it is no longer needed
        octree[id].first_child_id = static_cast<int>(octree.size());
        delete octree[id].gravity_octant;
        octree[id].gravity_octant = nullptr;

        max_res--;
        auto min_box = util::get_box(min, edge/2.0f);
        for(int i = 0; i < 8; i++) {
            octree.push_back(build_node(min_box[i], edge/2.0f, gravity, space, resolution));
        }

        for(int i = 0; i < 8; i++) {
            build_octree(precision, octree, octree[id].first_child_id + i, max_res, min_box[i], edge/2.0f, gravity, space, resolution);
        }
    }
}

void gravity::build_octree_with_integral(
    float precision,
    std::vector<node>& octree,
    int id,
    int max_res,
    glm::vec3 min,
    float edge,
    const std::vector<tube>& tubes,
    float G, float R
    ) {
    if(should_divide_with_integral(precision, octree[id], min, edge, tubes, G, R) && max_res > 0) {
        // inspecting node became internal node -> it gets children, and first_child_id must be updated with first
        // child id; the children are all adjacent; then delete gravity_octant, since it is no longer needed
        octree[id].first_child_id = static_cast<int>(octree.size());
        delete octree[id].gravity_octant;
        octree[id].gravity_octant = nullptr;

        max_res--;
        auto min_box = util::get_box(min, edge/2.0f);
        for(int i = 0; i < 8; i++) {
            octree.push_back(build_node_with_integral(min_box[i], edge/2.0f, tubes, G, R));
        }

        for(int i = 0; i < 8; i++) {
            build_octree_with_integral(precision, octree, octree[id].first_child_id + i, max_res, min_box[i], edge/2.0f, tubes, G, R);
        }
    }
}

void gravity::build_octree_with_integral_optimized(
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
            gravity_values_map.emplace(int_location[i], tmp_gravity_values.size());
            //tmp_gravity_values.push_back(gravity_values[new_value_index]);
            //cached_values.emplace(int_location[i], new_value_index);
            octree.push_back(-new_value_index);
        }
    }

    if(max_res > 0 && should_divide_with_integral_optimized(precision, id, octree, max_res,min, edge, int_min, int_edge, gravity_values, tmp_gravity_values, cached_values, gravity_values_map, tubes, G, R)) {
        // inspecting node became internal node -> it gets children, and first_child_id must be updated with first
        // child id; the children are all adjacent; then delete gravity_octant, since it is no longer needed

        auto mins = util::get_box(min, edge/2.f);
        auto int_mins = util::get_int_box(int_min, int_edge/2);

        for(int i = 0; i < 8; i++) {
            octree[id + i] = (int)octree.size();
            build_octree_with_integral_optimized(
                    precision, octree, octree[id + i], max_res - 1, mins[i],
                    edge/2.f, int_mins[i],
                    int_edge / 2, gravity_values, tmp_gravity_values, cached_values, gravity_values_map, tubes, G, R
                    );
        }
    }
}
