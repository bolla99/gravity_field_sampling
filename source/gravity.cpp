//
// Created by Giovanni Bollati on 27/10/23.
//

#include "gravity.hpp"
#include "util.hpp"
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <SDL.h>

std::vector<gravity::tube> gravity::get_tubes(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};

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
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(cube_edge, resolution, ray_dir, tubes, min, vertices, faces)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, r.origin, r.dir);

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


std::vector<gravity::mass> gravity::get_masses(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution,
        float* mass_R
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
    *mass_R = (float)std::cbrt(std::pow(cube_edge, 3) * (3.f / 4.f) / M_PI);

    auto cube_volume = (float)std::pow(cube_edge, 3);
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(resolution, cube_edge, cube_volume, ray_dir, tubes, min, volumes, vertices, faces, mass_R)
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
                *mass_R = std::min(*mass_R, (float)std::cbrt(std::pow(tube_unit, 3) * (3.f / 4.f) / M_PI));
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

glm::vec3 gravity::get_gravity_from_tubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point) {
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
    std::cout << "gravity tubes time elapsed " << (float)(SDL_GetTicks64() - start) / 1000.f << std::endl;
    glm::vec3 gravity = {0, 0, 0};
    for(int i = 0; i < omp_get_max_threads(); i++) {
        gravity = gravity + thread_gravity[i];
    }
    return gravity;
}

glm::vec3 gravity::get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, float R, glm::vec3 point) {
    Uint64 start = SDL_GetTicks64();
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }
#pragma omp parallel for default(none) shared(masses, point, G, thread_gravity, R, std::cout)
    for(auto & mass : masses) {
        glm::vec3 dir = mass.p - point;
        auto r3 = (float)pow(glm::length(dir), 3);
        if(glm::length(dir) < R) {
#pragma omp critical
            std::cout << "mass too near " << mass.p.x << " " << mass.p.y << " " << mass.p.z << std::endl;
            r3 = (float)std::pow(R, 3);
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

    /*
    std::cout << "indices debug" << std::endl;
    for(int i = 0; i < 8; i++) {
        std::cout << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << std::endl;
    }*/

    // translate 3d indices to 1d indices
    std::array<int, 8> indices_1d{};
    for(int i = 0; i < 8; i++) {
        indices_1d[i] = util::from_3d_indices_to_1d(indices[i], resolution);
    }

    /*
    std::cout << "1d indices debug" << std::endl;
    for(int i = 0; i < 8; i++) {
        std::cout << indices_1d[i] << std::endl;
    }
    */

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
        output_gravity += gravity[indices_1d[i]] * trilinear_coordinates[i];
    }
    return output_gravity;
}

std::vector<glm::vec3> gravity::get_discrete_space(glm::vec3 min, glm::vec3 max, int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 center = (max + min) * 0.5f;
    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float z_width = max.z - min.z;
    float max_extent = x_width;
    if (y_width > max_extent) max_extent = y_width;
    if (z_width > max_extent) max_extent = z_width;

    min = center - (max_extent / 2.0f);

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

glm::vec3 gravity::get_gravity_RT(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution, glm::vec3 point) {
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

#pragma omp parallel for default(none), shared(resolution, cube_edge, ray_dir, min, cube_volume, point, thread_gravity, G, vertices, faces)
    for(int i = 0; i < resolution + 1; i++) {
        // private vars for thread
        float units_per_tubef, tube_unit;
        float new_cube_volume, r3;
        int units_per_tube;
        gravity::tube t{};

        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            gravity::ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, r.origin, r.dir);

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

int gravity::build_octree(
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
    ) {
    if(!util::is_box_inside_mesh(util::get_box(min, edge), vertices, faces)
        && should_divide(precision, octree[id], min, edge, gravity, space, resolution)
        && max_res > 0
        ) {
        max_res--;
        auto min_box = util::get_box(min, edge/2.0f);
        for(int i = 0; i < 8; i++) {
            octree.push_back(gravity::build_node(min_box[i], edge/2.0f, gravity, space, resolution));
        }
        int new_id = next_id;
        int new_next_id = next_id + 8;
        for(int i = 0; i < 8; i++) {
            new_next_id = build_octree(precision, octree, new_id + i, new_next_id, max_res, min_box[i], edge/2.0f, gravity, space, resolution, vertices, faces);
        }
    }
    return next_id;
}
