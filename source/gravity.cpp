//
// Created by Giovanni Bollati on 27/10/23.
//

#include "gravity.hpp"
#include "util.hpp"
#include <iostream>
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
        int resolution) {
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
    auto cube_volume = (float)std::pow(cube_edge, 3);
    glm::vec3 ray_dir = {0.f, 0.f, 1.f};

#pragma omp parallel for default(none), shared(resolution, cube_edge, cube_volume, ray_dir, tubes, min, volumes, vertices, faces)
    for(int i = 0; i < resolution + 1; i++) {
        for(int j = 0; j < resolution + 1; j++) {
            // RAY
            ray r = {{min.x + (float)i * cube_edge, min.y + (float)j * cube_edge, min.z}, ray_dir};
            // FIND INTERSECTIONS
            std::vector<glm::vec3> intersections = util::ray_mesh_intersections_optimized(vertices, faces, r.origin, r.dir);

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

glm::vec3 gravity::get_gravity_from_masses(const std::vector<gravity::mass>& masses, float G, glm::vec3 point) {
    Uint64 start = SDL_GetTicks64();
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }
#pragma omp parallel for default(none) shared(masses, point, G, thread_gravity)
    for(auto & mass : masses) {
        glm::vec3 dir = mass.p - point;
        auto r3 = (float)pow(glm::length(dir), 3);
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

glm::vec3 gravity::get_gravity_from_1D_precomputed_vector(glm::vec3 point, const std::vector<glm::vec3>& gravity, const std::vector<glm::vec3>& space, glm::vec3 min, float range, int resolution) {
    // get indices of bounding box
    auto indices = util::get_box_indices(min, range, resolution, point);
    std::cout << "INDICES:" << std::endl;
    for(int i = 0; i < 8; i++) {
        std::cout << indices[i][0] << " " << indices[i][1] << " " << indices[i][2] << std::endl;
    }

    // translate 3d indices to 1d indices
    std::array<int, 8> indices_1d{};
    for(int i = 0; i < 8; i++) {
        indices_1d[i] = util::from_3d_indices_to_1d(indices[i], resolution);
    }
    std::cout << "indices 1d debug" << std::endl;
    for(int i = 0; i < 8; i++) {
        std::cout << indices_1d[i] << std::endl;
    }

    // get space coordinates of bounding box
    std::array<glm::vec3, 8> space_cube{};
    for(int i = 0; i < 8; i++) {
        space_cube[i] = space[indices_1d[i]];
    }
    for(int i = 0; i < 8; i++) {
        std::cout << "debug cube positions" << std::endl;
        std::cout << space_cube[i].x << " " << space_cube[i].y << " " << space_cube[i].z << " " << std::endl;
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

octree<gravity::gravity_cube>* gravity::get_gravity_octree_from_masses(
        glm::vec3 min, glm::vec3 max, int resolution, const std::vector<gravity::mass>& masses) {
    // get bounding box
    auto _cube = util::get_box(min, max);
    glm::vec3 center{_cube[0], _cube[1], _cube[2]};
    auto edge_length = _cube[3];
    min = center - (edge_length / 2.0f);

    // set octree root element
    gravity::gravity_cube gc = {{center, edge_length / 2.0f }, std::array<glm::vec3, 8>{}};
    gc.g[0] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y, min.z});
    gc.g[1] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y, min.z});
    gc.g[2] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y + edge_length, min.z});
    gc.g[3] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y + edge_length, min.z});
    gc.g[4] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y, min.z + edge_length});
    gc.g[5] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y, min.z + edge_length});
    gc.g[6] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y + edge_length, min.z + edge_length});
    gc.g[7] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y + edge_length, min.z + edge_length});

    // f function
    auto f = [masses](gravity::gravity_cube gc)->std::array<gravity::gravity_cube, 8> {
        float new_edge = gc.c.extent / 2.0f;
        glm::vec3 min = {gc.c.center.x - new_edge, gc.c.center.y - new_edge, gc.c.center.z - new_edge};
        std::array<gravity::cube, 8> new_cubes{};
        new_cubes[0] = {{min.x, min.y, min.z}, new_edge};
        new_cubes[1] = {{min.x + gc.c.extent, min.y, min.z}, new_edge};
        new_cubes[2] = {{min.x, min.y + gc.c.extent, min.z}, new_edge};
        new_cubes[3] = {{min.x + gc.c.extent, min.y + gc.c.extent, min.z}, new_edge};
        new_cubes[4] = {{min.x, min.y, min.z + gc.c.extent}, new_edge};
        new_cubes[5] = {{min.x + gc.c.extent, min.y, min.z + gc.c.extent}, new_edge};
        new_cubes[6] = {{min.x, min.y + gc.c.extent, min.z + gc.c.extent}, new_edge};
        new_cubes[7] = {{min.x + gc.c.extent, min.y + gc.c.extent, min.z + gc.c.extent}, new_edge};
        std::array<gravity::gravity_cube, 8> new_gravity_cubes{};
        for(int i = 0; i < 8; i++) {
            new_gravity_cubes[i].c = new_cubes[i];
            auto edge_length = new_gravity_cubes[i].c.extent;
            min = new_gravity_cubes[i].c.center - edge_length / 2.0f;
            new_gravity_cubes[i].g[0] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y, min.z});
            new_gravity_cubes[i].g[1] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y, min.z});
            new_gravity_cubes[i].g[2] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y + edge_length, min.z});
            new_gravity_cubes[i].g[3] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y + edge_length, min.z});
            new_gravity_cubes[i].g[4] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y, min.z + edge_length});
            new_gravity_cubes[i].g[5] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y, min.z + edge_length});
            new_gravity_cubes[i].g[6] = gravity::get_gravity_from_masses(masses, 10.0, {min.x, min.y + edge_length, min.z + edge_length});
            new_gravity_cubes[i].g[7] = gravity::get_gravity_from_masses(masses, 10.0, {min.x + edge_length, min.y + edge_length, min.z + edge_length});
        }
        return new_gravity_cubes;
    };

    // condition
    auto condition = [](std::array<gravity::gravity_cube, 8>)->bool {
        return true;
    };

    // create octree and call execute to construct it
    auto root = new octree<gravity::gravity_cube>(gc);
    root->execute(resolution, gc, f, condition);
    return root;
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

octree<gravity::cube>* gravity::get_discrete_space_as_octree(glm::vec3 min, glm::vec3 max, int resolution) {
    auto _cube = util::get_box(min, max);

    gravity::cube cube = {{_cube[0], _cube[1], _cube[2]}, _cube[3] / 2.0f };
    auto f = [](gravity::cube c)->std::array<gravity::cube, 8>{
        float new_extent = c.extent / 2.0f;
        glm::vec3 min = {c.center.x - new_extent, c.center.y - new_extent, c.center.z - new_extent};
        std::array<gravity::cube, 8> new_cubes{};
        new_cubes[0] = {{min.x, min.y, min.z}, new_extent};
        new_cubes[1] = {{min.x + c.extent, min.y, min.z}, new_extent};
        new_cubes[2] = {{min.x, min.y + c.extent, min.z}, new_extent};
        new_cubes[3] = {{min.x + c.extent, min.y + c.extent, min.z}, new_extent};
        new_cubes[4] = {{min.x, min.y, min.z + c.extent}, new_extent};
        new_cubes[5] = {{min.x + c.extent, min.y, min.z + c.extent}, new_extent};
        new_cubes[6] = {{min.x, min.y + c.extent, min.z + c.extent}, new_extent};
        new_cubes[7] = {{min.x + c.extent, min.y + c.extent, min.z + c.extent}, new_extent};
        return new_cubes;
    };
    auto condition = [](std::array<gravity::cube, 8>)->bool {
        return true;
    };
    auto root = new octree<gravity::cube>(cube);
    root->execute(resolution, cube, f, condition);
    return root;
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

glm::vec3 gravity::get_gravity_from_tetrahedrons(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) {
    glm::vec3 gravity = {0.f, 0.f, 0.f};
    for(auto & face : faces) {
        tetrahedron t = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                tetrahedrons_vertex};
        float mass = util::tetrahedron_volume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedron_barycentre(t.b1, t.b2, t.b3, t.v);
        float distance = glm::length(p - barycentre);
        gravity = gravity - ((p - barycentre)*mass) / (float)pow(distance, 3);
    }
    return gravity;
}

glm::vec3 gravity::get_gravity_from_tetrahedrons_corrected(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        const glm::vec3& p, const glm::vec3& tetrahedrons_vertex) {
    glm::vec3 gravity = {0.f, 0.f, 0.f};
    for(auto & face : faces) {
        tetrahedron t = {
                vertices[face.x - 1],
                vertices[face.y - 1],
                vertices[face.z - 1],
                tetrahedrons_vertex};
        float mass = util::tetrahedron_volume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedron_barycentre(t.b1, t.b2, t.b3, t.v);
        gravity = gravity + ((barycentre - p)*mass) / (float)pow(glm::length(barycentre - p), 3);
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
