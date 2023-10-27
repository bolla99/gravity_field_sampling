//
// Created by Giovanni Bollati on 27/10/23.
//

#include "gravity.hpp";
#include "util.hpp"
#include <iostream>
#include <omp.h>
#include <SDL.h>


std::vector<gravity::tube> gravity::getTubes(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = util::getMin(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::getMax(vertices);

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
            std::vector<glm::vec3> intersections = util::rayMeshIntersectionsOptimized(vertices, faces, r.origin, r.dir);

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


std::vector<gravity::mass> gravity::getMasses(
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::vec<3, unsigned int>>& faces,
        int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    std::vector<tube> tubes = {};
    std::vector<mass> volumes = {};

    Uint64 start = SDL_GetTicks64();

    // FIND XY PLANE
    glm::vec3 min = util::getMin(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::getMax(vertices);

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
            std::vector<glm::vec3> intersections = util::rayMeshIntersectionsOptimized(vertices, faces, r.origin, r.dir);

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

glm::vec3 gravity::getGravityFromTubes(const std::vector<glm::vec3>& vertices, int resolution, const std::vector<gravity::tube>& tubes, glm::vec3 point) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 thread_gravity[omp_get_max_threads()];
    for(int i = 0; i < omp_get_max_threads(); i++) {
        thread_gravity[i] = {0, 0, 0};
    }

    Uint64 start = SDL_GetTicks64();
    // FIND XY PLANE
    glm::vec3 min = util::getMin(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::getMax(vertices);

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

glm::vec3 gravity::getGravityFromMasses(const std::vector<gravity::mass>& masses, float G, glm::vec3 point) {
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

std::vector<glm::vec3> gravity::getDiscreteSpace(glm::vec3 min, glm::vec3 max, int resolution) {
    omp_set_num_threads(omp_get_max_threads());
    glm::vec3 center = (max + min) * 0.5f;
    float x_width = max.x - min.x;
    float y_width = max.y - min.y;
    float z_width = max.z - min.z;
    float max_extent = x_width;
    if (y_width > max_extent) max_extent = y_width;
    if (z_width > max_extent) max_extent = z_width;

    min = center - (max_extent / 2.0f);
    max = center + (max_extent / 2.0f);

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

glm::vec3 gravity::getGravityRT(
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
    glm::vec3 min = util::getMin(vertices);
    min.x += 0.0001; min.y += 0.0001;
    glm::vec3 max = util::getMax(vertices);

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
            std::vector<glm::vec3> intersections = util::rayMeshIntersectionsOptimized(vertices, faces, r.origin, r.dir);

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

glm::vec3 gravity::getGravityFromTetrahedrons(
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
        float mass = util::tetrahedronVolume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedronBarycentre(t.b1, t.b2, t.b3, t.v);
        float distance = glm::length(p - barycentre);
        gravity = gravity - ((p - barycentre)*mass) / (float)pow(distance, 3);
    }
    return gravity;
}

glm::vec3 gravity::getGravityFromTetrahedronsCorrected(
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
        float mass = util::tetrahedronVolume(t.b1, t.b2, t.b3, t.v);
        glm::vec3 barycentre = util::tetrahedronBarycentre(t.b1, t.b2, t.b3, t.v);
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
        volume += util::tetrahedronVolume(tetrahedron.b1, tetrahedron.b2, tetrahedron.b3, tetrahedron.v) * 10;
    }
    return volume;
}
