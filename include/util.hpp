//
// Created by Giovanni Bollati on 15/05/23.
//

#ifndef GL_TEST_PROJECT_UTIL_HPP
#define GL_TEST_PROJECT_UTIL_HPP

// std
#include <vector>
#include <array>

// glm
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

namespace util {
    glm::quat rotation_between_vectors(const glm::vec3& start, const glm::vec3& dest);

    // cramer rule used (Ax = b, x_i = det(A_i) / det(A), A_i = A with i-column replaced with b
    glm::vec3 barycentric_coords(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& p);

    float point_edge_distance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2);
    float point_triangle_distance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2, const glm::vec3& e3);

    float tetrahedron_volume(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);
    [[deprecated]] glm::vec3 tetrahedron_barycentre(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);

    bool ray_triangle_intersection(glm::vec3 ray_origin, glm::vec3 ray_dir, glm::vec3 t1, glm::vec3 t2, glm::vec3 t3, float* parameter);
    [[maybe_unused]] std::vector<glm::vec3> ray_mesh_intersections(const std::vector<glm::vec3>& vertices,
                                                                                     const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                                     glm::vec3 ray_origin, glm::vec3 ray_dir
    );
    std::vector<glm::vec3> ray_mesh_intersections_optimized(const std::vector<glm::vec3>& vertices,
                                                                              const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                              glm::vec3 ray_origin, glm::vec3 ray_dir
    );

    // returns a point which has the minimum coordinate along all the three axis
    glm::vec3 get_min(const std::vector<glm::vec3>& vertices);
    // returns a point which has the maximum coordinate along the three axis
    glm::vec3 get_max(const std::vector<glm::vec3>& vertices);

    // returns mesh max extent on xy plane
    inline float xy_mesh_max_extent(const std::vector<glm::vec3>& vertices) {
        return std::max(
            util::get_max(vertices).x - util::get_min(vertices).x,
            util::get_max(vertices).y - util::get_min(vertices).y
            );
    }

    // return an array of 4 floats; 3 first values are coordinates of the center of the cube
    // the fourth float is the edge length
    inline std::array<float, 4> get_box(const glm::vec3& min, const glm::vec3& max) {
        auto center = (max + min) * 0.5f;
        float edge = std::max(max.x - min.x, max.y - min.y);
        edge = std::max(edge, max.z - min.z);
        return {center.x, center.y, center.z, edge};
    }

    // implicit definition of box as std::array<glm::vec3, 8>
    inline std::array<glm::vec3, 8> get_box(const glm::vec3& min, float edge) {
        return {
            glm::vec3{min.x, min.y, min.z},
            glm::vec3{min.x + edge, min.y, min.z},
            glm::vec3{min.x, min.y + edge, min.z},
            glm::vec3{min.x + edge, min.y + edge, min.z},
            glm::vec3{min.x, min.y, min.z + edge},
            glm::vec3{min.x + edge, min.y, min.z + edge},
            glm::vec3{min.x, min.y + edge, min.z + edge},
            glm::vec3{min.x + edge, min.y + edge, min.z + edge}
        };
    }

    // box: first 3 floats are the center of the box, fourth is the edge
    inline std::array<glm::vec3, 8> get_box(const std::array<float, 4>& box) {
        return get_box(
            glm::vec3{box[0] - box[3] / 2.f, box[1] - box[3] / 2.f, box[2] - box[3] / 2.f},
            box[3]
            );
    }

    // allows to select a bounding box of a point from a 3d space vector
    std::array<std::array<int, 3>, 8> get_box_indices(glm::vec3 min, float range, int resolution, glm::vec3 point);

    inline int from_3d_indices_to_1d(std::array<int, 3> indices, int resolution) {
        return indices[0] * (resolution + 1) * (resolution + 1) + indices[1] * (resolution + 1) + indices[2];
    }

    // cube geometry
    // pre requisite: p is inside cube
    std::array<float, 8> trilinear_coordinates(const glm::vec3& p, const std::array<glm::vec3, 8>& box);

    inline glm::vec3 interpolate(const glm::vec3& p, const std::array<glm::vec3, 8>& box, const std::array<glm::vec3, 8> values) {
        auto weights = trilinear_coordinates(p, box);
        auto output = glm::vec3{};
        for(int i = 0; i < 8; i++) { output += values[i]*weights[i]; }
        return output;
    }

    // absolute comparison
    // precision -> max length of (v1 - v2 vector) allowed
    inline bool vectors_are_p_equal(glm::vec3 v1, glm::vec3 v2, float precision) {
        return glm::length(v1 - v2) < precision;
    }

    bool is_inside_box(const glm::vec3& p, const std::array<glm::vec3, 8>& box);
    bool is_box_inside_mesh(std::array<glm::vec3, 8> box, const std::vector<glm::vec3>& vertices, const std::vector<glm::vec<3, unsigned int>>& faces);
};

#endif //GL_TEST_PROJECT_UTIL_HPP
