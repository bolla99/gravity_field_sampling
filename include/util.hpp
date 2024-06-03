//
// Created by Giovanni Bollati on 15/05/23.
//

#ifndef GL_TEST_PROJECT_UTIL_HPP
#define GL_TEST_PROJECT_UTIL_HPP

// std
#include <vector>
#include <array>
#include <iostream>

// glm
#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>

namespace util {
    // returns quaternion that represents rotations from start vector to dest vector
    glm::quat rotation_between_vectors(const glm::vec3& start, const glm::vec3& dest);

    // cramer rule used (Ax = b, x_i = det(A_i) / det(A), A_i = A with i-column replaced with b
    glm::vec3 barycentric_coords(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& p);

    // returns distance between a point and an edge {e1, e2}
    float point_edge_distance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2);

    // returns the distance between a point and a triangle {e1, e2, e3}
    // DEPENDS on point_edge_distance
    float point_triangle_distance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2, const glm::vec3& e3);

    // returns volume of a tetrahedron {b1, b2, b3, v} assuming it has homogeneous density
    float tetrahedron_volume(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);

    // returns barycentre of a tetrahedron {b1, b2, b3, v}
    [[deprecated]] glm::vec3 tetrahedron_barycentre(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);

    // returns true if ray {ray_origin, ray_dir} intersects the triangle {t1, t2, t3}
    // if yes the distance between ray_origin and the intersection is stored in parameter
    // möller-trumbore + cramer's rule
    bool ray_triangle_intersection(glm::vec3 ray_origin, glm::vec3 ray_dir, glm::vec3 t1, glm::vec3 t2, glm::vec3 t3, float* parameter);

    // returns a vector of intersections between the ray {ray_origin, ray_dir} and multiple triangles
    // triangles are retrieved from vertices and faces collection: faces is a collection of triple of integers, which are
    // the indexes for the vertices collection, which contains actual triangles points.
    // DEPENDS on ray_triangle_intersection
    [[maybe_unused]] std::vector<glm::vec3> ray_mesh_intersections(const std::vector<glm::vec3>& vertices,
                                                                                     const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                                     glm::vec3 ray_origin, glm::vec3 ray_dir
    );

    // same as ray_mesh_intersections, but optimizations are done regarding how many and which triangles
    // are chosen for testing; instead base function tests all the triangles.
    // DEPENDS on ray_triangle_intersection
    std::vector<glm::vec3> ray_mesh_intersections_optimized(const std::vector<glm::vec3>& vertices,
                                                                              const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                              glm::vec3 ray_origin, glm::vec3 ray_dir
    );

    // returns a point which has the minimum coordinate along all the three axis
    glm::vec3 get_min(const std::vector<glm::vec3>& vertices);

    // returns a point which has the maximum coordinate along the three axis
    glm::vec3 get_max(const std::vector<glm::vec3>& vertices);

    // returns mesh max between x and y extent
    inline float xy_mesh_max_extent(const std::vector<glm::vec3>& vertices) {
        return std::max(
            util::get_max(vertices).x - util::get_min(vertices).x,
            util::get_max(vertices).y - util::get_min(vertices).y
            );
    }

    // return {c.x, c.y, c.z, edge} where c is the center of the cube
    // given mix and max vertices of the cube
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

    // integer version of get_box function
    inline std::array<glm::vec<3, int>, 8> get_int_box(const glm::vec<3, int>& min, int edge) {
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

    // returns a box in 8 vertices representation from the {center, edge} representation
    // DEPENDS on get_box(min, edge)
    inline std::array<glm::vec3, 8> get_box(const std::array<float, 4>& box) {
        return get_box(
            glm::vec3{box[0] - box[3] / 2.f, box[1] - box[3] / 2.f, box[2] - box[3] / 2.f},
            box[3]
            );
    }

    // returns the integer coordinated of the bounding box
    // of the point relative to the box {min, range} and given resolution;
    // resolution stands for number of edges in which the range is subdivided
    // (range / resolution is performed to get sub boxes edge length)
    std::array<std::array<int, 3>, 8> get_box_indices(glm::vec3 min, float range, int resolution, glm::vec3 point);

    // returns 1d coordinates given integer indices and resolution;
    // resolution stands for number of edges
    // the function implicitly defines the 1d vector representing a 3d grid
    inline int from_3d_indices_to_1d(std::array<int, 3> indices, int resolution) {
        return indices[0] * (resolution + 1) * (resolution + 1) + indices[1] * (resolution + 1) + indices[2];
    }

    // if p is inside the box, it RETURNS the barycentric coordinates
    std::array<float, 8> barycentric_coords(const glm::vec3& p, const std::array<glm::vec3, 8>& box);

    // RETURNS the interpolation of p given box and values
    // DEPENDS on util::barycentric_coords(box)
    inline glm::vec3 interpolate(const glm::vec3& p, const std::array<glm::vec3, 8>& box, const std::array<glm::vec3, 8> values) {
        auto weights = barycentric_coords(p, box);
        auto output = glm::vec3{};
        for(int i = 0; i < 8; i++) { output += values[i]*weights[i]; }
        return output;
    }

    // get gradient in p given box and values; derivative is computed for each of
    // the four edge of a given axis and then the axis component of gradient is retrieved
    // from these four derivative
    // DEPENDS on get_x_derivative_from_cube, get_y.., get_z..
    inline glm::vec3 get_gradient_from_box(const glm::vec3& p, const std::array<glm::vec3, 8>& locations, const std::array<float, 8>& values, float edge) {
        //auto edge = glm::distance(locations[0], locations[1]);

        auto gradient = glm::vec3{0.f, 0.f, 0.f};
        // x axis; interpolazione fatta sui valori y e z
        // 1 - 0;  3 - 2; 5 - 4; 7 - 6
        std::array<float, 4> x_values = {values[1] - values[0], values[3] - values[2], values[5] - values[4], values[7] - values[6]};
        float first_c = p.y - locations[0].y;
        float second_c = p.z - locations[0].z;
        std::array<float, 4> weights = {
                (edge - first_c) * (edge - second_c),
                first_c * (edge - second_c),
                (edge - first_c) * second_c,
                first_c * second_c
        };

        for(int i = 0; i < 4; i++) {
            gradient.x -= float(double(x_values[i]) * (double(weights[i]) / std::pow(edge, 3)));
        }

        // y axis x e z
        // 2 - 0; 3 - 1; 6 - 4; 7 - 5
        std::array<float, 4> y_values = {values[2] - values[0], values[3] - values[1], values[6] - values[4], values[7] - values[5]};
        first_c = p.x - locations[0].x;
        second_c = p.z - locations[0].z;
        weights = {
                (edge - first_c) * (edge - second_c),
                first_c * (edge - second_c),
                (edge - first_c) * second_c,
                first_c * second_c
        };

        for(int i = 0; i < 4; i++) {
            gradient.y -= float(double(y_values[i]) * (double(weights[i]) / std::pow(edge, 3)));
        }
        // z axis
        // 4 - 0; 5 - 1; 6 - 2; 7 - 3
        std::array<float, 4> z_values = {values[4] - values[0], values[5] - values[1], values[6] - values[2], values[7] - values[3]};
        first_c = p.x - locations[0].x;
        second_c = p.y - locations[0].y;
        weights = {
                (edge - first_c) * (edge - second_c),
                first_c * (edge - second_c),
                (edge - first_c) * second_c,
                first_c * second_c
        };

        for(int i = 0; i < 4; i++) {
            gradient.z -= float(double(z_values[i]) * (double(weights[i]) / std::pow(edge, 3)));;
        }
        return gradient;
    }

    // RETURNS true if (v1 - v2) vector length is less then precision
    inline bool vectors_are_p_equal(const glm::vec3& v1, const glm::vec3& v2, float precision) {
        auto v3 = v1;
        // avoid comparison with 0
        /*if(glm::length(v1) < 1) {
            if(glm::length(v2) < 1) return true;
            v3 = v2;
        } else if(glm::length(v2) < 1) {
            v3 = v1;
        }*/
        auto p = glm::length2(v1 - v2)/glm::length2(v3)*100;
        return p < precision;
    }

    // RETURNS edge derivatives for x axis given values
    inline std::array<float, 4> get_x_derivative_from_cube(const std::array<float, 8>& values) {
        return std::array<float, 4>{
            values[1] - values[0],
            values[3] - values[2],
            values[5] - values[4],
            values[7] - values[6]
            };
    }

    // RETURNS edge derivatives for y axis given values
    inline std::array<float, 4> get_y_derivative_from_cube(const std::array<float, 8>& values) {
        return std::array<float, 4>{
                values[2] - values[0],
                values[3] - values[1],
                values[6] - values[4],
                values[7] - values[5]
        };
    }

    // RETURNS edge derivatives for z axis given values
    inline std::array<float, 4> get_z_derivative_from_cube(const std::array<float, 8>& values) {
        return std::array<float, 4>{
                values[4] - values[0],
                values[5] - values[1],
                values[6] - values[2],
                values[7] - values[3]
        };
    }

    // RETURNS true is p is inside box
    inline bool is_inside_box(const glm::vec3& p, const std::array<glm::vec3, 8>& box) {
        return p.x >= box[0].x && p.y >= box[0].y && p.z >= box[0].z &&
               p.x <= box[7].x && p.y <= box[7].y && p.z <= box[7].z;
    }

    // RETURNS true if box is inside a mesh (vertices + faces)
    // DEPENDS on ray_mesh_intersections
    bool is_box_inside_mesh(std::array<glm::vec3, 8> box, const std::vector<glm::vec3>& vertices, const std::vector<glm::vec<3, unsigned int>>& faces);

    // computes mesh volume with tetrahedrons method
    float volume(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            const glm::vec3& t);

    inline glm::vec3 norm(const std::vector<glm::vec3>& vertices, glm::vec<3, int> face) {
        auto e1 = vertices[face.y] - vertices[face.x];
        auto e2 = vertices[face.z] - vertices[face.x];
        return glm::normalize(glm::cross(e1, e2));
    }

    std::vector<glm::vec3> random_locations(glm::vec3 min, float edge, int n);

    std::vector<glm::vec3> near_mesh_locations(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces
            );
    std::vector<glm::vec3> outside_mesh_locations(
            const std::vector<glm::vec3>& vertices,
            const std::vector<glm::vec<3, unsigned int>>& faces,
            glm::vec3 min, float edge
            );
};

#endif //GL_TEST_PROJECT_UTIL_HPP
