//
// Created by Giovanni Bollati on 15/05/23.
//

#ifndef GL_TEST_PROJECT_UTIL_HPP
#define GL_TEST_PROJECT_UTIL_HPP

#include <vector>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/norm.hpp>

namespace util {
    glm::quat rotationBetweenVectors(const glm::vec3& start, const glm::vec3& dest);

    // cramer rule used (Ax = b, x_i = det(A_i) / det(A), A_i = A with i-column replaced with b
    glm::vec3 barycentric_coords(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& p);

    float pointEdgeDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2);
    float pointTriangleDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2, const glm::vec3& e3);
    float tetrahedronVolume(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);
    glm::vec3 tetrahedronBarycentre(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v);

    bool rayTriangleIntersection(glm::vec3 ray_origin, glm::vec3 ray_dir, glm::vec3 t1, glm::vec3 t2, glm::vec3 t3, float* parameter);
    [[maybe_unused]] std::vector<glm::vec3> rayMeshIntersections(const std::vector<glm::vec3>& vertices,
                                                                                     const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                                     glm::vec3 ray_origin, glm::vec3 ray_dir
    );
    std::vector<glm::vec3> rayMeshIntersectionsOptimized(const std::vector<glm::vec3>& vertices,
                                                                              const std::vector<glm::vec<3, unsigned int>>& faces,
                                                                              glm::vec3 ray_origin, glm::vec3 ray_dir
    );

    // returns a point which has the minimum coordinate along all the three axis
    glm::vec3 getMin(const std::vector<glm::vec3>& vertices);
    // returns a point which has the maximum coordinate along the three axis
    glm::vec3 getMax(const std::vector<glm::vec3>& vertices);

    // cube geometry

    // pre requisite: p is inside cube
    std::array<float, 8> trilinearCoordinates(const glm::vec3& p, const std::array<glm::vec3, 8>& cube);
    bool isInsideCube(const glm::vec3& p, const std::array<glm::vec3, 8>& cube);
};

#endif //GL_TEST_PROJECT_UTIL_HPP
