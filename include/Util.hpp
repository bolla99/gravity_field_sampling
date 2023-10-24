//
// Created by Giovanni Bollati on 15/05/23.
//

#ifndef GL_TEST_PROJECT_UTIL_HPP
#define GL_TEST_PROJECT_UTIL_HPP

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/norm.hpp>
namespace util {
    glm::quat rotationBetweenVectors(const glm::vec3& start, const glm::vec3& dest);

    // cramer rule used (Ax = b, x_i = det(A_i) / det(A), A_i = A with i-column replaced with b
    glm::vec3 barycentric_coords(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& p);
};

#endif //GL_TEST_PROJECT_UTIL_HPP
