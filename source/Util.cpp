//
// Created by Giovanni Bollati on 15/05/23.
//

#include <Util.hpp>

glm::quat util::rotationBetweenVectors(const glm::vec3& start, const glm::vec3& dest) {
    // normalize vectors
    glm::vec3 start_n = glm::normalize(start);
    glm::vec3 dest_n = glm::normalize(dest);

    float cosTheta = dot(start_n, dest_n);
    glm::vec3 rotationAxis;

    // vectors are in opposite directions
    if (cosTheta < -1 + 0.001f){
        // arbitrary roation axis perpendicular to start
        rotationAxis = cross(glm::vec3(0.0f, 0.0f, 1.0f), start);

        // arbitrary was parallel to start, retry with another arbitrary axis
        if (glm::length2(rotationAxis) < 0.01 )
            rotationAxis = cross(glm::vec3(1.0f, 0.0f, 0.0f), start);

        // normalize rotation axis
        rotationAxis = normalize(rotationAxis);
        return glm::angleAxis(glm::radians(180.0f), rotationAxis);
    }

    rotationAxis = cross(start, dest);

    float s = sqrt( (1+cosTheta)*2 );

    // formula justified by rotationAxis length being sin(theta)
    // q[w, tx, ty, tz], w = cos(theta/2), t = sin(theta/2)
    float invs = 1 / s;

    return {
            s * 0.5f,
            rotationAxis.x * invs,
            rotationAxis.y * invs,
            rotationAxis.z * invs
    };
}

glm::vec3 util::barycentric_coords(const glm::vec3& a, const glm::vec3& b, const glm::vec3& c, const glm::vec3& p) {
    glm::vec3 v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = glm::dot(v0, v0);
    float d01 = glm::dot(v0, v1);
    float d11 = glm::dot(v1, v1);
    float d20 = glm::dot(v2, v0);
    float d21 = glm::dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return {u, v, w};
}