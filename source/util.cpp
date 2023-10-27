//
// Created by Giovanni Bollati on 15/05/23.
//

#include <util.hpp>

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

float util::pointEdgeDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2) {
    //std::cout << "Running point edge distance for point: " << point.toString() << " and edge: " << e1.toString() << " - " << e2.toString() << "\n";
    glm::vec3 e1e2 = e2 - e1;
    glm::vec3 e1point = point - e1;
    float e1e2mag2 = glm::length2(e1e2);
    float pos = glm::dot(e1e2, e1point) / e1e2mag2;

    // CLAMP
    if(pos < 0) pos = 0;
    if(pos > 1) pos = 1;
    //std::cout << "point edge distance " << (point - (e2*pos + e1*(1-pos))).mag() << "\n";
    return glm::length(point - (e2*pos + e1*(1-pos)));
}

float util::pointTriangleDistance(const glm::vec3& point, const glm::vec3& e1, const glm::vec3& e2, const glm::vec3& e3) {
    //std::cout << "Running point triangle distance for face with vertices: " << e1.toString() << e2.toString() << e3.toString() << "\n";
    glm::vec3 e12 = e2 - e1;
    glm::vec3 e13 = e3 - e1;
    glm::vec3 e1point = point - e1;
    glm::vec3 n = glm::cross(e12, e13);
    float nMag2 = glm::length2(n);

    // barycentric coordinates
    float alpha = glm::dot(glm::cross(e12, e1point), n) / nMag2;
    float beta = glm::dot(glm::cross(e1point, e13), n) / nMag2;
    float gamma = 1 - alpha - beta;

    // point projected on trinangle plane
    glm ::vec3 pointProjected = e1*alpha + e2*beta + e3*gamma;

    float semispace = glm::dot((point - pointProjected), n);
    float semispace_sign;
    if(semispace >= 0) semispace_sign = 1;
    if(semispace < 0) semispace_sign = -1;
    if(alpha >= 0 && beta >= 0 && gamma >= 0) return glm::length(point - pointProjected) * semispace_sign;
    else if(alpha < 0) return pointEdgeDistance(point, e2, e3) * semispace_sign;
    else if(beta < 0) return pointEdgeDistance(point, e1, e3) * semispace_sign;
    else return pointEdgeDistance(point, e1, e2) * semispace_sign;
}

float util::tetrahedronVolume(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v) {
    return glm::dot(glm::cross(b1 - b2, b3 - b2), v - b2) / 6.f;
}

glm::vec3 util::tetrahedronBarycentre(const glm::vec3& b1, const glm::vec3& b2, const glm::vec3& b3, const glm::vec3& v) {
    return (b1 * (1.f/4.f))
           + (b2 * (1.f/4.f))
           + (b3 * (1.f/4.f))
           + (v * (1.f/4.f));
}

bool util::rayTriangleIntersection(glm::vec3 ray_origin, glm::vec3 ray_dir, glm::vec3 t1, glm::vec3 t2, glm::vec3 t3, float* parameter) {
    const float EPSILON = 0.0000001;
    glm::vec3 e1 = t2 - t1;
    glm::vec3 e2 = t3 - t1;

    glm::mat3 m = glm::mat3(-ray_dir.x, -ray_dir.y, -ray_dir.z,
                            e1.x, e1.y, e1.z,
                            e2.x, e2.y, e2.z);
    glm::vec3 o = glm::vec3(ray_origin.x - t1.x, ray_origin.y - t1.y, ray_origin.z - t1.z);

    float d = glm::determinant(m);

    float u = glm::determinant(glm::mat3(-ray_dir.x, -ray_dir.y, -ray_dir.z, o.x, o.y, o.z, e2.x, e2.y, e2.z)) / d;
    if(u < 0.0 || u > 1.0) return false;
    float v = glm::determinant(glm::mat3(-ray_dir.x, -ray_dir.y, -ray_dir.z, e1.x, e1.y, e1.z, o.x, o.y, o.z)) / d;
    if(v < 0.0 || v > 1.0 || u + v > 1.0) return false;
    float t = glm::determinant(glm::mat3(o.x, o.y, o.z, e1.x, e1.y, e1.z, e2.x, e2.y, e2.z)) / d;
    if(t < EPSILON) return false;

    glm::vec3 tNorm = glm::cross(e1, e2);
    float raydotnorm = glm::dot(ray_dir, tNorm);

    if (raydotnorm > -EPSILON && raydotnorm < EPSILON) {
        return false;    // This ray is parallel to this triangle.
    }

    *parameter = t;
    return true;
}

std::vector<glm::vec3> util::rayMeshIntersections(const std::vector<glm::vec3>& vertices,
                                                  const std::vector<glm::vec<3, unsigned int>>& faces,
                                                  glm::vec3 ray_origin, glm::vec3 ray_dir
) {
    // parameter: intersection = origin + direction * parameter
    std::vector<glm::vec3> intersections = {};
    std::vector<float> parameters = {};
    float parameter;
    for(auto & face : faces) {
        if(util::rayTriangleIntersection(ray_origin, ray_dir, vertices[face.x - 1], vertices[face.y - 1], vertices[face.z - 1], &parameter)) {
            parameters.push_back(parameter);
        }
    }
    if(parameters.empty()) return intersections;
    // sort to return intersection in order of distance from ray origin
    std::sort(parameters.begin(), parameters.end());

    // discard unique values in order to avoid double intersection through edges
    auto new_last = std::unique(parameters.begin(), parameters.end(),
                                [](float f1, float f2) { return std::abs(f1 - f2) < std::numeric_limits<float>::epsilon(); });
    parameters.resize(std::distance(parameters.begin(), new_last));
    intersections.reserve(parameters.size());
    for(float p : parameters) {
        intersections.push_back(ray_origin + ray_dir * p);
    }
    return intersections;
}
std::vector<glm::vec3> util::rayMeshIntersectionsOptimized(const std::vector<glm::vec3>& vertices,
                                                           const std::vector<glm::vec<3, unsigned int>>& faces,
                                                           glm::vec3 ray_origin, glm::vec3 ray_dir
) {
    // parameter: intersection = origin + direction * parameter
    std::vector<glm::vec3> intersections = {};
    std::vector<float> parameters = {};
    float parameter;
    for(auto & face : faces) {
        glm::vec3 t0 = vertices[face.x - 1];
        glm::vec3 t1 = vertices[face.y - 1];
        glm::vec3 t2 = vertices[face.z - 1];
        glm::vec3 o = ray_origin;
        // CUSTOM RAY MESH INTERSECTION: since ray is parallel to z axis, its possible
        // to do a prececk and exclude triangles that are out of ray vision;

        if(t0.x > o.x && t1.x > o.x && t2.x > o.x) continue;
        if(t0.x < o.x && t1.x < o.x && t2.x < o.x) continue;
        if(t0.y > o.y && t1.y > o.y && t2.y > o.y) continue;
        if(t0.y < o.y && t1.y < o.y && t2.y < o.y) continue;
        float Do = -o.x + o.y;
        if(t0.y < t0.x + Do && t1.y < t1.x + Do && t2.y < t2.x + Do) continue;
        if(t0.y > t0.x + Do && t1.y > t1.x + Do && t2.y > t2.x + Do) continue;

        if(util::rayTriangleIntersection(ray_origin, ray_dir, t0, t1, t2, &parameter)) {
            parameters.push_back(parameter);
        }
    }
    if(parameters.empty()) return intersections;
    // sort to return intersection in order of distance from ray origin
    std::sort(parameters.begin(), parameters.end());

    // discard unique values in order to avoid double intersection through edges
    auto new_last = std::unique(parameters.begin(), parameters.end(),
                                [](float f1, float f2) { return std::abs(f1 - f2) < std::numeric_limits<float>::epsilon(); });
    parameters.resize(std::distance(parameters.begin(), new_last));
    intersections.reserve(parameters.size());
    for(float p : parameters) {
        intersections.push_back(ray_origin + ray_dir * p);
    }
    return intersections;
}

glm::vec3 util::getMin(const std::vector<glm::vec3>& vertices) {
    glm::vec3 min = vertices[0];
    for(int i = 1; i < vertices.size(); i++) {
        if(vertices[i].x < min.x) min.x = vertices[i].x;
        if(vertices[i].y < min.y) min.y = vertices[i].y;
        if(vertices[i].z < min.z) min.z = vertices[i].z;
    }
    return min;
}
glm::vec3 util::getMax(const std::vector<glm::vec3>& vertices) {
    glm::vec3 max = vertices[0];
    for(int i = 1; i < vertices.size(); i++) {
        if(vertices[i].x > max.x) max.x = vertices[i].x;
        if(vertices[i].y > max.y) max.y = vertices[i].y;
        if(vertices[i].z > max.z) max.z = vertices[i].z;
    }
    return max;
}