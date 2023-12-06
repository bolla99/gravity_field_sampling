//
// Created by bolla on 12/05/2023.
//

#ifndef GL_TEST_PROJECT_GPUCOMPUTING_HPP
#define GL_TEST_PROJECT_GPUCOMPUTING_HPP


// implementation will differ based on OS and available GPGPU APIs; implementation file is chosen during build phase;
namespace GPUComputing {
    // massesSize and discreteSize refer to previous argumnet length in bytes
    float* get_gravity_from_point_masses_and_discrete_space(const float* pointMasses, int massesSize, const float* discreteSpace, int spaceSize, float mass_radius);
};

#endif //GL_TEST_PROJECT_GPUCOMPUTING_HPP


