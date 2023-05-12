//
// Created by bolla on 12/05/2023.
//

#ifndef GL_TEST_PROJECT_GPUCOMPUTING_HPP
#define GL_TEST_PROJECT_GPUCOMPUTING_HPP


// implementation will differ based on OS and available GPGPU APIs; implementation file is chosen during build phase;
class GPUComputing {
public:
    static float* getGravityFromPointMassesAndDiscreteSpace(const float* pointMasses, const float* DiscreteSpace);
};

#endif //GL_TEST_PROJECT_GPUCOMPUTING_HPP


