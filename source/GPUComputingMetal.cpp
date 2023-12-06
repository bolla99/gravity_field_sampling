//
// Created by bolla on 12/05/2023.
//

#include <GPUComputing.hpp>
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>

float* GPUComputing::get_gravity_from_point_masses_and_discrete_space(const float* pointMasses, int massesSize, const float* discreteSpace, int spaceSize, float mass_radius) {
    std::ifstream file;
    file.open("../shaders/add.metal", std::ios::in);
    std::stringstream ss;
    ss << file.rdbuf();
    std::string metal_shader = ss.str();
    const char* metal_shader_c = metal_shader.c_str();

    auto source = NS::String::string(metal_shader_c, NS::UTF8StringEncoding);
    NS::Error* error;
    auto device = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
    auto library = NS::TransferPtr(device->newLibrary(source, nullptr, &error));

    auto str = NS::String::string("add_arrays", NS::ASCIIStringEncoding);
    auto function = NS::TransferPtr(library->newFunction(str));
    if(function.get() == nullptr) return nullptr;

    NS::Error* e;
    auto computePipelineState = NS::TransferPtr(device->newComputePipelineState(function.get(), &e));
    auto commandQueue = NS::TransferPtr(device->newCommandQueue());


    auto space_buffer = NS::TransferPtr(device->newBuffer(discreteSpace, spaceSize, MTL::ResourceStorageModeShared));
    auto masses_buffer = NS::TransferPtr(device->newBuffer(pointMasses, massesSize, MTL::ResourceStorageModeShared));

    // output buffer
    auto output_buffer = NS::TransferPtr(device->newBuffer(spaceSize, MTL::ResourceStorageModeShared));
    massesSize /= sizeof(float);
    auto masses_size_buffer = NS::TransferPtr(device->newBuffer(&massesSize, sizeof(int), MTL::ResourceStorageModePrivate));
    auto mass_radius_buffer = NS::TransferPtr(device->newBuffer(&mass_radius, sizeof(float), MTL::ResourceStorageModePrivate));

    auto commandBuffer = commandQueue->commandBuffer();
    auto computeCommandEncoder = commandBuffer->computeCommandEncoder();

    computeCommandEncoder->setComputePipelineState(computePipelineState.get());
    computeCommandEncoder->setBuffer(space_buffer.get(), 0, 0);
    computeCommandEncoder->setBuffer(masses_buffer.get(), 0, 1);
    computeCommandEncoder->setBuffer(output_buffer.get(), 0, 2);
    computeCommandEncoder->setBuffer(masses_size_buffer.get(), 0, 3);
    computeCommandEncoder->setBuffer(mass_radius_buffer.get(), 0, 4);

    for(int i = 0; i < 30; i+=3) {
        std::cout << "space buffer" << ((float *)space_buffer->contents())[i] << " "
                  << ((float *)space_buffer->contents())[i + 1] << " "
                  << ((float *)space_buffer->contents())[i + 2] << std::endl;
    }

    for(int i = 0; i < 30; i+=4) {
        std::cout << "masses output" << ((float *)masses_buffer->contents())[i] << " "
                  << ((float *)masses_buffer->contents())[i + 1] << " "
                  << ((float *)masses_buffer->contents())[i + 2]
                  << ((float *)masses_buffer->contents())[i + 3]<< std::endl;
    }

    auto size = spaceSize / (3 * sizeof(float));
    MTL::Size gridsize = MTL::Size::Make(size, 1, 1);
    NS::UInteger threadGroupSize = computePipelineState->maxTotalThreadsPerThreadgroup();
    if(threadGroupSize > size) {
        threadGroupSize = size;
    }
    MTL::Size threadsize = MTL::Size::Make(threadGroupSize, 1, 1);
    computeCommandEncoder->dispatchThreads(gridsize, threadsize);

    computeCommandEncoder->endEncoding();
    commandBuffer->commit();
    commandBuffer->waitUntilCompleted();

    auto output = (float*)malloc(spaceSize);
    memcpy(output, output_buffer->contents(), spaceSize);

    for(int i = 0; i < 30; i+=3) {
        std::cout << "direct output" << ((float *)output_buffer->contents())[i] << " "
                  << ((float *)output_buffer->contents())[i + 1] << " "
                  << ((float *)output_buffer->contents())[i + 2] << std::endl;
    }
    return output;
}