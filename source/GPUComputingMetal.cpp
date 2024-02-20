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
#include <algorithm>
#include <timer.hpp>

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

    /*for(int i = 0; i < 30; i+=3) {
        std::cout << "space buffer" << ((float *)space_buffer->contents())[i] << " "
                  << ((float *)space_buffer->contents())[i + 1] << " "
                  << ((float *)space_buffer->contents())[i + 2] << std::endl;
    }*/

    /*for(int i = 0; i < 30; i+=4) {
        std::cout << "masses output" << ((float *)masses_buffer->contents())[i] << " "
                  << ((float *)masses_buffer->contents())[i + 1] << " "
                  << ((float *)masses_buffer->contents())[i + 2]
                  << ((float *)masses_buffer->contents())[i + 3]<< std::endl;
    }*/

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

float* GPUComputing::get_gravity_from_tubes_with_integral(const float* tubes, int tubes_size, float* p, float R, float G) {
    // create autorelease pool
    auto pool = NS::TransferPtr(NS::AutoreleasePool::alloc()->init());
    static NS::String* source = nullptr;
    if(source == nullptr) {
        std::cout << "READING METAL SHADER";
        // read shader
        std::ifstream file;
        file.open("../shaders/get_gravity_from_tubes.metal", std::ios::in);
        std::stringstream ss;
        ss << file.rdbuf();

        source = NS::String::string(ss.str().c_str(), NS::UTF8StringEncoding);
    }
    static auto device = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
    NS::Error* e_1 = nullptr;
    static auto library = NS::TransferPtr(device->newLibrary(source, nullptr, &e_1));
    // shader compilation output
    if(e_1) std::cout << e_1->localizedDescription()->utf8String();

    static auto function = NS::TransferPtr(library->newFunction(
        NS::String::string("get_gravity_from_tubes", NS::ASCIIStringEncoding)
        ));
    NS::Error* e_2 = nullptr;
    static auto computePipelineState = NS::TransferPtr(device->newComputePipelineState(function.get(), &e_2));
    if(e_2) std::cout << e_2->localizedDescription()->utf8String();
    static auto commandQueue = NS::TransferPtr(device->newCommandQueue());

    // create buffers
    auto tubes_buffer = NS::TransferPtr(device->newBuffer(tubes, tubes_size * 6 * sizeof(float), MTL::ResourceStorageModeShared));
    auto tubes_as_float_size_buffer = NS::TransferPtr(device->newBuffer(&tubes_size, sizeof(int), MTL::ResourceStorageModeShared));
    auto point_buffer = NS::TransferPtr(device->newBuffer(p, sizeof(float)*3, MTL::ResourceStorageModeShared));
    auto G_buffer = NS::TransferPtr(device->newBuffer(&G, sizeof(float), MTL::ResourceStorageModeShared));
    auto R_buffer = NS::TransferPtr(device->newBuffer(&R, sizeof(float), MTL::ResourceStorageModeShared));
    auto output_buffer = NS::TransferPtr(device->newBuffer(sizeof(float)* 3 * tubes_size, MTL::ResourceStorageModeShared));

    auto commandBuffer = commandQueue->commandBuffer();
    auto computeCommandEncoder = commandBuffer->computeCommandEncoder();

    computeCommandEncoder->setComputePipelineState(computePipelineState.get());
    computeCommandEncoder->setBuffer(tubes_buffer.get(), 0, 0);
    computeCommandEncoder->setBuffer(tubes_as_float_size_buffer.get(), 0, 1);
    computeCommandEncoder->setBuffer(point_buffer.get(), 0, 2);
    computeCommandEncoder->setBuffer(G_buffer.get(), 0, 3);
    computeCommandEncoder->setBuffer(R_buffer.get(), 0, 4);
    computeCommandEncoder->setBuffer(output_buffer.get(), 0, 5);

    auto gridSize = MTL::Size::Make(tubes_size, 1, 1);
    auto threadGroupSize = computePipelineState->maxTotalThreadsPerThreadgroup();
    threadGroupSize = std::min(tubes_size, (int)threadGroupSize);
    auto threadSize = MTL::Size::Make(threadGroupSize, 1, 1);

    computeCommandEncoder->dispatchThreads(gridSize, threadSize);
    computeCommandEncoder->endEncoding();
    commandBuffer->commit();
    commandBuffer->waitUntilCompleted();

    auto output = (float*)malloc(tubes_size * 3 * sizeof(float));
    memcpy(output, output_buffer->contents(), tubes_size * 3 * sizeof(float));

    return output;
}

float* GPUComputing::get_gravities_from_tubes_with_integral(const float* tubes, int tubes_size, const float* points, int points_size, float R, float G) {
    // create autorelease pool
    auto pool = NS::TransferPtr(NS::AutoreleasePool::alloc()->init());
    static NS::String* source = nullptr;
    if(source == nullptr) {
        std::cout << "READING METAL SHADER";
        // read shader
        std::ifstream file;
        file.open("../shaders/get_gravities_from_tubes.metal", std::ios::in);
        std::stringstream ss;
        ss << file.rdbuf();

        source = NS::String::string(ss.str().c_str(), NS::UTF8StringEncoding);
    }
    static auto device = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
    NS::Error* e_1 = nullptr;
    static auto library = NS::TransferPtr(device->newLibrary(source, nullptr, &e_1));
    // shader compilation output
    if(e_1) std::cout << e_1->localizedDescription()->utf8String();

    static auto function = NS::TransferPtr(library->newFunction(
            NS::String::string("get_gravities_from_tubes", NS::ASCIIStringEncoding)
    ));
    NS::Error* e_2 = nullptr;
    static auto computePipelineState = NS::TransferPtr(device->newComputePipelineState(function.get(), &e_2));
    if(e_2) std::cout << e_2->localizedDescription()->utf8String();
    static auto commandQueue = NS::TransferPtr(device->newCommandQueue());

    // create buffers
    auto tubes_buffer = NS::TransferPtr(device->newBuffer(tubes, tubes_size * 6 * sizeof(float), MTL::ResourceStorageModeShared));
    auto tubes_size_buffer = NS::TransferPtr(device->newBuffer(&tubes_size, sizeof(int), MTL::ResourceStorageModeShared));
    auto points_buffer = NS::TransferPtr(device->newBuffer(points, sizeof(float)*3*points_size, MTL::ResourceStorageModeShared));
    auto G_buffer = NS::TransferPtr(device->newBuffer(&G, sizeof(float), MTL::ResourceStorageModeShared));
    auto R_buffer = NS::TransferPtr(device->newBuffer(&R, sizeof(float), MTL::ResourceStorageModeShared));
    auto output_buffer = NS::TransferPtr(device->newBuffer(sizeof(float) * 3 * points_size, MTL::ResourceStorageModeShared));

    if(output_buffer->contents() == nullptr) std::cout << "before computation output buffer is null" << std::endl;

    auto commandBuffer = commandQueue->commandBuffer();
    auto computeCommandEncoder = commandBuffer->computeCommandEncoder();

    computeCommandEncoder->setComputePipelineState(computePipelineState.get());
    computeCommandEncoder->setBuffer(tubes_buffer.get(), 0, 0);
    computeCommandEncoder->setBuffer(tubes_size_buffer.get(), 0, 1);
    computeCommandEncoder->setBuffer(points_buffer.get(),0, 2);
    computeCommandEncoder->setBuffer(G_buffer.get(), 0, 3);
    computeCommandEncoder->setBuffer(R_buffer.get(), 0, 4);
    computeCommandEncoder->setBuffer(output_buffer.get(), 0, 5);

    auto gridSize = MTL::Size::Make(points_size, 1, 1);
    auto threadGroupSize = computePipelineState->maxTotalThreadsPerThreadgroup();
    threadGroupSize = std::min(points_size, (int)threadGroupSize);
    auto threadSize = MTL::Size::Make(threadGroupSize, 1, 1);

    computeCommandEncoder->dispatchThreads(gridSize, threadSize);
    computeCommandEncoder->endEncoding();
    commandBuffer->commit();
    commandBuffer->waitUntilCompleted();

    auto output = (float*)malloc(points_size * 3 * sizeof(float));
    if(output == nullptr) {
        std::cout << "failed to allocate memory for gpu output array" << std::endl;
        if(output_buffer->contents() == nullptr) std::cout << "even gpu buffer output is null" << std::endl;
    }

    memcpy(output, output_buffer->contents(), points_size * 3 * sizeof(float));

    return output;
}