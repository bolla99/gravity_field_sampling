#include <metal_stdlib>
using namespace metal;

kernel void add_arrays(device float* space,
    device float* masses,
    device float* gravity,
    device int* n_of_masses,
    uint index [[thread_position_in_grid]]
                       ) {
        gravity[3*index] = 9.f;
        gravity[3*index + 1] = 0.f;
        gravity[3*index + 2] = 0.f;
        for(int i = 0; i < n_of_masses[0]; i+=4) {
            float3 mass = float3(masses[i], masses[i + 1], masses[i + 2]);
            float3 space_point = float3(space[3*index], space[3*index + 1], space[3*index + 2]);
            float3 dir = mass - space_point;
            float3 amount = float3(0.f, 0.f, 0.f);
            float r3 = pow(length(dir), 3);
            if(r3 > 0.000000000001 || r3 < -0.000000000001) {
                amount = (dir * 10 * masses[i + 3]) / r3;
            }
            gravity[3*index] += amount.x;
            gravity[3*index + 1] += amount.y;
            gravity[3*index + 2] += amount.z;
        }
}
