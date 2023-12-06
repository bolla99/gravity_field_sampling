#include <metal_stdlib>
using namespace metal;

kernel void add_arrays(device float* space,
    device float* masses,
    device float* gravity,
    device int* n_of_masses,
    device float* mass_radius,
    uint index [[thread_position_in_grid]]
                       ) {
        gravity[3*index] = 0.f;
        gravity[3*index + 1] = 0.f;
        gravity[3*index + 2] = 0.f;
        float3 space_point = float3(space[3*index], space[3*index + 1], space[3*index + 2]);
        for(int i = 0; i < n_of_masses[0]; i+=4) {
            float3 mass = float3(masses[i], masses[i + 1], masses[i + 2]);
            float3 dir = mass - space_point;
            // elevare a 3/2
            float r3 = pow(dot(dir, dir), 1.5);
            if(length(dir) < mass_radius[0]) {
                r3 = pow(mass_radius[0], 3);
            }
            if(r3 > -0.000001 && r3 < 0.000001) { continue; }
            gravity[3*index] += (dir.x * 10.f * masses[i + 3]) / r3;
            gravity[3*index + 1] += (dir.y * 10.f * masses[i + 3]) / r3;
            gravity[3*index + 2] += (dir.z * 10.f * masses[i + 3]) / r3;
            //}
        }
}
