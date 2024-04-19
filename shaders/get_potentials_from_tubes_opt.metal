#include <metal_stdlib>
#include <metal_relational>
using namespace metal;

kernel void get_potentials_from_tubes(
    device float* tubes,
    device int* n_of_tubes,
    device float* points,
    device float* G,
    device float* R,
    device float* potential,
    uint index [[thread_position_in_grid]]
                       ) {
        potential[index] = 0.f;
        float r = R[0];
        float PI = 3.14159265359;
        float3 p = float3(points[index*3], points[index*3 + 1], points[index*3 + 2]);

        for(int i = 0; i < n_of_tubes[0]; i++) {
            float3 t1 = float3(tubes[6*i], tubes[6*i + 1], tubes[6*i + 2]);
            float3 t2 = float3(tubes[6*i + 3], tubes[6*i + 4], tubes[6*i + 5]);
            float3 d = normalize(t2 - t1);
            float l = distance(t1, t2);
            float pot = -PI*G[0]*pow(r,2)*(log(sqrt(pow(l,2) + 2*l*dot(d, t1) -2*l*dot(d, p) + length_squared(t1) -2*dot(t1, p) + length_squared(p)) + l + dot(d, t1) - dot(d, p)) - log(sqrt(length_squared(t1) -2*dot(t1, p) + length_squared(p)) + dot(d, t1) - dot(d, p)));
            if(isnan(pot)) continue;
            potential[index] += pot;
        }
}
