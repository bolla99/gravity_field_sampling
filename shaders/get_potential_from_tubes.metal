#include <metal_stdlib>
#include <metal_relational>
using namespace metal;

kernel void get_potential_from_tubes(
    device float* tubes,
    device int* n_of_tubes,
    device float* point,
    device float* G,
    device float* R,
    device float* contributes,
    uint index [[thread_position_in_grid]]
                       ) {
        contributes[index] = 0.f;
        float r = R[0];
        float PI = 3.14159265359;
        float3 p = float3(point[0], point[1], point[2]);
        float3 t1 = float3(tubes[6*index], tubes[6*index + 1], tubes[6*index + 2]);
        float3 t2 = float3(tubes[6*index + 3], tubes[6*index + 4], tubes[6*index + 5]);

        float3 d = normalize(t2 - t1);
        float l = distance(t1, t2);
        float pot = -PI*G[0]*pow(r,2)*(log(sqrt(pow(l,2) + 2*l*dot(d, t1) -2*l*dot(d, p) + length_squared(t1) -2*dot(t1, p) + length_squared(p)) + l + dot(d, t1) - dot(d, p)) - log(sqrt(length_squared(t1) -2*dot(t1, p) + length_squared(p)) + dot(d, t1) - dot(d, p)));
        if(!isnan(pot)) contributes[index] = pot;


}