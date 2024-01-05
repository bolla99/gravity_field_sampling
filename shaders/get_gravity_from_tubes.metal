#include <metal_stdlib>
#include <metal_relational>
using namespace metal;

kernel void get_gravity_from_tubes(
    device float* tubes,
    device int* n_of_tubes,
    device float* point,
    device float* G,
    device float* R,
    device float* contributes,
    uint index [[thread_position_in_grid]]
                       ) {
        contributes[3*index] = 0.f;
        contributes[3*index + 1] = 0.f;
        contributes[3*index + 2] = 0.f;
        float3 t1 = float3(tubes[6*index], tubes[6*index + 1], tubes[6*index + 2]);
        float3 t2 = float3(tubes[6*index + 3], tubes[6*index + 4], tubes[6*index + 5]);
        float3 p = float3(point[0], point[1], point[2]);
        float3 t1_p = p - t1;
        float3 x = normalize(t2 - t1);
        if(distance_squared(t2, p) < distance_squared(t1, p)) {
            x = normalize(t1 - t2);
        }
        float3 o = t1 + x * dot(x, t1_p);
        float3 y = normalize(p - o);
        float3 z = normalize(cross(x, y));

        float a = distance(p, o);
        if(a < R[0]) return;

        float b = distance(t1, o);
        float c = distance(t2, o);

        if(b > c) {
            float tmp = b;
            b = c;
            c = tmp;
        }

        float PI = 3.14159265359;

        float fy = ((G[0]*PI*pow(R[0],2)) *
                   ((c / sqrt(pow(a, 2) + pow(c, 2))) - (b / sqrt(pow(a, 2) + pow(b, 2))))
                   ) / a;
        float fx = (G[0]*PI*pow(R[0],2)) *
                   ((1 / sqrt(pow(a, 2) + pow(b, 2))) - (1 / sqrt(pow(a, 2) + pow(c, 2))));

       if(isnan(fx) || isnan(fy)) return;

        fy = -fy;

        if(p.z > t1.z && p.z < t2.z) {
                    fy -= 2.f * ((G[0]*PI*pow(R[0],2)) *
                            ((b / sqrt(pow(a, 2) + pow(b, 2))))
                            ) / a;
                }

        float3 output_value = float3x3(x, y, z) * float3(fx, fy, 0.f);

        contributes[3*index] = output_value.x;
        contributes[3*index + 1] = output_value.y;
        contributes[3*index + 2] = output_value.z;
}
