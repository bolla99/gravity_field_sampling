#version 410 core

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 col;

// this attribute must have a default value, which is used in case a mesh doesn't have vertex color nor uvs, and
// a default color is given with a 1x1 texture and a constant (0, 0) uv for each vertex
layout (location = 2) in vec2 uv;
layout (location = 3) in vec3 normal;

out vec3 vertexColor;
out vec3 position;
out vec2 uvs;
out vec3 norm;

uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;
uniform mat4 modelMatrix;

void main() {
    gl_Position = projectionMatrix * viewMatrix * modelMatrix * vec4(pos.x, pos.y, pos.z, 1.0);

    vec4 worldPosition = modelMatrix * vec4(pos, 1);
    position = vec3(worldPosition.x, worldPosition.y, worldPosition.z);
    vertexColor = col;
    uvs = uv;
    norm = mat3(modelMatrix) * normal;
}