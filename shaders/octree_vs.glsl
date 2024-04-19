#version 410 core

layout (location = 0) in vec3 position;
layout (location = 1) in float _depth;

uniform mat4 projectionMatrix;
uniform mat4 viewMatrix;

out float depth;

void main() {
    gl_Position = projectionMatrix * viewMatrix * vec4(position, 1);
    depth = _depth;
}