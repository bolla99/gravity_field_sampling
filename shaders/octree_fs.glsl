#version 410 core

out vec4 FragColor;

in float depth;

uniform vec4 color;
uniform int max_depth;

void main() {
    FragColor = vec4(1 - (depth)/float(max_depth), 0, (depth)/float(max_depth), 1);
    //FragColor = vec4(0, depth, 0, 1);
}