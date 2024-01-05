#version 410 core

// INS
in vec3 position;
// vertex color; default color if mesh has no vertex color data
in vec3 vertexColor;
// uvs; (0, 0) if mesh has no uvs data
in vec2 uvs;
in vec3 norm;

// OUTS
out vec4 FragColor;

// UNIFORMS
// colorMode = 0 -> vertex color; colorMode = 1 -> texture mode (with default texture if no texture is loaded)
// default value -> 0 -> which means default color until a texture is loaded.
// switch to 1 if you want to display vertex color
uniform int colorMode;
uniform sampler2D sampler;

uniform float ambientLightIntensity;
uniform vec3 ambientLightColor;

uniform vec3 diffuseLightColor;
uniform vec3 diffuseLightPosition;

uniform float specularStrenght;
uniform int shininess;

uniform vec3 viewerPosition;

void main()
{
    vec3 diffuseLightDirection = normalize(diffuseLightPosition - position);
    vec3 normalized_normal = normalize(norm);
    float diffuseStrenght = max(dot(normalized_normal, diffuseLightDirection), 0.0);

    // SPECULAR LIGHT
    vec3 incidentLight = normalize(position - diffuseLightPosition);
    vec3 reflectionDirection = reflect(incidentLight, normalized_normal);
    vec3 viewDirection = normalize(viewerPosition - position);
    vec3 specularLight = specularStrenght * pow(max(dot(viewDirection, reflectionDirection), 0.0), shininess) * diffuseLightColor;

    vec3 diffuseLight = diffuseStrenght * diffuseLightColor;

    vec4 ObjectColor = (1 - colorMode) * texture(sampler, vec2(uvs.x, 1-uvs.y)) + colorMode * vec4(vertexColor, 1);

    // LIGHTNING
    vec3 ambientLight = ambientLightIntensity * ambientLightColor;

    vec4 light = vec4(ambientLight + diffuseLight + specularLight, 1);
    ObjectColor = light * ObjectColor;
    FragColor = ObjectColor;
}