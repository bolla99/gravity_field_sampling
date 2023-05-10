//
// Created by Giovanni Bollati on 09/04/23.
//

#ifndef GL_TEST_PROJECT_SHADER_HPP
#define GL_TEST_PROJECT_SHADER_HPP

#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include <iostream>
#include <fstream>
#include <SDL_opengl.h>

class Shader {
public:
    unsigned int programID;

    // constructor
    Shader(const std::string& vertexPath, const std::string& fragmentPath);

   void use() const;

    // get uniform location
    int getUniformLocation(const std::string& name) const;

    // load shader source and checkers for compilation and linkage
    static std::string loadShaderSource(const std::string& path);
    static bool checkShaderCompilation(unsigned int shader);
    static bool checkProgramLinkage(unsigned int program);
};

#endif //GL_TEST_PROJECT_SHADER_HPP
