//
// Created by Giovanni Bollati on 09/04/23.
//

#include "Shader.hpp"

Shader::Shader(const std::string& vertexPath, const std::string& fragmentPath) {
    std::string vertexSource;
    try {
        vertexSource = Shader::loadShaderSource(vertexPath);
    } catch(std::ifstream::failure &e) { std::cout << "vertex shader reading at" << vertexPath << " path failed: " << e.what(); }
    const char* vertexSource_cstr = vertexSource.c_str();

    std::string fragmentSource;
    try {
        fragmentSource = Shader::loadShaderSource(fragmentPath);
    } catch(std::ifstream::failure &e) { std::cout << "fragment shader reading at" << fragmentPath << " path failed: " << e.what(); }

    const char* fragmentSource_cstr = fragmentSource.c_str();

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertexShader, 1, &vertexSource_cstr, nullptr);
    glShaderSource(fragmentShader, 1, &fragmentSource_cstr, nullptr);

    // POINT OF FAILURE
    glCompileShader(vertexShader);
    glCompileShader(fragmentShader);

    checkShaderCompilation(vertexShader);
    checkShaderCompilation(fragmentShader);

    programID = glCreateProgram();
    glAttachShader(programID, vertexShader);
    glAttachShader(programID, fragmentShader);

    // POINT OF FAILURE
    glLinkProgram(programID);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}

void Shader::use() const {
    glUseProgram(programID);
}

int Shader::getUniformLocation(const std::string& name) const {
    const char* name_cstr = name.c_str();
    return glGetUniformLocation(this->programID, name_cstr);
}

std::string Shader::loadShaderSource(const std::string& path) {
    std::ifstream fstream;
    fstream.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    fstream.open(path, std::ios::in);
    std::string shader_source_code;
    std::string line;
    while(fstream.good()) {
        std::getline(fstream, line);
        shader_source_code.append(line + "\n");
    }
    return shader_source_code;
}

bool Shader::checkShaderCompilation(unsigned int shader) {
    int  success;
    char infoLog[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED -> " << infoLog;
        return false;
    }
    return true;
}
bool Shader::checkProgramLinkage(unsigned int program) {
    int success;
    char infoLog[512];
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(program, 512, nullptr, infoLog);
        std::cout << "ERROR::PROGRAM::LINKAGE FAILED -> %s\n" << infoLog;
        return false;
    }
    return true;
}
