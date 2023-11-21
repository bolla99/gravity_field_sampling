//
// Created by Giovanni Bollati on 09/04/23.
//

#include "Shader.hpp"

shader::shader(const std::string& vertex_path, const std::string& fragment_path) {
    std::string vertex_source;
    try {
        vertex_source = shader::load_shader_source(vertex_path);
    } catch(std::ifstream::failure &e) { std::cout << "vertex shader reading at" << vertex_path << " path failed: " << e.what(); }
    const char* vertexSource_cstr = vertex_source.c_str();

    std::string fragment_source;
    try {
        fragment_source = shader::load_shader_source(fragment_path);
    } catch(std::ifstream::failure &e) { std::cout << "fragment shader reading at" << fragment_path << " path failed: " << e.what(); }

    const char* fragmentSource_cstr = fragment_source.c_str();

    unsigned int vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    unsigned int fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertex_shader, 1, &vertexSource_cstr, nullptr);
    glShaderSource(fragment_shader, 1, &fragmentSource_cstr, nullptr);

    // POINT OF FAILURE
    glCompileShader(vertex_shader);
    glCompileShader(fragment_shader);

    check_shader_compilation(vertex_shader);
    check_shader_compilation(fragment_shader);

    programID = glCreateProgram();
    glAttachShader(programID, vertex_shader);
    glAttachShader(programID, fragment_shader);

    // POINT OF FAILURE
    glLinkProgram(programID);

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
}

void shader::use() const {
    glUseProgram(programID);
}

int shader::get_uniform_location(const std::string& name) const {
    const char* name_cstr = name.c_str();
    return glGetUniformLocation(this->programID, name_cstr);
}

std::string shader::load_shader_source(const std::string& path) {
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

bool shader::check_shader_compilation(unsigned int shader) {
    int  success;
    char info_log[512];
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if(!success)
    {
        glGetShaderInfoLog(shader, 512, nullptr, info_log);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED -> " << info_log;
        return false;
    }
    return true;
}
bool shader::check_program_linkage(unsigned int program) {
    int success;
    char info_log[512];
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(program, 512, nullptr, info_log);
        std::cout << "ERROR::PROGRAM::LINKAGE FAILED -> %s\n" << info_log;
        return false;
    }
    return true;
}
