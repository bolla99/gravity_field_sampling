//
// Created by Giovanni Bollati on 09/04/23.
//

#ifndef GL_TEST_PROJECT_SHADER_HPP
#define GL_TEST_PROJECT_SHADER_HPP

#ifndef GL_GLEXT_PROTOTYPES
#define GL_GLEXT_PROTOTYPES 1
#endif

#include <fstream>

#ifdef WIN32
#include <GL/glew.h>
#endif

class shader {
public:
    unsigned int programID;

    // constructor
    shader(const std::string& vertex_path, const std::string& fragment_path);

   void use() const;

    // get uniform location
    int get_uniform_location(const std::string& name) const;

    // load shader source and checkers for compilation and linkage
    static std::string load_shader_source(const std::string& path);
    static bool check_shader_compilation(unsigned int shader);
    static bool check_program_linkage(unsigned int program);
};

#endif //GL_TEST_PROJECT_SHADER_HPP
