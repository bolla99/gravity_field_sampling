## Efficient computation of gravity field induced by an arbitrarily shaped object
This project implements a method for the computation of gravity induced
by an arbitrarily shaped body which I developed as my bachelor's degree
thesis.

### Implementation notes
Dependencies:
- *GLM*: https://github.com/g-truc/glm
- *SDL2*: https://www.libsdl.org/
- *Metal*: https://developer.apple.com/metal/
- *OpenGL*: https://opengl.org/
- *OpenMP*: https://www.openmp.org/
- *Dear ImGUI*: https://github.com/ocornut/imgui

Only *GLM* and *Metal* are required for the basic functionality, 
which is provided by three namespaces, *gravity*, *util* and *GPUComputing*.
*util* and *GPUComputing* have no dependency other than the
above-mentioned ones; *util* provides general purpose 
funcionality, while *GPUComputing* provides GPU accelerated 
computations, which are implemented with *Metal*, which means 
the whole project can run only on *macOS*.
*util* and *GPUComputing* provide services for the  
*gravity* namespace, which exposes the interface for the core functionality.
Documentation can be found in the corresponding header file
(gravity.hpp).
The rest of the project consists of a mesh abstraction, a 
glsl shader utility class and a main, which provides an environment
which allow to try and test functionalities.