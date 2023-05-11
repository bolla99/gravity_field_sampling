cmake_minimum_required(VERSION 3.9)

project(gl_test_project)

# required for imgui
cmake_policy(SET CMP0079 NEW)

set(RUNTIME_OUTPUT_DIRECTORY, "${PROJECT_SOURCE_DIR}/win_install/")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")

include_directories("${PROJECT_SOURCE_DIR}/include/")

if(APPLE)
link_directories("/opt/homebrew/opt/libomp/lib")
include_directories("/opt/homebrew/opt/libomp/include")
endif()

add_executable(gl_test_project
		source/main.cpp
        source/Mesh.cpp
		source/Shader.cpp
		source/stb_image.cpp
		shaders/add.metal
        shaders/vs.glsl
        shaders/fs.glsl
		)

if(WIN32)
	message("so: windows")
endif()
if(APPLE)
	message("so: apple")
endif()

# SDL
find_package(SDL2 REQUIRED CONFIG REQUIRED COMPONENTS SDL2)

if(SDL2_FOUND)
	message("sdl2 founded")
endif()

# OPENGL
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
	message("opengl founded")
endif()

# IMGUI
add_subdirectory(imgui)

# GLM
add_subdirectory(glm)

# OPENMP
message("looking for openmp")

if(APPLE)
	message("setting openmp on macos")
	set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
	set(OpenMP_C_LIB_NAMES "omp")
	set(OpenMP_omp_LIBRARY omp)

	set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
	set(OpenMP_CXX_LIB_NAMES "omp")
	set(OpenMP_omp_LIBRARY omp)
endif()

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
	message("open mp founded")
endif()

# LINK LIBRARIES
target_link_libraries(gl_test_project PRIVATE
		SDL2::SDL2
		OpenGL::GL
		OpenGL::GLU
		imgui
		glm
		OpenMP::OpenMP_CXX
		)

if(APPLE)
	message("add metal")
	add_compile_definitions(METAL_READY)
	target_link_libraries(gl_test_project PRIVATE
		"-framework Metal"
		"-framework Foundation"
		"-framework QuartzCore"
		"-framework MetalKit"
		"-framework AppKit"
	)
endif()

target_include_directories(gl_test_project PRIVATE ./imgui/)
target_link_libraries(imgui PRIVATE SDL2::SDL2)
