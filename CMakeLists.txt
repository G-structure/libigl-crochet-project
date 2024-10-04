cmake_minimum_required(VERSION 3.16)
project(crochet)

set(CMAKE_BUILD_TYPE Debug)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  add_compile_options(-Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable)
endif()

# Libigl
include(libigl)

# Download tutorial data
include(libigl_tutorial_data)

# Eigen
include(eigen)

# Imgui
#include(imgui)

# Enable the target igl::glfw
igl_include(glfw)
# Other modules you could enable
#igl_include(embree)
igl_include(imgui)
#igl_include(opengl)
#igl_include(stb)
#igl_include(predicates)
#igl_include(xml)
#igl_include(copyleft cgal)
#igl_include(copyleft comiso)
#igl_include(copyleft core)
#igl_include(copyleft cork)
#igl_include(copyleft tetgen)
#igl_include(restricted matlab)
#igl_include(restricted mosek)
#igl_include(restricted triangle)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g -O0)
endif()

# Additional warnings
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  add_compile_options(-Wall -Wextra -Wpedantic -Wconversion)
endif()

# Address sanitizer (optional)
if(CMAKE_CXX_COMPILER_ID MATCHES GNU|Clang)
  add_compile_options(-fsanitize=address)
  add_link_options(-fsanitize=address)
endif()

# Add your project files
file(GLOB SRC_FILES src/crochet/*.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})

# Link igl (and the glfw module) to your project
target_link_libraries(${PROJECT_NAME} PUBLIC
  igl::glfw
  igl::tutorial_data
  Eigen3::Eigen
  ## Other modules you could link to
  # igl::embree
  igl::imgui
  # igl::opengl
  # igl::stb
  # igl::predicates
  # igl::xml
  # igl_copyleft::cgal
  # igl_copyleft::comiso
  # igl_copyleft::core
  # igl_copyleft::cork
  # igl_copyleft::tetgen
  # igl_restricted::matlab
  # igl_restricted::mosek
  # igl_restricted::triangle
  )
