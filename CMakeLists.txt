# cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

# Add a CMake option to enable or disable Tracy Profiler
option(USE_TRACY "Use Tracy Profiler" ON)

if (USE_TRACY)
    # Enable Tracy profiling by setting the preprocessor directive
    add_compile_definitions(TRACY_ENABLE)
    # Include Tracy's headers for the entire project
    include_directories(External/TracyProfiler/tracy)

    # Add the Tracy source code as part of the project
    add_library(tracyClient STATIC External/TracyProfiler/TracyClient.cpp)
endif()

# Common
file(GLOB_RECURSE COMMON_FILES Common/include/*.h Common/src/*.cpp)
add_library(Common ${COMMON_FILES})
target_include_directories(Common PUBLIC Common/include/)

# Physics
file(GLOB_RECURSE PHYSICS_FILES Physics/include/*.h Physics/src/*.cpp)
add_library(Physics ${PHYSICS_FILES})
target_include_directories(Physics PUBLIC Physics/include/)
target_include_directories(Physics PUBLIC Common/include/)

if (USE_TRACY)
    target_compile_definitions(Physics PUBLIC TRACY_ENABLE)
    # Link the TracyClient library
    target_link_libraries(Physics PRIVATE tracyClient)
endif()

# Sample library
file(GLOB_RECURSE SAMPLE_FILES samples/include/*.h samples/src/*.cpp)
add_library(PhysicsSamples ${SAMPLE_FILES})
target_include_directories(PhysicsSamples PUBLIC samples/include)
target_include_directories(PhysicsSamples PUBLIC Physics/include/)
target_include_directories(PhysicsSamples PUBLIC Common/include/)

if (USE_TRACY)
    target_compile_definitions(PhysicsSamples PUBLIC TRACY_ENABLE)
    # Link the TracyClient library
    target_link_libraries(PhysicsSamples PRIVATE tracyClient)
endif()

target_link_libraries(PhysicsSamples PUBLIC Physics Common imgui)

# Renderer
file(GLOB_RECURSE RENDERER_FILES Renderer/include/*.h Renderer/src/*.cpp)
add_library(Renderer ${RENDERER_FILES})
target_include_directories(Renderer PUBLIC Renderer/include/)
target_link_libraries(Renderer PUBLIC Falcor)

add_falcor_executable(Raytracing)

target_sources(Raytracing PRIVATE
    main.cpp
    Raytracing.cpp
    Raytracing.h
    Raytracing.3d.slang
    Raytracing.rt.slang
    triangle.slang
    SDF_Functions.slang
)

target_include_directories(Raytracing PUBLIC samples/include)
target_include_directories(Raytracing PUBLIC Physics/include/)
target_include_directories(Raytracing PUBLIC Common/include/)
target_link_libraries(Raytracing PUBLIC PhysicsSamples Renderer)

target_copy_shaders(Raytracing Samples/Raytracing)

target_source_group(Raytracing "Samples")


# Renderer for future direct x (actually opengl using glut testing)
# file(GLOB_RECURSE RENDERER_FILES Renderer/include/*.h Renderer/src/*.cpp)
# add_library(Renderer ${RENDERER_FILES})
# set_target_properties(Renderer PROPERTIES LINKER_LANGUAGE CXX)
# target_include_directories(Renderer PUBLIC Renderer/include/)
# target_include_directories(Renderer PUBLIC Common/include/)
# target_include_directories(Renderer PUBLIC External/imgui/include/)

# target_link_libraries(Renderer PUBLIC Samples GLUT::GLUT GLEW::GLEW OpenGL::GL)

# target_link_libraries(Samples PUBLIC Renderer)

# if (USE_TRACY)
#     target_compile_definitions(Renderer PUBLIC TRACY_ENABLE)
#     # Link the TracyClient library
#     target_link_libraries(Renderer PRIVATE tracyClient)
# endif()

# # Main
# add_executable(Main Main.cpp)
# target_link_libraries(Main PUBLIC Renderer)
