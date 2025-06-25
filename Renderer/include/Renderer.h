#pragma once

#include "Metrics.h"
#include "Core/Pass/RasterPass.h"
#include "RenderGraph/RenderGraph.h"
#include "Scene/SceneBuilder.h"

using namespace Falcor;

struct MarchingCubeVertex
{
    float3 position = float3(0.f, 0.f, 0.f);
    //float pad0;
    float3 normal = float3(0.f, 0.f, 1.f);
    //float pad1;
    //float3 tangent = float3(0.f, 0.f, 1.f);
    //float pad2;
    //float2 texcrd = float2(0.f, 0.f);
    //float pad3;
    //float pad4;
};

struct MarchingCubesTriangle
{
    MarchingCubeVertex vertexA;
    MarchingCubeVertex vertexB;
    MarchingCubeVertex vertexC;
};

class Renderer
{
public:
    Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept;

    void Init(RenderContext* render_context, bool rebuildBvh = false) noexcept;
    void RenderFrame(
        RenderContext* pRenderContext,
        const ref<Fbo>& pTargetFbo, const double& currentTime,
        const ref<Buffer>& bodies,
        const ref<Buffer>& SpatialIndices,
        const ref<Buffer>& SpatialOffsets) noexcept;
    void RenderUI(Gui* pGui, Gui::Window* app_gui_window, RenderContext* render_context) noexcept;
    void OnResize(uint32_t width, uint32_t height) noexcept;
    [[nodiscard]] bool onKeyEvent(const KeyboardEvent& keyEvent) const noexcept;
    [[nodiscard]] bool onMouseEvent(const MouseEvent& mouseEvent) const noexcept;
    void Deinit() noexcept;

    /**
     * \brief CreateRaytracingProgram is a method which creates the raytracing
     * program and links it to the scene's BVH. It should be called after all meshes were added to the scene.
     */
    void CreateRaytracingProgram(RenderContext* render_context) noexcept;
    void CreateRasterizationProgram() noexcept;

    void LaunchMarchingCubeComputePasses(RenderContext* render_context) noexcept;

    [[nodiscard]] NodeID AddSphereToScene(float3 pos, float radius) noexcept;
    [[nodiscard]] NodeID AddCubeToScene(float3 pos) noexcept;
    void UpdateSceneNodeTransform(NodeID nodeID, const Transform& transform) const noexcept;

    void RegisterParticleDensities(std::vector<float>* particle_densities) noexcept
    {
        particle_densities_ = particle_densities;
        particle_density_buffer_ = make_ref<Buffer>(
            device_,
            sizeof(float),
            particle_densities->size(),
            ResourceBindFlags::ShaderResource,
            MemoryType::DeviceLocal,
            particle_densities->data(),
            false
        );
    }

    uint32_t fluid_AABB_ID = 1;
    Transform fluid_transform;
    ref<RasterPass> raster_pass_;

private:
    void setPerFrameVariables(const double& currentTime) const noexcept;


    const ref<Device>& device_;
    const ref<Fbo>& target_fbo_;

    RenderGraph render_graph_;

    //ref<RasterPass> raster_pass_;
    //ref<ComputePass> density_map_pass_;

    ref<Program> program_;
    ref<ProgramVars> program_vars_;

    ref<Program> rt_program_;
    ref<RtProgramVars> rt_program_vars_;

    SceneBuilder* scene_builder_ = nullptr;
    ref<Scene> scene_;
    ref<Camera> camera_;

    ref<Texture> rt_output_tex_;
    ref<Texture> density_3d_tex_;
    ref<Texture> marching_cube_dens_tex;
    ref<Sampler> linearClampSampler_;
    ref<Buffer> bodies_buffer_ = nullptr;
    std::vector<float>* particle_densities_ = nullptr;
    ref<Buffer> particle_density_buffer_;

    MeshID sphere_mesh_id;
    MeshID cube_mesh_id;
    MeshID marching_cubes_mesh_id_;
    NodeID sphere_node_id_;
    NodeID raymarching_node_id;
    ref<TriangleMesh> sphere;

    ref<Buffer> b_pos;
    ref<Buffer> b_pos_readback;
    ref<Buffer> b_norm_readback;
    ref<Buffer> b_tang_readback;
    ref<Buffer> b_uv_readback;
    ref<Buffer> b_normal;
    ref<Buffer> b_tang;
    ref<Buffer> b_uv;

    uint32_t oldTriangleCount = 0;
    uint32_t triangleCount = 0;

    float var_ = 0.f;

    ref<ComputePass> compute_density_map_pass_ = nullptr;
    ref<ComputePass> marching_cubes_pass_ = nullptr;
    ref<ComputePass> compute_marching_cube_density_map_ = nullptr;
    ref<Buffer> marching_cubes_triangle_buffer_ = nullptr;
    ref<Buffer> read_back_triangle_buffer_ = nullptr;


    //ref<Material> dielectric_blue;
public:
    float DensityRayMarchMultiplier = 0.05f;
    float volumeValueOffset = 0.1f;
    float normalOffset = 0.1f;

    bool useMarchingCubes = false;
    bool draw_fluid_ = false;
    bool lightScattering = false;


    size_t MaxTriangleCount = 0;
    size_t MaxVertexCount = 0;
    unsigned marching_cube_vertex_count = 0;
    float3 march_mesh_scale = float3(Metrics::sim_bounds);
    //ref<Texture> density_3d_tex_ = nullptr;

    TriangleMesh::VertexList v;
    MeshID tri_id;

    std::map<std::string, ref<Buffer>> vertices;

    //bool mUseDOF = false;
    //uint32_t mSampleIndex = 0;

    // ===============================================================================
    //                                  Constants.                                   
    // ===============================================================================

    //int density_map_size = 64;
    //float3 sim_bounds = float3(Metrics::MetersToPixels(1.0f)) * 2.f;
    float bounds_size = Metrics::MetersToPixels(1.0f) * 2.f;
    int numPointsPerAxis = 64;
    float IsoLevel = 0.f;
    float SphereRadius = 15.f;

    float3 bg_clear_color = float3(.2, 1, .1);

    uint kMaxRayBounce = 3;

    float3 absorptionCoeff = float3(1.0, 0.4, 0.05);
    float3 scatteringCoeff = float3(2.19, 0.75, 0.55);
    float phaseG = 0.8f;

    float water_turbulence_ = 2.5f;
    float maxRayMarchingDistance = 5.f;
    float kMarchSize = 0.1f; //*Metrics::voxelSize;
    float maxLighMarchingDistance = 3.f;
    float sunLightMarchSize = 0.5f;

    Transform transfrom{};
    float3 translation{0.f, 0.f, 0.f};
    float3 rotation{0.f, 0.f, 0.f};
    float3 scale{Metrics::WALLDIST, Metrics::WALLDIST, Metrics::WALLDIST};
    bool useTransformations = false;
    bool useRecursiveRaytracing = true;
    bool debugNormals = false;

    float3 lightColor = float3(1, 1, 1);
    float3 lightDir = normalize(float3(1, -1, -1));

    float IoR = 1.33f;
};
