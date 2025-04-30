#pragma once

#include "Core/Pass/RasterPass.h"
#include "Scene/SceneBuilder.h"

using namespace Falcor;

class Renderer
{
public:
    Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept;

    void Init() noexcept;
    void RenderFrame(RenderContext* pRenderContext, const double& currentTime) const noexcept;
    void RenderUI(Gui* pGui, Gui::Window* app_gui_window) noexcept;
    void OnResize(uint32_t width, uint32_t height) noexcept;
    bool onKeyEvent(const KeyboardEvent& keyEvent) const noexcept;
    bool onMouseEvent(const MouseEvent& mouseEvent) const noexcept;
    void Deinit() noexcept;

    void SynchronizeSceneWithProgram() noexcept;
    [[nodiscard]] NodeID AddSphereToScene(float3 pos, float radius) noexcept;
    void UpdateSceneNodeTransform(NodeID nodeID, const Transform& transform) noexcept;

private:
    void setPerFrameVariables(const double& currentTime) const noexcept;

    const ref<Device>& device_;
    const ref<Fbo>& target_fbo_;

    ref<RasterPass> raster_pass_;

    ref<Program> program_;
    ref<ProgramVars> program_vars_;

    ref<Program> rt_program_;
    ref<RtProgramVars> rt_program_vars_;

    SceneBuilder* scene_builder_ = nullptr;
    ref<Scene> scene_;
    ref<Camera> camera_;

    ref<Texture> rt_output_tex_;
    ref<Texture> density_3d_tex_;

    MeshID sphere_mesh_id;

    bool draw_fluid_ = false;

    //bool mUseDOF = false;
    //uint32_t mSampleIndex = 0;

    // ===============================================================================
    //                                  Constants.                                   
    // ===============================================================================

    float3 bg_clear_color = float3(.2, 1, .1);

    uint kMaxRayBounce = 3;

    float3 absorptionCoeff = float3(1.0, 0.4, 0.05);
    float3 scatteringCoeff = float3(0.1, 0.2, 0.8);
    float phaseG = 0.8f;

    float water_turbulence_ = 2.5f;
    float maxRayMarchingDistance = 5.f;
    float kMarchSize = 0.1f;
    float maxLighMarchingDistance = 3.f;
    float sunLightMarchSize = 0.2f;

    float3 lightColor = float3(1, 1, 1);
    float3 lightDir = normalize(float3(1, -1, -1));

    float IoR = 1.33f;
};
