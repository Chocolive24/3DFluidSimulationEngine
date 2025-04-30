#pragma once

#include "Core/Pass/RasterPass.h"
#include "Scene/SceneBuilder.h"

using namespace Falcor;

class Renderer
{
public:
    Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept;

    void Init() noexcept;
    void RenderFrame(RenderContext* pRenderContext, const double& currentTime) noexcept;
    void RenderUI(Gui* pGui, Gui::Window* app_gui_window) noexcept;
    void OnResize(uint32_t width, uint32_t height) noexcept;
    bool onKeyEvent(const KeyboardEvent& keyEvent) const noexcept;
    bool onMouseEvent(const MouseEvent& mouseEvent) const noexcept;
    void Deinit() noexcept;

private:
    const ref<Device>& device_;
    const ref<Fbo>& target_fbo_;

    ref<RasterPass> raster_pass_;

    ref<Program> program_;
    ref<ProgramVars> program_vars_;
    ref<Buffer> vertex_buffer_;
    ref<VertexBufferLayout> vertex_buffer_layout_;
    ref<VertexLayout> vertex_layout_;
    ref<Vao> vao_;

    ref<Program> rt_program_;
    ref<RtProgramVars> rt_program_vars_;

    SceneBuilder* scene_builder_ = nullptr;
    ref<Scene> mpScene;
    ref<Camera> mpCamera;

    bool mUseDOF = false;
    uint32_t mSampleIndex = 0;

    ref<Texture> mpRtOut;
    ref<Texture> mpTexture3D;

    std::vector<NodeID> sphereNodeIDs;

    // ===============================================================================
    //                                  Constants.                                   
    // ===============================================================================

    uint32_t mSampleGuiWidth = 250;
    uint32_t mSampleGuiHeight = 200;
    uint32_t mSampleGuiPositionX = 20;
    uint32_t mSampleGuiPositionY = 40;

    float3 kClearColor = float3(.2, 1, .1);

    float waterTurbulence = 2.5f;

    uint kMaxRayBounce = 4;

    float3 absorptionCoeff = float3(1.0, 0.4, 0.05);
    float3 scatteringCoeff = float3(0.1, 0.2, 0.8);
    float phaseG = 0.8f;

    float maxRayMarchingDistance = 5.f;
    float kMarchSize = 0.1f;
    float maxLighMarchingDistance = 3.f;
    float sunLightMarchSize = 0.2f;

    float3 lightColor = float3(1, 1, 1);
    float3 lightDir = normalize(float3(1, -1, -1));

    float IoR = 1.33f;

    const std::string kEnvMapPath = "hallstatt4_hd.hdr";
};
