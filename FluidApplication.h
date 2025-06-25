/***************************************************************************
 # Copyright (c) 2015-23, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions
 # are met:
 #  * Redistributions of source code must retain the above copyright
 #    notice, this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright
 #    notice, this list of conditions and the following disclaimer in the
 #    documentation and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its
 #    contributors may be used to endorse or promote products derived
 #    from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
#pragma once
#include "Renderer.h"
#include "SampleManager.h"

#include <Core/SampleApp.h>

using namespace Falcor;

struct ParticleBody
{
    float3 Position = float3(0.f, 0.f, 0.f);
    float Density = 0.0f;
    // 16 bytes
    float3 Velocity = float3(0.f, 0.f, 0.f);
    float NearDensity = 1.0f;
    // 16 bytes
    float3 PredictedPosition = float3(0.f, 0.f, 0.f);
    float Pressure = 1.0f;
    // 16 bytes
    float3 Force = float3(0.f, 0.f, 0.f);
    float Mass = 0.f;
    // 16 bytes
    float SmoothingLength = 1.0f;
    float Viscosity = 0.1f;
    float pad1;
    float pad2;
    // 16 bytes
};

class FluidApplication : public SampleApp
{
public:
    FluidApplication(const SampleAppConfig& config);
    ~FluidApplication() override = default;

    void onLoad(RenderContext* pRenderContext) override;
    void onResize(uint32_t width, uint32_t height) override;

    void SimulationStep(RenderContext* pRenderContext);
    void onFrameRender(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo) override;
  
    void onGuiRender(Gui* pGui) override;
    bool onKeyEvent(const KeyboardEvent& keyEvent) override;
    bool onMouseEvent(const MouseEvent& mouseEvent) override;

private:
    void executeParticleComputePass(const ref<ComputePass>& compute_pass,
        RenderContext* pRenderContext,
        uint32_t total_threads_x,
        const uint32_t total_threads_y = 1,
        const uint32_t total_threads_z = 1) const noexcept;
    void renderPhysicsSampleGui();

    std::unique_ptr<Renderer> renderer_ = nullptr;
    SampleManager sample_manager_;

    std::vector<NodeID> sphereNodeIDs;

    float deltaTime = 0.f;

    //int density_map_size = 64;
    //float3 sim_bounds = float3(Metrics::WALLDIST) * 2.f;
    //float bounds_size = sim_bounds.x;

    World* world_ = nullptr;

    bool start_simul_ = true;
    bool regenrate_particles_ = false;

    float fixed_timer_ = kFixedDeltaTime;
    float time_since_last_fixed_update_ = 0.f;

    ref<ComputePass> spawn_particle_pass_ = nullptr;
    ref<ComputePass> compute_external_forces_pass_ = nullptr;
    ref<ComputePass> update_spatial_hash_pass_ = nullptr;
    ref<ComputePass> bitonic_sort_pass_ = nullptr;
    ref<ComputePass> calculate_offsets_pass_ = nullptr;
    ref<ComputePass> compute_neighbors_density_pass_ = nullptr;
    ref<ComputePass> compute_neighbors_pressure_pass_ = nullptr;
    ref<ComputePass> compute_neighbors_viscosity_pass_ = nullptr;
    ref<ComputePass> compute_bodies_positions_pass_ = nullptr;

    ref<Buffer> SpatialIndices = nullptr;
    ref<Buffer> SpatialOffsets = nullptr;
    ref<Buffer> bodies_buffer_ = nullptr;
    ref<Buffer> readback_spatial_indices = nullptr;
    ref<Buffer> readback_bodies_buffer_ = nullptr;
    ref<Buffer> regenrated_particles_ = nullptr;

    //ref<ComputePass> compute_density_map_pass_ = nullptr;
    //ref<Texture> density_map_;

    std::vector<ParticleBody> particle_bodies_{};

    uint32_t groupSize = 64;
    uint32_t totalThreadsX = ((Metrics::NbParticles + groupSize - 1) / groupSize) * groupSize;
};
