/***************************************************************************
 # Copyright (c) 2015-24, NVIDIA CORPORATION. All rights reserved.
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

#include "FluidApplication.h"

#include "Utils/UI/TextRenderer.h"
#include "Core/Program/Program.h"

#ifdef TRACY_ENABLE
#include <Tracy.hpp>
#include <TracyC.h>
#endif

FALCOR_EXPORT_D3D12_AGILITY_SDK

FluidApplication::FluidApplication(const SampleAppConfig& config) : SampleApp(config) {}

void FluidApplication::onLoad(RenderContext* pRenderContext)
{
    renderer_ = std::make_unique<Renderer>(getDevice(), getTargetFbo());

    sample_manager_.SetUp();
    world_ = &sample_manager_.GetWorldRef();

    renderer_->Init(getRenderContext());

    // particle_bodies_.resize(NbParticles);

    for (const auto& body : world_->_bodies)
    {
        if (!body.IsEnabled())
        {
            continue;
        }

        const float3 position{XMVectorGetX(body.Position), XMVectorGetY(body.Position), XMVectorGetZ(body.Position)};
        const float3 velocity{XMVectorGetX(body.Velocity), XMVectorGetY(body.Velocity), XMVectorGetZ(body.Velocity)};
        const float3 predictedPosition{
            XMVectorGetX(body.PredictedPosition), XMVectorGetY(body.PredictedPosition), XMVectorGetZ(body.PredictedPosition)
        };
        const float3 force{XMVectorGetX(body._force), XMVectorGetY(body._force), XMVectorGetZ(body._force)};

        if (!renderer_->useMarchingCubes)
        {
            const auto sphere_node_id = renderer_->AddSphereToScene(position, 1);
            sphereNodeIDs.push_back(sphere_node_id);
        }

        ParticleBody pb{};
        // pb.Position = XMFLOAT3{position.x, position.y, position.z};
        pb.Position = position;
        pb.Velocity = velocity;
        pb.PredictedPosition = predictedPosition;
        pb.Mass = body.Mass;
        pb.Force = force;
        particle_bodies_.push_back(pb);
    }

    spawn_particle_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "spawnParticles");

    compute_external_forces_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "computeExternalForces");

    update_spatial_hash_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "UpdateSpatialHash");

    bitonic_sort_pass_ = ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "BitonicSort");

    calculate_offsets_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "CalculateOffsets");

    compute_neighbors_density_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "computeNeighborsDensity");

    compute_neighbors_pressure_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "computeNeighborsPressure");

    compute_neighbors_viscosity_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "computeNeighborsViscosity");

    compute_bodies_positions_pass_ =
        ComputePass::create(getDevice(), "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang", "computeBodyPositions");

    bodies_buffer_ = make_ref<Buffer>(
        getDevice(),
        sizeof(ParticleBody),
        NbParticles,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
        MemoryType::DeviceLocal,
        particle_bodies_.data(),
        false
    );

    SpatialIndices = make_ref<Buffer>(
        getDevice(),
        sizeof(uint32_t) * 3,
        NbParticles,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
        MemoryType::DeviceLocal,
        nullptr,
        false
    );

    SpatialOffsets = make_ref<Buffer>(
        getDevice(),
        sizeof(uint32_t),
        NbParticles,
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
        MemoryType::DeviceLocal,
        nullptr,
        false
    );

    readback_bodies_buffer_ = make_ref<Buffer>(
        getDevice(),
        sizeof(ParticleBody),
        NbParticles,
        ResourceBindFlags::None, // No need for shader access
        MemoryType::ReadBack,    // CPU-readable
        nullptr,                 // No initial data
        false
    );

    readback_spatial_indices = make_ref<Buffer>(
        getDevice(),
        sizeof(uint3),
        NbParticles,
        ResourceBindFlags::None, // No need for shader access
        MemoryType::ReadBack,    // CPU-readable
        nullptr,                 // No initial data
        false
    );

    regenrated_particles_ = make_ref<Buffer>(
        getDevice(),
        sizeof(ParticleBody),
        NbParticles,
        ResourceBindFlags::None,
        MemoryType::Upload, // Must be Upload to map from CPU
        particle_bodies_.data(),
        false
    );

    auto compute_var = spawn_particle_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    // executeParticleComputePass(spawn_particle_pass_, pRenderContext, totalThreadsX);

    compute_var = compute_external_forces_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = update_spatial_hash_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = bitonic_sort_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = calculate_offsets_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = compute_neighbors_density_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = compute_neighbors_pressure_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = compute_neighbors_viscosity_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var = compute_bodies_positions_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    renderer_->CreateRaytracingProgram(getRenderContext());

    //renderer_->CreateRasterizationProgram();

    // renderer_->LaunchMarchingCubeComputePasses(getRenderContext());
}

void FluidApplication::onResize(uint32_t width, uint32_t height)
{
    renderer_->OnResize(width, height);
}

void FluidApplication::executeParticleComputePass(
    const ref<ComputePass>& compute_pass,
    RenderContext* pRenderContext,
    const uint32_t total_threads_x,
    const uint32_t total_threads_y,
    const uint32_t total_threads_z
) const noexcept
{
    const auto compute_var = compute_pass->getRootVar();
    // compute_var["gTexture3D"] = density_map_;
    compute_var["bodies"] = bodies_buffer_;
    compute_var["SpatialIndices"] = SpatialIndices;
    compute_var["SpatialOffsets"] = SpatialOffsets;

    compute_var["PerFrameCB"]["densityMapSize"] = Metrics::density_map_size;
    compute_var["PerFrameCB"]["simBounds"] = float3(Metrics::sim_bounds);
    compute_var["PerFrameCB"]["wallDist"] = Metrics::WALLDIST;
    compute_var["PerFrameCB"]["fixedDeltaTime"] = kFixedDeltaTime;
    compute_var["PerFrameCB"]["nbParticles"] = Metrics::NbParticles;
    compute_var["PerFrameCB"]["gravity"] = world_->Gravity;
    compute_var["PerFrameCB"]["smoothingRadius"] = SPH::SmoothingRadius;
    compute_var["PerFrameCB"]["targetDensity"] = SPH::TargetDensity;
    compute_var["PerFrameCB"]["pressureMultiplier"] = SPH::PressureMultiplier;
    compute_var["PerFrameCB"]["viscosityStrength"] = SPH::ViscosityStrength;
    compute_var["PerFrameCB"]["densityGraphicsMultiplier"] = densityGraphicsMultiplier;


    // Remove scaling for physics collision test
    Transform unscaledLocalToWorld = renderer_->fluid_transform;
    unscaledLocalToWorld.setScaling(renderer_->scale / WALLDIST); // Or normalize the axes manually
    float4x4 unscaledWorldToLocal = inverse(unscaledLocalToWorld.getMatrix());

    const float4x4 localToWorld = renderer_->fluid_transform.getMatrix();
    const float4x4 worldToLocal = inverse(localToWorld);

    compute_var["PerFrameCB"]["unscaledLocalToWorld"] = unscaledLocalToWorld.getMatrix();
    compute_var["PerFrameCB"]["unscaledWorldToLocal"] = unscaledWorldToLocal;
   

    compute_var["PerFrameCB"]["localToWorld"] = localToWorld;
    compute_var["PerFrameCB"]["worldToLocal"] = worldToLocal;

    //compute_var["PerFrameCB"]["ScaledSimBounds"] = renderer_->fluid_transform.getScaling() * 2.f;

    const float r = SPH::SmoothingRadius;
    const float spikyPow2 = 15.f / (2 * PI * Pow(r, 5));
    const float spikyPow3 = 15.f / (PI * Pow(r, 6));
    const float spikyPow2Grad = 15.f / (PI * Pow(r, 5));
    const float spikyPow3Grad = 45.f / (PI * Pow(r, 6));

    compute_var["PerFrameCB"]["K_SpikyPow2"] = spikyPow2;
    compute_var["PerFrameCB"]["K_SpikyPow3"] = spikyPow3;
    compute_var["PerFrameCB"]["K_SpikyPow2Grad"] = spikyPow2Grad;
    compute_var["PerFrameCB"]["K_SpikyPow3Grad"] = spikyPow3Grad;

    compute_var["PerFrameCB"]["collisionDamping"] = SPH::collisionDamping;

    compute_pass->execute(pRenderContext, total_threads_x, total_threads_y, total_threads_z);
}

uint32_t NextPowerOfTwo(int n)
{
    if (n == 0)
        return 1;

    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return n + 1;
}

void FluidApplication::onFrameRender(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
{
    if (regenrate_particles_)
    {
        pRenderContext->copyResource(bodies_buffer_.get(), regenrated_particles_.get());
        regenrate_particles_ = false;
    }

    if (start_simul_)
    {
        const auto delta_time = getGlobalClock().getDelta();
        fixed_timer_ += delta_time;
        time_since_last_fixed_update_ += delta_time;

        while (fixed_timer_ >= kFixedDeltaTime)
        {
            executeParticleComputePass(compute_external_forces_pass_, pRenderContext, totalThreadsX);
            executeParticleComputePass(update_spatial_hash_pass_, pRenderContext, totalThreadsX);

            const auto compute_var = bitonic_sort_pass_->getRootVar();
            // compute_var["gTexture3D"] = density_map_;
            compute_var["bodies"] = bodies_buffer_;
            compute_var["SpatialIndices"] = SpatialIndices;
            compute_var["SpatialOffsets"] = SpatialOffsets;

            compute_var["PerFrameCB"]["densityMapSize"] = Metrics::density_map_size;
            compute_var["PerFrameCB"]["simBounds"] = float3(Metrics::sim_bounds);
            compute_var["PerFrameCB"]["wallDist"] = Metrics::WALLDIST;
            compute_var["PerFrameCB"]["fixedDeltaTime"] = kFixedDeltaTime;
            compute_var["PerFrameCB"]["nbParticles"] = Metrics::NbParticles;
            compute_var["PerFrameCB"]["gravity"] = world_->Gravity;
            compute_var["PerFrameCB"]["smoothingRadius"] = SPH::SmoothingRadius;
            compute_var["PerFrameCB"]["targetDensity"] = SPH::TargetDensity;
            compute_var["PerFrameCB"]["pressureMultiplier"] = SPH::PressureMultiplier;
            compute_var["PerFrameCB"]["viscosityStrength"] = SPH::ViscosityStrength;
            compute_var["PerFrameCB"]["densityGraphicsMultiplier"] = densityGraphicsMultiplier;

            // Launch each step of the sorting algorithm (once the previous step is complete)
            // Number of steps = [log2(n) * (log2(n) + 1)] / 2
            // where n = nearest power of 2 that is greater or equal to the number of inputs
            int numStages = log2(NextPowerOfTwo(NbParticles));

            for (int stageIndex = 0; stageIndex < numStages; stageIndex++)
            {
                for (int stepIndex = 0; stepIndex < stageIndex + 1; stepIndex++)
                {
                    // Calculate some pattern stuff
                    int groupWidth = 1 << (stageIndex - stepIndex);
                    int groupHeight = 2 * groupWidth - 1;
                    compute_var["PerFrameCB"]["groupWidth"] = groupWidth;
                    compute_var["PerFrameCB"]["groupHeight"] = groupHeight;
                    compute_var["PerFrameCB"]["stepIndex"] = stepIndex;

                    bitonic_sort_pass_->execute(pRenderContext, NextPowerOfTwo(NbParticles) / 2, 1, 1);
                }
            }

            executeParticleComputePass(calculate_offsets_pass_, pRenderContext, totalThreadsX);

            executeParticleComputePass(compute_neighbors_density_pass_, pRenderContext, totalThreadsX);
            executeParticleComputePass(compute_neighbors_pressure_pass_, pRenderContext, totalThreadsX);
            executeParticleComputePass(compute_neighbors_viscosity_pass_, pRenderContext, totalThreadsX);

            executeParticleComputePass(compute_bodies_positions_pass_, pRenderContext, totalThreadsX);

            fixed_timer_ -= kFixedDeltaTime;
            time_since_last_fixed_update_ = 0.f;
        }
    }

    if (!renderer_->useMarchingCubes)
    {
        if (!renderer_->draw_fluid_)
        {
            pRenderContext->copyResource(readback_bodies_buffer_.get(), bodies_buffer_.get());

            const ParticleBody* body_data = static_cast<const ParticleBody*>(readback_bodies_buffer_->map());

            int sphere_iterator = 0;
            for (uint32_t i = 0; i < particle_bodies_.size(); i++)
            {
                const auto pos = body_data[i].Position;

                Transform transform;
                transform.setTranslation(pos);
                transform.setRotationEuler(float3(0.f, 0.f, 0.f));
                transform.setScaling(float3(1.f, 1.f, 1.f));

                // Update node transform
                renderer_->UpdateSceneNodeTransform(sphereNodeIDs[sphere_iterator], transform);
                sphere_iterator++;
            }

            readback_bodies_buffer_->unmap();
        }
        else
        {
            int sphere_iterator = 0;
            for (uint32_t i = 0; i < particle_bodies_.size(); i++)
            {
                Transform transform;
                transform.setTranslation(float3(0, -10'000, 0));
                transform.setRotationEuler(float3(0.f, 0.f, 0.f));
                transform.setScaling(float3(1.f, 1.f, 1.f));

                // Update node transform
                renderer_->UpdateSceneNodeTransform(sphereNodeIDs[sphere_iterator], transform);
                sphere_iterator++;
            }
        }
    }

    renderer_->RenderFrame(pRenderContext, pTargetFbo, getGlobalClock().getTime(), bodies_buffer_, SpatialIndices, SpatialOffsets);

    getTextRenderer().render(pRenderContext, getFrameRate().getMsg(), pTargetFbo, {20, 20});

#ifdef TRACY_ENABLE
    FrameMark;
#endif
}


void FluidApplication::onGuiRender(Gui* pGui)
{
    Gui::Window w(pGui, "Raytracing Fluid Rendering", {250, 200});
    renderGlobalUI(pGui);

    renderer_->RenderUI(pGui, &w, getRenderContext());

    renderPhysicsSampleGui();
}

bool FluidApplication::onKeyEvent(const KeyboardEvent& keyEvent)
{
    if (keyEvent.type == KeyboardEvent::Type::KeyPressed)
    {
        if (keyEvent.key == Input::Key::F)
        {
            start_simul_ = !start_simul_;
            sample_manager_.StopSample();
        }
        else if (keyEvent.key == Input::Key::R)
        {
            sample_manager_.RegenerateSample();
        }
    }

    return renderer_->onKeyEvent(keyEvent);
}

bool FluidApplication::onMouseEvent(const MouseEvent& mouseEvent)
{
    XMVECTOR mouse_pos = XMVectorSet(mouseEvent.screenPos.x, mouseEvent.screenPos.y, 0.f, 0.f);
    sample_manager_.GiveMousePositionToSample(mouse_pos);

    return renderer_->onMouseEvent(mouseEvent);
}

void FluidApplication::renderPhysicsSampleGui()
{
    static bool adjustWindow = true;

    if (adjustWindow)
    {
        ImGui::SetNextWindowSize(ImVec2(Metrics::Width / 2.f, static_cast<float>(Metrics::Height) / 3.f));
        ImGui::SetNextWindowPos(ImVec2(Metrics::Width / 2.f, 0.f));
        adjustWindow = false;
    }

    ImGui::Begin("Sample Manager");

    if (ImGui::BeginCombo("Select a Sample", sample_manager_.GetSampleName(sample_manager_.GetCurrentIndex()).c_str()))
    {
        for (std::size_t index = 0; index < sample_manager_.GetSampleNbr(); index++)
        {
            if (ImGui::Selectable(sample_manager_.GetSampleName(index).c_str(), sample_manager_.GetCurrentIndex() == index))
            {
                sample_manager_.ChangeSample(index);
            }
        }
        ImGui::EndCombo();
    }

    if (ImGui::Button("Regenerate"))
    {
        // sample_manager_.RegenerateSample();

        // std::cout << "regen\n";

        regenrate_particles_ = true;
        std::vector<XMVECTOR> particlePositions;

        for (size_t i = 0; i < NbParticles;)
        {
            XMVECTOR pos = XMVectorSet(
                Random::Range(-WALLDIST, WALLDIST * 0.2f),
                Random::Range(-WALLDIST * 0.8f, WALLDIST * 0.8f),
                Random::Range(-WALLDIST, WALLDIST * 0.2f),
                0.0f
            );

            bool overlaps = false;
            for (const auto& existing : particlePositions)
            {
                if (XMVectorGetX(XMVector3LengthSq(pos - existing)) < (PARTICLESIZE * PARTICLESIZE * 4))
                {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps)
            {
                particlePositions.push_back(pos);
                float x = XMVectorGetX(pos);
                float y = XMVectorGetY(pos);
                float z = XMVectorGetZ(pos);
                /*particle_bodies_[i] = ParticleBody{};
                particle_bodies_[i].Position = float3(x, y, z);*/
                ++i;
            }
        }

        // std::cout << "End regen\n";
    }


    ImGui::Spacing();

    ImGui::TextWrapped(sample_manager_.GetSampleDescription(sample_manager_.GetCurrentIndex()).c_str());

    ImGui::Spacing();

    sample_manager_.DrawImgui(sample_manager_.GetCurrentIndex());

    ImGui::SetCursorPosY(ImGui::GetWindowHeight() - (ImGui::GetFrameHeightWithSpacing()));

    if (ImGui::ArrowButton("PreviousSample", ImGuiDir_Left))
    {
        sample_manager_.PreviousSample();
    }

    ImGui::SameLine();

    if (ImGui::ArrowButton("NextSample", ImGuiDir_Right))
    {
        sample_manager_.NextSample();
    }

    ImGui::End();
}
