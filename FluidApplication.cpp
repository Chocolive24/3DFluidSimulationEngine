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

    particle_densities.resize(density_map_size * density_map_size * density_map_size);

    // for (int z = 0; z < density_map_size; ++z)
    //{
    //     for (int y = 0; y < density_map_size; ++y)
    //     {
    //         for (int x = 0; x < density_map_size; ++x)
    //         {
    //             // float3 id = float3(x, y, z);
    //             // float3 texturePos = id / float3(density_map_size - 1);      // normalized [0,1]
    //             // float3 worldPos = -(sim_bounds * 0.5f) + texturePos * sim_bounds; // map to [-100,100]

    //            //// TODO: peut être ? mettre world pos en Meters (PixelsTOMEters).

    //            // XMVECTOR xm_pos{worldPos.x, worldPos.y, worldPos.z};

    //            // float density = 1.f; // world_->CalculateDensityAtPosition(xm_pos) * 100;
    //            size_t index = x + y * density_map_size + z * density_map_size * density_map_size;
    //            particle_densities[index] = Random::Range(0.f, 1.f);
    //        }
    //    }
    //}

    renderer_->RegisterParticleDensities(&particle_densities);
    renderer_->Init();

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

        const auto sphere_node_id = renderer_->AddSphereToScene(position, 1);
        sphereNodeIDs.push_back(sphere_node_id);

        ParticleBody pb{};
        //pb.Position = XMFLOAT3{position.x, position.y, position.z};
        pb.Position = position;
        pb.Velocity = velocity;
        pb.PredictedPosition = predictedPosition;
        pb.Mass = body.Mass;
        pb.Force = force;
        particle_bodies_.push_back(pb);
    }
    
    compute_pass_ =
        ComputePass::create(getDevice(),
            "Samples/3DFluidSimulationEngine/Renderer/shaders/DensityMap.cs.slang", "updateBodies");

    bodies_buffer_ = make_ref<Buffer>(
        getDevice(),
        sizeof(particle_bodies_[0]),
        particle_bodies_.size(),
        ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
        MemoryType::DeviceLocal,
        particle_bodies_.data(),
        false
    );

    readback_bodies_buffer_ = make_ref<Buffer>(
        getDevice(),
        sizeof(particle_bodies_[0]),
        particle_bodies_.size(),
        ResourceBindFlags::None, // No need for shader access
        MemoryType::ReadBack,    // CPU-readable
        nullptr,                 // No initial data
        false
    );

    const auto compute_var = compute_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;

    // for (const auto& gd : sample_manager_.GetSampleData())
    //{
    //     if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
    //     {
    //         auto& sphere_gd = std::get<SphereF>(gd.Shape);
    //         const auto positionX = XMVectorGetX(sphere_gd.Center());
    //         const auto positionY = XMVectorGetY(sphere_gd.Center());
    //         const auto positionZ = XMVectorGetZ(sphere_gd.Center());

    //        const float3 pos = float3(positionX, positionY, positionZ);

    //        const auto sphere_node_id = renderer_->AddSphereToScene(pos, 1);
    //        sphereNodeIDs.push_back(sphere_node_id);
    //    }
    //    /*else if (gd.Shape.index() == static_cast<int>(ShapeType::Cuboid))
    //    {
    //        auto& cube_gd = std::get<CuboidF>(gd.Shape);
    //        const auto positionX = XMVectorGetX(cube_gd.Center());
    //        const auto positionY = XMVectorGetY(cube_gd.Center());
    //        const auto positionZ = XMVectorGetZ(cube_gd.Center());

    //        const float3 pos = float3(positionX, positionY, positionZ);

    //        const auto cube_node_id = renderer_->AddCubeToScene(pos);
    //    }*/
    //}

    renderer_->CreateRaytracingProgram();
}

void FluidApplication::onResize(uint32_t width, uint32_t height)
{
    renderer_->OnResize(width, height);
}

void FluidApplication::onFrameRender(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
{
    sample_manager_.UpdateSample();

    // #ifdef TRACY_ENABLE
    //     ZoneScoped;
    //
    //     for (int z = 0; z < density_map_size; ++z)
    //     {
    //         for (int y = 0; y < density_map_size; ++y)
    //         {
    //             for (int x = 0; x < density_map_size; ++x)
    //             {
    //                 float3 id = float3(x, y, z);
    //                 float3 texturePos = id / float3(density_map_size - 1);      // normalized [0,1]
    //                 float3 worldPos = -(sim_bounds * 0.5f) + texturePos * sim_bounds; // map to [-100,100]
    //
    //                 // TODO: peut être ? mettre world pos en Meters (PixelsTOMEters).
    //
    //                 XMVECTOR xm_pos{worldPos.x, worldPos.y, worldPos.z};
    //
    //                 float density = 1.f; // world_->CalculateDensityAtPosition(xm_pos) * 100;
    //                 size_t index = x + y * density_map_size + z * density_map_size * density_map_size;
    //                 particle_densities[index] = density;
    //             }
    //         }
    //     }
    // #endif

       
    constexpr uint32_t numBodies = 1024;

    const auto compute_var = compute_pass_->getRootVar();
    compute_var["bodies"] = bodies_buffer_;
    compute_var["PerFrameCB"]["deltaTime"] = 1.f / 60.f;
    compute_var["PerFrameCB"]["nbParticles"] = numBodies;
    compute_var["PerFrameCB"]["smoothingRadius"] = SPH::SmoothingRadius;

    constexpr uint32_t groupSize = 64;

    // Round up to the next multiple of 64
    constexpr uint32_t totalThreadsX = ((numBodies + groupSize - 1) / groupSize) * groupSize;

    compute_pass_->execute(pRenderContext, totalThreadsX, 1, 1);

    pRenderContext->copyResource(readback_bodies_buffer_.get(), bodies_buffer_.get());

    const ParticleBody* body_data = static_cast<const ParticleBody*>(readback_bodies_buffer_->map());

    int sphere_iterator = 0;
    for (uint32_t i = 0; i < particle_bodies_.size(); i++)
    {
        const auto pos = body_data[i].Position;

        /*const auto positionX = XMVectorGetX(pos);
        const auto positionY = XMVectorGetY(pos);
        const auto positionZ = XMVectorGetZ(pos);*/

        Transform transform;
        transform.setTranslation(pos);
        transform.setRotationEuler(float3(0.f, 0.f, 0.f));
        transform.setScaling(float3(1.f, 1.f, 1.f));

        // Update node transform
        renderer_->UpdateSceneNodeTransform(sphereNodeIDs[sphere_iterator], transform);
        sphere_iterator++;
    }

    readback_bodies_buffer_->unmap();

    // int sphere_iterator = 0;
    // for (const auto& gd : sample_manager_.GetSampleData())
    //{
    //     if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
    //     {
    //         const auto& sphere_gd = std::get<SphereF>(gd.Shape);
    //         const auto positionX = XMVectorGetX(sphere_gd.Center());
    //         const auto positionY = XMVectorGetY(sphere_gd.Center());
    //         const auto positionZ = XMVectorGetZ(sphere_gd.Center());

    //        Transform transform;
    //        transform.setTranslation(float3(positionX, positionY, positionZ));
    //        transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //        transform.setScaling(float3(1.f, 1.f, 1.f));

    //        // Update node transform
    //        renderer_->UpdateSceneNodeTransform(sphereNodeIDs[sphere_iterator], transform);
    //        sphere_iterator++;
    //    }
    //}

    renderer_->RenderFrame(pRenderContext, getGlobalClock().getTime());

    getTextRenderer().render(pRenderContext, getFrameRate().getMsg(), pTargetFbo, {20, 20});

#ifdef TRACY_ENABLE
    FrameMark;
#endif
}


void FluidApplication::onGuiRender(Gui* pGui)
{
    Gui::Window w(pGui, "Raytracing Fluid Rendering", {250, 200});
    renderGlobalUI(pGui);

    renderer_->RenderUI(pGui, &w);

    renderPhysicsSampleGui();
}

bool FluidApplication::onKeyEvent(const KeyboardEvent& keyEvent)
{
    if (keyEvent.type == KeyboardEvent::Type::KeyPressed)
    {
        if (keyEvent.key == Input::Key::F)
        {
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

    if (ImGui::Button("Regenerate"))
    {
        sample_manager_.RegenerateSample();
    }

    ImGui::SameLine();

    if (ImGui::ArrowButton("NextSample", ImGuiDir_Right))
    {
        sample_manager_.NextSample();
    }

    ImGui::End();
}
