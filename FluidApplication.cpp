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

FluidApplication::FluidApplication(const SampleAppConfig& config) : SampleApp(config)
{
    
}

FluidApplication::~FluidApplication() {}

void FluidApplication::onLoad(RenderContext* pRenderContext)
{
    sample_manager_.SetUp();
    renderer_ = std::make_unique<Renderer>(getDevice(), getTargetFbo());
    renderer_->Init();

    for (auto& gd : sample_manager_.GetSampleData())
    {
        if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
        {
            auto& sphere_gd = std::get<SphereF>(gd.Shape);
            const auto positionX = XMVectorGetX(sphere_gd.Center());
            const auto positionY = XMVectorGetY(sphere_gd.Center());
            const auto positionZ = XMVectorGetZ(sphere_gd.Center());

            const float3 pos = float3(positionX, positionY, positionZ);

            const auto sphere_node_id = renderer_->AddSphereToScene(pos, 3);
            sphereNodeIDs.push_back(sphere_node_id);
        }
    }

    renderer_->SynchronizeSceneWithProgram();
}

void FluidApplication::onResize(uint32_t width, uint32_t height)
{
    renderer_->OnResize(width, height);
}

void FluidApplication::onFrameRender(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
{
    sample_manager_.UpdateSample();

    for (size_t i = 0; i < sample_manager_.GetSampleData().size(); i++)
    {
        const auto& gd = sample_manager_.GetSampleData()[i];

        if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
        {
            const auto& sphere_gd = std::get<SphereF>(gd.Shape);
            const auto positionX = XMVectorGetX(sphere_gd.Center());
            const auto positionY = XMVectorGetY(sphere_gd.Center());
            const auto positionZ = XMVectorGetZ(sphere_gd.Center());

            Transform transform;
            transform.setTranslation(float3(positionX, positionY, positionZ));
            transform.setRotationEuler(float3(0.f, 0.f, 0.f));
            transform.setScaling(float3(1.f, 1.f, 1.f));

            // Update node transform
            renderer_->UpdateSceneNodeTransform(sphereNodeIDs[i], transform);
        }
    }

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

    static bool adjustWindow = true;

    if (adjustWindow)
    {
        ImGui::SetNextWindowSize(ImVec2(Metrics::Width / 2, static_cast<float>(Metrics::Height) / 2.5f));
        ImGui::SetNextWindowPos(ImVec2(0, 500));
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

bool FluidApplication::onKeyEvent(const KeyboardEvent& keyEvent)
{
    return renderer_->onKeyEvent(keyEvent);
}

bool FluidApplication::onMouseEvent(const MouseEvent& mouseEvent)
{
    return renderer_->onMouseEvent(mouseEvent);
}
