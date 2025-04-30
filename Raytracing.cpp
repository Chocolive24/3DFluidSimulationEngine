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
#include "Raytracing.h"

#include "SampleManager.h"

#include "Rendering/Lights/EnvMapSampler.h"
#include "Scene/SceneBuilder.h"
#include "Utils/Math/FalcorMath.h"
#include "Utils/UI/TextRenderer.h"
#include "Core/Program/Program.h"

#ifdef TRACY_ENABLE
#include <Tracy.hpp>
#include <TracyC.h>
#endif 

FALCOR_EXPORT_D3D12_AGILITY_SDK

SampleManager sample_manager;
static SceneBuilder* scene_builder;

uint32_t mSampleGuiWidth = 250;
uint32_t mSampleGuiHeight = 200;
uint32_t mSampleGuiPositionX = 20;
uint32_t mSampleGuiPositionY = 40;

static float3 kClearColor(.2, 1, .1);

static float waterTurbulence = 2.5f;

static uint kMaxRayBounce = 5;

static float3 absorptionCoeff(1.0, 0.4, 0.05);
static float3 scatteringCoeff(0.1, 0.2, 0.8);
static float phaseG = 0.8f;

static float maxRayMarchingDistance = 5.f;
static float kMarchSize = 0.1f;

static float maxLighMarchingDistance = 3.f;
static float sunLightMarchSize = 0.2f;

float3 lightColor = float3(1, 1, 1);
float3 lightDir = normalize(float3(1, -1, -1));

static float IoR = 1.33f;

static const std::string kEnvMapPath = "hallstatt4_hd.hdr";

Raytracing::Raytracing(const SampleAppConfig& config) : SampleApp(config)
{
    
}

Raytracing::~Raytracing() {}

void Raytracing::onLoad(RenderContext* pRenderContext)
{
    sample_manager.SetUp();
    renderer_ = std::make_unique<Renderer>(getDevice(), getTargetFbo());
    renderer_->Init();

//    static constexpr std::array<Vertex, 3> vertices
//    {
//                Vertex{float3(0.0f, 0.5f, 0.0f)},  // Top vertex (red)
//                Vertex{float3(0.5f, -0.5f, 0.0f)}, // Bottom right (green)
//                Vertex{float3(-0.5f, -0.5f, 0.0f)} // Bottom left (blue)
//    };
//        
//    vertex_buffer_ = getDevice()->createBuffer(
//                sizeof(Vertex) * vertices.size(),
//                ResourceBindFlags::Vertex,
//                MemoryType::Upload,
//                vertices.data());
//
//    ref<VertexLayout> pLayout = VertexLayout::create();
//    ref<VertexBufferLayout> pBufLayout = VertexBufferLayout::create();
//    pBufLayout->addElement("POSITION", 0, ResourceFormat::RGB32Float, 1, 0);
//    pLayout->addBufferLayout(0, pBufLayout);
//    Vao::BufferVec buffers{vertex_buffer_};
//
//    vao_ = Vao::create(Vao::Topology::TriangleList, pLayout, buffers);
//
//    // Create the RenderState
//    raster_pass_ = RasterPass::create(getDevice(),
//        "Samples/Raytracing/triangle.slang", "vsMain", "psMain");
//    auto& pState = raster_pass_->getState();
//
//    // create the depth-state
//    DepthStencilState::Desc dsDesc;
//    dsDesc.setDepthEnabled(false);
//    pState->setDepthStencilState(DepthStencilState::create(dsDesc));
//
//    // Rasterizer state
//    RasterizerState::Desc rsState;
//    rsState.setCullMode(RasterizerState::CullMode::None);
//    pState->setRasterizerState(RasterizerState::create(rsState));

    // Blend state
    //BlendState::Desc blendDesc;
    //blendDesc.setRtBlend(0, true).setRtParams(
    //    0,
    //    BlendState::BlendOp::Add,
    //    BlendState::BlendOp::Add,
    //    BlendState::BlendFunc::SrcAlpha,
    //    BlendState::BlendFunc::OneMinusSrcAlpha,
    //    BlendState::BlendFunc::One,
    //    BlendState::BlendFunc::One
    //);
    //pState->setBlendState(BlendState::create(blendDesc));

    

    //if (getDevice()->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    //{
    //    FALCOR_THROW("Device does not support raytracing!");
    //}

    //Settings settings{};

    //// Create the SceneBuilder
    //SceneBuilder::Flags flags =
    //    SceneBuilder::Flags::RTDontMergeStatic | SceneBuilder::Flags::RTDontMergeDynamic |
    //                            SceneBuilder::Flags::RTDontMergeInstanced |
    //                            SceneBuilder::Flags::DontOptimizeGraph;
    //scene_builder = new SceneBuilder(getDevice(), settings, flags);

    //// Define triangle vertices (position, normal, and texture coordinates)
    //TriangleMesh::VertexList rt_vertices = {
    //    {float3(-0.5f, -0.5f, 0.0f), float3(0, 0, 1), float2(0, 0)}, // v0
    //    {float3(0.5f, -0.5f, 0.0f), float3(0, 0, 1), float2(1, 0)},  // v1
    //    {float3(0.0f, 0.5f, 0.0f), float3(0, 0, 1), float2(0.5, 1)}  // v2
    //};

    //// Define triangle indices (single triangle)
    //std::vector<uint32_t> indices = {0, 1, 2};

    //// Create a TriangleMesh object
    ////auto pTriangleMesh = TriangleMesh::create(rt_vertices, indices);

    //auto sphere = TriangleMesh::createSphere(3.f);
    //auto cube = TriangleMesh::createCube(float3(1.f, 1.f, 1.f));

    //// Create a lambertian material
    //ref<Material> lambertian = StandardMaterial::create(getDevice(), "Lambertian");
    ////ref<Texture> checkerTile = Texture::createFromFile(getDevice(), "CheckerTile_BaseColor.png", true, false);
    ////lambertian->setTexture(Material::TextureSlot::BaseColor, checkerTile);
    //lambertian->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
    //lambertian->setRoughnessMollification(1.f);
    //lambertian->setIndexOfRefraction(0.f);

    //ref<Material> water_particle_mat = StandardMaterial::create(getDevice(), "water particle");
    ////ref<Texture> checkerTile = Texture::createFromFile(getDevice(), "CheckerTile_BaseColor.png", true, false);
    ////lambertian->setTexture(Material::TextureSlot::BaseColor, checkerTile);
    //water_particle_mat->toBasicMaterial()->setBaseColor3(float3(0.1f, 0.2f, 1.0f));
    //water_particle_mat->setRoughnessMollification(1.f);
    //water_particle_mat->setIndexOfRefraction(0.f);

    ////lambertian->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));

    ////ref<Material> dielectric_red = StandardMaterial::create(getDevice(), "DielecRed");
    ////dielectric_red->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));
    ////dielectric_red->setDoubleSided(true);
    ////dielectric_red->setIndexOfRefraction(1.f);
    //////dielectric_red->toBasicMaterial()->setDiffuseTransmission(1.f);

    //ref<Material> dielectric_blue = StandardMaterial::create(getDevice(), "DielecBlue");
    //dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    //dielectric_blue->setDoubleSided(true);
    //dielectric_blue->setIndexOfRefraction(1.f);
    //dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    ////auto triangle_mesh_id_1 = scene_builder_.addTriangleMesh(sphere, dielectric_red);
    ////auto triangle_mesh_id_2 = scene_builder_.addTriangleMesh(cube, dielectric_blue);
    ////auto triangle_mesh_id_3 = scene_builder_.addTriangleMesh(cube, dielectric_blue);
    //auto triangle_mesh_id_4 = scene_builder->addTriangleMesh(cube, lambertian);
    //auto sphere_mesh_id = scene_builder->addTriangleMesh(sphere, dielectric_blue);

    ////AABB raymarching_AABB = AABB(float3(-50, -5, -5) + float3(0, 10, 0),
    ////    float3(50, 5, 5) + float3(0, 10, 0));
    ////uint32_t raymarching_AABB_ID = 1;
    ////scene_builder_.addCustomPrimitive(raymarching_AABB_ID, raymarching_AABB);

    ////auto raymarching_node = SceneBuilder::Node();
    ////raymarching_node.name = "RaymarchingNode";
    ////auto raymarching_transform = Transform();
    ////raymarching_transform.setTranslation(float3(0.f, 10.f, 0.f));
    ////raymarching_transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    ////raymarching_transform.setScaling(float3(1.f, 1.f, 1.f));
    ////raymarching_node.transform = raymarching_transform.getMatrix();
    ////auto raymarching_node_id = scene_builder_.addNode(raymarching_node);

    ////auto node = SceneBuilder::Node();
    ////node.name = "Sphere1";
    ////auto transform = Transform();
    ////transform.setTranslation(float3(0.f, 0.f, 0.f));
    ////transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    ////transform.setScaling(float3(1.f, 1.f, 1.f));
    ////node.transform = transform.getMatrix();
    ////auto node_id = scene_builder->addNode(node);

    ////// Add Mesh Instances
    ////scene_builder->addMeshInstance(node_id, triangle_mesh_id_1);

    ////auto node_2 = SceneBuilder::Node();
    ////node_2.name = "Sphere1";
    ////auto transform_2 = Transform();
    ////transform_2.setTranslation(float3(10.f, 0.f, 0.f));
    ////transform_2.setRotationEuler(float3(0.f, 0.f, 0.f));
    ////transform_2.setScaling(float3(1.f));
    ////node_2.transform = transform_2.getMatrix();
    ////auto node_id_2 = scene_builder_.addNode(node_2);

    ////// Add Mesh Instances
    ////scene_builder_.addMeshInstance(node_id_2, triangle_mesh_id_2);

    ////auto node_3 = SceneBuilder::Node();
    ////node_3.name = "cube";
    ////auto transform_3 = Transform();
    ////transform_3.setTranslation(float3(10.f, 0.f, 10.f));
    ////transform_3.setRotationEuler(float3(0.f, 0.f, 0.f));
    ////transform_3.setScaling(float3(10.f, 1.f, 10.f));
    ////node_3.transform = transform_3.getMatrix();
    ////auto node_id_3 = scene_builder_.addNode(node_3);

    ////// Add Mesh Instances
    ////scene_builder_.addMeshInstance(node_id_3, triangle_mesh_id_3);

    ////auto node_4 = SceneBuilder::Node();
    ////node_4.name = "cube";
    ////auto transform_4 = Transform();
    ////transform_4.setTranslation(float3(0.f, -1.2f, 0.f));
    ////transform_4.setRotationEuler(float3(0.f, 0.f, 0.f));
    ////transform_4.setScaling(float3(1000.f, 1.f, 1000.f));
    ////node_4.transform = transform_4.getMatrix();
    ////auto node_id_4 = scene_builder->addNode(node_4);

    ////// Add Mesh Instances
    ////scene_builder->addMeshInstance(node_id_4, triangle_mesh_id_4);

    //int i = 0;
    //for (auto& gd : sample_manager.GetSampleData())
    //{
    //    
    //    if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
    //    {
    //        i++;
    //        auto& sphere_gd = std::get<SphereF>(gd.Shape);
    //        const auto positionX = XMVectorGetX(sphere_gd.Center());
    //        const auto positionY = XMVectorGetY(sphere_gd.Center());
    //        const auto positionZ = XMVectorGetZ(sphere_gd.Center());

    //        auto node = SceneBuilder::Node();
    //        std::string name = "Sphere " + std::to_string(i);
    //        node.name = name;
    //        auto transform = Transform();
    //        transform.setTranslation(float3(positionX, positionY, positionZ));
    //        transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //        transform.setScaling(float3(1.f, 1.f, 1.f));
    //        node.transform = transform.getMatrix();
    //        auto node_id = scene_builder->addNode(node);

    //        sphereNodeIDs.push_back(node_id);

    //        // Add Mesh Instances
    //        scene_builder->addMeshInstance(node_id, sphere_mesh_id);
    //    }
    //}

    //for (int j = 0; j < 10; j++)
    //{
    //    auto node = SceneBuilder::Node();
    //    std::string name = "Sphere " + std::to_string(i);
    //    node.name = name;
    //    auto transform = Transform();
    //    transform.setTranslation(float3(1, 0, 1));
    //    transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //    transform.setScaling(float3(1.f, 1.f, 1.f));
    //    node.transform = transform.getMatrix();
    //    auto node_id = scene_builder->addNode(node);

    //    // Add Mesh Instances
    //    scene_builder->addMeshInstance(node_id, sphere_mesh_id);
    //}

  //  auto envMap = EnvMap::createFromFile(getDevice(), "hallstatt4_hd.hdr");
  //  envMap->setIntensity(1.0);
  //  scene_builder->setEnvMap(envMap);

  //  ref<Camera> camera = ref<Camera>(new Camera("Camera"));
  //  camera->setPosition(float3(0, 0.0, -250));
  //  camera->setTarget(float3(0, 0.0, 0));
  //  camera->setUpVector(float3(0, 1, 0));
  //  camera->setFocalLength(35);
  //  camera->setDepthRange(0.1f, 10000.f);

  //  scene_builder->addCamera(camera);

  //  mpScene = scene_builder->getScene();

  //  //mpScene = Scene::create(getDevice(), kDefaultScene);
  //  
  //  mpCamera = mpScene->getCamera();

  //  // Update the controllers
  //  float radius = mpScene->getSceneBounds().radius();
  //  mpScene->setCameraSpeed(50.f);
  ///*  float nearZ = std::max(0.1f, radius / 750.0f);
  //  float farZ = radius * 10;*/
  //  /*mpCamera->setDepthRange(nearZ, farZ);*/
  //  auto pTargetFbo = getTargetFbo().get();
  //  mpCamera->setAspectRatio(static_cast<float>(pTargetFbo->getWidth()) / static_cast<float>(pTargetFbo->getHeight()));

  //  //mpScene->setBlasUpdateMode(Scene::UpdateMode::Rebuild);

  //  // Get shader modules and type conformances for types used by the scene.
  //  // These need to be set on the program in order to use Falcor's material system.
  //  auto shaderModules = mpScene->getShaderModules();
  //  auto typeConformances = mpScene->getTypeConformances();

  //  // Get scene defines. These need to be set on any program using the scene.
  //  auto defines = mpScene->getSceneDefines();

  // // Create a triangle using a Scene object
  //  // Create a raytracing program description
  //  ProgramDesc rtProgDesc;
  //  ProgramDesc::ShaderModule s_module = ProgramDesc::ShaderModule("Samples/Raytracing/SDF_Functions.slang");
  //  s_module.addFile("Samples/Raytracing/SDF_Functions.slang");
  //  shaderModules.emplace_back(s_module);
  //  rtProgDesc.addShaderModules(shaderModules);
  //  rtProgDesc.addShaderLibrary("Samples/Raytracing/Raytracing.rt.slang");
  //  rtProgDesc.addTypeConformances(typeConformances);
  //  rtProgDesc.setMaxTraceRecursionDepth(kMaxRayBounce + 1); 
  //  rtProgDesc.setMaxPayloadSize(48);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
  //                                           // should be set as small as possible for maximum performance.

  //  ref<RtBindingTable> sbt = RtBindingTable::create(1, 1, mpScene->getGeometryCount());
  //  sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
  //  sbt->setMiss(0, rtProgDesc.addMiss("miss"));

  //  auto primary = rtProgDesc.addHitGroup("closestHit", "");
  //  sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);

  //  auto raymarching_hit_group = rtProgDesc.addHitGroup(
  //      "RaymarchingClosestHit", "", "RaymarchingIntersection");
  //  sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::Custom), raymarching_hit_group);

  //  rt_program_ = Program::create(getDevice(), rtProgDesc, defines);
  //  rt_program_vars_ = RtProgramVars::create(getDevice(), rt_program_, sbt);



    //// Create the raytracing program
    //auto rtProgram = Program::create(getDevice(), rtProgDesc);

    //// Create shader variables (for passing data to shaders)
    //auto rtVars = RtProgramVars::create(getDevice(), rtProgram);
 /*   if (getDevice()->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    {
        FALCOR_THROW("Device does not support raytracing!");
    }

    loadScene(kDefaultScene, getTargetFbo().get());*/
}

void Raytracing::onResize(uint32_t width, uint32_t height)
{
    renderer_->OnResize(width, height);

    //const float h = static_cast<float>(height);
    //const float w = static_cast<float>(width);

    //if (mpCamera)
    //{
    //    mpCamera->setFocalLength(18);
    //    const float aspectRatio = (w / h);
    //    mpCamera->setAspectRatio(aspectRatio);
    //}

    //mpRtOut = getDevice()->createTexture2D(
    //    width, height, ResourceFormat::RGBA16Float, 1, 1, nullptr, ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    //);
}

void Raytracing::onFrameRender(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
{
    //pRenderContext->clearFbo(pTargetFbo.get(), float4(kClearColor, 1), 1.0f, 0, FboAttachmentType::All);

    sample_manager.UpdateSample();

    //for (size_t i = 0; i < sample_manager.GetSampleData().size(); i++)
    //{
    //    const auto& gd = sample_manager.GetSampleData()[i];

    //    if (gd.Shape.index() == static_cast<int>(ShapeType::Sphere))
    //    {
    //        const auto& sphere_gd = std::get<SphereF>(gd.Shape);
    //        const auto positionX = XMVectorGetX(sphere_gd.Center());
    //        const auto positionY = XMVectorGetY(sphere_gd.Center());
    //        const auto positionZ = XMVectorGetZ(sphere_gd.Center());

    //        Transform transform;
    //        transform.setTranslation(float3(positionX, positionY, positionZ));
    //        transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //        transform.setScaling(float3(1.f, 1.f, 1.f));

    //        // Update node transform
    //        mpScene->updateNodeTransform(sphereNodeIDs[i].get(), transform.getMatrix());
    //    }
    //}

    renderer_->RenderFrame(pRenderContext, getGlobalClock().getTime());

    //raster_pass_->getState()->setVao(vao_);
    //raster_pass_->getState()->setFbo(pTargetFbo);
    //raster_pass_->draw(pRenderContext, 3, 0);

    //IScene::UpdateFlags updates = mpScene->update(pRenderContext, getGlobalClock().getTime());
    //if (is_set(updates, IScene::UpdateFlags::GeometryChanged))
    //    FALCOR_THROW("This sample does not support scene geometry changes.");
    //if (is_set(updates, IScene::UpdateFlags::RecompileNeeded))
    //    FALCOR_THROW("This sample does not support scene changes that require shader recompilation.");

    //FALCOR_ASSERT(mpScene)
    ////FALCOR_PROFILE(pRenderContext, "renderRT");

    //auto var = rt_program_vars_->getRootVar();
    //var["PerFrameCB"]["invView"] = inverse(mpCamera->getViewMatrix());
    //var["PerFrameCB"]["viewportDims"] = float2(pTargetFbo->getWidth(), pTargetFbo->getHeight());
    //float fovY = focalLengthToFovY(mpCamera->getFocalLength(), Camera::kDefaultFrameHeight);
    //var["PerFrameCB"]["tanHalfFovY"] = std::tan(fovY * 0.5f);
    //var["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
    //var["PerFrameCB"]["useDOF"] = mUseDOF;

    //var["PerFrameCB"]["backgroundColor"] = kClearColor;

    //var["PerFrameCB"]["waterTurbulence"] = waterTurbulence;

    //var["PerFrameCB"]["maxRayBounce"] = kMaxRayBounce;

    //var["PerFrameCB"]["absorptionCoeff"] = absorptionCoeff;
    //var["PerFrameCB"]["scatteringCoeff"] = scatteringCoeff;
    //var["PerFrameCB"]["phaseG"] = phaseG;
    //
    //var["PerFrameCB"]["maxRaymarchingDistance"] = maxRayMarchingDistance;
    //var["PerFrameCB"]["marchSize"] = kMarchSize;

    //var["PerFrameCB"]["maxLighMarchingDistance"] = maxLighMarchingDistance;
    //var["PerFrameCB"]["sunLightMarchSize"] = sunLightMarchSize;

    //var["PerFrameCB"]["lightColor"] = lightColor;
    //var["PerFrameCB"]["lightDir"] = lightDir;

    //var["PerFrameCB"]["IoR"] = IoR;

    //var["PerFrameCB"]["time"] = static_cast<float>(getGlobalClock().getTime());
    //static int frame = 0;
    //var["PerFrameCB"]["iFrame"] = frame++;

    //var["gOutput"] = mpRtOut;
    //var["gTexture3D"] = mpTexture3D;
    //

    //pRenderContext->clearUAV(mpRtOut->getUAV().get(), float4(kClearColor, 1));
    //mpScene->raytrace(pRenderContext, rt_program_.get(), rt_program_vars_, uint3(pTargetFbo->getWidth(), pTargetFbo->getHeight(), 1));
    //pRenderContext->blit(mpRtOut->getSRV(), pTargetFbo->getRenderTargetView(0));

    getTextRenderer().render(pRenderContext, getFrameRate().getMsg(), pTargetFbo, {20, 20});

#ifdef TRACY_ENABLE
    FrameMark;
#endif
}

void Raytracing::onGuiRender(Gui* pGui)
{
    ////Gui::Window w(pGui, "Hello DXR Settings", {300, 400}, {10, 80});

    ////w.checkbox("Ray Trace", mRayTrace);
    ////w.checkbox("Use Depth of Field", mUseDOF);
    ////if (w.button("Load Scene"))
    ////{
    ////    std::filesystem::path path;
    ////    if (openFileDialog(Scene::getFileExtensionFilters(), path))
    ////    {
    ////        loadScene(path, getTargetFbo().get());
    ////    }
    ////}

    ////mpScene->renderUI(w);

    Gui::Window w(pGui, "Raytracing Fluid Rendering", {250, 200});
    renderGlobalUI(pGui);

    renderer_->RenderUI(pGui, &w);

    //w.rgbColor("Background color", kClearColor);

    //w.var("Water Turbulance", waterTurbulence);

    //w.var("MaxRayBounce", kMaxRayBounce);

    //w.var("absorptionCoeff", absorptionCoeff);
    //w.var("scatteringCoeff", scatteringCoeff);
    //w.var("Phase G ", phaseG);
    //
    //w.var("maxRaymarchingDistance", maxRayMarchingDistance);
    //w.var("MarchSize", kMarchSize);
    //w.var("maxLighMarchingDistance", maxLighMarchingDistance);
    //w.var("sunLightMarchSize", sunLightMarchSize);

    //w.rgbColor("Light color", lightColor);
    //static float3 ImGUI_LightDir = lightDir;
    //w.var("Light Direction", ImGUI_LightDir);
    //lightDir = math::normalize(ImGUI_LightDir);

    //w.var("IoR", IoR);

    //w.checkbox("Use Depth of Field", mUseDOF);

    //mpScene->renderUI(w);

    static bool adjustWindow = true;

    if (adjustWindow)
    {
        ImGui::SetNextWindowSize(ImVec2(Metrics::Width / 2, static_cast<float>(Metrics::Height) / 2.5f));
        ImGui::SetNextWindowPos(ImVec2(0, 500));
        adjustWindow = false;
    }

    ImGui::Begin("Sample Manager");

    if (ImGui::BeginCombo("Select a Sample", sample_manager.GetSampleName(sample_manager.GetCurrentIndex()).c_str()))
    {
        for (std::size_t index = 0; index < sample_manager.GetSampleNbr(); index++)
        {
            if (ImGui::Selectable(sample_manager.GetSampleName(index).c_str(), sample_manager.GetCurrentIndex() == index))
            {
                sample_manager.ChangeSample(index);
            }
        }
        ImGui::EndCombo();
    }

    ImGui::Spacing();

    ImGui::TextWrapped(sample_manager.GetSampleDescription(sample_manager.GetCurrentIndex()).c_str());

    ImGui::Spacing();

    sample_manager.DrawImgui(sample_manager.GetCurrentIndex());

    ImGui::SetCursorPosY(ImGui::GetWindowHeight() - (ImGui::GetFrameHeightWithSpacing()));

    if (ImGui::ArrowButton("PreviousSample", ImGuiDir_Left))
    {
        sample_manager.PreviousSample();
    }

    ImGui::SameLine();

    if (ImGui::Button("Regenerate"))
    {
        sample_manager.RegenerateSample();
    }

    ImGui::SameLine();

    if (ImGui::ArrowButton("NextSample", ImGuiDir_Right))
    {
        sample_manager.NextSample();
    }

    ImGui::End();
}

bool Raytracing::onKeyEvent(const KeyboardEvent& keyEvent)
{
    //if (keyEvent.key == Input::Key::Space && keyEvent.type == KeyboardEvent::Type::KeyPressed)
    //{
    //    mRayTrace = !mRayTrace;
    //    return true;
    //}

    //if (mpScene && mpScene->onKeyEvent(keyEvent))
    //    return true;

    //return false;

    return renderer_->onKeyEvent(keyEvent);
}

bool Raytracing::onMouseEvent(const MouseEvent& mouseEvent)
{
    return renderer_->onMouseEvent(mouseEvent);
   /* return mpScene && mpScene->onMouseEvent(mouseEvent);*/
}
//
//void Raytracing::loadScene(const std::filesystem::path& path, const Fbo* pTargetFbo)
//{
//    mpScene = Scene::create(getDevice(), path);
//    mpCamera = mpScene->getCamera();
//
//    // Update the controllers
//    float radius = mpScene->getSceneBounds().radius();
//    mpScene->setCameraSpeed(radius * 0.25f);
//    float nearZ = std::max(0.1f, radius / 750.0f);
//    float farZ = radius * 10;
//    mpCamera->setDepthRange(nearZ, farZ);
//    mpCamera->setAspectRatio((float)pTargetFbo->getWidth() / (float)pTargetFbo->getHeight());
//
//    // Get shader modules and type conformances for types used by the scene.
//    // These need to be set on the program in order to use Falcor's material system.
//    auto shaderModules = mpScene->getShaderModules();
//    auto typeConformances = mpScene->getTypeConformances();
//
//    // Get scene defines. These need to be set on any program using the scene.
//    auto defines = mpScene->getSceneDefines();
//
//    // Create raster pass.
//    // This utility wraps the creation of the program and vars, and sets the necessary scene defines.
//   /* ProgramDesc rasterProgDesc;
//    rasterProgDesc.addShaderModules(shaderModules);
//    rasterProgDesc.addShaderLibrary("Samples/HelloDXR/HelloDXR.3d.slang").vsEntry("vsMain").psEntry("psMain");
//    rasterProgDesc.addTypeConformances(typeConformances);*/
//
//    //mpRasterPass = RasterPass::create(getDevice(), rasterProgDesc, defines);
//
//    // We'll now create a raytracing program. To do that we need to setup two things:
//    // - A program description (ProgramDesc). This holds all shader entry points, compiler flags, macro defintions,
//    // etc.
//    // - A binding table (RtBindingTable). This maps shaders to geometries in the scene, and sets the ray generation and
//    // miss shaders.
//    //
//    // After setting up these, we can create the Program and associated RtProgramVars that holds the variable/resource
//    // bindings. The Program can be reused for different scenes, but RtProgramVars needs to binding table which is
//    // Scene-specific and needs to be re-created when switching scene. In this example, we re-create both the program
//    // and vars when a scene is loaded.
//
//    ProgramDesc rtProgDesc;
//    rtProgDesc.addShaderModules(shaderModules);
//    rtProgDesc.addShaderLibrary("Samples/Raytracing/Raytracing.rt.slang");
//    rtProgDesc.addTypeConformances(typeConformances);
//    rtProgDesc.setMaxTraceRecursionDepth(3); // 1 for calling TraceRay from RayGen, 1 for calling it from the
//                                             // primary-ray ClosestHit shader for reflections, 1 for reflection ray
//                                             // tracing a shadow ray
//    rtProgDesc.setMaxPayloadSize(24);        // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
//                                             // should be set as small as possible for maximum performance.
//
//    ref<RtBindingTable> sbt = RtBindingTable::create(2, 2, mpScene->getGeometryCount());
//    sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
//    sbt->setMiss(0, rtProgDesc.addMiss("primaryMiss"));
//    sbt->setMiss(1, rtProgDesc.addMiss("shadowMiss"));
//    auto primary = rtProgDesc.addHitGroup("primaryClosestHit", "primaryAnyHit");
//    auto shadow = rtProgDesc.addHitGroup("", "shadowAnyHit");
//    sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);
//    sbt->setHitGroup(1, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), shadow);
//
//    mpRaytraceProgram = Program::create(getDevice(), rtProgDesc, defines);
//    mpRtVars = RtProgramVars::create(getDevice(), mpRaytraceProgram, sbt);
//}

//void Raytracing::setPerFrameVars(const Fbo* pTargetFbo)
//{
//    auto var = mpRtVars->getRootVar();
//    var["PerFrameCB"]["invView"] = inverse(mpCamera->getViewMatrix());
//    var["PerFrameCB"]["viewportDims"] = float2(pTargetFbo->getWidth(), pTargetFbo->getHeight());
//    float fovY = focalLengthToFovY(mpCamera->getFocalLength(), Camera::kDefaultFrameHeight);
//    var["PerFrameCB"]["tanHalfFovY"] = std::tan(fovY * 0.5f);
//    var["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
//    var["PerFrameCB"]["useDOF"] = mUseDOF;
//    var["gOutput"] = mpRtOut;
//}

//void Raytracing::renderRaster(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
//{
//    FALCOR_ASSERT(mpScene);
//    FALCOR_PROFILE(pRenderContext, "renderRaster");
//
//    mpRasterPass->getState()->setFbo(pTargetFbo);
//    mpScene->rasterize(pRenderContext, mpRasterPass->getState().get(), mpRasterPass->getVars().get());
//}

//void Raytracing::renderRT(RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo)
//{
//    FALCOR_ASSERT(mpScene);
//    FALCOR_PROFILE(pRenderContext, "renderRT");
//
//    setPerFrameVars(pTargetFbo.get());
//
//    pRenderContext->clearUAV(mpRtOut->getUAV().get(), kClearColor);
//    mpScene->raytrace(pRenderContext, mpRaytraceProgram.get(), mpRtVars, uint3(pTargetFbo->getWidth(), pTargetFbo->getHeight(), 1));
//    pRenderContext->blit(mpRtOut->getSRV(), pTargetFbo->getRenderTargetView(0));
//}

int runMain(int argc, char** argv)
{
    SampleAppConfig config;
    config.windowDesc.title = "Raytracing";
    config.windowDesc.resizableWindow = true;

    Raytracing helloDXR(config);
    return helloDXR.run();
}

int main(int argc, char** argv)
{
    return catchAndReportAllExceptions([&]() { return runMain(argc, argv); });
}
