#include "Renderer.h"

#include "Core/Program/Program.h"
#include "Scene/Material/StandardMaterial.h"
#include "Utils/Math/FalcorMath.h"

Renderer::Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept
    : device_(device), target_fbo_(target_fbo)
{
}

void Renderer::Init() noexcept
{
    if (device_->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    {
        FALCOR_THROW("Device does not support raytracing!");
    }

    //    static constexpr std::array<Vertex, 3> vertices
    //    {
    //                Vertex{float3(0.0f, 0.5f, 0.0f)},  // Top vertex (red)
    //                Vertex{float3(0.5f, -0.5f, 0.0f)}, // Bottom right (green)
    //                Vertex{float3(-0.5f, -0.5f, 0.0f)} // Bottom left (blue)
    //    };
    //
    //    vertex_buffer_ = device_->createBuffer(
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
    //    raster_pass_ = RasterPass::create(device_,
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
    // BlendState::Desc blendDesc;
    // blendDesc.setRtBlend(0, true).setRtParams(
    //    0,
    //    BlendState::BlendOp::Add,
    //    BlendState::BlendOp::Add,
    //    BlendState::BlendFunc::SrcAlpha,
    //    BlendState::BlendFunc::OneMinusSrcAlpha,
    //    BlendState::BlendFunc::One,
    //    BlendState::BlendFunc::One
    //);
    // pState->setBlendState(BlendState::create(blendDesc));

    Settings settings{};

    // Create the SceneBuilder
    SceneBuilder::Flags flags = SceneBuilder::Flags::RTDontMergeStatic | SceneBuilder::Flags::RTDontMergeDynamic |
                                SceneBuilder::Flags::RTDontMergeInstanced | SceneBuilder::Flags::DontOptimizeGraph;
    scene_builder_ = new SceneBuilder(device_, settings, flags);

    auto sphere = TriangleMesh::createSphere(3.f);
    auto cube = TriangleMesh::createCube(float3(1.f, 1.f, 1.f));

    // Create a lambertian material
    ref<Material> lambertian = StandardMaterial::create(device_, "Lambertian");
    lambertian->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
    lambertian->setRoughnessMollification(1.f);
    lambertian->setIndexOfRefraction(0.f);

    ref<Material> water_particle_mat = StandardMaterial::create(device_, "water particle");
    water_particle_mat->toBasicMaterial()->setBaseColor3(float3(0.1f, 0.2f, 1.0f));
    water_particle_mat->setRoughnessMollification(1.f);
    water_particle_mat->setIndexOfRefraction(0.f);

    // lambertian->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));

    // ref<Material> dielectric_red = StandardMaterial::create(device_, "DielecRed");
    // dielectric_red->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));
    // dielectric_red->setDoubleSided(true);
    // dielectric_red->setIndexOfRefraction(1.f);
    ////dielectric_red->toBasicMaterial()->setDiffuseTransmission(1.f);

    ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    dielectric_blue->setDoubleSided(true);
    dielectric_blue->setIndexOfRefraction(1.f);
    dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    // auto triangle_mesh_id_1 = scene_builder__.addTriangleMesh(sphere, dielectric_red);
    // auto triangle_mesh_id_2 = scene_builder__.addTriangleMesh(cube, dielectric_blue);
    // auto triangle_mesh_id_3 = scene_builder__.addTriangleMesh(cube, dielectric_blue);
    auto triangle_mesh_id_4 = scene_builder_->addTriangleMesh(cube, lambertian);
    auto sphere_mesh_id = scene_builder_->addTriangleMesh(sphere, dielectric_blue);

    // AABB raymarching_AABB = AABB(float3(-50, -5, -5) + float3(0, 10, 0),
    //     float3(50, 5, 5) + float3(0, 10, 0));
    // uint32_t raymarching_AABB_ID = 1;
    // scene_builder__.addCustomPrimitive(raymarching_AABB_ID, raymarching_AABB);

    // auto raymarching_node = SceneBuilder::Node();
    // raymarching_node.name = "RaymarchingNode";
    // auto raymarching_transform = Transform();
    // raymarching_transform.setTranslation(float3(0.f, 10.f, 0.f));
    // raymarching_transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    // raymarching_transform.setScaling(float3(1.f, 1.f, 1.f));
    // raymarching_node.transform = raymarching_transform.getMatrix();
    // auto raymarching_node_id = scene_builder__.addNode(raymarching_node);

    // auto node = SceneBuilder::Node();
    // node.name = "Sphere1";
    // auto transform = Transform();
    // transform.setTranslation(float3(0.f, 0.f, 0.f));
    // transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    // transform.setScaling(float3(1.f, 1.f, 1.f));
    // node.transform = transform.getMatrix();
    // auto node_id = scene_builder_->addNode(node);

    //// Add Mesh Instances
    // scene_builder_->addMeshInstance(node_id, triangle_mesh_id_1);

    // auto node_2 = SceneBuilder::Node();
    // node_2.name = "Sphere1";
    // auto transform_2 = Transform();
    // transform_2.setTranslation(float3(10.f, 0.f, 0.f));
    // transform_2.setRotationEuler(float3(0.f, 0.f, 0.f));
    // transform_2.setScaling(float3(1.f));
    // node_2.transform = transform_2.getMatrix();
    // auto node_id_2 = scene_builder__.addNode(node_2);

    //// Add Mesh Instances
    // scene_builder__.addMeshInstance(node_id_2, triangle_mesh_id_2);

    // auto node_3 = SceneBuilder::Node();
    // node_3.name = "cube";
    // auto transform_3 = Transform();
    // transform_3.setTranslation(float3(10.f, 0.f, 10.f));
    // transform_3.setRotationEuler(float3(0.f, 0.f, 0.f));
    // transform_3.setScaling(float3(10.f, 1.f, 10.f));
    // node_3.transform = transform_3.getMatrix();
    // auto node_id_3 = scene_builder__.addNode(node_3);

    //// Add Mesh Instances
    // scene_builder__.addMeshInstance(node_id_3, triangle_mesh_id_3);

    // auto node_4 = SceneBuilder::Node();
    // node_4.name = "cube";
    // auto transform_4 = Transform();
    // transform_4.setTranslation(float3(0.f, -1.2f, 0.f));
    // transform_4.setRotationEuler(float3(0.f, 0.f, 0.f));
    // transform_4.setScaling(float3(1000.f, 1.f, 1000.f));
    // node_4.transform = transform_4.getMatrix();
    // auto node_id_4 = scene_builder_->addNode(node_4);

    //// Add Mesh Instances
    // scene_builder_->addMeshInstance(node_id_4, triangle_mesh_id_4);

    auto envMap = EnvMap::createFromFile(device_, "hallstatt4_hd.hdr");
    envMap->setIntensity(1.0);
    scene_builder_->setEnvMap(envMap);

    ref<Camera> camera = ref<Camera>(new Camera("Camera"));
    camera->setPosition(float3(0, 0.0, -250));
    camera->setTarget(float3(0, 0.0, 0));
    camera->setUpVector(float3(0, 1, 0));
    camera->setFocalLength(35);
    camera->setDepthRange(0.1f, 10000.f);

    scene_builder_->addCamera(camera);

    mpScene = scene_builder_->getScene();

    mpCamera = mpScene->getCamera();

    // Update the controllers
    float radius = mpScene->getSceneBounds().radius();
    mpScene->setCameraSpeed(50.f);
    auto pTargetFbo = target_fbo_.get();
    mpCamera->setAspectRatio(static_cast<float>(pTargetFbo->getWidth()) / static_cast<float>(pTargetFbo->getHeight()));

    // Get shader modules and type conformances for types used by the scene.
    // These need to be set on the program in order to use Falcor's material system.
    auto shaderModules = mpScene->getShaderModules();
    auto typeConformances = mpScene->getTypeConformances();

    // Get scene defines. These need to be set on any program using the scene.
    auto defines = mpScene->getSceneDefines();

    // Create a triangle using a Scene object
    // Create a raytracing program description
    ProgramDesc rtProgDesc;
    ProgramDesc::ShaderModule s_module = ProgramDesc::ShaderModule("Samples/Raytracing/SDF_Functions.slang");
    s_module.addFile("Samples/Raytracing/SDF_Functions.slang");
    shaderModules.emplace_back(s_module);
    rtProgDesc.addShaderModules(shaderModules);
    rtProgDesc.addShaderLibrary("Samples/Raytracing/Raytracing.rt.slang");
    rtProgDesc.addTypeConformances(typeConformances);
    rtProgDesc.setMaxTraceRecursionDepth(kMaxRayBounce + 1);
    rtProgDesc.setMaxPayloadSize(48); // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
                                      // should be set as small as possible for maximum performance.

    ref<RtBindingTable> sbt = RtBindingTable::create(1, 1, mpScene->getGeometryCount());
    sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
    sbt->setMiss(0, rtProgDesc.addMiss("miss"));

    auto primary = rtProgDesc.addHitGroup("closestHit", "");
    sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);

    auto raymarching_hit_group = rtProgDesc.addHitGroup("RaymarchingClosestHit", "", "RaymarchingIntersection");
    sbt->setHitGroup(0, mpScene->getGeometryIDs(Scene::GeometryType::Custom), raymarching_hit_group);

    rt_program_ = Program::create(device_, rtProgDesc, defines);
    rt_program_vars_ = RtProgramVars::create(device_, rt_program_, sbt);

    //// Create the raytracing program
    // auto rtProgram = Program::create(device_, rtProgDesc);

    //// Create shader variables (for passing data to shaders)
    // auto rtVars = RtProgramVars::create(device_, rtProgram);
    /*   if (device_->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
       {
           FALCOR_THROW("Device does not support raytracing!");
       }

       loadScene(kDefaultScene, getTargetFbo().get());*/
}

void Renderer::RenderFrame(RenderContext* pRenderContext, const double& currentTime) noexcept
{
    pRenderContext->clearFbo(target_fbo_.get(), float4(kClearColor, 1),
        1.0f, 0, FboAttachmentType::All);

    // raster_pass_->getState()->setVao(vao_);
    // raster_pass_->getState()->setFbo(pTargetFbo);
    // raster_pass_->draw(pRenderContext, 3, 0);

    IScene::UpdateFlags updates = mpScene->update(pRenderContext, currentTime);
    if (is_set(updates, IScene::UpdateFlags::GeometryChanged))
        FALCOR_THROW("This sample does not support scene geometry changes.");
    if (is_set(updates, IScene::UpdateFlags::RecompileNeeded))
        FALCOR_THROW("This sample does not support scene changes that require shader recompilation.");

    FALCOR_ASSERT(mpScene)
    // FALCOR_PROFILE(pRenderContext, "renderRT");

    auto var = rt_program_vars_->getRootVar();
    var["PerFrameCB"]["invView"] = inverse(mpCamera->getViewMatrix());
    var["PerFrameCB"]["viewportDims"] = float2(target_fbo_->getWidth(), target_fbo_->getHeight());
    float fovY = focalLengthToFovY(mpCamera->getFocalLength(), Camera::kDefaultFrameHeight);
    var["PerFrameCB"]["tanHalfFovY"] = std::tan(fovY * 0.5f);
    var["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
    var["PerFrameCB"]["useDOF"] = mUseDOF;

    var["PerFrameCB"]["backgroundColor"] = kClearColor;

    var["PerFrameCB"]["waterTurbulence"] = waterTurbulence;

    var["PerFrameCB"]["maxRayBounce"] = kMaxRayBounce;

    var["PerFrameCB"]["absorptionCoeff"] = absorptionCoeff;
    var["PerFrameCB"]["scatteringCoeff"] = scatteringCoeff;
    var["PerFrameCB"]["phaseG"] = phaseG;

    var["PerFrameCB"]["maxRaymarchingDistance"] = maxRayMarchingDistance;
    var["PerFrameCB"]["marchSize"] = kMarchSize;

    var["PerFrameCB"]["maxLighMarchingDistance"] = maxLighMarchingDistance;
    var["PerFrameCB"]["sunLightMarchSize"] = sunLightMarchSize;

    var["PerFrameCB"]["lightColor"] = lightColor;
    var["PerFrameCB"]["lightDir"] = lightDir;

    var["PerFrameCB"]["IoR"] = IoR;

    var["PerFrameCB"]["time"] = static_cast<float>(currentTime);
    static int frame = 0;
    var["PerFrameCB"]["iFrame"] = frame++;

    var["gOutput"] = mpRtOut;
    var["gTexture3D"] = mpTexture3D;

    pRenderContext->clearUAV(mpRtOut->getUAV().get(), float4(kClearColor, 1));
    mpScene->raytrace(pRenderContext, rt_program_.get(), rt_program_vars_,
        uint3(target_fbo_->getWidth(), target_fbo_->getHeight(), 1));
    pRenderContext->blit(mpRtOut->getSRV(), target_fbo_->getRenderTargetView(0));
}

void Renderer::RenderUI(Gui* pGui, Gui::Window* app_gui_window) noexcept
{
    app_gui_window->rgbColor("Background color", kClearColor);

    app_gui_window->var("Water Turbulance", waterTurbulence);

    app_gui_window->var("MaxRayBounce", kMaxRayBounce);

    app_gui_window->var("absorptionCoeff", absorptionCoeff);
    app_gui_window->var("scatteringCoeff", scatteringCoeff);
    app_gui_window->var("Phase G ", phaseG);

    app_gui_window->var("maxRaymarchingDistance", maxRayMarchingDistance);
    app_gui_window->var("MarchSize", kMarchSize);
    app_gui_window->var("maxLighMarchingDistance", maxLighMarchingDistance);
    app_gui_window->var("sunLightMarchSize", sunLightMarchSize);

    app_gui_window->rgbColor("Light color", lightColor);
    static float3 ImGUI_LightDir = lightDir;
    app_gui_window->var("Light Direction", ImGUI_LightDir);
    lightDir = math::normalize(ImGUI_LightDir);

    app_gui_window->var("IoR", IoR);

    app_gui_window->checkbox("Use Depth of Field", mUseDOF);

    mpScene->renderUI(*app_gui_window);
}

void Renderer::OnResize(uint32_t width, uint32_t height) noexcept
{
    const float h = static_cast<float>(height);
    const float w = static_cast<float>(width);

    if (mpCamera)
    {
        mpCamera->setFocalLength(18);
        const float aspectRatio = (w / h);
        mpCamera->setAspectRatio(aspectRatio);
    }

    mpRtOut = device_->createTexture2D(width, height, ResourceFormat::RGBA16Float,
        1, 1, nullptr,
        ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    );
}

bool Renderer::onKeyEvent(const KeyboardEvent& keyEvent) const noexcept
{
    if (mpScene && mpScene->onKeyEvent(keyEvent))
        return true;

    return false;
}
bool Renderer::onMouseEvent(const MouseEvent& mouseEvent) const noexcept
{
    return mpScene && mpScene->onMouseEvent(mouseEvent);
}

void Renderer::Deinit() noexcept {}


