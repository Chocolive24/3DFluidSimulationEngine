#include "Renderer.h"
#include "Metrics.h"

#include "Core/Program/Program.h"
#include "Scene/Material/StandardMaterial.h"
#include "Utils/Threading.h"
#include "Utils/Math/FalcorMath.h"
#include "Utils/Timing/Profiler.h"
#include "../../Physics/include/SPH.h"
#include "Scene/SDFs/SparseVoxelSet/SDFSVS.h"

Renderer::Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept
    : device_(device), target_fbo_(target_fbo),
render_graph_(device, "FluidRenderGraph")
{}

void Renderer::Init(RenderContext* render_context) noexcept
{
    if (device_->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    {
        FALCOR_THROW("Device does not support raytracing!");
    }

    //createRasterizationProgram();

    Settings settings{};

    // Create the SceneBuilder
    SceneBuilder::Flags flags = SceneBuilder::Flags::RTDontMergeStatic | SceneBuilder::Flags::RTDontMergeDynamic |
                                SceneBuilder::Flags::RTDontMergeInstanced | SceneBuilder::Flags::DontOptimizeGraph;
    scene_builder_ = new SceneBuilder(device_, settings, flags);

    auto sphere_mesh = TriangleMesh::createSphere(Metrics::PARTICLESIZE);
    //auto cube_mesh = TriangleMesh::createCube(float3(Metrics::sim_bounds - 1.f));

    ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    dielectric_blue->setDoubleSided(true);
    dielectric_blue->setIndexOfRefraction(1.f);
    dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, dielectric_blue, true);
    //cube_mesh_id = scene_builder_->addTriangleMesh(cube_mesh, dielectric_blue);

    // Create a lambertian material
    ref<Material> lambertian = StandardMaterial::create(device_, "Lambertian");
    lambertian->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
    lambertian->setRoughnessMollification(1.f);
    lambertian->setIndexOfRefraction(0.f);

    //sphere = TriangleMesh::createQuad(float2(5.f));
    //sphere_mesh_id = scene_builder_->addTriangleMesh(sphere, dielectric_blue, true);

    //auto node = SceneBuilder::Node();
    //const std::string name = "Sphere " /* + std::to_string(i)*/;
    //node.name = name;
    //auto transform = Transform();
    //transform.setTranslation(float3(0.f, 0.f, 0.f));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(1.f, 1.f, 1.f));
    //node.transform = transform.getMatrix();
    //sphere_node_id_ = scene_builder_->addNode(node);

    //// Add Mesh Instances
    //scene_builder_->addMeshInstance(sphere_node_id_, sphere_mesh_id);

    AABB fluid_AABB = AABB(float3(-Metrics::WALLDIST), float3(Metrics::WALLDIST));
    uint32_t fluid_AABB_ID = 1;
    scene_builder_->addCustomPrimitive(fluid_AABB_ID, fluid_AABB);

    auto fluid_node = SceneBuilder::Node();
    fluid_node.name = "RaymarchingNode";
    fluid_transform = Transform();
    fluid_transform.setTranslation(float3(0.f, 0, 0.f));
    fluid_transform.setRotationEuler(float3(0.f));
    //fluid_transform.setRotationEuler(float3(3.14 / 4, 3.14 / 4, 3.14 / 4));
    fluid_transform.setScaling(float3(1.f, 1.f, 1.f));
    fluid_node.transform = fluid_transform.getMatrix();
    auto raymarching_node_id = scene_builder_->addNode(fluid_node);

    //auto sdf = SDFSVS::create(device_);
    //sdf->generateCheeseValues(64, 0);
    //scene_builder_->addSDFGrid(sdf, lambertian);

    //auto node = SceneBuilder::Node();
    //node.name = "Cube Density Map Size";
    //auto transform = Transform();
    //transform.setTranslation(float3(0.f, 0.f, 0));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(1, 1.f, 1.f));
    //node.transform = transform.getMatrix();
    //auto node_id = scene_builder_->addNode(node);

    auto envMap = EnvMap::createFromFile(device_, "data/images/hallstatt4_hd.hdr");
    envMap->setIntensity(1.0);
    scene_builder_->setEnvMap(envMap);

    ref<Camera> camera = ref<Camera>(new Camera("Camera"));
    camera->setPosition(float3(0, 0.0, -250));
    camera->setTarget(float3(0, 0.0, 0));
    camera->setUpVector(float3(0, 1, 0));
    camera->setFocalLength(35);
    camera->setDepthRange(0.1f, 10000.f);

    scene_builder_->addCamera(camera);

    compute_density_map_pass_ =
       ComputePass::create(device_,
           "Samples/3DFluidSimulationEngine/Renderer/shaders/SPH.cs.slang",
           "computeDensityMap");

    density_3d_tex_ = device_->createTexture3D(
       Metrics::density_map_size,
       Metrics::density_map_size,
       Metrics::density_map_size,
        ResourceFormat::R32Float,
        1, // mips
        nullptr,
        ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    );

    Sampler::Desc sampler_desc{};
    sampler_desc.setFilterMode(
        TextureFilteringMode::Linear,
        TextureFilteringMode::Linear,
        TextureFilteringMode::Linear
    );
    sampler_desc.setAddressingMode(
        TextureAddressingMode::Clamp,
        TextureAddressingMode::Clamp,
        TextureAddressingMode::Clamp
    );
    linearClampSampler_ = make_ref<Sampler>(device_, sampler_desc);


    //    std::vector<float3> positions;
    //std::vector<float3> normals;
    //std::vector<float3> tangents;
    //std::vector<float2> uv;
    //for (const auto& vertex : sphere->getVertices())
    //{
    //    // std::cout << vertex.position.x << " " << vertex.position.y << " " << vertex.position.z << "\n";
    //    const auto new_Pos = vertex.position + float3(0, 1, 0);
    //    positions.push_back(new_Pos);
    //    normals.push_back(vertex.normal);
    //    tangents.push_back(vertex.normal);
    //    uv.push_back(float2(vertex.texCoord.x, vertex.texCoord.y));
    //}

    //b_pos = make_ref<Buffer>(
    //    device_,
    //    sizeof(float) * 3,
    //    positions.size(),
    //    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
    //    MemoryType::DeviceLocal,
    //    positions.data(),
    //    false
    //);

    //b_pos_readback =
    //    make_ref<Buffer>(device_,
    //        sizeof(positions[0]),
    //        positions.size(),
    //        ResourceBindFlags::None,
    //        MemoryType::Upload,
    //        nullptr,
    //        false);

    //b_normal = make_ref<Buffer>(
    //    device_,
    //    sizeof(normals[0]),
    //    normals.size(),
    //    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
    //    MemoryType::Upload,
    //    normals.data(),
    //    false
    //);

    //b_tang = make_ref<Buffer>(
    //    device_,
    //    sizeof(tangents[0]),
    //    tangents.size(),
    //    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
    //    MemoryType::Upload,
    //    tangents.data(),
    //    false
    //);

    //b_uv = make_ref<Buffer>(
    //    device_,
    //    sizeof(uv[0]),
    //    uv.size(),
    //    ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
    //    MemoryType::Upload,
    //    uv.data(),
    //    false
    //);

    //render_context->copyResource(b_pos_readback.get(), b_pos.get());

    //const float3* poses = static_cast<const float3*>(b_pos_readback->map());

    //  for (int i = 0; i < 4; i++)
    //  {
    //      const auto p = poses[i];

    //      std::cout << p.x << " " << p.y << " " << p.z << '\n';
    //  }

    //b_pos_readback->unmap();

    //auto r = RtAccelerationStructure::create(device_, {});
    //RtAccelerationStructure::BuildDesc b;
   /* b.
    render_context->buildAccelerationStructure()*/

  /*  scene_->getMeshVerticesAndIndices()

    render_context->in*/
}

void Renderer::RenderFrame(RenderContext* pRenderContext, const double& currentTime,
    const ref<Buffer>& bodies,
    const ref<Buffer>& SpatialIndices,
    const ref<Buffer>& SpatialOffsets) noexcept
{
    pRenderContext->clearFbo(target_fbo_.get(),
        float4(bg_clear_color, 1), 1.0f, 0,
        FboAttachmentType::All);

    //auto transform = Transform();
    //transform.setTranslation(float3(10, 0, 0));
    //scene_->updateNodeTransform(sphere_node_id_.get(), transform.getMatrix());

    //static bool hasUpdated = false;
    //if (true)
    //{
    //    //hasUpdated = true;
    //    //std::cout << "\nUPDATES VERTICES\n";

    //    const std::map < std::string, ref<Buffer>> vertices{
    //        {"positions", b_pos},
    //        {"normals", b_normal},
    //        {"tangents", b_tang},
    //        {"texcrds", b_uv},
    //    };

    //    scene_->setMeshVertices(sphere_mesh_id, vertices);
    //}

     const auto compute_var = compute_density_map_pass_->getRootVar();
     compute_var["bodies"] = bodies;
     compute_var["SpatialIndices"] = SpatialIndices;
     compute_var["SpatialOffsets"] = SpatialOffsets;
     compute_var["gTexture3D"] = density_3d_tex_;

     compute_var["PerFrameCB"]["densityMapSize"] = Metrics::density_map_size;
     compute_var["PerFrameCB"]["simBounds"] = float3(Metrics::sim_bounds);
     compute_var["PerFrameCB"]["wallDist"] = Metrics::WALLDIST;
     compute_var["PerFrameCB"]["fixedDeltaTime"] = Metrics::kFixedDeltaTime;
     compute_var["PerFrameCB"]["nbParticles"] = Metrics::NbParticles;
     //compute_var["PerFrameCB"]["gravity"] = world_->Gravity;
     compute_var["PerFrameCB"]["smoothingRadius"] = SPH::SmoothingRadius;
     compute_var["PerFrameCB"]["targetDensity"] = SPH::TargetDensity;
     compute_var["PerFrameCB"]["pressureMultiplier"] = SPH::PressureMultiplier;
     compute_var["PerFrameCB"]["viscosityStrength"] = SPH::ViscosityStrength;
     compute_var["PerFrameCB"]["densityGraphicsMultiplier"] = Metrics::densityGraphicsMultiplier;

     const float4x4 localToWorld = fluid_transform.getMatrix();
     const float4x4 worldToLocal = inverse(localToWorld);

     compute_var["PerFrameCB"]["localToWorld"] = localToWorld;
     compute_var["PerFrameCB"]["worldToLocal"] = worldToLocal;

     const float r = SPH::SmoothingRadius;
     const float spikyPow2 = 15.f / (2 * PI * pow(r, 5));
     const float spikyPow3 = 15.f / (PI * pow(r, 6));
     const float spikyPow2Grad = 15.f / (PI * pow(r, 5));
     const float spikyPow3Grad = 45.f / (PI * pow(r, 6));

     compute_var["PerFrameCB"]["K_SpikyPow2"] = spikyPow2;
     compute_var["PerFrameCB"]["K_SpikyPow3"] = spikyPow3;
     compute_var["PerFrameCB"]["K_SpikyPow2Grad"] = spikyPow2Grad;
     compute_var["PerFrameCB"]["K_SpikyPow3Grad"] = spikyPow3Grad;

    // const uint32_t thread_groups = (density_map_size + 7) / 8;
     compute_density_map_pass_->execute(pRenderContext,
         Metrics::density_map_size, Metrics::density_map_size, Metrics::density_map_size);

    //std::cout << "Before scene update\n";
    IScene::UpdateFlags updates = scene_->update(pRenderContext, currentTime);

    //if (is_set(updates, IScene::UpdateFlags::GeometryChanged))
    //{
    //    std::cout << "GeometryChanged\n";
    //}
    //if (is_set(updates, IScene::UpdateFlags::GeometryMoved))
    //{
    //    std::cout << "GeometryMoved\n";
    //}
    //if (is_set(updates, IScene::UpdateFlags::MeshesChanged))
    //{
    //    std::cout << "MeshesChanged\n";
    //}
    //std::cout << "After scene update\n";

    if (is_set(updates, IScene::UpdateFlags::GeometryChanged))
        FALCOR_THROW("This sample does not support scene geometry changes.");
    if (is_set(updates, IScene::UpdateFlags::RecompileNeeded))
        FALCOR_THROW("This sample does not support scene changes that require shader recompilation.");

 /*   FALCOR_ASSERT(scene_);
    FALCOR_PROFILE(pRenderContext, "renderRT");*/

    setPerFrameVariables(currentTime);

    pRenderContext->clearUAV(rt_output_tex_->getUAV().get(), float4(bg_clear_color, 1));
    scene_->raytrace(pRenderContext, rt_program_.get(), rt_program_vars_, uint3(target_fbo_->getWidth(), target_fbo_->getHeight(), 1));
    pRenderContext->blit(rt_output_tex_->getSRV(), target_fbo_->getRenderTargetView(0));
}

void Renderer::RenderUI(Gui* pGui, Gui::Window* app_gui_window) noexcept
{
    app_gui_window->rgbColor("Background color", bg_clear_color);

    app_gui_window->slider("densityGraphicsMultiplier",
        Metrics::densityGraphicsMultiplier, 0.f, 200.f);
    app_gui_window->slider("volumeValueOffset", volumeValueOffset, 0.f, 1.f);
    app_gui_window->slider("DensityDepth", DensityRayMarchMultiplier, 0.f, 1.f);
    app_gui_window->slider("SphereRadius", SphereRadius, 0.f, 200.f);

    app_gui_window->var("ISO Level", IsoLevel);
    app_gui_window->var("normalOffset", normalOffset);

    app_gui_window->checkbox("Draw Fluid ?", draw_fluid_);
    app_gui_window->checkbox("Light Scattering ?", lightScattering);
    //if (draw_fluid_)
    //{
        app_gui_window->var("Water Turbulance", water_turbulence_);

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

        // app_gui_window->checkbox("Use Depth of Field", mUseDOF);
    //}

    scene_->renderUI(*app_gui_window);
}

void Renderer::OnResize(uint32_t width, uint32_t height) noexcept
{
    const float h = static_cast<float>(height);
    const float w = static_cast<float>(width);

    if (camera_)
    {
        camera_->setFocalLength(18);
        const float aspectRatio = (w / h);
        camera_->setAspectRatio(aspectRatio);
    }

    rt_output_tex_ = device_->createTexture2D(
        width, height, ResourceFormat::RGBA16Float, 1, 1, nullptr, ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    );
}

bool Renderer::onKeyEvent(const KeyboardEvent& keyEvent) const noexcept
{
    if (scene_ && scene_->onKeyEvent(keyEvent))
        return true;

    return false;
}
bool Renderer::onMouseEvent(const MouseEvent& mouseEvent) const noexcept
{
    return scene_ && scene_->onMouseEvent(mouseEvent);
}

void Renderer::Deinit() noexcept {}

void Renderer::CreateRaytracingProgram() noexcept
{
    scene_ = scene_builder_->getScene();
    camera_ = scene_->getCamera();

    // Update the controllers
    scene_->setCameraSpeed(200.f);
    const auto pTargetFbo = target_fbo_.get();
    camera_->setAspectRatio(static_cast<float>(pTargetFbo->getWidth()) / static_cast<float>(pTargetFbo->getHeight()));

    // Get shader modules and type conformances for types used by the scene.
    // These need to be set on the program in order to use Falcor's material system.
    auto shaderModules = scene_->getShaderModules();
    const auto typeConformances = scene_->getTypeConformances();

    // Get scene defines. These need to be set on any program using the scene.
    const auto defines = scene_->getSceneDefines();

    // Create a raytracing program description
    ProgramDesc rtProgDesc;
    ProgramDesc::ShaderModule s_module = ProgramDesc::ShaderModule("Samples/3DFluidSimulationEngine/Renderer/shaders/SDF_Functions.slang");
    s_module.addFile("Samples/3DFluidSimulationEngine/Renderer/shaders/SDF_Functions.slang");
    shaderModules.emplace_back(s_module);
    rtProgDesc.addShaderModules(shaderModules);
    rtProgDesc.addShaderLibrary("Samples/3DFluidSimulationEngine/Renderer/shaders/Raytracing.rt.slang");
    rtProgDesc.addTypeConformances(typeConformances);
    rtProgDesc.setMaxTraceRecursionDepth(kMaxRayBounce + 1);
    rtProgDesc.setMaxPayloadSize(48); // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
                                      // should be set as small as possible for maximum performance.

    const ref<RtBindingTable> sbt = RtBindingTable::create(1, 1, scene_->getGeometryCount());
    sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
    sbt->setMiss(0, rtProgDesc.addMiss("miss"));

    const auto primary = rtProgDesc.addHitGroup("closestHit", "");
    sbt->setHitGroup(0, scene_->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);

    const auto raymarching_hit_group = rtProgDesc.addHitGroup("RaymarchingClosestHit", "", "RaymarchingIntersection");
    sbt->setHitGroup(0, scene_->getGeometryIDs(Scene::GeometryType::Custom), raymarching_hit_group);

    rt_program_ = Program::create(device_, rtProgDesc, defines);
    rt_program_vars_ = RtProgramVars::create(device_, rt_program_, sbt);
}

void Renderer::LaunchMarchingCubeComputePasses(RenderContext* render_context) noexcept
{
    density_3d_tex_ = device_->createTexture3D(
        Metrics::density_map_size,
        Metrics::density_map_size,
        Metrics::density_map_size,
        ResourceFormat::R32Float,
        1, // mips
        nullptr,
        ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    );

    // 1) Buffer allocation
    const auto maxTriangleCount = 5 * (numPointsPerAxis - 1) * (numPointsPerAxis - 1) * (numPointsPerAxis - 1);
    const size_t triangleStructSize = sizeof(MarchingCubesTriangle);

    // The structured buffer that your HLSL AppendStructuredBuffer<Triangle> will write into:
    marching_cubes_triangle_buffer_ = make_ref<Buffer>(
        device_,                               // Falcor device
        triangleStructSize,                    // structSize (bytes per element)
        maxTriangleCount,                      // elementCount
        ResourceBindFlags::UnorderedAccess |   // UAV for compute
            ResourceBindFlags::ShaderResource, // SRV for rendering or readback
        MemoryType::DeviceLocal,               // GPU-only (faster)
        nullptr,                               // no init data
        true                                   // create hidden counter
    );

    marching_cubes_pass_ =
        ComputePass::create(device_, "Samples/3DFluidSimulationEngine/Renderer/shaders/MarchingCubes.cs.slang", "ProcessCube");

    render_context->clearUAVCounter(marching_cubes_triangle_buffer_, 0);

    const auto compute_var = marching_cubes_pass_->getRootVar();
    compute_var["DensityTexture"] = density_3d_tex_;
    compute_var["triangles"] = marching_cubes_triangle_buffer_;

    compute_var["PerFrameCB"]["numPointsPerAxis"] = numPointsPerAxis;
    compute_var["PerFrameCB"]["isoLevel"] = IsoLevel;
    compute_var["PerFrameCB"]["textureSize"] = Metrics::density_map_size;
    compute_var["PerFrameCB"]["boundSize"] = Metrics::sim_bounds;
    compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;

    marching_cubes_pass_->execute(render_context, 64, 64, 64);

    render_context->uavBarrier(marching_cubes_triangle_buffer_.get());
    marching_cubes_triangle_count_ = marching_cubes_triangle_buffer_->getUAVCounter()->getElement<uint>(0);

    // A small readback buffer to fetch the append counter:
    read_back_triangle_buffer_ = make_ref<Buffer>(
        device_, triangleStructSize, marching_cubes_triangle_count_, ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false
    );

    render_context->copyBufferRegion(
        read_back_triangle_buffer_.get(), 0, marching_cubes_triangle_buffer_.get(), 0, marching_cubes_triangle_count_ * triangleStructSize
    );

    const MarchingCubesTriangle* triangles = static_cast<const MarchingCubesTriangle*>(read_back_triangle_buffer_->map());

    for (uint i = 0; i < 12; i++)
    {
        std::cout << "Init " << '\n';
        std::cout << "UAV counter aka triangle count: " << marching_cubes_triangle_count_ << '\n';

        std::cout << triangles[i].vertexA.position.x << " " << triangles[i].vertexA.position.y << " " << triangles[i].vertexA.position.z
                  << '\n';
    }

    read_back_triangle_buffer_->unmap();
}

NodeID Renderer::AddSphereToScene(const float3 pos, const float radius) noexcept
{
    auto node = SceneBuilder::Node();
    const std::string name = "Sphere " /* + std::to_string(i)*/;
    node.name = name;
    auto transform = Transform();
    transform.setTranslation(pos);
    transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    transform.setScaling(float3(radius, radius, radius));
    node.transform = transform.getMatrix();
    const auto node_id = scene_builder_->addNode(node);

    // Add Mesh Instances
    scene_builder_->addMeshInstance(node_id, sphere_mesh_id);

    return node_id;
}

NodeID Renderer::AddCubeToScene(const float3 pos) noexcept
{
    // ref<TriangleMesh> sphere_mesh = TriangleMesh::createSphere(radius);

    // ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    // dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    // dielectric_blue->setDoubleSided(true);
    // dielectric_blue->setIndexOfRefraction(1.f);
    // dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    // MeshID sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, dielectric_blue);

    auto node = SceneBuilder::Node();
    const std::string name = "Cube " /* + std::to_string(i)*/;
    node.name = name;
    auto transform = Transform();
    transform.setTranslation(pos);
    transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    transform.setScaling(float3(1, 1, 1));
    node.transform = transform.getMatrix();
    const auto node_id = scene_builder_->addNode(node);

    // sphereNodeIDs.push_back(node_id);

    // Add Mesh Instances
    scene_builder_->addMeshInstance(node_id, cube_mesh_id);

    return node_id;
}

void Renderer::UpdateSceneNodeTransform(const NodeID nodeID, const Transform& transform) const noexcept
{
    scene_->updateNodeTransform(nodeID.get(), transform.getMatrix());
}

// void Renderer::CreateDensityMap(const std::vector<float>& particle_densities) noexcept
//{
//
// }

// void Renderer::ResetDensityMap(const std::vector<float>& particle_densities) noexcept
//{
//
//     density_3d_tex_.reset(particle_densities);
// }

void Renderer::setPerFrameVariables(const double& currentTime) const noexcept
{
    const auto var = rt_program_vars_->getRootVar();

    var["PerFrameCB"]["invView"] = inverse(camera_->getViewMatrix());
    var["PerFrameCB"]["viewportDims"] = float2(target_fbo_->getWidth(), target_fbo_->getHeight());
    const float fovY = focalLengthToFovY(camera_->getFocalLength(), Camera::kDefaultFrameHeight);
    var["PerFrameCB"]["tanHalfFovY"] = std::tan(fovY * 0.5f);

    /*var["PerFrameCB"]["sampleIndex"] = mSampleIndex++;
    var["PerFrameCB"]["useDOF"] = mUseDOF;*/

    var["PerFrameCB"]["drawFluid"] = draw_fluid_;
    var["PerFrameCB"]["lightScattering"] = lightScattering;

    var["PerFrameCB"]["backgroundColor"] = bg_clear_color;

    var["PerFrameCB"]["waterTurbulence"] = water_turbulence_;

    var["PerFrameCB"]["maxRayBounce"] = kMaxRayBounce;

    var["PerFrameCB"]["absorptionCoeff"] = absorptionCoeff;
    var["PerFrameCB"]["scatteringCoeff"] = scatteringCoeff;
    var["PerFrameCB"]["phaseG"] = phaseG;

    var["PerFrameCB"]["normalOffset"] = normalOffset;
    var["PerFrameCB"]["isoLevel"] = IsoLevel;
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

    var["PerFrameCB"]["DensityRayMarchMultiplier"] = DensityRayMarchMultiplier;
    var["PerFrameCB"]["densityMapSize"] = Metrics::density_map_size;
    var["PerFrameCB"]["simBounds"] = float3(Metrics::sim_bounds);
    var["PerFrameCB"]["volumeValueOffset"] = volumeValueOffset;

    var["gOutput"] = rt_output_tex_;
    var["gTexture3D"] = density_3d_tex_;
    var["linearClampSampler"] = linearClampSampler_;
}

void Renderer::createRasterizationProgram() const noexcept
{
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
}
