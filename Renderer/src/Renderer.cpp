#include "Renderer.h"
#include "Metrics.h"

#include "Core/Program/Program.h"
#include "Scene/Material/StandardMaterial.h"
#include "Utils/Math/FalcorMath.h"

Renderer::Renderer(const ref<Device>& device, const ref<Fbo>& target_fbo) noexcept
    : device_(device), target_fbo_(target_fbo), render_graph_(device, "FluidRenderGraph")
{}

void Renderer::Init(RenderContext* render_context) noexcept
{
    if (device_->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    {
        FALCOR_THROW("Device does not support raytracing!");
    }

    createRasterizationProgram();

    Settings settings{};

    // Create the SceneBuilder
    SceneBuilder::Flags flags = SceneBuilder::Flags::RTDontMergeStatic | SceneBuilder::Flags::RTDontMergeDynamic |
                                SceneBuilder::Flags::RTDontMergeInstanced | SceneBuilder::Flags::DontOptimizeGraph;
    scene_builder_ = new SceneBuilder(device_, settings, flags);

    auto sphere_mesh = TriangleMesh::createSphere(Metrics::MetersToPixels(0.05f));
    // auto cube_mesh = TriangleMesh::createCube(float3(Metrics::MetersToPixels(1.0f)) * 2.f);
    auto cube_mesh = TriangleMesh::createCube(float3(density_map_size));

    ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    dielectric_blue->setDoubleSided(true);
    dielectric_blue->setIndexOfRefraction(1.f);
    dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, dielectric_blue);
    cube_mesh_id = scene_builder_->addTriangleMesh(cube_mesh, dielectric_blue);

    // Create a lambertian material
    ref<Material> lambertian = StandardMaterial::create(device_, "Lambertian");
    lambertian->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
    lambertian->setRoughnessMollification(1.f);
    lambertian->setIndexOfRefraction(0.f);

    // lambertian->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));

    // ref<Material> dielectric_red = StandardMaterial::create(device_, "DielecRed");
    // dielectric_red->toBasicMaterial()->setBaseColor3(float3(1.f, 0.05f, 0.05f));
    // dielectric_red->setDoubleSided(true);
    // dielectric_red->setIndexOfRefraction(1.f);
    ////dielectric_red->toBasicMaterial()->setDiffuseTransmission(1.f);

    // ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    // dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    // dielectric_blue->setDoubleSided(true);
    // dielectric_blue->setIndexOfRefraction(1.f);
    // dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    // auto triangle_mesh_id_1 = scene_builder__.addTriangleMesh(sphere, dielectric_red);
    // auto triangle_mesh_id_2 = scene_builder__.addTriangleMesh(cube, dielectric_blue);
    // auto triangle_mesh_id_3 = scene_builder__.addTriangleMesh(cube, dielectric_blue);
    // auto triangle_mesh_id_4 = scene_builder_->addTriangleMesh(cube, lambertian);
    // auto sphere_mesh_id = scene_builder_->addTriangleMesh(sphere, dielectric_blue);

    // AABB raymarching_AABB = AABB(float3(-50, -5, -5) + float3(0, 10, 0),
    //     float3(50, 5, 5) + float3(0, 10, 0));
    // uint32_t raymarching_AABB_ID = 1;
    // scene_builder_->addCustomPrimitive(raymarching_AABB_ID, raymarching_AABB);

    // auto raymarching_node = SceneBuilder::Node();
    // raymarching_node.name = "RaymarchingNode";
    // auto raymarching_transform = Transform();
    // raymarching_transform.setTranslation(float3(0.f, 10.f, 0.f));
    // raymarching_transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    // raymarching_transform.setScaling(float3(1.f, 1.f, 1.f));
    // raymarching_node.transform = raymarching_transform.getMatrix();
    // auto raymarching_node_id = scene_builder_->addNode(raymarching_node);

    //auto node = SceneBuilder::Node();
    //node.name = "Cube simul bounds";
    //auto transform = Transform();
    //transform.setTranslation(float3(0.f, 0.f, 0));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(1, 1.f, 1.f));
    //node.transform = transform.getMatrix();
    //auto node_id = scene_builder_->addNode(node);

    //// Add Mesh Instances
    //scene_builder_->addMeshInstance(node_id, cube_mesh_id);

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

    // const auto half_density_map_size = 0.5f * sim_bounds;
    // AABB raymarching_AABB = AABB(float3(-half_density_map_size), float3(half_density_map_size));
    // uint32_t raymarching_AABB_ID = 1;
    // scene_builder_->addCustomPrimitive(raymarching_AABB_ID, raymarching_AABB);

    // auto raymarching_node = SceneBuilder::Node();
    // raymarching_node.name = "RaymarchingNode";
    // auto raymarching_transform = Transform();
    // raymarching_transform.setTranslation(float3(0.f, 0, 0.f));
    // raymarching_transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    // raymarching_transform.setScaling(float3(1.f, 1.f, 1.f));
    // raymarching_node.transform = raymarching_transform.getMatrix();
    // auto raymarching_node_id = scene_builder_->addNode(raymarching_node);

    // auto aabb_mesh = TriangleMesh::createCube(float3(density_map_size));

    // auto aabb_mesh_id = scene_builder_->addTriangleMesh(aabb_mesh, dielectric_blue);

    // auto raymarching_node = SceneBuilder::Node();
    // raymarching_node.name = "RaymarchingNode";
    // auto raymarching_transform = Transform();
    // raymarching_transform.setTranslation(float3(0.f, 0, 0.f));
    // raymarching_transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    // raymarching_transform.setScaling(float3(1.f, 1.f, 1.f));
    // raymarching_node.transform = raymarching_transform.getMatrix();
    // auto raymarching_node_id = scene_builder_->addNode(raymarching_node);

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

    std::vector<float> data(density_map_size * density_map_size * density_map_size, 0.0f);

    // for (int z = 0; z < density_map_size; ++z)
    //{
    //     for (int y = 0; y < density_map_size; ++y)
    //     {
    //         for (int x = 0; x < density_map_size; ++x)
    //         {
    //             int idx = z * density_map_size * density_map_size + y * density_map_size + x;

    //            // Option 1: Linear ramp along Z
    //            data[idx] = float(x) / float(density_map_size);

    //             // Option 2: Checker pattern
    //             // data[idx] = ((x + y + z) % 2 == 0) ? 1.0f : 0.0f;

    //             // Option 3: Spherical gradient
    //             /*float cx = density_map_size / 2.0f, cy = density_map_size / 2.0f, cz = density_map_size / 2.0f;
    //             float dx = x - cx, dy = y - cy, dz = z - cz;
    //             float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
    //             data[idx] = 1.0f - std::min(dist / (density_map_size / 2.0f), 1.0f);*/
    //             //data[idx] = 1;
    //        }
    //    }
    //}

    // auto quad = TriangleMesh::createQuad(float2(density_map_size, density_map_size));

    // auto id = scene_builder_->addTriangleMesh(quad, dielectric_blue);

    // auto node_4 = SceneBuilder::Node();
    // node_4.name = "cube";
    // auto transform_4 = Transform();
    // transform_4.setTranslation(float3(0.f, 0, 0.f));
    // transform_4.setRotationEuler(float3(1.5708, 0.f, 0));
    // transform_4.setScaling(float3(1, 1.f, 1));
    // node_4.transform = transform_4.getMatrix();
    // auto node_id_4 = scene_builder_->addNode(node_4);

    //// Add Mesh Instances
    // scene_builder_->addMeshInstance(node_id_4, id);

    // density_3d_tex_ = device_->createTexture3D(
    //     density_map_size,
    //     density_map_size,
    //     density_map_size,
    //     ResourceFormat::R32Float,
    //     1, // mips
    //     data.data(),
    //     ResourceBindFlags::ShaderResource
    //);

    // density_3d_tex_ = device_->createTexture3D(
    //     density_map_size,
    //     density_map_size,
    //     density_map_size,
    //     ResourceFormat::R32Float,
    //     1, // mips
    //     particle_densities_->data(),
    //     ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    //);

    // density_map_pass_ = ComputePass::create(device_,
    //     "Samples/3DFluidSimulationEngine/Renderer/shaders/DensityMap.cs.slang", "createDensityMap");

    // render_graph_.addPass(density_map_pass_, "ComputeDensityMap");

        // 1) Buffer allocation
    const auto maxTriangleCount = 5 * (numPointsPerAxis - 1) * (numPointsPerAxis - 1) * (numPointsPerAxis - 1);
    const size_t triangleStructSize = sizeof(MarchingCubesTriangle); // 3 verts × float3 = 36 bytes

    std::vector<MarchingCubesTriangle> test(maxTriangleCount);

    // The structured buffer that your HLSL AppendStructuredBuffer<Triangle> will write into:
    marching_cubes_triangle_buffer_ = make_ref<Buffer>(
        device_,                               // Falcor device
        triangleStructSize,                    // structSize (bytes per element)
        maxTriangleCount,                      // elementCount
        ResourceBindFlags::UnorderedAccess |   // UAV for compute
            ResourceBindFlags::ShaderResource, // SRV for rendering or readback
        MemoryType::DeviceLocal,               // GPU-only (faster)
        test.data(),                                  // no init data
        true                                   // create hidden counter
    );

    // A small readback buffer to fetch the append counter:
    read_back_triangle_buffer_ = make_ref<Buffer>(
        device_,
        triangleStructSize,
        maxTriangleCount, 
        ResourceBindFlags::None,
        MemoryType::ReadBack,
        nullptr,
        false
    );

    compute_density_map_pass_ =
        ComputePass::create(device_,
            "Samples/3DFluidSimulationEngine/Renderer/shaders/MarchingCubes.cs.slang",
            "ComputeDensityTexture");

    density_3d_tex_ = device_->createTexture3D(
        density_map_size,
        density_map_size,
        density_map_size,
        ResourceFormat::R32Float,
        1, // mips
        nullptr,
        ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
    );

     marching_cubes_pass_ =
        ComputePass::create(device_,
            "Samples/3DFluidSimulationEngine/Renderer/shaders/MarchingCubes.cs.slang",
            "ProcessCube");

    auto compute_var = compute_density_map_pass_->getRootVar();
    compute_var["DensityTexture"] = density_3d_tex_;
    compute_var["triangles"] = marching_cubes_triangle_buffer_;

    compute_var = marching_cubes_pass_->getRootVar();
    compute_var["DensityTexture"] = density_3d_tex_;
    compute_var["triangles"] = marching_cubes_triangle_buffer_;

    compute_var = compute_density_map_pass_->getRootVar();
    compute_var["PerFrameCB"]["numPointsPerAxis"] = 64;
    compute_var["PerFrameCB"]["isoLevel"] = 0;
    compute_var["PerFrameCB"]["textureSize"] = density_map_size;
    compute_var["PerFrameCB"]["boundSize"] = sim_bounds.x;
    compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;

    compute_density_map_pass_->execute(render_context, 64, 64, 64);

    compute_var = marching_cubes_pass_->getRootVar();
    compute_var["PerFrameCB"]["numPointsPerAxis"] = 64;
    compute_var["PerFrameCB"]["isoLevel"] = 0;
    compute_var["PerFrameCB"]["textureSize"] = density_map_size;
    compute_var["PerFrameCB"]["boundSize"] = sim_bounds.x;
    compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;

    marching_cubes_pass_->execute(render_context, 64, 64, 64);

    render_context->copyResource(read_back_triangle_buffer_.get(), marching_cubes_triangle_buffer_.get());

    const MarchingCubesTriangle* triangles = static_cast<const MarchingCubesTriangle*>(read_back_triangle_buffer_->map());

 /*   for (int i = 0; i < maxTriangleCount; i++)
    {
        std::cout << triangles[i].vertexA.position.x << " "
        << triangles[i].vertexA.position.y << " " << triangles[i].vertexA.position.z
                  << '\n';
    }*/

    // 3) Extract positions and build a linear index list
    //TriangleMesh::VertexList positions;
    //TriangleMesh::IndexList indices;
    //positions.reserve(maxTriangleCount * 3);
    //indices.reserve(maxTriangleCount * 3);

    //for (uint32_t i = 0; i < maxTriangleCount; i++)
    //{
    //    /*if (triangles[i].vertexA.position.x == -20000.f)
    //    {
    //        std::cout << "STOOOOP\n";
    //        break;
    //    }*/

    //    // read the three vertices
    //    positions.push_back({triangles[i].vertexA.position, {0, 0, 1}, {0, 0}});
    //    positions.push_back({triangles[i].vertexB.position, {0, 0, 1}, {0, 0}});
    //    positions.push_back({triangles[i].vertexC.position, {0, 0, 1}, {0, 0}});
    //    // line up indices as 0,1,2, 3,4,5, ...
    //    uint32_t base = i * 3;
    //    indices.push_back(base + 0);
    //    indices.push_back(base + 1);
    //    indices.push_back(base + 2);
    //}

    //// 4) Unmap when done reading
    //read_back_triangle_buffer_->unmap();

    //auto marching_cubes_mesh = TriangleMesh::create(positions, indices);
    //auto marchingCubesMeshID = scene_builder_->addTriangleMesh(marching_cubes_mesh, dielectric_blue);

    //
    //auto node = SceneBuilder::Node();
    //node.name = "Marching Cubes Node";
    //auto transform = Transform();
    //transform.setTranslation(float3(0.f, 0.f, 0));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(1, 1.f, 1.f));
    //node.transform = transform.getMatrix();
    //auto node_id = scene_builder_->addNode(node);

    //// Add Mesh Instances
    //scene_builder_->addMeshInstance(node_id, marchingCubesMeshID);
}

void Renderer::RenderFrame(RenderContext* pRenderContext, const double& currentTime) noexcept
{
    pRenderContext->clearFbo(target_fbo_.get(), float4(bg_clear_color, 1), 1.0f, 0, FboAttachmentType::All);

    //static bool marchingCubes = false;

    //if (!marchingCubes)
    //{
    //    //marchingCubes = true;

    //     // 1) Clear the hidden append-counter.  Must happen *before* you dispatch.
    //    pRenderContext->clearUAVCounter(marching_cubes_triangle_buffer_, 0);

    //    auto compute_var = compute_density_map_pass_->getRootVar();
    //    compute_var["DensityTexture"] = density_3d_tex_;
    //    //compute_var["triangles"] = marching_cubes_triangle_buffer_;

    //    /*compute_var["PerFrameCB"]["numPointsPerAxis"] = numPointsPerAxis;
    //    compute_var["PerFrameCB"]["isoLevel"] = IsoLevel;*/
    //    compute_var["PerFrameCB"]["textureSize"] = density_map_size;
    //    compute_var["PerFrameCB"]["boundSize"] = sim_bounds.x;
    //    compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;

    //    compute_density_map_pass_->execute(pRenderContext, 64, 64, 64);

    //    compute_var = marching_cubes_pass_->getRootVar();
    //    compute_var["DensityTexture"] = density_3d_tex_;
    //    compute_var["triangles"] = marching_cubes_triangle_buffer_;

    //    compute_var["PerFrameCB"]["numPointsPerAxis"] = numPointsPerAxis;
    //    compute_var["PerFrameCB"]["isoLevel"] = IsoLevel;
    //    compute_var["PerFrameCB"]["textureSize"] = density_map_size;
    //    compute_var["PerFrameCB"]["boundSize"] = sim_bounds.x;
    //    compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;

    //    marching_cubes_pass_->execute(pRenderContext, 64, 64, 64);

    //    //TODO: suivre chatgpt -> UAVCounter est un buffer, il faut donc lire dedans, pas simplement le print
    //    //TODO: capter pk frame 0 ne marche pas pour l'algo marching cubes.
    //    //TODO: voir script GenTest.cs de seb lague, il fait un buffer pour copier le nbr de triangles
    //    //TODO: ATTENTION au memory layout de MarchCUbeTriangle. ca doit être 3x16 = 48. actuellemtn c'est 36 sur CPU.
    //    //TODO: https://chatgpt.com/c/681b50ec-3ed0-8003-913f-bb3f1dec326c

    //    std::cout << "UAV counter: " << marching_cubes_triangle_buffer_->getUAVCounter()->getElement<uint>(0)
    //    << '\n';

    //    pRenderContext->copyResource(read_back_triangle_buffer_.get(), marching_cubes_triangle_buffer_.get());

    //    const MarchingCubesTriangle* triangles = static_cast<const MarchingCubesTriangle*>(read_back_triangle_buffer_->map());

    //    std::cout << triangles[100].vertexA.position.x << '\n';
    //    std::cout << triangles[500].vertexA.position.x << '\n';
    //    std::cout << triangles[324].vertexA.position.x << '\n';
    //}

  

    // const auto compute_var = density_map_pass_->getRootVar();
    // compute_var["gTexture3D"] = density_3d_tex_;
    // compute_var["PerFrameCB"]["densityMapSize"] = density_map_size;
    // compute_var["PerFrameCB"]["simBounds"] = sim_bounds;
    // compute_var["particleDensities"] = particle_density_buffer_;

    // const uint32_t thread_groups = (density_map_size + 7) / 8;
    // density_map_pass_->execute(pRenderContext, 64, 64, 64);

    // pRenderContext->updateTextureData(density_3d_tex_.get(), particle_densities_->data());

    IScene::UpdateFlags updates = scene_->update(pRenderContext, currentTime);
    if (is_set(updates, IScene::UpdateFlags::GeometryChanged))
        FALCOR_THROW("This sample does not support scene geometry changes.");
    if (is_set(updates, IScene::UpdateFlags::RecompileNeeded))
        FALCOR_THROW("This sample does not support scene changes that require shader recompilation.");

    // ALCOR_ASSERT(mpScene)
    //  FALCOR_PROFILE(pRenderContext, "renderRT");

    setPerFrameVariables(currentTime);

    pRenderContext->clearUAV(rt_output_tex_->getUAV().get(), float4(bg_clear_color, 1));
    scene_->raytrace(pRenderContext, rt_program_.get(), rt_program_vars_, uint3(target_fbo_->getWidth(), target_fbo_->getHeight(), 1));
    pRenderContext->blit(rt_output_tex_->getSRV(), target_fbo_->getRenderTargetView(0));
}

void Renderer::RenderUI(Gui* pGui, Gui::Window* app_gui_window) noexcept
{
    app_gui_window->rgbColor("Background color", bg_clear_color);

    app_gui_window->slider("DensityDepth", DensityDepth, 0.f, 200.f);
    app_gui_window->slider("SphereRadius", SphereRadius, 0.f, 200.f);

    app_gui_window->checkbox("Draw Fluid ?", draw_fluid_);
    if (draw_fluid_)
    {
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
    }

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

NodeID Renderer::AddSphereToScene(const float3 pos, const float radius) noexcept
{
    // ref<TriangleMesh> sphere_mesh = TriangleMesh::createSphere(radius);

    // ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    // dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    // dielectric_blue->setDoubleSided(true);
    // dielectric_blue->setIndexOfRefraction(1.f);
    // dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    // MeshID sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, dielectric_blue);

    auto node = SceneBuilder::Node();
    const std::string name = "Sphere " /* + std::to_string(i)*/;
    node.name = name;
    auto transform = Transform();
    transform.setTranslation(pos);
    transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    transform.setScaling(float3(radius, radius, radius));
    node.transform = transform.getMatrix();
    const auto node_id = scene_builder_->addNode(node);

    // sphereNodeIDs.push_back(node_id);

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

    var["PerFrameCB"]["backgroundColor"] = bg_clear_color;

    var["PerFrameCB"]["waterTurbulence"] = water_turbulence_;

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

    var["PerFrameCB"]["DensityDepth"] = DensityDepth;
    var["PerFrameCB"]["densityMapSize"] = density_map_size;

    var["gOutput"] = rt_output_tex_;
    var["gTexture3D"] = density_3d_tex_;
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
