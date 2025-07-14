#include "Renderer.h"

#include <DirectXMath.h>
#include <Tracy.hpp>

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

void Renderer::Init(RenderContext* render_context, bool rebuildBvh) noexcept
{
    if (device_->isFeatureSupported(Device::SupportedFeatures::Raytracing) == false)
    {
        FALCOR_THROW("Device does not support raytracing!");
    }

    Settings settings{};

    // Create the SceneBuilder
    SceneBuilder::Flags flags = SceneBuilder::Flags::RTDontMergeStatic | SceneBuilder::Flags::RTDontMergeDynamic |
                                SceneBuilder::Flags::RTDontMergeInstanced | SceneBuilder::Flags::DontOptimizeGraph;
    scene_builder_ = new SceneBuilder(device_, settings, flags);

    ref<Material> dielectric_blue = StandardMaterial::create(device_, "DielecBlue");
    dielectric_blue->toBasicMaterial()->setBaseColor3(float3(0.05f, 0.05f, 1.0f));
    dielectric_blue->setDoubleSided(true);
    dielectric_blue->setIndexOfRefraction(1.f);
    dielectric_blue->toBasicMaterial()->setDiffuseTransmission(1.f);

    // cube_mesh_id = scene_builder_->addTriangleMesh(cube_mesh, dielectric_blue);

    // Create a lambertian material
    lambertianSphere = StandardMaterial::create(device_, "LambertianSphere");
    lambertianSphere->toBasicMaterial()->setBaseColor3(float3(0.1f, 0.33f, 0.9f));
    lambertianSphere->setRoughnessMollification(1.f);
    lambertianSphere->setIndexOfRefraction(0.f);

    lambertianCube = StandardMaterial::create(device_, "LambertianCube");
    lambertianCube->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
    lambertianCube->setRoughnessMollification(1.f);
    lambertianCube->setIndexOfRefraction(0.f);

    ref<Material> lambertianTexture = StandardMaterial::create(device_, "LambertianTexture");
    auto texture = Texture::createFromFile(device_,
        "data/images/CheckerTile_BaseColor.png", true, false);
 /*   auto texture = Texture::createFromFile(device_,
        "Samples/3DFluidSimulationEngine/data/images/CheckerTile_BaseColor.png",
        true, false);*/
    lambertianTexture->toBasicMaterial()->setBaseColorTexture(texture);
    lambertianTexture->setRoughnessMollification(1.f);
    lambertianTexture->setIndexOfRefraction(0.f);

    if (useMarchingCubes)
    {
        // 1) Buffer allocation
        MaxTriangleCount =
            static_cast<size_t>(5) * (Metrics::density_map_size - 1) * (Metrics::density_map_size - 1) * (Metrics::density_map_size - 1);
        MaxVertexCount = MaxTriangleCount * 3;
        // constexpr size_t triangleStructSize = sizeof(MarchingCubesTriangle);

        // v = {
        //    {float3(0.0f, 1.0f, -10), float3(0.0f, 0.0f, 1.0f), float2(0.5f, 1.0f)},   // Top
        //    {float3(-1.0f, -1.0f, -10), float3(0.0f, 0.0f, 1.0f), float2(0.0f, 0.0f)}, // Left
        //    {float3(1.0f, -1.0f, -10), float3(0.0f, 0.0f, 1.0f), float2(1.0f, 0.0f)}   // Right
        //};

        //    v = {
        //        //        position               normal             texcoord
        //        {{-0.5f, -0.5f, -0.5f}, {0, 0, -1}, {0, 0}}, // 0 - Back face
        //        {{0.5f, -0.5f, -0.5f}, {0, 0, -1}, {1, 0}},  // 1
        //        {{0.5f, 0.5f, -0.5f}, {0, 0, -1}, {1, 1}},   // 2
        //        {{-0.5f, 0.5f, -0.5f}, {0, 0, -1}, {0, 1}},  // 3
        //
        //        {{-0.5f, -0.5f, 0.5f}, {0, 0, 1}, {0, 0}}, // 4 - Front face
        //        {{0.5f, -0.5f, 0.5f}, {0, 0, 1}, {1, 0}},  // 5
        //        {{0.5f, 0.5f, 0.5f}, {0, 0, 1}, {1, 1}},   // 6
        //        {{-0.5f, 0.5f, 0.5f}, {0, 0, 1}, {0, 1}},  // 7
        //    };
        //
        //

        //
        //    std::vector<uint32_t> indices = {
        //    // Back face
        //    0, 1, 2,
        //    0, 2, 3,
        //
        //    // Front face
        //    4, 6, 5,
        //    4, 7, 6,
        //
        //    // Left face
        //    4, 5, 1,
        //    4, 1, 0,
        //
        //    // Right face
        //    3, 2, 6,
        //    3, 6, 7,
        //
        //    // Bottom face
        //    4, 0, 3,
        //    4, 3, 7,
        //
        //    // Top face
        //    1, 5, 6,
        //    1, 6, 2,
        //};

        TriangleMesh::Vertex vert{float3(10, 10, 10), float3(0, 0, 1), float2(0, 1)};
        v.resize(MaxVertexCount, vert);

        TriangleMesh::IndexList indices;
        for (unsigned idx = 0; idx < MaxVertexCount; idx++)
        {
            indices.push_back(idx);
        }

        auto tri_mesh = TriangleMesh::create(v, indices);
        tri_id = scene_builder_->addTriangleMesh(tri_mesh, dielectric_blue, true);

        auto node_tri = SceneBuilder::Node();
        auto tri_name = "Marhcing cube mesh " /* + std::to_string(i)*/;
        node_tri.name = tri_name;
        auto transform_tri = Transform();
        transform_tri.setTranslation(float3(0, 0.f, 0.f));
        transform_tri.setRotationEuler(float3(0.f, 0.f, 0.f));
        transform_tri.setScaling(float3(1, 1, 1));
        node_tri.transform = transform_tri.getMatrix();
        auto tri_node_id = scene_builder_->addNode(node_tri);

        scene_builder_->addMeshInstance(tri_node_id, tri_id);
    }

    if (!useMarchingCubes)
    {
        auto sphere_mesh = TriangleMesh::createSphere(SPH::SmoothingRadius);
        sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, lambertianSphere);

        auto cube_mesh = TriangleMesh::createCube(float3(Metrics::WALLDIST / 4.f,
            Metrics::WALLDIST * 2.5f, Metrics::WALLDIST / 4.f));
        cube_mesh_id = scene_builder_->addTriangleMesh(cube_mesh, lambertianCube);

        auto node = SceneBuilder::Node();
        const std::string name = "Cube ";
        node.name = name;
        auto transform = Transform();
        transform.setTranslation(float3(0, 0, 0));
        transform.setRotationEuler(float3(35.f, 0.f, 0.f));
        transform.setScaling(float3(1, 1, 1));
        node.transform = transform.getMatrix();
        const auto node_id = scene_builder_->addNode(node);

         // Add Mesh Instances
        scene_builder_->addMeshInstance(node_id, cube_mesh_id);

        auto plan_mesh = TriangleMesh::createQuad(float2(500, 500));
        plane_mesh_id = scene_builder_->addTriangleMesh(plan_mesh, lambertianTexture);

        auto node_p = SceneBuilder::Node();
        const std::string name_p = "Plane";
        node_p.name = name;
        auto transform_p = Transform();
        transform_p.setTranslation(float3(0, -Metrics::WALLDIST + 1, 0));
        transform_p.setRotationEuler(float3(0, 0.f, 0.f));
        transform_p.setScaling(float3(1, 1, 1));
        node_p.transform = transform_p.getMatrix();
        const auto node_id_p = scene_builder_->addNode(node_p);

        // Add Mesh Instances
        scene_builder_->addMeshInstance(node_id_p, plane_mesh_id);
    }


    // auto cube_mesh = TriangleMesh::createCube(float3(Metrics::sim_bounds -1.f));

    // sphere = TriangleMesh::createQuad(float2(5.f));
    // sphere_mesh_id = scene_builder_->addTriangleMesh(sphere, dielectric_blue, true);

    //auto node = SceneBuilder::Node();
    //std::string name = "Sphere " /* + std::to_string(i)*/;
    //node.name = name;
    //auto transform = Transform();
    //transform.setTranslation(float3(50.f, 0.f, 0.f));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(20, 20, 20));
    //node.transform = transform.getMatrix();
    //sphere_node_id_ = scene_builder_->addNode(node);

    //scene_builder_->addMeshInstance(sphere_node_id_, sphere_mesh_id);

    //cube_mesh_id = scene_builder_->addTriangleMesh(cube_mesh, dielectric_blue, true);

    //node = SceneBuilder::Node();
    //auto name2 = "Cube " /* + std::to_string(i)*/;
    //node.name = name;
    //transform = Transform();
    //transform.setTranslation(float3(0, 0.f, 0.f));
    //transform.setRotationEuler(float3(0.f, 0.f, 0.f));
    //transform.setScaling(float3(1, 1, 1));
    //node.transform = transform.getMatrix();
    //auto cube_node = scene_builder_->addNode(node);

    ////// Add Mesh Instances
    //scene_builder_->addMeshInstance(cube_node, cube_mesh_id);

    if (!useMarchingCubes)
    {
        float3 voxelGridRes{Metrics::voxelGridResolution[0], Metrics::voxelGridResolution[1], Metrics::voxelGridResolution[2]};
        float3 cellSize = float3(Metrics::sim_bounds) / voxelGridRes;
        float3 halfSize = cellSize * 0.5f;

        int primIndex = 0;

        //for (int x = 0; x < Metrics::voxelGridResolution[0]; ++x)
        //{
        //    for (int y = 0; y < Metrics::voxelGridResolution[1]; ++y)
        //    {
        //        for (int z = 0; z < Metrics::voxelGridResolution[2]; ++z)
        //        {
        //            float3 center = float3(
        //                -Metrics::WALLDIST + (x + 0.5f) * cellSize.x,
        //                -Metrics::WALLDIST + (y + 0.5f) * cellSize.y,
        //                -Metrics::WALLDIST + (z + 0.5f) * cellSize.z
        //            );

        //            float3 min = center - float3(halfSize);
        //            float3 max = center + float3(halfSize);
        //            AABB box(min, max);

        //            // Add the AABB with a unique ID (must not conflict with other primitives)
        //            scene_builder_->addCustomPrimitive(primIndex, box);

        //            //// Step 2 — Create and add node
        //            //auto node = SceneBuilder::Node();
        //            //node.name = "AABB_Cell_" + std::to_string(primIndex);

        //            //// Set transform (here: identity, since AABB is already in world space)
        //            //auto transform = Transform();
        //            //transform.setTranslation(center);
        //            //transform.setRotationEulerDeg(float3(0));
        //            //transform.setScaling(float3(1));

        //            //node.transform = transform.getMatrix();

        //            //// Associate the node with the primitive
        //            //scene_builder_->addNode(node);

        //            primIndex++;
        //        }
        //    }
        //}

        cutomPrimitveMasks = make_ref<Buffer>(
            device_,
            sizeof(uint32_t),
            masks.size(),
            ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
            MemoryType::DeviceLocal,
            masks.data(),
            false
        );

        const AABB fluid_AABB = AABB(float3(-1), float3(1));
        //AABB fluid_AABB = AABB(float3(-Metrics::WALLDIST), float3(Metrics::WALLDIST));

        fluid_transform = Transform();
        fluid_transform.setTranslation(translation);
        fluid_transform.setRotationEulerDeg(rotation);
        fluid_transform.setScaling(scale);

        const AABB transformed_aabb = fluid_AABB.transform(fluid_transform.getMatrix());
        scene_builder_->addCustomPrimitive(fluid_AABB_ID, transformed_aabb);

        auto fluid_node = SceneBuilder::Node();
        fluid_node.name = "RaymarchingNode";

        fluid_node.transform = fluid_transform.getMatrix();
        raymarching_node_id = scene_builder_->addNode(fluid_node);
    }

    //auto envMap = EnvMap::createFromFile(device_,
    //    "Samples/3DFluidSimulationEngine/data/images/hallstatt4_hd.hdr");
    auto envMap = EnvMap::createFromFile(device_, "data/images/hallstatt4_hd.hdr");
    envMap->setIntensity(1.0);
    scene_builder_->setEnvMap(envMap);

    ref<Camera> camera = ref<Camera>(new Camera("Camera"));
    camera->setPosition(float3(0, 0.0, -100));
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
    sampler_desc.setFilterMode(TextureFilteringMode::Linear, TextureFilteringMode::Linear, TextureFilteringMode::Linear);
    sampler_desc.setAddressingMode(TextureAddressingMode::Clamp, TextureAddressingMode::Clamp, TextureAddressingMode::Clamp);
    linearClampSampler_ = make_ref<Sampler>(device_, sampler_desc);

    if (!rebuildBvh)
    {
        if (useMarchingCubes)
        {
            marching_cube_dens_tex = device_->createTexture3D(
                Metrics::density_map_size,
                Metrics::density_map_size,
                Metrics::density_map_size,
                ResourceFormat::R32Float,
                1, // mips
                nullptr,
                ResourceBindFlags::UnorderedAccess | ResourceBindFlags::ShaderResource
            );

            // The structured buffer that your HLSL AppendStructuredBuffer<Triangle> will write into:
            marching_cubes_triangle_buffer_ = make_ref<Buffer>(
                device_,                               // Falcor device
                sizeof(MarchingCubesTriangle),         // structSize (bytes per element)
                MaxTriangleCount,                      // elementCount
                ResourceBindFlags::UnorderedAccess |   // UAV for compute
                    ResourceBindFlags::ShaderResource, // SRV for rendering or readback
                MemoryType::DeviceLocal,               // GPU-only (faster)
                nullptr,                               // no init data
                true                                   // create hidden counter
            );

            marching_cubes_pass_ =
                ComputePass::create(device_, "Samples/3DFluidSimulationEngine/Renderer/shaders/MarchingCubes.cs.slang", "ProcessCube");

            /*compute_marching_cube_density_map_ = ComputePass::create(
                device_, "Samples/3DFluidSimulationEngine/Renderer/shaders/MarchingCubes.cs.slang", "ComputeDensityTexture"
            );*/

            // A small readback buffer to fetch the append counter:
            read_back_triangle_buffer_ = make_ref<Buffer>(
                device_, sizeof(MarchingCubesTriangle), MaxTriangleCount, ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false
            );

            std::vector<float3> positions;
            std::vector<float3> normals;
            std::vector<float3> tangents;
            std::vector<float2> uv;
            for (const auto& vertex : v)
            {
                // std::cout << vertex.position.x << " " << vertex.position.y << " " << vertex.position.z << "\n";
                const auto new_Pos = vertex.position; // + float3(5, 1, 0);
                positions.push_back(new_Pos);
                normals.push_back(vertex.normal);
                tangents.push_back(vertex.normal);
                uv.push_back(float2(vertex.texCoord.x, vertex.texCoord.y));
            }
            // for (const auto& vertex : tri_mesh->getVertices())
            //{
            //     // std::cout << vertex.position.x << " " << vertex.position.y << " " << vertex.position.z << "\n";
            //     const auto new_Pos = vertex.position; // + float3(5, 1, 0);
            //     positions.push_back(new_Pos);
            //     normals.push_back(vertex.normal);
            //     tangents.push_back(vertex.normal);
            //     uv.push_back(float2(vertex.texCoord.x, vertex.texCoord.y));
            // }

            b_pos = make_ref<Buffer>(
                device_,
                sizeof(float) * 3,
                positions.size(),
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::DeviceLocal,
                positions.data(),
                false
            );

            b_pos_readback = make_ref<Buffer>(
                device_, sizeof(float) * 3, positions.size(), ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false
            );

            b_normal = make_ref<Buffer>(
                device_,
                sizeof(normals[0]),
                normals.size(),
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::Upload,
                normals.data(),
                false
            );

            b_norm_readback = make_ref<Buffer>(
                device_, sizeof(normals[0]), normals.size(), ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false
            );

            b_tang = make_ref<Buffer>(
                device_,
                sizeof(tangents[0]),
                tangents.size(),
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::Upload,
                tangents.data(),
                false
            );

            b_tang_readback = make_ref<Buffer>(
                device_, sizeof(tangents[0]), tangents.size(), ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false
            );

            b_uv = make_ref<Buffer>(
                device_,
                sizeof(uv[0]),
                uv.size(),
                ResourceBindFlags::ShaderResource | ResourceBindFlags::UnorderedAccess,
                MemoryType::Upload,
                uv.data(),
                false
            );

            b_uv_readback =
                make_ref<Buffer>(device_, sizeof(uv[0]), uv.size(), ResourceBindFlags::None, MemoryType::ReadBack, nullptr, false);

            vertices = {
                {"positions", b_pos},
                {"normals", b_normal},
                {"tangents", b_tang},
                {"texcrds", b_uv},
            };
        }
    }
    else
    {
        if (useMarchingCubes)
        {
            scene_->setMeshVertices(
                tri_id,
                {
                    {"positions", b_pos},
                    {"normals", b_normal},
                    {"tangents", b_tang},
                    {"texcrds", b_uv},
                }
            );
        }
    }
}

void Renderer::RenderFrame(
    RenderContext* pRenderContext, const ref<Fbo>& pTargetFbo,  const double& currentTime,
    const ref<Buffer>& bodies,
    const ref<Buffer>& SpatialIndices,
    const ref<Buffer>& SpatialOffsets) noexcept
{
    pRenderContext->clearFbo(target_fbo_.get(),
        float4(bg_clear_color, 1), 1.0f, 0,
        FboAttachmentType::All);

    if (!useMarchingCubes)
    {
        {
            FALCOR_PROFILE(pRenderContext, "Reset raymarching custom primitve masks");

            cutomPrimitveMasks->setBlob(masks.data(), 0, masks.size() * sizeof(uint32_t));
        }

        {
            FALCOR_PROFILE(pRenderContext, "Update Raymarching custom primitve transformation");

            const AABB fluid_AABB = AABB(float3(-1), float3(1));
            //const AABB fluid_AABB = AABB(float3(-Metrics::WALLDIST), float3(Metrics::WALLDIST));
            const AABB transformed_aabb = fluid_AABB.transform(fluid_transform.getMatrix());

            scene_->updateCustomPrimitive(0, transformed_aabb);
        }
    }

    {
        FALCOR_PROFILE(pRenderContext, "Compute DensityMap Pass");

        const auto compute_var = compute_density_map_pass_->getRootVar();
        compute_var["bodies"] = bodies;
        compute_var["SpatialIndices"] = SpatialIndices;
        compute_var["SpatialOffsets"] = SpatialOffsets;
        compute_var["cutomPrimitveMasks"] = cutomPrimitveMasks;
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

        //compute_var["PerFrameCB"]["ScaledSimBounds"] = fluid_transform.getScaling() * 2.f;

        compute_var["PerFrameCB"]["useTransformations"] = useTransformations;
        uint3 voxelGridRes{Metrics::voxelGridResolution[0], Metrics::voxelGridResolution[1], Metrics::voxelGridResolution[2]};
        compute_var["PerFrameCB"]["voxelGridResolution"] = voxelGridRes;

        const float r = SPH::SmoothingRadius;
        const float spikyPow2 = 15.f / (2 * PI * Pow(r, 5));
        const float spikyPow3 = 15.f / (PI * Pow(r, 6));
        const float spikyPow2Grad = 15.f / (PI * Pow(r, 5));
        const float spikyPow3Grad = 45.f / (PI * Pow(r, 6));

        compute_var["PerFrameCB"]["K_SpikyPow2"] = spikyPow2;
        compute_var["PerFrameCB"]["K_SpikyPow3"] = spikyPow3;
        compute_var["PerFrameCB"]["K_SpikyPow2Grad"] = spikyPow2Grad;
        compute_var["PerFrameCB"]["K_SpikyPow3Grad"] = spikyPow3Grad;

        // const uint32_t thread_groups = (density_map_size + 7) / 8;
        compute_density_map_pass_->execute(pRenderContext,
             Metrics::density_map_size, Metrics::density_map_size, Metrics::density_map_size);
    }


    if (draw_fluid_ && useMarchingCubes)
    {
        LaunchMarchingCubeComputePasses(pRenderContext);

        //pRenderContext->copyResource(b_pos_readback.get(), b_pos.get());
        //pRenderContext->copyResource(b_norm_readback.get(), b_normal.get());
        // pRenderContext->copyResource(b_tang_readback.get(), b_tang.get());
        // pRenderContext->copyResource(b_uv_readback.get(), b_uv.get());

        //const float3* poses = static_cast<const float3*>(b_pos_readback->map());
        // const float3* normals = static_cast<const float3*>(b_norm_readback->map());
        // const float3* tangents = static_cast<const float3*>(b_tang_readback->map());
        // const float2* uvs = static_cast<const float2*>(b_uv_readback->map());

        // std::cout << marching_cube_vertex_count << '\n';
        /*  for (int i = 0; i < v.size(); i++)
          {
              std::cout << "POS[" << i << "] " << poses[i].x << " " << poses[i].y << " " << poses[i].z << '\n';
              std::cout << "NORMAL[" << i << "] " << normals[i].x << " " << normals[i].y << " " << normals[i].z << '\n';
              std::cout << "TANGENT[" << i << "] " << tangents[i].x << " " << tangents[i].y << " " << tangents[i].z << '\n';
              std::cout << "UV[" << i << "] " << uvs[i].x << " " << uvs[i].y << '\n';
          }*/

        // std::cout << "POS 0: " << poses[0].x << " " << poses[0].y << " " << poses[0].z << '\n';
        // std::cout << "POS " << 1 << ": " << poses[1].x << " "
        //           << poses[1].y << " " << poses[1].z
        //           << '\n';
        // std::cout << "POS " << 10 << ": " << poses[10].x << " " << poses[10].y << " " << poses[10].z << '\n';
        // std::cout << "POS " << marching_cube_vertex_count << ": " << poses[marching_cube_vertex_count].x << " "
        //           << poses[marching_cube_vertex_count].y << " " << poses[marching_cube_vertex_count].z << '\n';

        /*   b_pos_readback->unmap();
           b_normal->unmap();
           b_tang->unmap();
           b_uv->unmap();*/

        //if (b_pos->getElementCount() != scene_->getMesh(tri_id).getVertexCount())
        //{
        //    std::cout << b_pos->getElementCount() << '\n';
        //    std::cout << scene_->getMesh(tri_id).getVertexCount() << '\n';
        //    std::cout << scene_->getMesh(sphere_mesh_id).getVertexCount() << '\n';
        //    std::cout << "BUG VERTEX COUNT AND B_POS\n";
        //    std::exit(666);
        //}

        /*if (triangleCount != oldTriangleCount)
        {
            scene_->setMeshVertices(tri_id, vertices);
        }*/

        /*static bool tamere = false;
        if (tamere)
        {
            scene_->setMeshVertices(tri_id, vertices);
            tamere = true;
        }*/

         
        {
            FALCOR_PROFILE(pRenderContext, "Set Mesh Vertices");
            scene_->setMeshVertices(tri_id, vertices);
        }
    }

    {
        FALCOR_PROFILE(pRenderContext, "Update Scene");

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
    }

 /*   FALCOR_ASSERT(scene_);
    FALCOR_PROFILE(pRenderContext, "renderRT");*/

    {
        //FALCOR_PROFILE(pRenderContext, "Render Raytracing Pass");

        setPerFrameVariables(currentTime);

        pRenderContext->clearUAV(rt_output_tex_->getUAV().get(), float4(bg_clear_color, 1));

        /* raster_pass_->getState()->setFbo(pTargetFbo);
         scene_->rasterize(pRenderContext, raster_pass_->getState().get(), raster_pass_->getVars().get());*/

        scene_->raytrace(pRenderContext, rt_program_.get(), rt_program_vars_, uint3(target_fbo_->getWidth(), target_fbo_->getHeight(), 1));
        pRenderContext->blit(rt_output_tex_->getSRV(), target_fbo_->getRenderTargetView(0));
    }

    oldTriangleCount = triangleCount;
}

void Renderer::RenderUI(Gui* pGui, Gui::Window* app_gui_window, RenderContext* render_context) noexcept
{
    app_gui_window->rgbColor("Background color", bg_clear_color);
    lambertianCube->toBasicMaterial()->setBaseColor3(bg_clear_color);

     if (app_gui_window->button("Reload Scene"))
     {
         Init(render_context, true);
         auto sphere_mesh = TriangleMesh::createSphere(Metrics::PARTICLESIZE);

         lambertianCube = StandardMaterial::create(device_, "Lambertian");
         lambertianCube->toBasicMaterial()->setBaseColor3(float3(0.2f, 0.9f, 0.1f));
         lambertianCube->setRoughnessMollification(1.f);
         lambertianCube->setIndexOfRefraction(0.f);

         sphere_mesh_id = scene_builder_->addTriangleMesh(sphere_mesh, lambertianCube, true);

         auto node = SceneBuilder::Node();
         std::string name = "Sphere " /* + std::to_string(i)*/;
         node.name = name;
         auto transform = Transform();
         transform.setTranslation(float3(50.f, 0.f, 0.f));
         transform.setRotationEuler(float3(0.f, 0.f, 0.f));
         transform.setScaling(float3(20, 20, 20));
         node.transform = transform.getMatrix();
         sphere_node_id_ = scene_builder_->addNode(node);

         scene_builder_->addMeshInstance(sphere_node_id_, sphere_mesh_id);

         /*scene_->setMeshVertices(
             tri_id,
             {
                 {"positions", b_pos},
                 {"normals", b_normal},
                 {"tangents", b_tang},
                 {"texcrds", b_uv},
             }
         );*/

         //scene_ = scene_builder_->getScene();
         CreateRaytracingProgram(render_context);
     }

    app_gui_window->var("densityGraphicsMultiplier",
        Metrics::densityGraphicsMultiplier);
    app_gui_window->var("volumeValueOffset", volumeValueOffset, 0.f, 100.f);
    app_gui_window->var("DensityRayMarchMultiplier", DensityRayMarchMultiplier);
    app_gui_window->checkbox("Draw Shadow ?", drawShadow);
    if (drawShadow)
    {
        app_gui_window->var("shadowDensityMultiplier", shadowDensityMultiplier);
    }
    
    //app_gui_window->var("SphereRadius", SphereRadius, 0.f, 200.f);

    app_gui_window->var("ISO Level", IsoLevel);
    app_gui_window->var("normalOffset", normalOffset);

    app_gui_window->checkbox("useVoxelOpti", useVoxelOpti);
    app_gui_window->checkbox("debugVoxelGrid", debugVoxelGrid);

    app_gui_window->checkbox("Draw Fluid ?", draw_fluid_);
    app_gui_window->checkbox("Use Transformations ?", useTransformations);
    app_gui_window->checkbox("Use Recursive Raytracing ?", useRecursiveRaytracing);
    app_gui_window->checkbox("Use Debug Normals ?", debugNormals);
    app_gui_window->checkbox("Light Scattering ?", lightScattering);
    

    const auto oldDensityMapSize = Metrics::density_map_size;

    app_gui_window->var("Density Texture Size", Metrics::density_map_size, 1, 500);

    if (oldDensityMapSize != Metrics::density_map_size)
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
    }

    app_gui_window->slider("translation", translation, -500.f, 500.f);
    app_gui_window->slider("rotation", rotation, 0.f, 360.f);
    app_gui_window->slider("scale", scale, 1.f, 500.f);

    app_gui_window->slider("march_mesh_scale", march_mesh_scale, 0.1f, 100.f);

    fluid_transform.setTranslation(translation);
    fluid_transform.setRotationEulerDeg(rotation);
    fluid_transform.setScaling(scale);

    //scene_->updateNodeTransform(raymarching_node_id.get(), fluid_transform.getMatrix());

    //if (draw_fluid_)
    //{
        //app_gui_window->var("Water Turbulance", water_turbulence_);

        app_gui_window->var("MaxRayBounce", kMaxRayBounce);

        //app_gui_window->var("absorptionCoeff", absorptionCoeff);
        app_gui_window->var("scatteringCoeff", scatteringCoeff);
        //app_gui_window->var("Phase G ", phaseG);

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

void Renderer::CreateRaytracingProgram(RenderContext* render_context) noexcept
{
    //Init(render_context);

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
    rtProgDesc.setCompilerFlags(SlangCompilerFlags::GenerateDebugInfo); // Ajoute le debug Slang uniquement en Debug
    rtProgDesc.addTypeConformances(typeConformances);
    rtProgDesc.setMaxTraceRecursionDepth(6);
    rtProgDesc.setMaxPayloadSize(36); // The largest ray payload struct (PrimaryRayData) is 24 bytes. The payload size
                                      // should be set as small as possible for maximum performance.

    const ref<RtBindingTable> sbt = RtBindingTable::create(2, 2, scene_->getGeometryCount());
    sbt->setRayGen(rtProgDesc.addRayGen("rayGen"));
    sbt->setMiss(0, rtProgDesc.addMiss("miss"));
    sbt->setMiss(1, rtProgDesc.addMiss("shadowMiss"));

    const auto primary = rtProgDesc.addHitGroup("closestHit", "");
    sbt->setHitGroup(0, scene_->getGeometryIDs(Scene::GeometryType::TriangleMesh), primary);

    const auto shadow = rtProgDesc.addHitGroup("", "shadowAnyHit");
    sbt->setHitGroup(1, scene_->getGeometryIDs(Scene::GeometryType::TriangleMesh), shadow);

    const auto raymarching_hit_group = rtProgDesc.addHitGroup(
        "RaymarchingClosestHit", "", "RaymarchingIntersection");
    sbt->setHitGroup(0, scene_->getGeometryIDs(Scene::GeometryType::Custom), raymarching_hit_group);

    const auto shadowFluid = rtProgDesc.addHitGroup("", "shadowAnyHitFluid", "AABB_SimulationIntersection");
    sbt->setHitGroup(1, scene_->getGeometryIDs(Scene::GeometryType::Custom), shadowFluid);

    rt_program_ = Program::create(device_, rtProgDesc, defines);
    rt_program_vars_ = RtProgramVars::create(device_, rt_program_, sbt);
}

void Renderer::LaunchMarchingCubeComputePasses(RenderContext* render_context) noexcept
{
    /*const auto compute_var2 = compute_marching_cube_density_map_->getRootVar();
    compute_var2["DensityTexture"] = marching_cube_dens_tex;
    compute_var2["triangles"] = marching_cubes_triangle_buffer_;

    compute_var2["PerFrameCB"]["numPointsPerAxis"] = numPointsPerAxis;
    compute_var2["PerFrameCB"]["isoLevel"] = IsoLevel;
    compute_var2["PerFrameCB"]["textureSize"] = Metrics::density_map_size;
    compute_var2["PerFrameCB"]["boundSize"] = Metrics::sim_bounds;
    compute_var2["PerFrameCB"]["SphereRadius"] = SphereRadius;
    compute_var2["PerFrameCB"]["var"] = var_;

    compute_marching_cube_density_map_->execute(
        render_context, Metrics::density_map_size, Metrics::density_map_size, Metrics::density_map_size
    );*/

    {
        FALCOR_PROFILE(render_context, "Compute Marching Cube Pass");

        render_context->clearUAVCounter(marching_cubes_triangle_buffer_, 0);

        const auto compute_var = marching_cubes_pass_->getRootVar();
        compute_var["DensityTexture"] = density_3d_tex_;
        compute_var["triangles"] = marching_cubes_triangle_buffer_;

        compute_var["PerFrameCB"]["numPointsPerAxis"] = numPointsPerAxis;
        compute_var["PerFrameCB"]["isoLevel"] = IsoLevel;
        compute_var["PerFrameCB"]["textureSize"] = Metrics::density_map_size;
        compute_var["PerFrameCB"]["boundSize"] = Metrics::sim_bounds;
        compute_var["PerFrameCB"]["SphereRadius"] = SphereRadius;
        compute_var["PerFrameCB"]["var"] = var_;

        compute_var["PerFrameCB"]["simBounds"] = float3(Metrics::sim_bounds);
        compute_var["PerFrameCB"]["normalOffset"] = normalOffset;
        compute_var["PerFrameCB"]["volumeValueOffset"] = volumeValueOffset;
        compute_var["PerFrameCB"]["scale"] = march_mesh_scale;

        compute_var["linearClampSampler"] = linearClampSampler_;

        render_context->clearUAVCounter(marching_cubes_triangle_buffer_, 0);

        int numVoxelsPerX = Metrics::density_map_size - 1;
        int numVoxelsPerY = Metrics::density_map_size - 1;
        int numVoxelsPerZ = Metrics::density_map_size - 1;
        marching_cubes_pass_->execute(render_context, numVoxelsPerX, numVoxelsPerY, numVoxelsPerZ);
    }

    render_context->uavBarrier(marching_cubes_triangle_buffer_.get());
    triangleCount = marching_cubes_triangle_buffer_->getUAVCounter()->getElement<uint>(0);
   
    uint32_t vertexCount = triangleCount * 3;
    //size_t maxVertexCount = static_cast<size_t>(MaxTriangleCount * 3);

    /*render_context->copyBufferRegion(
        read_back_triangle_buffer_.get(), 0,
        marching_cubes_triangle_buffer_.get(), 0,
        triangleCount * sizeof(MarchingCubesTriangle)
    );*/

    //if (triangleCount != oldTriangleCount)
    {
        // std::cout << "New triangle count: " << triangleCount << '\n';

        {
            FALCOR_PROFILE(render_context, "Copy Triangle Data on CPU");

            render_context->copyResource(read_back_triangle_buffer_.get(), marching_cubes_triangle_buffer_.get());
        }

        {
            FALCOR_PROFILE(render_context, "Read Triangle Data on CPU");

            const MarchingCubesTriangle* triangles = static_cast<const MarchingCubesTriangle*>(read_back_triangle_buffer_->map());

            std::vector<float3> new_pos;
            new_pos.resize(MaxVertexCount, float3(100, 100, 100));
            std::vector<float3> new_normals;
            new_normals.resize(MaxVertexCount, float3(100, 100, 100));
            for (uint32_t i = 0; i < triangleCount; ++i)
            {
                new_pos[i * 3 + 0] = triangles[i].vertexA.position;
                new_pos[i * 3 + 1] = triangles[i].vertexB.position;
                new_pos[i * 3 + 2] = triangles[i].vertexC.position;

                new_normals[i * 3 + 0] = triangles[i].vertexA.normal;
                new_normals[i * 3 + 1] = triangles[i].vertexB.normal;
                new_normals[i * 3 + 2] = triangles[i].vertexC.normal;
            }

            // for (uint32_t i = vertexCount; i < MaxVertexCount; ++i)
            //     new_pos[i] = float3(0, 0, 0);

            b_pos->setBlob(new_pos.data(), 0, MaxVertexCount * sizeof(float3));
            b_normal->setBlob(new_normals.data(), 0, MaxVertexCount * sizeof(float3));

            read_back_triangle_buffer_->unmap();
        }
    }
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
    var["PerFrameCB"]["shadowDensityMultiplier"] = shadowDensityMultiplier;
    var["PerFrameCB"]["drawShadow"] = drawShadow;

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
    var["cutomPrimitveMasks"] = cutomPrimitveMasks;

    var["PerFrameCB"]["useVoxelOpti"] = useVoxelOpti;
    var["PerFrameCB"]["debugVoxelGrid"] = debugVoxelGrid;
    var["PerFrameCB"]["voxelGridTotalResolution"] = Metrics::voxelGridTotalResolution;

    const float4x4 localToWorld = fluid_transform.getMatrix();

    //float3 renderScale = float3(0.5f); // render at half-size
    //const float4x4 scaledLocalToWorld = math::scale(localToWorld, renderScale);

    const float4x4 worldToLocal = inverse(localToWorld);

    var["PerFrameCB"]["localToWorld"] = localToWorld;
    var["PerFrameCB"]["worldToLocal"] = worldToLocal;

    var["PerFrameCB"]["ScaledSimBounds"] = fluid_transform.getScaling() * 2.f;
    var["PerFrameCB"]["useTransformations"] = useTransformations;
    var["PerFrameCB"]["useRecursiveRaytracing"] = useRecursiveRaytracing;
    var["PerFrameCB"]["debugNormals"] = debugNormals;
}

void Renderer::CreateRasterizationProgram() noexcept
{
    // Create raster pass.
    // This utility wraps the creation of the program and vars, and sets the necessary scene defines.
    ProgramDesc rasterProgDesc;
    rasterProgDesc.addShaderModules(scene_->getShaderModules());
    rasterProgDesc.addShaderLibrary("Samples/Raytracing/HelloDXR.3d.slang").vsEntry("vsMain").psEntry("psMain");
    rasterProgDesc.addTypeConformances(scene_->getTypeConformances());

    raster_pass_ = RasterPass::create(device_, rasterProgDesc, scene_->getSceneDefines());

    //// Create the RenderState
    //raster_pass_ = RasterPass::create(device_,
    //        "Samples/Raytracing/triangle.slang", "vsMain", "psMain");

    auto& pState = raster_pass_->getState();
    
    // create the depth-state
    DepthStencilState::Desc dsDesc;
    dsDesc.setDepthEnabled(false);
    pState->setDepthStencilState(DepthStencilState::create(dsDesc));
    
    // Rasterizer state
    RasterizerState::Desc rsState;
    rsState.setCullMode(RasterizerState::CullMode::None);
    pState->setRasterizerState(RasterizerState::create(rsState));

     // Blend state
     BlendState::Desc blendDesc;
     blendDesc.setRtBlend(0, true).setRtParams(
        0,
        BlendState::BlendOp::Add,
        BlendState::BlendOp::Add,
        BlendState::BlendFunc::SrcAlpha,
        BlendState::BlendFunc::OneMinusSrcAlpha,
        BlendState::BlendFunc::One,
        BlendState::BlendFunc::One
    );
     pState->setBlendState(BlendState::create(blendDesc));
}
