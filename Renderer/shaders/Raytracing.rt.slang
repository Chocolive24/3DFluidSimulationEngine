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
//#pragma add - search - path "Samples/Raytracing" // or wherever your SDF module is
import SDF_Functions;

import Scene.Raytracing;
import Utils.Sampling.TinyUniformSampleGenerator;
import Rendering.Materials.Fresnel;
import Rendering.Materials.TexLODHelpers;
import Rendering.Lights.LightHelpers;

RWTexture2D<float4> gOutput;
Texture3D<float> gTexture3D;
SamplerState linearClampSampler;

cbuffer PerFrameCB
{
    float4x4 invView;
    float2 viewportDims;
    float tanHalfFovY;
    uint sampleIndex;
    bool useDOF;

    bool drawFluid;
    bool lightScattering;
    float3 backgroundColor;
    float waterTurbulence;
    uint maxRayBounce;
    float3 absorptionCoeff;
    float3 scatteringCoeff;
    float phaseG;

    float isoLevel;
    float normalOffset;
    float maxRaymarchingDistance;
    float marchSize;
    float maxLighMarchingDistance;
    float sunLightMarchSize;

    float3 lightColor;
    float3 lightDir;
    float IoR;
    float time;
    int iFrame;

    float DensityRayMarchMultiplier;
    int densityMapSize;
    float3 simBounds;
    float3 ScaledSimBounds;
    float volumeValueOffset;

    float4x4 worldToLocal;
    float4x4 localToWorld;

    bool useTransformations;
    bool useRecursiveRaytracing;
    bool debugNormals;
};



struct SDF_GeometryAttributes
{
    float3 pos = float3(0, 0, 0);
    float3 normal = float3(0, 0, 0);
    float densityAlongRay = 0.f;
    bool isInsideFluid = false;
}

struct PrimaryRayData
{
    float4 color = float4(0, 0, 0, 0);
    float3 transmittance = float3(1, 1, 1);
    float3 debugNormal;
    float2 uv = float2(0, 0);
    uint depth = 0;
    float hitT = -1.f;
    float absorptionDistance = 0.f;
};

struct CubeHitInfo
{
    bool didHit;
    bool isInside;
    float dst;
    float3 hitPoint;
    float3 normal;
}

struct SurfaceInfo
{
    float3 pos;
    float3 normal;
    float densityAlongRay;
    bool foundSurface;
};

struct LightResponse
{
    float3 reflectDir;
    float3 refractDir;
    float reflectWeight;
    float refractWeight;
};

static const float TinyNudge = 0.01;

// PCG (permuted congruential generator). Thanks to:
// www.pcg-random.org and www.shadertoy.com/view/XlGcRh
uint NextRandom(inout uint state)
{
    state = state * 747796405 + 2891336453;
    uint result = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
    result = (result >> 22) ^ result;
    return result;
}

float RandomValue(inout uint state)
{
    return NextRandom(state) / 4294967295.0; // 2^32 - 1
}

/**
 * Returns the shading normal flipped for backfacing hits on double-sided materials.
 * This reproduces the legacy behavior before the conventions were updated.
 */
float3 getFlippedShadingNormal(const ShadingData sd)
{
    return !sd.frontFacing && sd.mtl.isDoubleSided() ? -sd.frame.N : sd.frame.N;
}

/**
 * Helper for computing relatice index of refraction (eta) at a hit.
 * This is computed based on the IoRs of the media outside and inside the material.
 */
float computeEta(const ShadingData sd)
{
    float insideIoR = gScene.materials.evalIoR(sd.materialID);
    float eta = sd.frontFacing ? (sd.IoR / insideIoR) : (insideIoR / sd.IoR);
    return eta;
}


float3 BeerLambert(float3 absorptionCoefficient, float distanceTraveled)
{
    return exp(-absorptionCoefficient * distanceTraveled);
}

float PhaseHG(float cosTheta, float g) {
    float denom = 1.0 + g * g - 2.0 * g * cosTheta;
    return (1.0 - g * g) / (4.0 * 3.1415926538 * pow(denom, 1.5));
}

float3 CalculateTransmittedSunLightAtPos(float3 position)
{
    float3 transmittedSunLight = lightColor;
    for (float sun_t = 0.f; sun_t < maxLighMarchingDistance; sun_t += sunLightMarchSize)
    {
        float3 posAlongSunRay = position + (-lightDir) * sun_t;
        bool isInVolume = map(posAlongSunRay, waterTurbulence, time, iFrame) < 0.001f;

        if (isInVolume)
        {
            transmittedSunLight *= exp(-absorptionCoeff * sunLightMarchSize);
        }
    }

    return transmittedSunLight;
}

float3 CalculateTransmittanceRaymarching(SDF_GeometryAttributes attribs, float3 rayDir)
{
    float3 totalLight = float3(0.f, 0.f, 0.f);
    float3 transmittance = float3(1.f, 1.f, 1.f);

    for (float t = 0.f; t < maxRaymarchingDistance; t += marchSize)
    {
        float3 position = attribs.pos + rayDir * t;
        bool isInVolume = map(position, waterTurbulence, time, iFrame) < 0.001f;

        if (isInVolume)
        {
            float3 transmittedSunLight = CalculateTransmittedSunLightAtPos(position);

            float cosTheta = dot(lightDir, -rayDir);
            float phase = PhaseHG(cosTheta, phaseG);
            float3 scattered = transmittedSunLight * phase * scatteringCoeff;

            totalLight += scattered * transmittance * marchSize;

            transmittance *= exp(-absorptionCoeff * marchSize);

            if (max(transmittance.r, max(transmittance.g, transmittance.b)) < 0.01f)
                break;
        }
    }

    return totalLight;
}

// Miss shader
[shader("miss")]
void miss(inout PrimaryRayData hitData)
{
    // Use the ray direction to sample the environment map
    float3 envColor = gScene.envMap.eval(WorldRayDirection());
    // hitData.color = hasHit ? float4(0.1, 0.3, 0.8, 1.0) : float4(envColor, 1.f);
    hitData.color = float4(envColor, 1.f);
    hitData.hitT = -1.f;
    hitData.absorptionDistance = 0;
}

[shader("closesthit")]
void closestHit(inout PrimaryRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
{
    // // Get the hit-point data.
    float3 rayDir = WorldRayDirection();
    float currentRayDistance = RayTCurrent();
    uint triangleIndex = PrimitiveIndex();

    // if (!drawFluid)
    // {
    //     hitData.color = float4(0.1f, 0.2f, 1.f, 1.f);
    //     hitData.hitT = currentRayDistance;
    //     hitData.absorptionDistance = 0;
    //     return;
    // }

    // if (hitData.depth >= maxRayBounce)
    // {
    //     float3 envColor = gScene.envMap.eval(WorldRayDirection(), 0);
    //     hitData.color = float4(envColor, 1);
    //     hitData.hitT = RayTCurrent();
    //     hitData.absorptionDistance = hitData.hitT;
    //     return;
    // }

    // Retrieve hit information
    const GeometryInstanceID instanceID = getGeometryInstanceID();
    VertexData v = getVertexData(instanceID, triangleIndex, attribs);
    float3 world_pos = v.posW;
    float3 normal = v.normalW;

    uint materialID = gScene.getMaterialID(instanceID);
    ShadingData sd = gScene.materials.prepareShadingData(v, materialID, -rayDir);
    float3 orientedFaceNormal = getFlippedShadingNormal(sd);

    // Create material instance and query its properties.
    let lod = ExplicitLodTextureSampler(0.f);
    let mi = gScene.materials.getMaterialInstance(sd, lod);
    let bsdfProperties = mi.getProperties(sd);

    if (gScene.materials.evalIoR(sd.materialID) < 1.f)
    {
        // Étape 1 : Récupération des UVs
        float2 uv = sd.uv;

        BasicMaterialData matData = gScene.materials.getBasicMaterialData(materialID);

        TextureHandle textureHandle = matData.texBaseColor;
        float4 baseColorFactor = matData.baseColor;

        float4 baseColor = gScene.materials.sampleTexture(textureHandle, linearClampSampler, sd.uv, baseColorFactor, lod);

        hitData.color = float4(baseColor.rgb, 1.0);
        
        hitData.color = baseColor;
        hitData.hitT = RayTCurrent();
        hitData.absorptionDistance = 0;
        return;
    }

    // 1. Compute Fresnel reflectance
    float etaI = sd.IoR; // Incident IoR (Air by default).
    float etaT = gScene.materials.evalIoR(sd.materialID); // Tranmission IoR (the one from the medium).
    float eta = sd.frontFacing ? etaI / etaT : etaT / etaI;

    float F0 = pow((etaI - etaT) / (etaI + etaT), 2.0);
    float F90 = 1.f;
    float cosTheta = max(0.00001f, dot(getFlippedShadingNormal(sd), -rayDir));
    float reflectionRatio = evalFresnelSchlick(F0, F90, cosTheta);
    float refractionRatio = 1.f - reflectionRatio;

    // 2. Prepare color storage
    float4 reflectedColor = float4(0, 0, 0, 0);
    float4 refractedColor = float4(0, 0, 0, 0);
    float absorptionTotalDistance = 0.f;

    // 3. Check for total internal reflection
    bool tir = false;
    float3 refractDir = rayDir; // Default if eta == 1.0

    if (eta != 1.0)
    {
        if (!refractWithTIR(rayDir, orientedFaceNormal, eta, refractDir))
            tir = true;
    }

    // 4. Trace reflection ray if needed
    if (reflectionRatio > 0.001f || tir)
    {
        // Handle total internal reflection
        float3 reflectDir = reflect(rayDir, orientedFaceNormal); // Compute the reflection direction

        // Define reflection ray
        RayDesc reflectRay;
        reflectRay.Origin = sd.computeRayOrigin(true); // Offset to avoid self-intersection
        reflectRay.Direction = reflectDir;
        reflectRay.TMin = 0.001f;
        reflectRay.TMax = 100000.f;

        // Trace reflection ray
        PrimaryRayData reflectedRayData = {};
        reflectedRayData.depth = hitData.depth + 1;
        reflectedRayData.hitT = 0;
        reflectedRayData.absorptionDistance = hitData.absorptionDistance;
        TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, reflectRay, reflectedRayData);

        reflectedColor = reflectedRayData.color;
        absorptionTotalDistance = reflectedRayData.absorptionDistance;

        // If we are in a TiR case, the next reflection ray will be inside the medium. 
        // That's why we add its distance to the total absorption distance.
        if (tir)
        {
            absorptionTotalDistance += reflectedRayData.hitT;
        }
    }

    // 5. Trace refraction ray if not TIR
    if (!tir && refractionRatio > 0.001f)
    {
        // Trace a refracted ray to get the color after the transmission.
        RayDesc refractRay;
        refractRay.Origin = sd.computeRayOrigin(false);
        refractRay.Direction = refractDir;
        refractRay.TMin = 0.001f;
        refractRay.TMax = 100000.f;

        PrimaryRayData refractedRayData = {};
        refractedRayData.depth = hitData.depth + 1;
        refractedRayData.hitT = 0.f;
        refractedRayData.absorptionDistance = hitData.absorptionDistance;
        TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, refractRay, refractedRayData);

        refractedColor = refractedRayData.color;
        absorptionTotalDistance = refractedRayData.absorptionDistance;
    }

    // If the current ray is inside the medium, add its distance to the total absorb distance.
    if (!sd.frontFacing)
    {
        absorptionTotalDistance += currentRayDistance;
    }

    hitData.absorptionDistance = absorptionTotalDistance;

    // Apply Beer's Law attenuation
    float3 absorb = exp(-absorptionCoeff * absorptionTotalDistance);
    refractedColor *= float4(absorb, 1.f);

    // 6. Mix results
    float4 finalColor = tir ? reflectedColor :
            reflectionRatio * reflectedColor + refractionRatio * refractedColor;

    hitData.color = finalColor;
}

SurfaceInfo FindNextSurface(float3 origin, float3 rayDir, bool findNextFluidEntryPoint, 
    uint rngState, float rngWeight, float maxDst)
{
    SurfaceInfo info = (SurfaceInfo)0;
    // Safety check
    if (dot(rayDir, rayDir) < 0.5) return info;

    //float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, origin, rayDir);
    //  Use your new matrix-aware ray-box distance function
    float2 boundsDstInfo = RayBoxDst_Matrix(origin, rayDir, localToWorld, worldToLocal);
    //float r = (RandomValue(rngState) - 0.5) * marchSize * 0.4 * 1;
    float r = (1 - 0.5) * marchSize * 0.4 * 1;
    bool hasExittedFluid = !IsInsideFluid(origin, rayDir);
    origin = origin + rayDir * (boundsDstInfo.x + r);

    float stepSize = marchSize;
    bool hasEnteredFluid = false;
    float3 lastPosInFluid = origin;

    float dstToTest = boundsDstInfo[1] - (TinyNudge) * 2;

    for (float dst = 0; dst < dstToTest; dst += stepSize)
    {
        bool isLastStep = dst + stepSize >= dstToTest;
        float3 samplePos = origin + rayDir * dst;
        float thickness = SampleDensityMap(samplePos) * DensityRayMarchMultiplier * stepSize;
        bool insideFluid = thickness > 0;
        if (insideFluid)
        {
            hasEnteredFluid = true;
            lastPosInFluid = samplePos;
            if (dst <= maxDst)
            {
                info.densityAlongRay += thickness;
            }
        }

        if (!insideFluid) hasExittedFluid = true;

        bool found;
        if (findNextFluidEntryPoint) found = insideFluid && hasExittedFluid;
        else found = hasEnteredFluid && (!insideFluid || isLastStep);

        if (found)
        {
            info.pos = lastPosInFluid;
            info.foundSurface = true;
            break;
        }
    }

    return info;
}

// SurfaceInfo FindNextSurface(float3 origin, float3 rayDir, bool findNextFluidEntryPoint, float maxDst)
// {
//     SurfaceInfo info;
//     info.pos = float3(0, 0, 0);
//     info.densityAlongRay = 0;
//     info.foundSurface = false;

//     // Safety check
//     if (dot(rayDir, rayDir) < 0.5) return info;

//     // Compute AABB intersection
//     float3 aabbMin = -0.5f * simBounds;
//     float3 aabbMax = 0.5f * simBounds;
//     float2 bounds = RayBoxDst(aabbMin, aabbMax, origin, rayDir);

//     float entryT = max(bounds.x, 0.0);
//     float travelT = bounds.y;
//     float dstToTest = travelT - 1e-4;

//     origin = origin + rayDir * entryT;

//     float stepSize = marchSize;

//     bool hasEnteredFluid = false;
//     bool hasExitedFluid = false;
//     float3 lastPosInFluid = origin;

//     for (float dst = 0; dst < dstToTest; dst += stepSize)
//     {
//         float3 samplePos = origin + rayDir * dst;
//         float density = SampleDensityMap(samplePos) * DensityRayMarchMultiplier * stepSize;
//         bool insideFluid = (density > 0);

//         if (insideFluid)
//         {
//             hasEnteredFluid = true;
//             lastPosInFluid = samplePos;

//             if (dst <= maxDst)
//             {
//                 info.densityAlongRay += density;
//             }
//         }
//         else
//         {
//             hasExitedFluid = true;
//         }

//         bool found = false;
//         if (findNextFluidEntryPoint)
//         {
//             found = insideFluid && hasExitedFluid;
//         }
//         else
//         {
//             bool isLastStep = (dst + stepSize >= dstToTest);
//             found = hasEnteredFluid && (!insideFluid || isLastStep);
//         }

//         if (found)
//         {
//             info.pos = lastPosInFluid;
//             info.foundSurface = true;
//             break;
//         }
//     }

//     return info;
// }

// bool IsInsideFluid(float3 worldPos)
// {
//     // Use a fake ray pointing in any direction — here +Z (doesn't matter much)
//     float3 worldRayDir = float3(0, 0, 1);

//     // Call RayBoxWithMatrix which internally transforms to local space
//     CubeHitInfo hit = RayBoxWithMatrix(worldPos, worldRayDir, localToWorld, worldToLocal);

//     // Check if the ray starts inside the box
//     bool insideBox = hit.didHit && hit.isInside;

//     // Optional: only consider "inside" if density is nonzero
//     return insideBox && SampleDensityMap(worldPos) > 0;
// }

bool IsInsideFluid(float3 pos, float3 rayDir)
{
    //float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, pos, float3(0, 0, 1));
    float2 boundsDstInfo = RayBoxDst_Matrix(pos, rayDir, localToWorld, worldToLocal);

    return (boundsDstInfo.x <= 0 && boundsDstInfo.y > 0) && SampleDensityMap(pos) > 0;
}

// [shader("closesthit")]
// void RaymarchingClosestHit(inout PrimaryRayData payload, SDF_GeometryAttributes attribs)
// {
//     Ray ray;
//     ray.dir = WorldRayDirection();

//     float densityStepSize = sunLightMarchSize * (payload.depth + 1); // increase step size with each iteration

//     // === Raymarch once ===
//     payload.transmittance *= exp(-attribs.densityAlongRay * scatteringCoeff);

//     // === Surface interaction ===
//     float3 normal = ComputeDensityNormal(attribs.pos);
//     bool localTravel = dot(normal, ray.dir) > 0;
//     if (localTravel) normal = -normal;

//     // Compute IoR ratio based on fluid state
//     float eta = localTravel ? (IoR / 1.0) : (1.0 / IoR);

//     // Fresnel reflectance using Schlick approximation
//     float cosTheta = abs(dot(normalize(ray.dir), normal));
//     float F0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
//     float reflectionRatio = evalFresnelSchlick(F0, 1.0, cosTheta);
//     float refractionRatio = 1.0 - reflectionRatio;

//     // Compute directions
//     float3 reflectDir = reflect(ray.dir, normal);
//     float3 refractDir = refract(ray.dir, normal, eta);

//     float densityRefrac = CalculateDensityAlongRay(attribs.pos, refractDir, densityStepSize);
//     float densityReflect = CalculateDensityAlongRay(attribs.pos, reflectDir, densityStepSize);

//     bool followRefrac = densityRefrac * refractionRatio > densityReflect * reflectionRatio;

//     // Accumulate light from the less dominant path
//     if (followRefrac)
//     {
//         payload.light += gScene.envMap.eval(reflectDir) * payload.transmittance * exp(-densityReflect * scatteringCoeff) * reflectionRatio;
//     }
//     else
//     {
//         payload.light += gScene.envMap.eval(refractDir) * payload.transmittance * exp(-densityRefrac * scatteringCoeff) * refractionRatio;
//     }

//     // Advance ray
//     float3 nextOrigin = attribs.pos + ray.dir * 1e-3f; // Avoid self-hit
//     float3 nextDir = followRefrac ? refractDir : reflectDir;

//     // Update transmittance
//     payload.transmittance *= followRefrac ? refractionRatio : reflectionRatio;

//     // === Spawn new recursive ray ===
//     if (payload.depth < maxRayBounce)
//     {
//         PrimaryRayData newPayload;
//         newPayload.transmittance = payload.transmittance;
//         newPayload.light = float3(0, 0, 0);
//         newPayload.depth = payload.depth + 1;

//         RayDesc newRay;
//         newRay.Origin = nextOrigin;
//         newRay.Direction = nextDir;
//         newRay.TMin = 0.001;
//         newRay.TMax = 1e+4;

//         TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, 1, 0, newRay, newPayload);

//         //payload.light += newPayload.light; // Recover new light data for final shading.
//         payload.transmittance = newPayload.transmittance;
//         //payload.light = newPayload.light; // Recover new light data for final shading.
//         payload.color = newPayload.color;
//     }

//     if (payload.depth >= maxRayBounce)
//     {
//         // Max depth: sample env map
//         // float densityRemainder = CalculateDensityAlongRay(nextOrigin, nextDir, sunLightMarchSize);
//         // payload.light += gScene.envMap.eval(nextDir) * payload.transmittance * exp(-densityRemainder * scatteringCoeff);
//         //payload.light += gScene.envMap.eval(WorldRayDirection());
//         payload.color.rgb += gScene.envMap.eval(WorldRayDirection());
//     }

//     //payload.color = float4(payload.light, 1.0);
//     payload.hitT = RayTCurrent();
// }

// static bool travellingThroughFluid = false;
// static bool firstIntersection = true;

float3 Refract(float3 inDir, float3 normal, float iorA, float iorB)
{
    float refractRatio = iorA / iorB;
    float cosAngleIn = -dot(inDir, normal);
    float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
    if (sinSqrAngleOfRefraction > 1) return 0; // Ray is fully reflected, no refraction occurs

    float3 refractDir = refractRatio * inDir + (refractRatio * cosAngleIn - sqrt(1 - sinSqrAngleOfRefraction)) * normal;
    return refractDir;
}

// Calculate the proportion of light that is reflected at the boundary between two media (via the fresnel equations)
// Note: the amount of light refracted can be calculated as 1 minus this value
float CalculateReflectance(float3 inDir, float3 normal, float iorA, float iorB)
{
    float refractRatio = iorA / iorB;
    float cosAngleIn = -dot(inDir, normal);
    float sinSqrAngleOfRefraction = refractRatio * refractRatio * (1 - cosAngleIn * cosAngleIn);
    if (sinSqrAngleOfRefraction >= 1) return 1; // Ray is fully reflected, no refraction occurs

    float cosAngleOfRefraction = sqrt(1 - sinSqrAngleOfRefraction);
    // Perpendicular polarization
    float rPerpendicular = (iorA * cosAngleIn - iorB * cosAngleOfRefraction) / (iorA * cosAngleIn + iorB * cosAngleOfRefraction);
    rPerpendicular *= rPerpendicular;
    // Parallel polarization
    float rParallel = (iorB * cosAngleIn - iorA * cosAngleOfRefraction) / (iorB * cosAngleIn + iorA * cosAngleOfRefraction);
    rParallel *= rParallel;

    // Return the average of the perpendicular and parallel polarizations
    return (rPerpendicular + rParallel) / 2;
}

LightResponse CalculateReflectionAndRefraction(float3 inDir, float3 normal, float iorA, float iorB)
{
    LightResponse result;

    result.reflectWeight = CalculateReflectance(inDir, normal, iorA, iorB);
    result.refractWeight = 1 - result.reflectWeight;

    result.reflectDir = reflect(inDir, normal);
    result.refractDir = Refract(inDir, normal, iorA, iorB);

    return result;
}

float3 Transmittance(float thickness)
{
    return exp(-thickness * scatteringCoeff);
}

[shader("closesthit")]
void RaymarchingClosestHit(inout PrimaryRayData hitData, SDF_GeometryAttributes attribs)
{
    float3 debugNormal;

    if (useRecursiveRaytracing)
    {
        float3 rayOrigin = attribs.pos;
        float3 normal = attribs.normal;

        float3 rayDir = normalize(WorldRayDirection());
        float3 light = float3(0, 0, 0);

        // Decide whether we're entering or exiting
        bool travellingThroughFluid = attribs.isInsideFluid;

        // Apply absorption from entry to this surface
        hitData.transmittance *= Transmittance(attribs.densityAlongRay);

        if (dot(normal, rayDir) > 0) normal = -normal;

        if (debugNormals) hitData.debugNormal = normal;

        float iorA = travellingThroughFluid ? IoR : 1.0f;
        float iorB = travellingThroughFluid ? 1.0f : IoR;

        LightResponse lightResponse = CalculateReflectionAndRefraction(rayDir, normal, iorA, iorB);

        float densityStepSize = sunLightMarchSize * (hitData.depth + 1); // increase step size with each iteration

        float densityAlongRefractRay = CalculateDensityAlongRay(attribs.pos, lightResponse.refractDir, densityStepSize);
        float densityAlongReflectRay = CalculateDensityAlongRay(attribs.pos, lightResponse.reflectDir, densityStepSize);

        bool traceRefractedRay = densityAlongRefractRay * lightResponse.refractWeight > densityAlongReflectRay * lightResponse.reflectWeight;

        // Accumulate environment lighting from less dominant path
        if (traceRefractedRay)
        {
            light += gScene.envMap.eval(lightResponse.reflectDir)
              * hitData.transmittance
              * Transmittance(densityAlongReflectRay)
              * lightResponse.reflectWeight;
        }
        else
        {
            light += gScene.envMap.eval(lightResponse.refractDir)
              * hitData.transmittance
              * Transmittance(densityAlongRefractRay)
              * lightResponse.refractWeight;
        }

        // Prepare bounce ray
        PrimaryRayData bounceResult = (PrimaryRayData)0;
        bounceResult.depth = hitData.depth + 1;

        // Update transmittance for main path (before bounce)
        float3 pathTransmittance = hitData.transmittance
                              * (traceRefractedRay ? lightResponse.refractWeight : lightResponse.reflectWeight);

        bounceResult.transmittance = pathTransmittance;
        bounceResult.absorptionDistance = attribs.densityAlongRay;

        float3 nextRayDir = traceRefractedRay ? lightResponse.refractDir : lightResponse.reflectDir;

        if (hitData.depth < maxRayBounce)
        {
            RayDesc bounceRay;
            bounceRay.Origin = rayOrigin;
            bounceRay.Direction = nextRayDir;
            bounceRay.TMin = 0.001f;
            bounceRay.TMax = 1e6;

            TraceRay(gScene.rtAccel, RAY_FLAG_NONE, ~0, 0, 1, 0, bounceRay, bounceResult);

            // Apply absorption for this segment after bounce
            float densityThisSegment = traceRefractedRay ? densityAlongRefractRay : densityAlongReflectRay;

            light += bounceResult.color.rgb
               * Transmittance(densityThisSegment); // Apply Beer's Law for the current segment

            debugNormal = bounceResult.debugNormal;
            hitData.transmittance = bounceResult.transmittance;
        }
        else
        {
            float densityRemainder = CalculateDensityAlongRay(rayOrigin, nextRayDir, sunLightMarchSize);
            light += gScene.envMap.eval(nextRayDir)
              * hitData.transmittance
              * Transmittance(densityRemainder);
        }

        debugNormal += float3(1);
        debugNormal *= 0.5f;

        float3 color = debugNormals ? debugNormal : light;
        hitData.color = float4(color, 1.0);
        hitData.hitT = RayTCurrent();
        return;
    }
    else 
    {
        uint rngState = (uint)(hitData.uv.x * 1243 + hitData.uv.y * 96456);

        Ray ray;
        ray.origin = WorldRayOrigin();
        ray.dir = normalize(WorldRayDirection());

        bool travellingThroughFluid = IsInsideFluid(ray.origin, ray.dir);

        float3 transmittance = float3(1, 1, 1);
        float3 light = float3(0, 0, 0);

        for (int i = 0; i < maxRayBounce; i++)
        {
            float densityStepSize = sunLightMarchSize * (i + 1); // increase step size with each iteration

            bool searchForNextFluidEntryPoint = !travellingThroughFluid;
            // bool searchForNextFluidEntryPoint = !IsInsideFluid(ray.origin, ray.dir);

            SurfaceInfo surfaceInfo = FindNextSurface(ray.origin, ray.dir, searchForNextFluidEntryPoint, rngState, i == 0 ? 1 : 0, 1000);
            
            if (!surfaceInfo.foundSurface) break;

            transmittance *= Transmittance(surfaceInfo.densityAlongRay);

            // If light hits the floor it will be scattered in all directions (in hemisphere)
            // Not sure how to handle this in real-time, so just break out of loop here
            if (surfaceInfo.pos.y < -ScaledSimBounds.y / 2 + 0.05)
            {
                break;
            }

            float3 normal = ComputeDensityNormal(surfaceInfo.pos);
            if (dot(normal, ray.dir) > 0) normal = -normal;

            if (debugNormals) debugNormal = normal;

            // Indicies of refraction
            float iorA = travellingThroughFluid ? IoR : 1.f;
            float iorB = travellingThroughFluid ? 1.f : IoR;

            // Calculate reflection and refraction, and choose which path to follow
            LightResponse lightResponse = CalculateReflectionAndRefraction(ray.dir, normal, iorA, iorB);

            float densityAlongRefractRay = CalculateDensityAlongRay(surfaceInfo.pos, lightResponse.refractDir, densityStepSize);
            float densityAlongReflectRay = CalculateDensityAlongRay(surfaceInfo.pos, lightResponse.reflectDir, densityStepSize);
            
            bool traceRefractedRay = densityAlongRefractRay * lightResponse.refractWeight > densityAlongReflectRay * lightResponse.reflectWeight;
            travellingThroughFluid = traceRefractedRay != travellingThroughFluid;

            // Accumulate light from the less dominant path
            if (traceRefractedRay)
            {
                light += gScene.envMap.eval(lightResponse.reflectDir) * transmittance * Transmittance(densityAlongReflectRay) * lightResponse.reflectWeight;
            }
            else
            {
                light += gScene.envMap.eval(lightResponse.refractDir) * transmittance * Transmittance(densityAlongRefractRay) * lightResponse.refractWeight;
            }

            // Set up ray for more interesting path
            ray.origin = surfaceInfo.pos;
            ray.dir = traceRefractedRay ? lightResponse.refractDir : lightResponse.reflectDir;
            transmittance *= (traceRefractedRay ? lightResponse.refractWeight : lightResponse.reflectWeight);
        }

        // Final raymarch if no more bounces
        float densityRemainder = CalculateDensityAlongRay(ray.origin, ray.dir, sunLightMarchSize);
        light += gScene.envMap.eval(ray.dir) * transmittance * Transmittance(densityRemainder);

        debugNormal += float3(1);
        debugNormal *= 0.5f;
        float3 color = debugNormals ? debugNormal : light;
        hitData.color = float4(color, 1.0);
        hitData.hitT = RayTCurrent();
        return;
    }
}

//     // float densityAlongViewRay = 0.f;
//     // densityAlongViewRay = CalculateDensityAlongRay(attribs.w_position, WorldRayDirection(), marchSize);

//     // float3 totalScatteredLight;
//     // for (float dst = 0.f; dst < attribs.tMax - attribs.tMin; dst += marchSize)
//     // {
//     //     float3 p = attribs.w_position + WorldRayDirection() * dst;
        
//     //     float densityAlongStep = SampleDensityMap(p) * DensityRayMarchMultiplier * marchSize;
//     //     densityAlongViewRay += densityAlongStep;

//     //     float3 transmittedSunLight;
//     //     if (lightScattering)
//     //     {
//     //         // TODO: Not sure if it works correctly
//     //         float densityAlongSunRay = CalculateDensityAlongRay(p, -lightDir, sunLightMarchSize);
//     //         transmittedSunLight = exp(-densityAlongSunRay * scatteringCoeff);
//     //     }
//     //     else
//     //     {
//     //         transmittedSunLight = float3(1, 1, 1);
//     //     }

//     //     float3 inScatteredLight = transmittedSunLight * densityAlongStep * scatteringCoeff;
//     //     float3 viewTransmittance = exp(-densityAlongViewRay * scatteringCoeff);
//     //     totalScatteredLight += inScatteredLight * viewTransmittance;
//     // }

//     // hitData.color = float4(densityAlongViewRay.xxx, 0);
//     // hitData.hitT = RayTCurrent();
//     // return;

//     // if (length(totalScatteredLight) <= isoLevel)
//     // {
//     //     totalScatteredLight = gScene.envMap.eval(WorldRayDirection());
//     // }

//     // hitData.color = float4(densityAlongViewRay.xxx, 0);
//     // hitData.hitT = RayTCurrent();
//     // return;

//     // float3 n = attribs.w_normal;
//     // n += float3(1, 1, 1);
//     // n *= 0.5f;
//     // hitData.color = float4(float3(n.x, n.y, n.z), 1);
//     // hitData.hitT = attribs.tMax - attribs.tMin;
//     // return;

//     // if (hitData.depth >= maxRayBounce)
//     // {
//     //     hitData.color = float4(backgroundColor, 0);
//     //     hitData.hitT = 0.f;
//     //     return;
//     // }

//     // float3 rayDir = WorldRayDirection();
//     // float currentRayT = RayTCurrent();

//     // float3 reflectDir = reflect(rayDir, attribs.w_normal); // Compute the reflection direction

//     // // Define reflection ray
//     // RayDesc reflectRay;
//     // reflectRay.Origin = attribs.w_position + reflectDir * 0.1f; // Offset to avoid self-intersection
//     // reflectRay.Direction = reflectDir;
//     // reflectRay.TMin = 0.001f;
//     // reflectRay.TMax = 100000.f;

//     // // Trace reflection ray
//     // PrimaryRayData reflectedRayData = {};
//     // reflectedRayData.depth = hitData.depth + 1;
//     // reflectedRayData.hitT = 0;
//     // TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, reflectRay, reflectedRayData);

//     // hitData.color = reflectedRayData.color;
//     // hitData.hitT = reflectedRayData.hitT;
//     // return;

//     // 1. Compute Fresnel reflectance
//     // float etaI = 1.f;                                // Incident IoR (Air by default).
//     // float etaT = IoR;                                // Tranmission IoR (the one from the medium).
//     // bool isEntering = dot(rayDir, attribs.w_normal) > 0.0;
//     // float eta = isEntering ? etaI / etaT : etaT / etaI;
//     // float3 flippedShadingNormal = isEntering ? attribs.w_normal : -attribs.w_normal;

//     // float F0 = pow((etaI - etaT) / (etaI + etaT), 2.0);
//     // float F90 = 1.f;
//     // float cosTheta = max(0.00001f, dot(flippedShadingNormal, rayDir));
//     // float reflectionRatio = evalFresnelSchlick(F0, F90, cosTheta);
//     // float refractionRatio = 1.f - reflectionRatio;

//     // // 2. Prepare color storage
//     // float3 reflectedColor = float3(0, 0, 0);
//     // float3 refractedColor = float3(0, 0, 0);

//     // // 3. Check for total internal reflection
//     // bool tir = false;
//     // float3 refractDir = rayDir; // Default if eta == 1.0

//     // if (eta != 1.0)
//     // {
//     //     if (!refractWithTIR(rayDir, flippedShadingNormal, eta, refractDir))
//     //         tir = true;
//     // }

//     // // 4. Trace reflection ray if needed
//     // if (reflectionRatio > 0.001f || tir)
//     // {
//     //     // Handle total internal reflection
//     //     // float3 reflectDir = reflect(rayDir, attribs.w_normal); // Compute the reflection direction

//     //     // // Define reflection ray
//     //     // RayDesc reflectRay;
//     //     // reflectRay.Origin = attribs.w_position + reflectDir * 0.1f; // Offset to avoid self-intersection
//     //     // reflectRay.Direction = reflectDir;
//     //     // reflectRay.TMin = 0.001f;
//     //     // reflectRay.TMax = 100000.f;

//     //     // // Trace reflection ray
//     //     // PrimaryRayData reflectedRayData = {};
//     //     // reflectedRayData.depth = hitData.depth + 1;
//     //     // reflectedRayData.hitT = 0;
//     //     // TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, reflectRay, reflectedRayData);

//     //     // reflectedColor = reflectedRayData.color.rgb;
//     // }

//     // // 5. Trace refraction ray if not TIR
//     // if (!tir && refractionRatio > 0.001f)
//     // {
//     //     refractedColor = CalculateTransmittanceRaymarching(attribs, refractDir);

//     //     // // Trace a refracted ray to get the color after the transmission.
//     //     // RayDesc refractRay;
//     //     // refractRay.Origin = attribs.w_position + refractDir * 0.001f;
//     //     // refractRay.Direction = refractDir;
//     //     // refractRay.TMin = 0.001f;
//     //     // refractRay.TMax = 100000.f;

//     //     // PrimaryRayData refractedRayData = {};
//     //     // refractedRayData.depth = hitData.depth + 1;
//     //     // refractedRayData.hitT = 0.f;
//     //     // refractedRayData.absorptionDistance = hitData.absorptionDistance;
//     //     // TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, refractRay, refractedRayData);

//     //     // refractedColor = refractedRayData.color.rgb;
//     // }

//     // //float3 extinctionCoeff = absorptionCoeff + scatteringCoeff;

//     // // 6. Mix results
//     // float3 totalLight = tir ? reflectedColor :
//     //         reflectionRatio * reflectedColor + refractionRatio * refractedColor;

//     // //float3 totalLight = CalculateTransmittanceRaymarching(attribs, rayDir);

//     // reflectedColor *= reflectionRatio;
//     // refractedColor = float3(0, 0, 1) * refractionRatio;

//     // float3 finalColor = reflectedColor + refractedColor;

//     // hitData.color = float4(totalLight, 1);
//     // hitData.hitT = currentRayT;
// }

float SampleDensityMap(float3 pos)
{
    if (useTransformations)
    {
        // // Apply render-scale *before* transforming to local space
        // float3 renderVolumeCenter = float3(0, 0, 0);
        // float3 renderScale = float3(1, 1, 1);
        // float3 scaledPos = (pos - renderVolumeCenter) / renderScale + renderVolumeCenter;

        // float3 local = mul(worldToLocal, float4(scaledPos, 1.0f)).xyz;
        // float3 uvw = local * 0.5f + 0.5f;

        // const float epsilon = 0.0001f;
        // if (any(uvw < epsilon || uvw > 1.0f - epsilon))
        //     return 0.0f;

        // return gTexture3D.SampleLevel(linearClampSampler, uvw, 0).r - volumeValueOffset;

        // Step 1 — Transform to local space
        float3 localPos = mul(worldToLocal, float4(pos, 1.0f)).xyz;

        // Step 2 — Convert from [-1,1] to [0,1]
        float3 uvw = localPos * 0.5f + 0.5f;

        const float epsilon = 0.0001f;
        if (any(uvw <= epsilon || uvw >= 1.0f - epsilon))
            return 0.0f; // or some safe value like -volumeValueOffset

        // Step 3 — Sample with clamp (no need for early out if clamp is used)
        return gTexture3D.SampleLevel(linearClampSampler, uvw, 0).r - volumeValueOffset;
    }
    else
    {
        float3 uvw = (pos + simBounds * 0.5f) / simBounds;

        const float epsilon = 0.0001f;
        if (any(uvw <= epsilon || uvw >= 1.0f - epsilon))
            return 0.0f; // or some safe value like -volumeValueOffset

        return gTexture3D.SampleLevel(linearClampSampler, uvw, 0).r - volumeValueOffset;
    }
}

// Returns the closest axis-aligned face normal of the local AABB box [-1, 1]^3
float3 CalculateClosestFaceNormal(float3 boxSize, float3 p)
{
    float3 halfSize = boxSize * 0.5;
    float3 o = (halfSize - abs(p));
    return (o.x < o.y && o.x < o.z) ? 
        float3(sign(p.x), 0, 0) : (o.y < o.z) ? float3(0, sign(p.y), 0) : float3(0, 0, sign(p.z));
}

float3 ComputeDensityNormal(float3 pos)
{
    if (useTransformations)
    {
        // Transform the world position to local fluid space
        float3 posLocal = mul(worldToLocal, float4(pos, 1.0f)).xyz;

        float3 offsetX = float3(1, 0, 0) * normalOffset;
        float3 offsetY = float3(0, 1, 0) * normalOffset;
        float3 offsetZ = float3(0, 0, 1) * normalOffset;

        // Compute sample positions in local space
        float3 posX1 = posLocal - offsetX;
        float3 posX2 = posLocal + offsetX;
        float3 posY1 = posLocal - offsetY;
        float3 posY2 = posLocal + offsetY;
        float3 posZ1 = posLocal - offsetZ;
        float3 posZ2 = posLocal + offsetZ;

        // Transform to world space before sampling
        float dx = SampleDensityMap(mul(localToWorld, float4(posX1, 1)).xyz)
                - SampleDensityMap(mul(localToWorld, float4(posX2, 1)).xyz);
        float dy = SampleDensityMap(mul(localToWorld, float4(posY1, 1)).xyz)
                - SampleDensityMap(mul(localToWorld, float4(posY2, 1)).xyz);
        float dz = SampleDensityMap(mul(localToWorld, float4(posZ1, 1)).xyz)
                - SampleDensityMap(mul(localToWorld, float4(posZ2, 1)).xyz);

        float3 volumeNormal = normalize(float3(dx, dy, dz));

        // Correct boundary flattening
        float3 localHalf = ScaledSimBounds * 0.5;
        float3 o = localHalf - abs(posLocal);
        float faceWeight = min(o.x, min(o.y, o.z));

        float3 faceNormalLocal = CalculateClosestFaceNormal(ScaledSimBounds, posLocal);
        float3 faceNormalWorld = normalize(mul((float3x3)localToWorld, faceNormalLocal));

        const float smoothDst = 0.3;
        const float smoothPow = 5;
        faceWeight = (1 - smoothstep(0, smoothDst, faceWeight)) * (1 - pow(saturate(volumeNormal.y), smoothPow));

        return normalize(volumeNormal * (1 - faceWeight) + faceNormalWorld * faceWeight);
    }
    else
    {
        float3 uvw = (pos + ScaledSimBounds * 0.5) / ScaledSimBounds;

        // const float s = 0.1;
        float3 offsetX = float3(1, 0, 0) * normalOffset;
        float3 offsetY = float3(0, 1, 0) * normalOffset;
        float3 offsetZ = float3(0, 0, 1) * normalOffset;

        float dx = SampleDensityMap(pos - offsetX) - SampleDensityMap(pos + offsetX);
        float dy = SampleDensityMap(pos - offsetY) - SampleDensityMap(pos + offsetY);
        float dz = SampleDensityMap(pos - offsetZ) - SampleDensityMap(pos + offsetZ);

        float3 volumeNormal = normalize(float3(dx, dy, dz));

        // Smoothly flatten normals out at boundary edges
        float3 o = ScaledSimBounds / 2 - abs(pos);
        float faceWeight = min(o.x, min(o.y, o.z));
        float3 faceNormal = CalculateClosestFaceNormal(ScaledSimBounds, pos);
        const float smoothDst = 0.3;
        const float smoothPow = 5;
        faceWeight = (1 - smoothstep(0, smoothDst, faceWeight)) * (1 - pow(saturate(volumeNormal.y), smoothPow));

        return normalize(volumeNormal * (1 - faceWeight) + faceNormal * (faceWeight));
    }

}

// Test if a ray segment <RayTMin(), RayTCurrent()> intersects an AABB.
// Limitation: this test does not take RayFlags into consideration and does not calculate a surface normal.
// Ref: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
bool RayAABBIntersectionTest(Ray ray, float3 aabb[2], out float tmin, out float tmax)
{
    float3 tmin3 = float3(0.f);
    float3 tmax3 = float3(0.f);
    int3 sign3 = ray.dir > 0;

    // Handle rays parallel to any x|y|z slabs of the AABB.
    // If a ray is within the parallel slabs,
    //  the tmin, tmax will get set to -inf and +inf
    //  which will get ignored on tmin/tmax = max/min.
    // If a ray is outside the parallel slabs, -inf/+inf will
    //  make tmax > tmin fail (i.e. no intersection).
    // TODO: handle cases where ray origin is within a slab
    //  that a ray direction is parallel to. In that case
    //  0 * INF => NaN
    const float FLT_INFINITY = 1e+32;
    float3 invRayDirection = select(ray.dir != 0, 1 / ray.dir, select(ray.dir > 0, FLT_INFINITY, -FLT_INFINITY));

    tmin3.x = (aabb[1 - sign3.x].x - ray.origin.x) * invRayDirection.x;
    tmax3.x = (aabb[sign3.x].x - ray.origin.x) * invRayDirection.x;

    tmin3.y = (aabb[1 - sign3.y].y - ray.origin.y) * invRayDirection.y;
    tmax3.y = (aabb[sign3.y].y - ray.origin.y) * invRayDirection.y;

    tmin3.z = (aabb[1 - sign3.z].z - ray.origin.z) * invRayDirection.z;
    tmax3.z = (aabb[sign3.z].z - ray.origin.z) * invRayDirection.z;

    tmin = max(max(tmin3.x, tmin3.y), tmin3.z);
    tmax = min(min(tmax3.x, tmax3.y), tmax3.z);

    return tmax > tmin && tmax >= RayTMin() && tmin <= RayTCurrent();
}

// Returns (dstToBox, dstInsideBox) in world space
// If ray misses the transformed box, dstInsideBox == 0
float2 RayBoxDst_Matrix(float3 worldRayOrigin, float3 worldRayDir, float4x4 localToWorld, float4x4 worldToLocal)
{
    if (useTransformations)
    {
        // Step 1 — Transform ray to local space
        float3 localOrigin = mul(worldToLocal, float4(worldRayOrigin, 1)).xyz;
        float3 localDir = normalize(mul(worldToLocal, float4(worldRayDir, 0)).xyz);

        // Step 2 — Classic slab intersection in local AABB space [-1, 1]
        float3 boundsMin = float3(-1.0, -1.0, -1.0);
        float3 boundsMax = float3(1.0, 1.0, 1.0);
        float3 invDir = 1.0 / localDir;

        float3 t0 = (boundsMin - localOrigin) * invDir;
        float3 t1 = (boundsMax - localOrigin) * invDir;

        float3 tmin = min(t0, t1);
        float3 tmax = max(t0, t1);

        float tEnter = max(max(tmin.x, tmin.y), tmin.z);
        float tExit = min(min(tmax.x, tmax.y), tmax.z);

        // Step 3 — If ray misses the box
        if (tEnter > tExit || tExit < 0)
            return float2(0, 0);

        // Step 4 — Transform distances back to world space length (by ray direction projection)
        float3 localEntry = localOrigin + tEnter * localDir;
        float3 localExit = localOrigin + tExit * localDir;

        float3 worldEntry = mul(localToWorld, float4(localEntry, 1)).xyz;
        float3 worldExit = mul(localToWorld, float4(localExit, 1)).xyz;

        float dstToBox = max(0, dot(worldEntry - worldRayOrigin, normalize(worldRayDir)));
        float dstInsideBox = max(0, dot(worldExit - worldEntry, normalize(worldRayDir)));

        return float2(dstToBox, dstInsideBox);
    }
    else 
    {
        return RayBoxDst(-simBounds * 0.5, simBounds * 0.5, worldRayOrigin, worldRayDir);
    }
  
}

// Returns (dstToBox, dstInsideBox). If ray misses box, dstInsideBox will be zero
float2 RayBoxDst(float3 boundsMin, float3 boundsMax, float3 rayOrigin, float3 rayDir)
{
    float3 invRayDir = 1 / rayDir;
    // Adapted from: http://jcgt.org/published/0007/03/04/
    float3 t0 = (boundsMin - rayOrigin) * invRayDir;
    float3 t1 = (boundsMax - rayOrigin) * invRayDir;
    float3 tmin = min(t0, t1);
    float3 tmax = max(t0, t1);

    float dstA = max(max(tmin.x, tmin.y), tmin.z);
    float dstB = min(tmax.x, min(tmax.y, tmax.z));

    // CASE 1: ray intersects box from outside (0 <= dstA <= dstB)
    // dstA is dst to nearest intersection, dstB dst to far intersection

    // CASE 2: ray intersects box from inside (dstA < 0 < dstB)
    // dstA is the dst to intersection behind the ray, dstB is dst to forward intersection

    // CASE 3: ray misses box (dstA > dstB)

    float dstToBox = max(0, dstA);
    float dstInsideBox = max(0, dstB - dstToBox);
    return float2(dstToBox, dstInsideBox);
}

float CalculateDensityAlongRay(float3 rayPos, float3 rayDir, float stepSize)
{
    // Test for non-normalize ray and return 0 in that case.
    // This happens when refract direction is calculated, but ray is totally reflected
    if (dot(rayDir, rayDir) < 0.9) return 0;

    //float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, rayPos, rayDir);
    float2 boundsDstInfo = RayBoxDst_Matrix(rayPos, rayDir, localToWorld, worldToLocal);
    float dstToBounds = boundsDstInfo[0];
    float dstThroughBounds = boundsDstInfo[1];
    if (dstThroughBounds <= 0) return 0;

    float dstTravelled = 0;
    float opticalDepth = 0;
    float nudge = stepSize * 0.5;
    float3 entryPoint = rayPos + rayDir * (dstToBounds + nudge);
    dstThroughBounds -= (nudge + TinyNudge);

    while (dstTravelled < dstThroughBounds)
    {
        rayPos = entryPoint + rayDir * dstTravelled;
        float density = SampleDensityMap(rayPos) * DensityRayMarchMultiplier * stepSize;
        if (density > 0)
        {
            opticalDepth += density;
        }
        dstTravelled += stepSize;
    }

    return opticalDepth;
}

// Test intersection of ray with unit box centered at origin
CubeHitInfo RayUnitBox(float3 pos, float3 dir)
{
    float3 boxMin = -1;
    float3 boxMax = 1;

    // if (!unitBox)
    // {
    //     boxMin = -0.5f * simBounds;
    //     boxMax = 0.5f * simBounds;
    // }

    float3 invDir = 1 / dir;

    // Thanks to https://tavianator.com/2011/ray_box.html
    float3 tMin = (boxMin - pos) * invDir;
    float3 tMax = (boxMax - pos) * invDir;
    float3 t1 = min(tMin, tMax);
    float3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);

    // Set hit info
    CubeHitInfo hitInfo = (CubeHitInfo)0;
    hitInfo.dst = asfloat(0x7F800000); // +Infinity.
    hitInfo.didHit = tFar >= tNear && tFar > 0;
    hitInfo.isInside = tFar > tNear && tNear <= 0;

    if (hitInfo.didHit)
    {
        float hitDst = hitInfo.isInside ? tFar : tNear;
        float3 hitPos = pos + dir * hitDst;

        hitInfo.dst = hitDst;
        hitInfo.hitPoint = hitPos;

        // Calculate normal
        float3 o = (1 - abs(hitPos));
        float3 absNormal = (o.x < o.y && o.x < o.z) ? float3(1, 0, 0) : (o.y < o.z) ? float3(0, 1, 0) : float3(0, 0, 1);
        hitInfo.normal = absNormal * sign(hitPos) * (hitInfo.isInside ? -1 : 1);
    }

    return hitInfo;
}

CubeHitInfo RayBoxWithMatrix(float3 rayPos, float3 rayDir, float4x4 localToWorld, float4x4 worldToLocal)
{
    float3 posLocal = mul(worldToLocal, float4(rayPos, 1)).xyz;
    float3 dirLocal = normalize(mul(worldToLocal, float4(rayDir, 0)).xyz);
    CubeHitInfo hitInfo = RayUnitBox(posLocal, dirLocal);
    hitInfo.normal = normalize(mul(localToWorld, float4(hitInfo.normal, 0.0)).xyz);
    hitInfo.hitPoint = mul(localToWorld, float4(hitInfo.hitPoint, 1.0)).xyz;

    if (hitInfo.didHit) hitInfo.dst = length(hitInfo.hitPoint - rayPos);
    return hitInfo;
}

[shader("intersection")]
void RaymarchingIntersection()
{
    if (!drawFluid)
        return;

    float3 rayOrigin = WorldRayOrigin();
    float3 rayDir = normalize(WorldRayDirection());

    // Step 1: Intersect with the transformed fluid bounding box (world space)
    CubeHitInfo cubeHitInfo = RayBoxWithMatrix(rayOrigin, rayDir, localToWorld, worldToLocal);

    if (!cubeHitInfo.didHit)
        return;

    if (useRecursiveRaytracing)
    {
        bool travellingThroughFluid = IsInsideFluid(rayOrigin, rayDir);

        uint dummyRNG = 0; // RNG state is not useful here unless noise affects surface detection

        // Run the raymarch to find the first fluid surface hit
        SurfaceInfo surfaceInfo = FindNextSurface(rayOrigin, rayDir, !travellingThroughFluid, dummyRNG, 1, cubeHitInfo.dst);

        if (!surfaceInfo.foundSurface)
            return;

        SDF_GeometryAttributes attr = {};
        attr.pos = surfaceInfo.pos;
        attr.normal = ComputeDensityNormal(surfaceInfo.pos); // Don't flip here — let closest hit handle that
        attr.densityAlongRay = surfaceInfo.densityAlongRay;
        attr.isInsideFluid = travellingThroughFluid;

        // Report the intersection to continue in closest hit
        float hitDistance = length(surfaceInfo.pos - rayOrigin); // Accurate distance from ray start
        ReportHit(hitDistance, 0, attr);

        return;
    }
    else
    {
        SDF_GeometryAttributes attr2 = {};
        attr2.pos = cubeHitInfo.hitPoint;
        ReportHit(1, 0, attr2);
        return;
    }

    // if (!drawFluid)
    //     return;

    // float3 rayOrigin = WorldRayOrigin();
    // float3 rayDir = WorldRayDirection();

    // // Test against the transformed fluid container
    // CubeHitInfo cubeHitInfo = RayBoxWithMatrix(rayOrigin, rayDir, localToWorld, worldToLocal);

    // if (cubeHitInfo.didHit)
    // {
    //     SDF_GeometryAttributes attr = {};
    //     attr.pos = cubeHitInfo.hitPoint;
    //     attr.normal = cubeHitInfo.normal;
    //     // Optionally encode hit.normal or hit.hitPoint into attr
    //     ReportHit(cubeHitInfo.dst, 0 /*geometry index*/, attr);
    // }

    // return;

    // Ray ray;
    // ray.origin = rayOrigin;
    // ray.dir = rayDir;

    // float3 aabb[2] = { float3(-0.5f) * simBounds, float3(0.5f) * simBounds };

    // float tMin = 0.f;
    // float tMax = 0.f;
    // SDF_GeometryAttributes attr = {};
    
    // if (RayAABBIntersectionTest(ray, aabb, tMin, tMax))
    // {
    //     // Ensure tMin is in the AABB when the camera is already inside it.
    //     tMin = max(tMin, RayTMin());
    //     tMax = min(tMax, RayTCurrent());

    //     float3 startPos = rayOrigin + rayDir * tMin;

    //     // float prevDensity = SampleDensityMap(startPos);

    //     ReportHit(tMax - tMin, 0, attr);
    //     return;
    // }

    // bool searchEntry = !IsInsideFluid(rayOrigin); // This can cause white artifacts on top of fluid surface.
    // SurfaceInfo hitInfo = FindNextSurface(rayOrigin, rayDir, searchEntry, RayTCurrent());

    // if (!hitInfo.foundSurface)
    //     return;


    // attr.pos = hitInfo.pos;
    // attr.normal = ComputeDensityNormal(hitInfo.pos);
    // attr.densityAlongRay = hitInfo.densityAlongRay;

    // // tHit computation
    // float tHit = dot(hitInfo.pos - rayOrigin, rayDir);

    // ReportHit(tHit, 0, attr);
}

// float3 aabb[2] = { float3(-0.5f) * simBounds, float3(0.5f) * simBounds };

// float tMin = 0.f;
// float tMax = 0.f;

// if (RayAABBIntersectionTest(ray, aabb, tMin, tMax))
// {
//     // Ensure tMin is in the AABB when the camera is already inside it.
//     tMin = max(tMin, RayTMin());
//     tMax = min(tMax, RayTCurrent());

//     float3 startPos = rayOrigin + rayDir * tMin;

//     // float prevDensity = SampleDensityMap(startPos);

//     // ReportHit(tMax - tMin, 0, attr);
//     // return;
// }

// ReportHit(tMin, /*hitKind*/ 0, attr);
// return;

//     ReportHit(tMin, /*hitKind*/ 0, attr);
//     return;

//     float prevDensity = SampleDensityMap(startPos);

//     // If the entry point intersection is already inside the fluid. ReportHit.
//     if (prevDensity >= isoLevel)
//     {
//         attr.w_normal = ComputeDensityNormal(rayOrigin + rayDir * tMin);
//         attr.w_position = startPos;
//         // attr.tMin = tMin;
//         // attr.tMax = tMax;

//         ReportHit(tMin, /*hitKind*/ 0, attr);
//         return;
//     }

//     for (float dst = tMin + marchSize; dst < tMax; dst += marchSize)
//     {
//         float3 p = rayOrigin + rayDir * dst;
//         float currDensity = SampleDensityMap(p);

//         // Detect a crossing from outside to inside fluid
//         if (prevDensity < isoLevel && currDensity >= isoLevel)
//         {
//             // Optional: refine hit position with linear interpolation
//             float alpha = (isoLevel - prevDensity) / (currDensity - prevDensity);
//             float hitT = dst - marchSize + alpha * marchSize;

//             attr.w_normal = ComputeDensityNormal(rayOrigin + rayDir * hitT);
//             attr.w_position = rayOrigin + rayDir * hitT;
//             // attr.tMin = hitT;
//             // attr.tMax = tMax;

//             ReportHit(hitT, /*hitKind*/ 0, attr);
//             return;
//         }

//         prevDensity = currDensity;
//     }
// }

// bool findEntry = prevDensity < 0;

// If the entry point intersection is already inside the fluid. ReportHit.
// if (findEntry)
// {
//     attr.w_normal = ComputeDensityNormal(rayOrigin + rayDir * tMin);
//     attr.w_position = startPos;
//     // attr.tMin = tMin;
//     // attr.tMax = tMax;

//     ReportHit(tMin, /*hitKind*/ 0, attr);
//     return;
// }

// for (float dst = tMin + marchSize; dst < tMax + marchSize; dst += marchSize)
// {
//     float3 p = rayOrigin + rayDir * dst;
//     float currDensity = SampleDensityMap(p);

//     // Detect a crossing from outside to inside fluid
//     if (prevDensity < 0 && currDensity >= 0 && findEntry)
//     {
//         // Optional: refine hit position with linear interpolation
//         float alpha = (0 - prevDensity) / (currDensity - prevDensity);
//         float hitT = dst - marchSize + alpha * marchSize;

//         attr.w_normal = ComputeDensityNormal(rayOrigin + rayDir * hitT);
//         attr.w_position = rayOrigin + rayDir * hitT;
//         // attr.tMin = hitT;
//         // attr.tMax = tMax;

//         ReportHit(hitT, /*hitKind*/ 0, attr);
//         return;
//     }

//     // Detect a crossing from inside to outside fluid
//     if (prevDensity >= 0 && currDensity < 0 && !findEntry)
//     {
//         // Optional: refine hit position with linear interpolation
//         float alpha = (0 - prevDensity) / (currDensity - prevDensity);
//         float hitT = dst - marchSize + alpha * marchSize;

//         attr.w_normal = ComputeDensityNormal(rayOrigin + rayDir * hitT);
//         attr.w_normal = -attr.w_normal;
//         attr.w_position = rayOrigin + rayDir * hitT;
//         // attr.tMin = hitT;
//         // attr.tMax = tMax;

//         ReportHit(hitT, /*hitKind*/ 0, attr);
//         return;
//     }

//     prevDensity = currDensity;
// }

// Ray generation (just a simple ray going through the scene)
[shader("raygeneration")]
void rayGen()
{
    uint3 launchIndex = DispatchRaysIndex();
    uint2 dim = DispatchRaysDimensions().xy;
    float2 uv = launchIndex.xy / (float2)dim;

    RayDesc ray;
    ray = gScene.camera.computeRayPinhole(launchIndex.xy, vector<uint, 2>(viewportDims), false).toRayDesc();

    PrimaryRayData hitData;
    hitData.color = float4(0.f, 0.f, 0.f, 0.f);
    hitData.depth = 0;
    hitData.uv = uv;
    TraceRay(gScene.rtAccel, RAY_FLAG_NONE /*rayFlags*/, 0xFF, /*Test all Geometry*/ 0 /* ray index*/, rayTypeCount, 0, ray, hitData);

    gOutput[launchIndex.xy] = hitData.color;

    // Use RayQuery to find bounds of the fluid AABB
    // RayQuery<RAY_FLAG_NONE> rq;
    // rq.TraceRayInline(Scene, RAY_FLAG_NONE, 0xFF, rayOrigin, 0.001f, rayDir, 1e38f);

    // float tMin = 0.0f;
    // float tMax = 0.0f;
    // bool fluidAABBHit = false;

    // while (rq.Proceed())
    // {
    //     if (rq.CandidateType() == CANDIDATE_PROCEDURAL_PRIMITIVE) // you must define this index
    //     {
    //         tMin = rq.(); // entry into AABB
    //         tMax = rq.CandidateTriangleRayTMax(); // exit from AABB
    //         fluidAABBHit = true;
    //         rq.CommitProceduralPrimitiveHit(); // Needed to finalize hit
    //     }
    // }

    // if (fluidAABBHit)
    // {
    //     SDF_GeometryAttributes attr = RaymarchWithRayQuery(rayOrigin, rayDir, tMin, tMax);

    //     if (attr.t > 0)
    //     {
    //         // Fluid surface was hit → shade
    //         gOutput[launchIndex] = float4(0.1, 0.3, 0.8, 1); // blue fluid
    //         return;
    //     }
    // }

    // // Otherwise trace scene normally
    // PrimaryRayData payload;
    // payload.color = float4(0, 0, 0, 0);
    // TraceRay(Scene, RAY_FLAG_NONE, 0xFF, rayOrigin, 0.001f, rayDir, 1e38f, 0, 1, payload);
    // gOutput[launchIndex] = payload.color;
}

