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
// RWTexture3D<float> gTexture3D;
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
    float volumeValueOffset;
};



struct SDF_GeometryAttributes
{
    float3 pos;
    float densityAlongRay;
    // float t;
    //float tMin;
    float3 normal;
    //float tMax;
}

struct PrimaryRayData
{
    float4 color = float4(0, 0, 0, 0);
    float3 transmittance = float3(1, 1, 1);
    float3 light = float3(0, 0, 0);
    uint depth = 0;
    float hitT = -1.f;
    float absorptionDistance = 0.f;
};

struct SurfaceInfo
{
    float3 pos;
    float densityAlongRay;
    bool foundSurface;
};

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
    // float3 rayDir = WorldRayDirection();
    // float3 rayOrigin = WorldRayOrigin();

    // SDF_GeometryAttributes attr;

    // float3 finalColor = float3(0, 0, 0);
    // bool hasHit = false;

    // float t = RayTMin(); // total distance travelled

    // // Raymarching.
    // for (int i = 0; i < 80 && t < RayTCurrent(); i++)
    // {
    //     float3 p = rayOrigin + rayDir * t; // Position along the ray.

    //     float d = map(p, waterTurbulence, time, iFrame); // Current distance to the scene.

    //     // Early stop if close enough or if too far.
    //     if (d < 0.001f)
    //     {
    //         hasHit = true;
    //         finalColor = float3(0.1f, 0.3f, 0.8f);
    //         // attr.w_position = p;
    //         // attr.w_normal = sdCalculateNormal(p, waterTurbulence, time, iFrame);
    //         // attr.t = max(t, 0.001f);

    //         // ReportHit(t, /*hitKind*/ 0, attr);
    //         // break;
    //     }

    //     t += max(d, 0.01f); // March the ray.

    //     if (t > 100.f)
    //     {
    //         break;
    //     }
    // }


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
    // Get the hit-point data.
    float3 rayDir = WorldRayDirection();
    float currentRayDistance = RayTCurrent();
    uint triangleIndex = PrimitiveIndex();

    if (!drawFluid)
    {
        hitData.color = float4(0.1f, 0.2f, 1.f, 1.f);
        hitData.hitT = currentRayDistance;
        hitData.absorptionDistance = 0;
        return;
    }

    if (hitData.depth >= maxRayBounce)
    {
        hitData.color = float4(1.f, 1.f, 0.f, 1);
        hitData.hitT = 0;
        hitData.absorptionDistance = 0;
        return;
    }

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
        hitData.color = float4(backgroundColor, 1.f);
        // hitData.color = float4(bsdfProperties.diffuseReflectionAlbedo, 1.f);
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

SurfaceInfo FindNextSurface(float3 origin, float3 rayDir, bool findNextFluidEntryPoint, float maxDst)
{
    SurfaceInfo info;
    info.pos = float3(0, 0, 0);
    info.densityAlongRay = 0;
    info.foundSurface = false;

    // Safety check
    if (dot(rayDir, rayDir) < 0.5) return info;

    // Compute AABB intersection
    float3 aabbMin = -0.5f * simBounds;
    float3 aabbMax = 0.5f * simBounds;
    float2 bounds = RayBoxDst(aabbMin, aabbMax, origin, rayDir);

    float entryT = max(bounds.x, 0.0);
    float travelT = bounds.y;
    float dstToTest = travelT - 1e-4;

    origin = origin + rayDir * entryT;

    float stepSize = marchSize;

    bool hasEnteredFluid = false;
    bool hasExitedFluid = false;
    float3 lastPosInFluid = origin;

    for (float dst = 0; dst < dstToTest; dst += stepSize)
    {
        float3 samplePos = origin + rayDir * dst;
        float density = SampleDensityMap(samplePos) * DensityRayMarchMultiplier * stepSize;
        bool insideFluid = (density > 0);

        if (insideFluid)
        {
            hasEnteredFluid = true;
            lastPosInFluid = samplePos;

            if (dst <= maxDst)
            {
                info.densityAlongRay += density;
            }
        }
        else
        {
            hasExitedFluid = true;
        }

        bool found = false;
        if (findNextFluidEntryPoint)
        {
            found = insideFluid && hasExitedFluid;
        }
        else
        {
            bool isLastStep = (dst + stepSize >= dstToTest);
            found = hasEnteredFluid && (!insideFluid || isLastStep);
        }

        if (found)
        {
            info.pos = lastPosInFluid;
            info.foundSurface = true;
            break;
        }
    }

    return info;
}

bool IsInsideFluid(float3 pos)
{
    float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, pos, float3(0, 0, 1));
    return (boundsDstInfo.x <= 0 && boundsDstInfo.y > 0) && SampleDensityMap(pos) > 0;
}

[shader("closesthit")]
void RaymarchingClosestHit(inout PrimaryRayData payload, SDF_GeometryAttributes attribs)
{
    Ray ray;
    ray.origin = WorldRayOrigin();
    ray.dir = WorldRayDirection();

    float densityStepSize = sunLightMarchSize * (payload.depth + 1); // increase step size with each iteration

    // === Raymarch once ===
    payload.transmittance *= exp(-attribs.densityAlongRay * scatteringCoeff);

    // === Surface interaction ===
    float3 normal = ComputeDensityNormal(attribs.pos);
    bool localTravel = payload.depth == 1 || payload.depth == 3;
    //bool localTravel = dot(normal, ray.dir) > 0;
    if (localTravel) normal = -normal;

    // Compute IoR ratio based on fluid state
    float eta = localTravel ? (IoR / 1.0) : (1.0 / IoR);

    // Fresnel reflectance using Schlick approximation
    float cosTheta = abs(dot(normalize(ray.dir), normal));
    float F0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
    float reflectionRatio = evalFresnelSchlick(F0, 1.0, cosTheta);
    float refractionRatio = 1.0 - reflectionRatio;

    // Compute directions
    float3 reflectDir = reflect(ray.dir, normal);
    float3 refractDir = refract(ray.dir, normal, eta);

    float densityRefrac = CalculateDensityAlongRay(attribs.pos, refractDir, densityStepSize);
    float densityReflect = CalculateDensityAlongRay(attribs.pos, reflectDir, densityStepSize);

    bool followRefrac = densityRefrac * refractionRatio > densityReflect * reflectionRatio;

    // Accumulate light from the less dominant path
    if (followRefrac)
    {
        payload.light += gScene.envMap.eval(reflectDir) * payload.transmittance * exp(-densityReflect * scatteringCoeff) * reflectionRatio;
    }
    else
    {
        payload.light += gScene.envMap.eval(refractDir) * payload.transmittance * exp(-densityRefrac * scatteringCoeff) * refractionRatio;
    }

    // Advance ray
    float3 nextOrigin = attribs.pos + ray.dir * 1e-3f; // Avoid self-hit
    float3 nextDir = followRefrac ? refractDir : reflectDir;

    // Update transmittance
    payload.transmittance *= followRefrac ? refractionRatio : reflectionRatio;

    // === Spawn new recursive ray ===
    if (payload.depth < maxRayBounce)
    {
        PrimaryRayData newPayload;
        newPayload.transmittance = payload.transmittance;
        newPayload.light = float3(0, 0, 0);
        newPayload.depth = payload.depth + 1;

        RayDesc newRay;
        newRay.Origin = nextOrigin;
        newRay.Direction = nextDir;
        newRay.TMin = 0.001;
        newRay.TMax = 1e+4;

        TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, 1, 0, newRay, newPayload);

        payload.light += newPayload.light; // Recover new light data for final shading.
    }

    if (payload.depth >= maxRayBounce)
    {
        // Max depth: sample env map
        float densityRemainder = CalculateDensityAlongRay(nextOrigin, nextDir, sunLightMarchSize);
        payload.light += gScene.envMap.eval(nextDir) * payload.transmittance * exp(-densityRemainder * scatteringCoeff);
    }

    payload.color = float4(payload.light, 1.0);
    payload.hitT = RayTCurrent();
}

// [shader("closesthit")]
// void RaymarchingClosestHit(inout PrimaryRayData hitData, SDF_GeometryAttributes attribs)
// {
//     Ray ray;
//     ray.origin = WorldRayOrigin();
//     ray.dir = WorldRayDirection();

//     bool travelThrougFluid = IsInsideFluid(ray.origin);

//     float3 transmittance = float3(1, 1, 1);
//     float3 light = float3(0, 0, 0);
//     // float3 transmittance = hitData.transmittance;
//     // float3 light = hitData.light;

//     for (int i = 0; i < maxRayBounce; i++)
//     //if (hitData.depth < maxRayBounce)
//     {
//         float densityStepSize = sunLightMarchSize * (i + 1); // increase step size with each iteration
//         //float densityStepSize = sunLightMarchSize * (hitData.depth + 1); // increase step size with each iteration

//         bool searchEntry = !travelThrougFluid;

//         SurfaceInfo hitInfo = FindNextSurface(ray.origin, ray.dir, searchEntry, RayTCurrent());

//         transmittance *= exp(-hitInfo.densityAlongRay * scatteringCoeff);
//         //transmittance *= exp(-attribs.densityAlongRay * scatteringCoeff);
//         //hitData.transmittance *= exp(-attribs.densityAlongRay * scatteringCoeff);

//         // If no surface found, exit loop but still accumulate final light
//         if (!hitInfo.foundSurface) break;
//         //if (!hitInfo.foundSurface) return;

//         float3 normal = ComputeDensityNormal(hitInfo.pos);
//         bool travelThrougFluid = (dot(normal, ray.dir) > 0);
//         if (travelThrougFluid) normal = -normal;
//         // float3 normal = attribs.normal;
//         // bool travelThrougFluid = (dot(attribs.normal, ray.dir) > 0);
//         // if (travelThrougFluid) attribs.normal = -attribs.normal;

//         // Compute IoR ratio based on fluid state
//         // float eta = travelThrougFluid ? 1.0 : IoR;
//         float eta = travelThrougFluid ? (IoR / 1.0) : (1.0 / IoR);

//         // Fresnel reflectance using Schlick approximation
//         float cosTheta = abs(dot(normalize(ray.dir), normal));
//         float F0 = pow((1.0 - eta) / (1.0 + eta), 2.0);
//         float reflectionRatio = evalFresnelSchlick(F0, 1.0, cosTheta);
//         float refractionRatio = 1.0 - reflectionRatio;

//         // Compute directions
//         float3 reflectDir = reflect(ray.dir, normal);
//         float3 refractDir = refract(ray.dir, normal, eta);
//         //bool tir = all(refractDir == float3(0.0, 0.0, 0.0)); // Total internal reflection

//         // Compute transmittance of both paths
//         //float densityRefrac = tir ? 0.0 : CalculateDensityAlongRay(hitInfo.pos, refractDir, densityStepSize);
//         float densityRefrac =  CalculateDensityAlongRay(hitInfo.pos, refractDir, densityStepSize);
//         float densityReflect = CalculateDensityAlongRay(hitInfo.pos, reflectDir, densityStepSize);
//         // float densityRefrac = tir ? 0.0 : CalculateDensityAlongRay(attribs.pos, refractDir, densityStepSize);
//         // float densityReflect = CalculateDensityAlongRay(attribs.pos, reflectDir, densityStepSize);

//         // float3 transRefract = exp(-densityRefrac * scatteringCoeff) * refractionRatio;
//         // float3 transReflect = exp(-densityReflect * scatteringCoeff) * reflectionRatio;

//         //bool followRefrac = !tir && length(transRefract) > length(transReflect);
//         bool followRefrac = densityRefrac * refractionRatio > densityReflect * reflectionRatio;

//         // Accumulate light from the less dominant path
//         if (followRefrac)
//         {
//             // light += 1.f * transmittance * transReflect;
//             // light += gScene.envMap.eval(reflectDir) * transmittance * transReflect;
//             light += gScene.envMap.eval(reflectDir) * transmittance * exp(-densityReflect * scatteringCoeff) * reflectionRatio;
//             //hitData.light += 1.f * hitData.transmittance * transReflect;
//         }
//         else
//         {
//             // light += 1.f * transmittance * transRefract;
//             // light += gScene.envMap.eval(refractDir) * transmittance * transReflect;
//             light += gScene.envMap.eval(refractDir) * transmittance * exp(-densityRefrac * scatteringCoeff) * refractionRatio;
//             //hitData.light += 1.f * hitData.transmittance * transRefract;
//         }

//         // Advance ray
//         ray.origin = hitInfo.pos + ray.dir * 1e-3f; // Avoid self-hit
//         //ray.origin = attribs.pos + ray.dir * 1e-3f; // Avoid self-hit
//         ray.dir = followRefrac ? refractDir : reflectDir;

//         // Update transmittance
//         transmittance *= followRefrac ? refractionRatio : reflectionRatio;
//         // hitData.transmittance *= followRefrac ? refractionRatio : reflectionRatio;
//     }

//     // Final raymarch if no more bounces
//     float densityRemainder = CalculateDensityAlongRay(ray.origin, ray.dir, sunLightMarchSize);
//     //light += 1.f * transmittance * exp(-densityRemainder * scatteringCoeff);
//     light += gScene.envMap.eval(ray.dir) * transmittance * exp(-densityRemainder * scatteringCoeff);

//     //hitData.color = float4(1 - light, 1.0);
//     hitData.color = float4(light, 1.0);
//     hitData.hitT = RayTCurrent();
//     return;
// }

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
    float3 uvw = (pos + simBounds * 0.5f) / simBounds;

    const float epsilon = 0.0001f;
    bool isEdge = any(uvw >= 1 - epsilon || uvw <= epsilon);
    if (isEdge) return -volumeValueOffset;

    return gTexture3D.SampleLevel(linearClampSampler, uvw, 0).r - volumeValueOffset;
}

float3 CalculateClosestFaceNormal(float3 boxSize, float3 p)
{
    float3 halfSize = boxSize * 0.5;
    float3 o = (halfSize - abs(p));
    return (o.x < o.y && o.x < o.z) ? float3(sign(p.x), 0, 0) : (o.y < o.z) ? float3(0, sign(p.y), 0) : float3(0, 0, sign(p.z));
}

float3 ComputeDensityNormal(float3 pos)
{
    float3 offsetX = float3(1, 0, 0) * normalOffset;
    float3 offsetY = float3(0, 1, 0) * normalOffset;
    float3 offsetZ = float3(0, 0, 1) * normalOffset;

    float dx = SampleDensityMap(pos - offsetX) - SampleDensityMap(pos + offsetX);
    float dy = SampleDensityMap(pos - offsetY) - SampleDensityMap(pos + offsetY);
    float dz = SampleDensityMap(pos - offsetZ) - SampleDensityMap(pos + offsetZ);

    float3 volumeNormal = normalize(float3(dx, dy, dz));

    // Smoothly flatten normals out at boundary edges
    float3 o = simBounds / 2 - abs(pos);
    float faceWeight = min(o.x, min(o.y, o.z));
    float3 faceNormal = CalculateClosestFaceNormal(simBounds, pos);
    const float smoothDst = 0.3;
    const float smoothPow = 5;
    faceWeight = (1 - smoothstep(0, smoothDst, faceWeight)) * (1 - pow(saturate(volumeNormal.y), smoothPow));

    return normalize(volumeNormal * (1 - faceWeight) + faceNormal * (faceWeight));
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

    float opticalDepth = 0.f;

    float3 aabbMin = -0.5f * simBounds;
    float3 aabbMax = 0.5f * simBounds;

    // Compute AABB intersection with full ray
    float2 boundsDst = RayBoxDst(aabbMin, aabbMax, rayPos, rayDir);
    float dstToBox = boundsDst.x;
    float dstInsideBox = boundsDst.y;

    if (dstInsideBox <= 0) return 0;

    float entryT = max(dstToBox, 0.0); // clamp if already inside box
    float3 startPos = rayPos + rayDir * (entryT + stepSize * 0.5);

    float maxTravel = dstInsideBox - (entryT + 1e-4);
    if (maxTravel <= 0) return 0;

    for (float dst = 0.f; dst < maxTravel; dst += stepSize)
    {
        float3 p = startPos + rayDir * dst;
        float density = SampleDensityMap(p);

        if (density > 0)
        {
            opticalDepth += density * DensityRayMarchMultiplier * stepSize;
        }        
    }

    return opticalDepth;
}

[shader("intersection")]
void RaymarchingIntersection()
{
    if (!drawFluid)
    {
        return;
    }

    SDF_GeometryAttributes attr = {};

    float3 rayDir = WorldRayDirection();
    float3 rayOrigin = WorldRayOrigin();

    Ray ray;
    ray.origin = rayOrigin;
    ray.dir = rayDir;

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

    bool searchEntry = !IsInsideFluid(rayOrigin);
    SurfaceInfo hitInfo = FindNextSurface(rayOrigin, ray.dir, searchEntry, RayTCurrent());

    if (hitInfo.foundSurface)
    {
        attr.pos = hitInfo.pos;
        attr.normal = ComputeDensityNormal(hitInfo.pos);
        attr.densityAlongRay = hitInfo.densityAlongRay;

        float tHit = length(hitInfo.pos - rayOrigin); // or dot(ray.dir, hitInfo.pos - ray.origin)

        ReportHit(tHit, /*hitKind*/ 0, attr);
        return;
    }
}

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

    RayDesc ray;
    ray = gScene.camera.computeRayPinhole(launchIndex.xy, vector<uint, 2>(viewportDims), true).toRayDesc();

    PrimaryRayData hitData;
    hitData.color = float4(0.f, 0.f, 0.f, 0.f);
    hitData.depth = 0;
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

