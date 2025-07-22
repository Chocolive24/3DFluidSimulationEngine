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
//import SDF_Functions;

import Scene.Raytracing;
import Utils.Sampling.TinyUniformSampleGenerator;
import Rendering.Materials.Fresnel;
import Rendering.Materials.TexLODHelpers;
import Rendering.Lights.LightHelpers;

RWTexture2D<float4> gOutput;
Texture3D<float> gTexture3D;
SamplerState linearClampSampler;

RWStructuredBuffer<uint> cutomPrimitveMasks;

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
    float smoothDst;
    float smoothPow;

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

    bool useVoxelOpti;
    bool debugVoxelGrid;
    uint voxelGridTotalResolution;

    float shadowDensityMultiplier;
    bool drawShadow;

    bool approximateSecondaryRayBounce;
    bool useRaymarching;
};



struct ImplicitGeometryAttributes
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
    // float3 debugNormal;
    // float2 uv = float2(0, 0);
    uint depth = 0;
    float hitT = -1.f;
    float absorptionDistance = 0.f;
};

struct ShadowRayData
{
    float3 color = float3(1);
    float absorptionDistance = 0.f;
    float lastT;
    bool hit = false;
    bool insideFluid;
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

// Miss shader
[shader("miss")]
void miss(inout PrimaryRayData hitData)
{
    // Use the ray direction to sample the environment map
    float3 envColor = gScene.envMap.eval(WorldRayDirection());
    // hitData.color = hasHit ? float4(0.1, 0.3, 0.8, 1.0) : float4(envColor, 1.f);
    hitData.color = float4(envColor, 1.f);
    hitData.hitT = -1.f;
    //hitData.absorptionDistance = 0;
}

[shader("closesthit")]
void closestHit(inout PrimaryRayData hitData, BuiltInTriangleIntersectionAttributes attribs)
{
    float3 rayDir = WorldRayDirection();
    float currentRayDistance = hitData.depth == 0 ? 0 : RayTCurrent();
    uint triangleIndex = PrimitiveIndex();

    if (!drawFluid)
    {
        hitData.color = float4(0.1f, 0.2f, 1.f, 1.f);
        hitData.hitT = currentRayDistance;
        hitData.absorptionDistance = 0;
        return;
    }

    // --- Get intersection and shading data ---
    const GeometryInstanceID instanceID = getGeometryInstanceID();
    VertexData v = getVertexData(instanceID, triangleIndex, attribs);

    uint materialID = gScene.getMaterialID(instanceID);

    ShadingData sd = gScene.materials.prepareShadingData(v, materialID, -rayDir);
    float3 normal = getFlippedShadingNormal(sd); // Always points against the incident ray

    // Handle non-dielectric materials (fallback to base color texture with optional shadows)
    float matIoR = gScene.materials.evalIoR(sd.materialID);

    if (matIoR <= 0.f) // Considered non-refractive
    {
        BasicMaterialData matData = gScene.materials.getBasicMaterialData(materialID);
        float4 baseColor = gScene.materials.sampleTexture(
            matData.texBaseColor, linearClampSampler, sd.uv, matData.baseColor, ExplicitLodTextureSampler(0) );

        float4 finalColor = baseColor;

        if (drawShadow)
        {
            RayDesc shadowRay;
            shadowRay.Origin = sd.computeRayOrigin(true);
            shadowRay.Direction = normalize(-lightDir);
            shadowRay.TMin = 0.001f;
            shadowRay.TMax = 1000000.f;

            ShadowRayData shadowData = traceShadowRay(shadowRay);

            if (!useRaymarching)
            {
                // --- Beer's Law (absorption) ---
                float3 absorb = exp(-absorptionCoeff * shadowData.absorptionDistance);
                finalColor *= float4(absorb, 1.0f);
            }
            else
            {
                finalColor *= float4(shadowData.color, 1.f);
            }
        }

        hitData.color = finalColor;
        hitData.hitT = currentRayDistance;
        hitData.absorptionDistance = hitData.absorptionDistance + currentRayDistance; ;
        return;
    }

    // --- Dielectric material: reflection + refraction ---
    float etaI = 1.0f; // air
    float etaT = matIoR;
    float eta = sd.frontFacing ? etaI / etaT : etaT / etaI;

    float cosTheta = clamp(dot(normal, -rayDir), 0.00001f, 1.0f);
    float F0 = pow((etaI - etaT) / (etaI + etaT), 2.0f);
    float reflectionRatio = evalFresnelSchlick(F0, 1.0f, cosTheta);
    float refractionRatio = 1.0f - reflectionRatio;

    float3 reflectDir = reflect(rayDir, normal);
    float3 refractDir;
    bool tir = !refractWithTIR(rayDir, normal, eta, refractDir);

    float4 reflectedColor = float4(0, 0, 0, 0);
    float4 refractedColor = float4(0, 0, 0, 0);
    float absorptionTotalDistance = hitData.absorptionDistance;

    if (hitData.depth < maxRayBounce)
    {
        // --- Reflection ray ---
        RayDesc reflectRay;
        reflectRay.Origin = sd.computeRayOrigin(true);
        reflectRay.Direction = reflectDir;
        reflectRay.TMin = 0.001f;
        reflectRay.TMax = 100000.f;

        PrimaryRayData reflectedData = {};
        reflectedData.depth = hitData.depth + 1;
        reflectedData.absorptionDistance = absorptionTotalDistance;
        TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, reflectRay, reflectedData);
        reflectedColor = reflectedData.color;

        if (tir)
        {
            // Total internal reflection, stay in the medium
            absorptionTotalDistance = reflectedData.absorptionDistance + currentRayDistance;
        }
    }

    // --- Refraction ray ---
    if (!tir && hitData.depth < maxRayBounce)
    {
        RayDesc refractRay;
        refractRay.Origin = sd.computeRayOrigin(false);
        refractRay.Direction = refractDir;
        refractRay.TMin = 0.001f;
        refractRay.TMax = 100000.f;

        PrimaryRayData refractedData = {};
        refractedData.depth = hitData.depth + 1;
        refractedData.absorptionDistance = absorptionTotalDistance;
        TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, refractRay, refractedData);
        refractedColor = refractedData.color;

        // Add the travel distance inside the medium
        absorptionTotalDistance = refractedData.absorptionDistance + currentRayDistance;
    }

    if (hitData.depth >= maxRayBounce)
    {
        float3 envColor = gScene.envMap.eval(rayDir, 0);
        hitData.color = float4(envColor, 1);
        hitData.hitT = currentRayDistance;
        //hitData.absorptionDistance = hitData.absorptionDistance + currentRayDistance;
        return;
    }

    // --- Beer's Law (absorption) ---
    float3 absorb = exp(-absorptionCoeff * absorptionTotalDistance);
    refractedColor *= float4(absorb, 1.0f);

    // --- Final color ---
    float4 finalColor = tir ? reflectedColor :
                        reflectionRatio * reflectedColor + refractionRatio * refractedColor;

    hitData.color = finalColor;
    hitData.hitT = currentRayDistance;
    hitData.absorptionDistance = absorptionTotalDistance;
}

[shader("closesthit")]
void RaymarchingClosestHit(inout PrimaryRayData hitData, ImplicitGeometryAttributes attribs)
{
    float3 debugNormal;

    // if (debugVoxelGrid)
    // {
    //     float3 color = float3(float(
    //         gScene.getCustomPrimitiveIndex(getGeometryInstanceID())) / float(voxelGridTotalResolution));
    //     hitData.color = float4(color, 1);
    //     return;
    // }

    if (useRecursiveRaytracing)
    {
        // Position et direction du rayon
        float3 rayOrigin = attribs.pos;
        float3 rayDir = normalize(WorldRayDirection());

        float3 normal = attribs.normal;
        if (dot(normal, rayDir) > 0) normal = -normal;

        // Appliquer l’absorption depuis l’entrée jusqu’à cette surface
        hitData.transmittance *= Transmittance(attribs.densityAlongRay);

        // Indices de réfraction selon sens du trajet
        bool insideFluid = attribs.isInsideFluid;
        float iorA = insideFluid ? IoR : 1.0f;
        float iorB = insideFluid ? 1.0f : IoR;

        LightResponse lightResponse = CalculateReflectionAndRefraction(rayDir, normal, iorA, iorB);

        float3 refractDir = lightResponse.refractDir;
        float3 reflectDir = lightResponse.reflectDir;
        float refractWeight = lightResponse.refractWeight;
        float reflectWeight = lightResponse.reflectWeight;

        // Densités d'absorption sur chaque direction
        float stepSize = sunLightMarchSize * (hitData.depth + 1);
        float densityRefract = CalculateDensityAlongRay(rayOrigin, refractDir, stepSize);
        float densityReflect = CalculateDensityAlongRay(rayOrigin, reflectDir, stepSize);

        float3 colorRefract = float3(0, 0, 0);
        float3 colorReflect = float3(0, 0, 0);

        float3 envReflect = gScene.envMap.eval(reflectDir, 0) * Transmittance(densityReflect);
        float3 envRefract = gScene.envMap.eval(refractDir, 0) * Transmittance(densityRefract);

        bool traceRefract = refractWeight > 0 && (!approximateSecondaryRayBounce || refractWeight >= reflectWeight);
        bool traceReflect = reflectWeight > 0 && (!approximateSecondaryRayBounce || reflectWeight > refractWeight);

        // Rayon réfracté
        if (hitData.depth < maxRayBounce)
        {
            RayDesc refrRay;
            refrRay.Origin = rayOrigin;
            refrRay.Direction = refractDir;
            refrRay.TMin = 0.001f;
            refrRay.TMax = 1e6f;

            PrimaryRayData refrHit = (PrimaryRayData)0;
            refrHit.depth = hitData.depth + 1;
            refrHit.transmittance = hitData.transmittance * refractWeight;

            TraceRay(gScene.rtAccel, RAY_FLAG_NONE, ~0, 0, 1, 0, refrRay, refrHit);
            colorRefract = refrHit.color.rgb * Transmittance(densityRefract);
        }
        // else
        // {
        //     colorRefract = envRefract;
        // }

        // Rayon réfléchi
        if (hitData.depth < maxRayBounce)
        {
            RayDesc reflRay;
            reflRay.Origin = rayOrigin;
            reflRay.Direction = reflectDir;
            reflRay.TMin = 0.001f;
            reflRay.TMax = 1e6f;

            PrimaryRayData reflHit = (PrimaryRayData)0;
            reflHit.depth = hitData.depth + 1;
            reflHit.transmittance = hitData.transmittance * reflectWeight;

            TraceRay(gScene.rtAccel, RAY_FLAG_NONE, ~0, 0, 1, 0, reflRay, reflHit);
            colorReflect = reflHit.color.rgb * Transmittance(densityReflect);
        }
        // else
        // {
        //     colorReflect = envReflect;
        // }

        // Contribution si max depth atteint
        if (hitData.depth >= maxRayBounce)
        {
            float3 rayDir = normalize(WorldRayDirection());
            float3 envColor = gScene.envMap.eval(rayDir, 0);
            hitData.color = float4(envColor, 1);
            hitData.hitT = RayTCurrent();
            return;
        }

        // Combinaison pondérée
        float3 light = reflectWeight * colorReflect + refractWeight * colorRefract;

        hitData.color = float4(debugNormals ? (normal * 0.5f + 0.5f) : light, 1.0);
        hitData.hitT = RayTCurrent();
        return;
    }
    else
    {
        // uint rngState = (uint)(hitData.uv.x * 1243 + hitData.uv.y * 96456);
        uint rngState = 1;

        Ray ray;
        ray.origin = WorldRayOrigin();
        ray.dir = normalize(WorldRayDirection());

        bool travellingThroughFluid = IsInsideFluid(ray.origin, ray.dir);

        float3 transmittance = float3(1, 1, 1);
        float3 light = float3(0, 0, 0);

        for (int i = 0; i < maxRayBounce + 1; i++)
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

float3 SecantSearch(float3 A, float3 B)
{
    float dA = SampleDensityMap(A) - volumeValueOffset;
    float dB = SampleDensityMap(B) - volumeValueOffset;
    float t = dA / (dA - dB); // interpolation linéaire
    return lerp(A, B, t);
}

SurfaceInfo FindNextSurface(float3 origin, float3 rayDir, bool findNextFluidEntryPoint, 
    uint rngState, float rngWeight, float maxDst)
{
    SurfaceInfo info = (SurfaceInfo)0;
    // Safety check
    if (dot(rayDir, rayDir) < 0.5) return info;

    // float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, origin, rayDir);
    //   Use your new matrix-aware ray-box distance function
    float2 boundsDstInfo = RayBoxDst_Matrix(origin, rayDir, localToWorld, worldToLocal);
    float r = (RandomValue(rngState) - 0.5) * marchSize * 0.4 * 1;
    //float r = (1 - 0.5) * marchSize * 0.4 * 1;
    bool hasExittedFluid = !IsInsideFluid(origin, rayDir);
    origin = origin + rayDir * (boundsDstInfo.x + r);

    float stepSize = marchSize;
    bool hasEnteredFluid = false;
    float3 lastPosInFluid = origin;

    float dstToTest = boundsDstInfo[1] - (TinyNudge)*2;
    // float dstToTest = 7.5f;

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
            //info.pos = SecantSearch(lastPosInFluid, samplePos);
            info.foundSurface = true;
            break;
        }
    }

    return info;
}

bool IsInsideFluid(float3 pos, float3 rayDir)
{
    //float2 boundsDstInfo = RayBoxDst(-simBounds * 0.5, simBounds * 0.5, pos, float3(0, 0, 1));
    float2 boundsDstInfo = RayBoxDst_Matrix(pos, rayDir, localToWorld, worldToLocal);

    return (boundsDstInfo.x <= 0 && boundsDstInfo.y > 0) && SampleDensityMap(pos) > 0;
}

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

float3 TraceLightRay(float3 pos, float3 dir)
{
    // 1. Test géométrie opaque
    PrimaryRayData geomHit = (PrimaryRayData)0;
    geomHit.depth = 1;

    RayDesc geometryRay;
    geometryRay.Origin = pos;
    geometryRay.Direction = dir;
    geometryRay.TMin = 0.001f;
    geometryRay.TMax = 1e6f;

    TraceRay(gScene.rtAccel, RAY_FLAG_SKIP_PROCEDURAL_PRIMITIVES, ~0u, 0, 1, 0, geometryRay, geomHit);

    float distToGeom = geomHit.hitT;

    // If there is some geometry, returns the closesthit color, otherwise return the miss color (the envMap).
    return geomHit.color.rgb;
}

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

float SampleSafeDensity(float3 coord) {
    float3 clamped = clamp(coord, float3(0, 0, 0), densityMapSize - 1);
    return SampleDensityMap(clamped);
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

        faceWeight = (1 - smoothstep(0, smoothDst, faceWeight)) * (1 - pow(saturate(volumeNormal.y), smoothPow));

        return normalize(volumeNormal * (1 - faceWeight) + faceNormalWorld * faceWeight);
    }
    else
    {
        float3 uvw = (pos + ScaledSimBounds * 0.5) / ScaledSimBounds;

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

    // if (useVoxelOpti)
    // {
    //     if (cutomPrimitveMasks[gScene.getCustomPrimitiveIndex(getGeometryInstanceID())] == 0)
    //     {
    //         return;
    //     }
    // }

    
    // ImplicitGeometryAttributes attr2 = {};
    // ReportHit(1, 0, attr2);
    // return;

    float3 rayOrigin = WorldRayOrigin();
    float3 rayDir = normalize(WorldRayDirection());

    // Step 1: Intersect with the transformed fluid bounding box (world space)
    CubeHitInfo cubeHitInfo = RayBoxWithMatrix(rayOrigin, rayDir, localToWorld, worldToLocal);

    if (!cubeHitInfo.didHit)
        return;

    if (useRecursiveRaytracing)
    {
        bool travellingThroughFluid = IsInsideFluid(rayOrigin + rayDir * normalOffset * 0.5f, rayDir);

        uint dummyRNG = 0; // RNG state is not useful here unless noise affects surface detection

        // Run the raymarch to find the first fluid surface hit
        SurfaceInfo surfaceInfo = FindNextSurface(rayOrigin, rayDir, !travellingThroughFluid, dummyRNG, 1, cubeHitInfo.dst);

        if (!surfaceInfo.foundSurface)
            return;

        ImplicitGeometryAttributes attr = {};
        attr.pos = surfaceInfo.pos;
        attr.normal = normalize(ComputeDensityNormal(surfaceInfo.pos)); // Don't flip here — let closest hit handle that
        attr.densityAlongRay = surfaceInfo.densityAlongRay;
        attr.isInsideFluid = travellingThroughFluid;

        // Report the intersection to continue in closest hit
        float hitDistance = length(surfaceInfo.pos - rayOrigin); // Accurate distance from ray start
        ReportHit(hitDistance, 0, attr);

        return;
    }
    else
    {
        ImplicitGeometryAttributes attr2 = {};
        attr2.pos = cubeHitInfo.hitPoint;
        ReportHit(1, 0, attr2);
        return;
    }
}

[shader("intersection")]
void AABB_SimulationIntersection()
{
    if (!drawFluid)
        return;

    ImplicitGeometryAttributes attr = {};

    float3 rayOrigin = WorldRayOrigin();
    float3 rayDir = normalize(WorldRayDirection());

    // Step 1: Intersect with the transformed fluid bounding box (world space)
    CubeHitInfo cubeHitInfo = RayBoxWithMatrix(rayOrigin, rayDir, localToWorld, worldToLocal);

    if (!cubeHitInfo.didHit)
        return;

    attr.pos = cubeHitInfo.hitPoint;
    ReportHit(cubeHitInfo.dst, 0, attr);
}

[shader("miss")]
void shadowMiss(inout ShadowRayData hitData)
{
    hitData.hit = false;
}

[shader("anyhit")]
void shadowAnyHit(inout ShadowRayData hitData, in BuiltInTriangleIntersectionAttributes attribs)
{
    float currentT = RayTCurrent();

    uint triangleIndex = PrimitiveIndex();
    GeometryInstanceID instanceID = getGeometryInstanceID();
    VertexData v = getVertexData(instanceID, triangleIndex, attribs);
    uint materialID = gScene.getMaterialID(instanceID);
    ShadingData sd = gScene.materials.prepareShadingData(v, materialID, -WorldRayDirection());

    if (hitData.insideFluid)
    {
        float deltaT = currentT - hitData.lastT;
        hitData.absorptionDistance += deltaT * shadowDensityMultiplier;
    }

    // Mettre à jour l'état
    hitData.insideFluid = sd.frontFacing; // On rentre si on touche face frontale, on sort si c’est une face arrière
    hitData.lastT = currentT;

    IgnoreHit();
}

[shader("closesthit")]
void ShadowClosestHit(inout ShadowRayData hitData, in BuiltInTriangleIntersectionAttributes attribs)
{
    hitData.hit = true; 
}

[shader("anyhit")]
void shadowAnyHitFluid(inout ShadowRayData hitData, in ImplicitGeometryAttributes attribs)
{
    hitData.hit = true;

    float density = CalculateDensityAlongRay(WorldRayOrigin(), WorldRayDirection(), sunLightMarchSize);
    hitData.absorptionDistance = density;
    hitData.color = Transmittance(density * shadowDensityMultiplier);
}

ShadowRayData traceShadowRay(RayDesc ray)
{
    ShadowRayData rayData;
    rayData.hit = false;
    rayData.color = float3(1, 1, 1);
    TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 1, rayTypeCount, 1, ray, rayData);
    return rayData;
}

// Ray generation
[shader("raygeneration")]
void rayGen()
{
    uint3 launchIndex = DispatchRaysIndex();
    uint2 dim = DispatchRaysDimensions().xy;
    float2 uv = launchIndex.xy / (float2)dim;

    RayDesc ray;
    ray = gScene.camera.computeRayPinhole(launchIndex.xy, vector<uint, 2>(viewportDims), false).toRayDesc();

    PrimaryRayData hitData = (PrimaryRayData)0;
    TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, 1, 0, ray, hitData);

    gOutput[launchIndex.xy] = hitData.color;
}