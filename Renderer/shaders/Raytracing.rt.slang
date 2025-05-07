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
SamplerState s;

cbuffer PerFrameCB
{
    float4x4 invView;
    float2 viewportDims;
    float tanHalfFovY;
    uint sampleIndex;
    bool useDOF;
    bool drawFluid;
    float3 backgroundColor;
    float waterTurbulence;
    uint maxRayBounce;
    float3 absorptionCoeff;
    float3 scatteringCoeff;
    float phaseG;
    float maxRaymarchingDistance;
    float marchSize;
    float maxLighMarchingDistance;
    float sunLightMarchSize;
    float3 lightColor;
    float3 lightDir;
    float IoR;
    float time;
    int iFrame;
    float DensityDepth;
    int densityMapSize;
};



struct SDF_GeometryAttributes
{
    float3 w_position;
    float3 w_normal;
    float t;
}

struct PrimaryRayData
{
    float4 color = float4(0, 0, 0, 0);
    uint depth = 0;
    float hitT = -1.f;
    float absorptionDistance = 0.f;
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
        float3 position = attribs.w_position + rayDir * t;
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
    //hitData.color = hasHit ? float4(0.1, 0.3, 0.8, 1.0) : float4(envColor, 1.f);
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

[shader("closesthit")]
void RaymarchingClosestHit(inout PrimaryRayData hitData, SDF_GeometryAttributes attribs)
{
    hitData.color = float4(float3(0.1f, 0.3f, 0.8f), 0);
    hitData.hitT = RayTCurrent();
    return;

    if (hitData.depth >= maxRayBounce)
    {
        hitData.color = float4(backgroundColor, 0);
        hitData.hitT = 0.f;
        return;
    }

    float3 rayDir = WorldRayDirection();
    float currentRayT = RayTCurrent();

    // 1. Compute Fresnel reflectance
    float etaI = 1.f;                                // Incident IoR (Air by default).
    float etaT = IoR;                                // Tranmission IoR (the one from the medium).
    bool isEntering = dot(rayDir, attribs.w_normal) > 0.0;
    float eta = isEntering ? etaI / etaT : etaT / etaI;
    float3 flippedShadingNormal = isEntering ? attribs.w_normal : -attribs.w_normal;

    float F0 = pow((etaI - etaT) / (etaI + etaT), 2.0);
    float F90 = 1.f;
    float cosTheta = max(0.00001f, dot(flippedShadingNormal, rayDir));
    float reflectionRatio = evalFresnelSchlick(F0, F90, cosTheta);
    float refractionRatio = 1.f - reflectionRatio;

    // 2. Prepare color storage
    float3 reflectedColor = float3(0, 0, 0);
    float3 refractedColor = float3(0, 0, 0);

    // 3. Check for total internal reflection
    bool tir = false;
    float3 refractDir = rayDir; // Default if eta == 1.0

    if (eta != 1.0)
    {
        if (!refractWithTIR(rayDir, flippedShadingNormal, eta, refractDir))
            tir = true;
    }

    // 4. Trace reflection ray if needed
    if (reflectionRatio > 0.001f || tir)
    {
        // Handle total internal reflection
        // float3 reflectDir = reflect(rayDir, attribs.w_normal); // Compute the reflection direction

        // // Define reflection ray
        // RayDesc reflectRay;
        // reflectRay.Origin = attribs.w_position + reflectDir * 0.1f; // Offset to avoid self-intersection
        // reflectRay.Direction = reflectDir;
        // reflectRay.TMin = 0.001f;
        // reflectRay.TMax = 100000.f;

        // // Trace reflection ray
        // PrimaryRayData reflectedRayData = {};
        // reflectedRayData.depth = hitData.depth + 1;
        // reflectedRayData.hitT = 0;
        // TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, reflectRay, reflectedRayData);

        // reflectedColor = reflectedRayData.color.rgb;
    }

    // 5. Trace refraction ray if not TIR
    if (!tir && refractionRatio > 0.001f)
    {
        refractedColor = CalculateTransmittanceRaymarching(attribs, refractDir);

        // // Trace a refracted ray to get the color after the transmission.
        // RayDesc refractRay;
        // refractRay.Origin = attribs.w_position + refractDir * 0.001f;
        // refractRay.Direction = refractDir;
        // refractRay.TMin = 0.001f;
        // refractRay.TMax = 100000.f;

        // PrimaryRayData refractedRayData = {};
        // refractedRayData.depth = hitData.depth + 1;
        // refractedRayData.hitT = 0.f;
        // refractedRayData.absorptionDistance = hitData.absorptionDistance;
        // TraceRay(gScene.rtAccel, RAY_FLAG_NONE, 0xFF, 0, rayTypeCount, 0, refractRay, refractedRayData);

        // refractedColor = refractedRayData.color.rgb;
    }

    //float3 extinctionCoeff = absorptionCoeff + scatteringCoeff;

    // 6. Mix results
    float3 totalLight = tir ? reflectedColor :
            reflectionRatio * reflectedColor + refractionRatio * refractedColor;

    //float3 totalLight = CalculateTransmittanceRaymarching(attribs, rayDir);

    reflectedColor *= reflectionRatio;
    refractedColor = float3(0, 0, 1) * refractionRatio;

    float3 finalColor = reflectedColor + refractedColor;

    hitData.color = float4(totalLight, 1);
    hitData.hitT = currentRayT;
}

[shader("intersection")]
void RaymarchingIntersection()
{
    float3 rayDir = WorldRayDirection();
    float3 rayOrigin = WorldRayOrigin();

    SDF_GeometryAttributes attr;

    float3 finalColor = float3(0, 0, 0);

    float t = RayTMin();
    float tMax = RayTCurrent();
    float distInsideAABB = abs(tMax - t);
    int nbrOfRaymarchingSteps = 50;
    float stepSize = distInsideAABB / nbrOfRaymarchingSteps;

    for (int i = 0; i < nbrOfRaymarchingSteps; i++)
    {
        float3 p = rayOrigin + rayDir * t;

        float halfDensityMapSize = densityMapSize * 0.5f;

        // Define the world-space bounds of the volume
        float3 volumeMin = float3(-halfDensityMapSize);
        float3 volumeMax = float3(halfDensityMapSize);

        // Normalize position `p` to [0,1] texture coordinates (UVW)
        float3 uvw = (p - volumeMin) / (volumeMax - volumeMin); // Don't saturate yet!

        // Get texture dimensions (assuming cube texture)
        int3 volumeDim = int3(densityMapSize, densityMapSize, densityMapSize);

        // Convert UVW to texel coordinates
        int3 texCoord = int3(uvw * (volumeDim - 1));

        // Clamp to avoid out-of-bounds access
        texCoord = clamp(texCoord, int3(0, 0, 0), volumeDim - 1);

        // Sample using Load()
        float density = gTexture3D.Load(int4(texCoord, 0));
        //float density = gTexture3D.Sample(s, p);

        // Example early exit
        if (density > DensityDepth)
        {
            attr.w_position = p;
            //attr.w_normal = sdCalculateNormal(p, waterTurbulence, time, iFrame);
            attr.t = max(t, 0.001f);

            ReportHit(t, /*hitKind*/ 0, attr);
            break;
        }

        t += stepSize;
    }

    // Raymarching.
    // for (int i = 0; i < 80 && t < RayTCurrent(); i++)
    // {
    //     float3 p = rayOrigin + rayDir * t; // Position along the ray.

    //     float d = map(p, waterTurbulence, time, iFrame); // Current distance to the scene.

    //     // Early stop if close enough or if too far.
    //     if (d < 0.001f)
    //     {
    //         attr.w_position = p;
    //         attr.w_normal = sdCalculateNormal(p, waterTurbulence, time, iFrame);
    //         attr.t = max(t, 0.001f);

    //         ReportHit(t, /*hitKind*/ 0, attr);
    //         break;
    //     }

    //     t += max(d, 0.01f); // March the ray.

    //     if (t > 100.f)
    //     {
    //         break;
    //     }
    // }
}

// Ray generation (just a simple ray going through the scene)
[shader("raygeneration")]
void rayGen()
{
    uint3 launchIndex = DispatchRaysIndex();

    RayDesc ray;
    ray = gScene.camera.computeRayPinhole(launchIndex.xy, vector<uint, 2>(viewportDims)).toRayDesc();

    PrimaryRayData hitData;
    hitData.color = float4(0.f, 0.f, 0.f, 0.f);
    hitData.depth = 0;
    TraceRay(gScene.rtAccel, RAY_FLAG_NONE /*rayFlags*/, 0xFF, /*Test all Geometry*/ 0 /* ray index*/, rayTypeCount, 0, ray, hitData);
    
    // if (hitData.hitT > 0)
    // {
    //     float3 worldPos = ray.Origin + hitData.hitT * ray.Direction;

    //     float halfDensityMapSize = densityMapSize * 0.5f;

    //     // Define your volume bounds (must match what you used in CPU to generate the texture)
    //     float3 volumeMin = float3(-halfDensityMapSize, -halfDensityMapSize, -halfDensityMapSize);
    //     float3 volumeMax = float3(halfDensityMapSize, halfDensityMapSize, halfDensityMapSize);

    //     // Normalize to [0,1]
    //     float3 uvw = saturate((worldPos - volumeMin) / (volumeMax - volumeMin));

    //     // Option 1: sample if allowed
    //     // float density = gTexture3D.SampleLevel(gSampler, uvw, 0.0);

    //     // Option 2: Load (lib_6_6 compatible)
    //     float3 volumeDim = float3(densityMapSize, densityMapSize, densityMapSize);
    //     int3 texCoord = int3(uvw * float3(volumeDim)); // volumeDim is (width, height, depth)
    //     texCoord = clamp(texCoord, int3(0), volumeDim - 1);

    //     //texCoord.z = DensityDepth;

    //     float density = gTexture3D.Load(int4(texCoord, 0));

    //     gOutput[launchIndex.xy] = float4(density, density, density, 1); // or use for absorption, etc.
    // }
    // else 
    // {
    //     gOutput[launchIndex.xy] = hitData.color;
    // }

    gOutput[launchIndex.xy] = hitData.color;
}

