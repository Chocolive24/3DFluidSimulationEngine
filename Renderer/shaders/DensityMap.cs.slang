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

RWTexture3D<float> gTexture3D;

StructuredBuffer<float> particleDensities;

cbuffer PerFrameCB
{
    int densityMapSize;
    float3 simBounds;
    float deltaTime;
}

enum class BodyType { DYNAMIC, STATIC, FLUID, NONE };

struct Body {
    float3 Position = float3(0, 0, 0);
    float3 Velocity = float3(0, 0, 0);

    float3 PredictedPosition = float3(0, 0, 0);

    float Mass = -1.f; // Body is disabled if mass is negative
    BodyType Type = BodyType::DYNAMIC;

    float3 Force = float3(0, 0, 0); // Total force acting on the body
};
struct ParticleData
{
    float Density = 1.0f;
    float NearDensity = 1.0f;
    float Pressure = 1.0f;
    float SmoothingLength = 1.0f;
    float Viscosity = 0.1f;
    float3 Position = {0, 0, 0};
    // XMVECTOR Acceleration = XMVectorZero();
};

RWStructuredBuffer<Body> bodies;

[numthreads(8, 8, 8)]
void updateBodies(uint3 threadID: SV_DispatchThreadID)
{
    uint index = threadID.x; // Or flatten (x + y * width + z * width * height) for 3D
    Body body = bodies[index];
    body.Force += ( 0,-500, 0 );

    float3 acceleration = body.Force * (1 / body.Mass);

    body.Velocity += acceleration * deltaTime;
    body.Position += body.Velocity * deltaTime;
    //body.PredictedPosition = XMVectorAdd(body.Position, XMVectorScale(body.Velocity, deltaTime));
    body.Force = (0, 0, 0);
}

[numthreads(8, 8, 8)]
void createDensityMap(uint3 threadID: SV_DispatchThreadID)
{
    if (any(threadID >= densityMapSize)) return;

    float3 texturePos = threadID / float3(densityMapSize - 1);      // normalized [0,1]
    float3 worldPos = (texturePos - 0.5f) * simBounds; // map to [-100,100]

    // float density = 1.f; // world_->CalculateDensityAtPosition(xm_pos) * 100;
    size_t index = threadID.x + threadID.y * densityMapSize + threadID.z * densityMapSize * densityMapSize;
    gTexture3D[threadID] = particleDensities[index];

    //gTexture3D[threadID] = float(threadID.x) / float(densityMapSize);
}