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
    float3 Position = float3(0.f, 0.f, 0.f); float pad0;
    float3 Velocity = float3(0.f, 0.f, 0.f); float pad1;
    float3 PredictedPosition = float3(0.f, 0.f, 0.f); float pad2;
    float Mass = 0.f;
    uint Type = 0; float pad3[2]; // Align to 16 bytes total
    float3 Force = float3(0.f, 0.f, 0.f); float pad4;
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

[numthreads(8, 8, 1)]
void updateBodies(uint3 threadID: SV_DispatchThreadID, uint3 groupID: SV_GroupID)
{
    uint index = threadID.x + threadID.y * 8; // Local thread index in the group (0-63)
    index += groupID.x * 64;                  // Global thread index across groups (0-1023)

    if (index >= 1024) return; // Ensure we do not go out of bounds

    Body body = bodies[index]; // Access the body at the calculated index

    body.Force += float3(0, -5, 0); // Apply gravity or force

    float3 acceleration = body.Force * (1 / body.Mass); // Calculate acceleration

    body.Velocity += acceleration * deltaTime;  // Update velocity
    body.Position += body.Velocity * deltaTime; // Update position

    body.Force = float3(0, 0, 0); // Reset force after applying

    bodies[index] = body; // Write the updated body back to the buffer
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