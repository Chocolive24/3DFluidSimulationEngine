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
    int nbParticles;
    float smoothingRadius;
}

enum class BodyType { DYNAMIC, STATIC, FLUID, NONE };

struct ParticleBody {
    float3 Position = float3(0.f, 0.f, 0.f);
    float Density = 1.0f;
    // 16 bytes
    float3 Velocity = float3(0.f, 0.f, 0.f);
    float NearDensity = 1.0f;
    // 16 bytes
    float3 PredictedPosition = float3(0.f, 0.f, 0.f);
    float Pressure = 1.0f;
    // 16 bytes
    float3 Force = float3(0.f, 0.f, 0.f);
    float Mass = 0.f;
    // 16 bytes
    float SmoothingLength = 1.0f;
    float Viscosity = 0.1f;
    float pad1;
    float pad2;
    // 16 bytes
};

RWStructuredBuffer<ParticleBody> bodies;

void ApplyForce(inout ParticleBody Particlebody, float3 force) {
    Particlebody.Force += force;
}

static const float WALLDIST = 100;
static const float PI = 3.14159265359;

void SimulationBounds(inout ParticleBody Particlebody) {

    if (Particlebody.Position.x <= -WALLDIST)
    {
        Particlebody.Velocity.x = 0;
        Particlebody.Position.x = -WALLDIST;
    }
    else if (Particlebody.Position.x >= WALLDIST)
    {
        Particlebody.Velocity.x = 0;
        Particlebody.Position.x = WALLDIST;
    }
    if (Particlebody.Position.y <= -WALLDIST)
    {
        Particlebody.Velocity.y = 0;
        Particlebody.Position.y = -WALLDIST;
    }
    else if (Particlebody.Position.y >= WALLDIST)
    {
        Particlebody.Velocity.y = 0;
        Particlebody.Position.y = WALLDIST;
    }
    if (Particlebody.Position.z <= -WALLDIST)
    {
        Particlebody.Velocity.z = 0;
        Particlebody.Position.z = -WALLDIST;
    }
    else if (Particlebody.Position.z >= WALLDIST)
    {
        Particlebody.Velocity.z = 0;
        Particlebody.Position.z = WALLDIST;
    }
}

[numthreads(64, 1, 1)]
void updateBodies(uint3 threadID: SV_DispatchThreadID)
{
    uint index = threadID.x;

    if (index >= nbParticles) return; // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[index]; // Access the Particlebody at the calculated index

    ApplyForce(particlebody, float3(0, -500, 0));                       // Apply gravity
    float3 acceleration = particlebody.Force * (1 / particlebody.Mass); // Calculate acceleration

    particlebody.Velocity += acceleration * deltaTime;          // Update velocity
    particlebody.Position += particlebody.Velocity * deltaTime; // Update position

    particlebody.Force = float3(0, 0, 0); // Reset force after applying
    SimulationBounds(particlebody);

    bodies[index] = particlebody; // Write the updated particlebody back to the buffer
}

float SmoothingKernel(float radius, float distance)
{
    if (distance >= radius)
    {
        return 0;
    }

    float volume = (PI * pow(radius, 4)) / 6;
    return (radius - distance) * (radius - distance) / volume;
}

[numthreads(64, 1, 1)]
void computeNeighborsDensity(uint3 threadID: SV_DispatchThreadID) 
{
    if (threadID.x >= nbParticles) return;  // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[threadID.x]; // Access the Particlebody at the calculated index

    // std::vector<ParticleBodyRef> neighbors = grid.findNeighbors(data.Position);

    float density = 0;
    for (int i = 0; i < nbParticles; i++) {
        
        if (i == threadID.x)
            continue;
        
        ParticleBody otherParticleBody = bodies[i];

        float dx = particlebody.Position.x - otherParticleBody.Position.x;
        float dy = particlebody.Position.y - otherParticleBody.Position.y;
        float dz = particlebody.Position.z - otherParticleBody.Position.z;

        float distance = sqrt(dx * dx + dy * dy + dz * dz);

        float influence = SmoothingKernel(smoothingRadius, distance);
        density += influence;
    }

    particlebody.Density = density;


    // for (auto &otherRef: neighbors)
    // {

    // }
    bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer
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