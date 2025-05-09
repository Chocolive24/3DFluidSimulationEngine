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
static const int3 offsets[27] = {
{ -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 }, { -1, 0, -1 }, { -1, 0, 0 },  { -1, 0, 1 }, { -1, 1, -1 }, { -1, 1, 0 },   { -1, 1, 1 },
{ 0, -1, -1 }, { 0, -1, 0 },  { 0, -1, 1 },  { 0, 0, -1 }, { 0, 0, 0 }, { 0, 0, 1 },    { 0, 1, -1 },  { 0, 1, 0 },   { 0, 1, 1 },
{ 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, -1 },   { 1, 0, 0 },   { 1, 0, 1 },   { 1, 1, -1 },  { 1, 1, 0 },   { 1, 1, 1 }
};

RWTexture3D<float> gTexture3D;

cbuffer PerFrameCB
{
    int densityMapSize;
    float3 simBounds;
    float wallDist;
    float fixedDeltaTime;
    int nbParticles;
    float gravity;
    float smoothingRadius;
    float targetDensity;
    float pressureMultiplier;
    float viscosityStrength;
    float densityGraphicsMultiplier;
}

struct ParticleBody 
{
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
    float pad1 = -1.f;
    float pad2 = -1.f;
    // 16 bytes
};

RWStructuredBuffer<ParticleBody> bodies;

static const float PI = 3.14159265359;

struct GridCell
{
    float3 position = float3(0, 0, 0);
    RWStructuredBuffer<ParticleBody> bodies;
}

void ApplyForce(inout ParticleBody Particlebody, float3 force) {
    Particlebody.Force += force;
}

void SimulationBounds(inout ParticleBody Particlebody) {
    if (Particlebody.Position.x <= -wallDist)
    {
        Particlebody.Velocity.x = abs(Particlebody.Velocity.x);
        //Particlebody.Position.x = -WALLDIST;
    }
    else if (Particlebody.Position.x >= wallDist)
    {
        Particlebody.Velocity.x = -abs(Particlebody.Velocity.x);
        //Particlebody.Position.x = WALLDIST;
    }
    if (Particlebody.Position.y <= -wallDist)
    {
        Particlebody.Velocity.y = abs(Particlebody.Velocity.y);
        Particlebody.Position.y = -wallDist;
    }
    else if (Particlebody.Position.y >= wallDist)
    {
        Particlebody.Velocity.y = -abs(Particlebody.Velocity.y);
        //Particlebody.Position.y = WALLDIST;
    }
    if (Particlebody.Position.z <= -wallDist)
    {
        Particlebody.Velocity.z = abs(Particlebody.Velocity.z);
        //Particlebody.Position.z = -WALLDIST;
    }
    else if (Particlebody.Position.z >= wallDist)
    {
        Particlebody.Velocity.z = -abs(Particlebody.Velocity.z);
        //Particlebody.Position.z = WALLDIST;
    }
}

[numthreads(64, 1, 1)]
void updateBodies(uint3 threadID: SV_DispatchThreadID)
{
    uint index = threadID.x;

    if (index >= nbParticles) return; // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[index]; // Access the Particlebody at the calculated index

    ApplyForce(particlebody, float3(0, -gravity, 0));                       // Apply gravity
    float3 acceleration = particlebody.Force * (1 / particlebody.Mass); // Calculate acceleration

    particlebody.Velocity += acceleration * fixedDeltaTime;          // Update velocity
    particlebody.Position += particlebody.Velocity * fixedDeltaTime; // Update position

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

float SmoothingKernelDerivative(float radius, float distance)
{
    if (distance >= radius)
    {
        return 0;
    }

    float scale = 12 / (pow(radius, 4) * PI);
    return (distance - radius) * scale;
}

[numthreads(64, 1, 1)]
void computeNeighborsDensity(uint3 threadID: SV_DispatchThreadID) 
{
    if (threadID.x >= nbParticles) return;  // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[threadID.x]; // Access the Particlebody at the calculated index

    // std::vector<ParticleBodyRef> neighbors = grid.findNeighbors(data.Position);

    float density = 0;
    for (int i = 0; i < nbParticles; i++) {
        
        // if (i == threadID.x)
        //     continue;
        
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

float ConvertDensityToPressure(float density)
{
    float densityError = density - targetDensity;
    float pressure = densityError * pressureMultiplier;

    return pressure;
}

float CalculateSharedPressure(float pressure1, float pressure2)
{
    return (ConvertDensityToPressure(pressure1) + ConvertDensityToPressure(pressure2)) * 0.5f;
}

[numthreads(64, 1, 1)]
void computeNeighborsPressure(uint3 threadID: SV_DispatchThreadID)
{
    if (threadID.x >= nbParticles) return; // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[threadID.x]; // Access the Particlebody at the calculated index

    // std::vector<BodyRef> neighbors = grid.findNeighbors(data.Position);

    float3 pressureForce = float3(0.f, 0.f, 0.f);

    for (int i = 0; i < nbParticles; i++)
    {
        if (i == threadID.x)
            continue;

        ParticleBody otherParticleBody = bodies[i];
        float3 offset = particlebody.Position - otherParticleBody.Position;
        float distance = length(offset);
        float3 direction = distance == 0 ? float3(0.f, 1.f, 0.f) : offset / distance;
        float slope = SmoothingKernelDerivative(smoothingRadius, distance);
        float otherDensity = otherParticleBody.Density;
        float sharedPressure = CalculateSharedPressure(otherDensity, particlebody.Density);
        pressureForce += sharedPressure * direction * slope * otherParticleBody.Mass / otherDensity;
    }

    ApplyForce(particlebody, pressureForce / particlebody.Density);
   
    bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer
}

[numthreads(64, 1, 1)]
void computeNeighborsViscosity(uint3 threadID: SV_DispatchThreadID)
{
    if (threadID.x >= nbParticles) return; // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[threadID.x]; // Access the Particlebody at the calculated index

    //std::vector<BodyRef> neighbors = grid.findNeighbors(data.Position);

    float3 viscosityForce = float3(0, 0, 0);

    for (int i = 0; i < nbParticles; i++)
    {
        if (i == threadID.x)
            continue;

        ParticleBody otherParticleBody = bodies[i];

        float3 r = particlebody.Position - otherParticleBody.Position;
        float distance = length(r);

        if (distance > smoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, distance);
        viscosityForce += (otherParticleBody.Velocity - particlebody.Velocity) * influence * particlebody.Mass;
    }

    ApplyForce(particlebody, viscosityForce * viscosityStrength);

    bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer
}

float CalculateDensityAtPosition(float3 posiiton)
{
    //return 0.5f;
    //const std::vector<BodyRef> neighbors = grid.findNeighbors(position);

    float density = 0;
    for (int i = 0; i < nbParticles; i++)
    {
        // if (i == threadID.x)
        //     continue;

        ParticleBody otherParticleBody = bodies[i];

        float dx = posiiton.x - otherParticleBody.Position.x;
        float dy = posiiton.y - otherParticleBody.Position.y;
        float dz = posiiton.z - otherParticleBody.Position.z;

        float distance = sqrt(dx * dx + dy * dy + dz * dz);

        if (distance > smoothingRadius) continue;

        float influence = SmoothingKernel(smoothingRadius, distance);
        density += influence;
    }

    density *= densityGraphicsMultiplier;

    return density;
}

[numthreads(8, 8, 8)]
void computeDensityMap(uint3 threadID: SV_DispatchThreadID)
{
    if (any(threadID >= densityMapSize)) return;

    float3 texturePos = threadID / float3(densityMapSize - 1);      // normalized [0,1]
    float3 worldPos = (texturePos - 0.5f) * simBounds; // map to [-100,100]

    float density = CalculateDensityAtPosition(worldPos);
    gTexture3D[threadID] = density;
}