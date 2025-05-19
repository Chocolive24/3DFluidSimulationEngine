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
static const int3 offsets3D[27] = {
{ -1, -1, -1 }, { -1, -1, 0 }, { -1, -1, 1 }, { -1, 0, -1 }, { -1, 0, 0 },  { -1, 0, 1 }, { -1, 1, -1 }, { -1, 1, 0 },   { -1, 1, 1 },
{ 0, -1, -1 }, { 0, -1, 0 },  { 0, -1, 1 },  { 0, 0, -1 }, { 0, 0, 0 }, { 0, 0, 1 },    { 0, 1, -1 },  { 0, 1, 0 },   { 0, 1, 1 },
{ 1, -1, -1 }, { 1, -1, 0 }, { 1, -1, 1 }, { 1, 0, -1 },   { 1, 0, 0 },   { 1, 0, 1 },   { 1, 1, -1 },  { 1, 1, 0 },   { 1, 1, 1 }
};

RWTexture3D<float> gTexture3D;

RWStructuredBuffer<uint3> SpatialIndices; // used for spatial hashing
RWStructuredBuffer<uint> SpatialOffsets;  // used for spatial hashing

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

    // Bitonic Sort.
    uint groupWidth;
    uint groupHeight;
    uint stepIndex;
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

[numthreads(64, 1, 1)]
void spawnParticles(uint3 threadID: SV_DispatchThreadID)
{
    uint index = threadID.x;
    if (index >= nbParticles) return;

    int particlesPerAxis = (int)pow(nbParticles, 1.0 / 3.0) + 1;
    float spacing = 5 * 2 + 5;

    int xIndex = index % particlesPerAxis;
    int yIndex = (index / particlesPerAxis) % particlesPerAxis;
    int zIndex = index / (particlesPerAxis * particlesPerAxis);

    float x = (xIndex - particlesPerAxis / 2.f + 0.5f) * spacing;
    float y = (yIndex - particlesPerAxis / 2.f + 0.5f) * spacing;
    float z = (zIndex - particlesPerAxis / 2.f + 0.5f) * spacing;

    ParticleBody pb;
    pb.Position = float3(x, y, z);

    bodies[index] = pb;
}

static const float PI = 3.14159265359;

// Constants used for hashing
static const uint hashK1 = 15823;
static const uint hashK2 = 9737333;
static const uint hashK3 = 440817757;

// Convert floating point position into an integer cell coordinate
int3 GetCell3D(float3 position, float radius)
{
    return (int3)floor(position / radius);
}

// Hash cell coordinate to a single unsigned integer
uint HashCell3D(int3 cell)
{
    cell = (uint3)cell;
    return (cell.x * hashK1) + (cell.y * hashK2) + (cell.z * hashK3);
}

uint KeyFromHash(uint hash, uint tableSize)
{
    return hash % tableSize;
}

[numthreads(64, 1, 1)]
void UpdateSpatialHash(uint3 id: SV_DispatchThreadID)
{
    if (id.x >= nbParticles) return;

    // Reset offsets
    SpatialOffsets[id.x] = nbParticles;
    // Update index buffer
    uint index = id.x;
    int3 cell = GetCell3D(bodies[index].Position, smoothingRadius);
    uint hash = HashCell3D(cell);
    uint key = KeyFromHash(hash, nbParticles);
    SpatialIndices[id.x] = uint3(index, hash, key);
}

// Sort the given entries by their keys (smallest to largest)
// This is done using bitonic merge sort, and takes multiple iterations
[numthreads(128, 1, 1)]
void BitonicSort(uint3 id: SV_DispatchThreadID)
{
    uint i = id.x;

    uint hIndex = i & (groupWidth - 1);
    uint indexLeft = hIndex + (groupHeight + 1) * (i / groupWidth);
    uint rightStepSize = stepIndex == 0 ? groupHeight - 2 * hIndex : (groupHeight + 1) / 2;
    uint indexRight = indexLeft + rightStepSize;

    // Exit if out of bounds (for non-power of 2 input sizes)
    if (indexRight >= nbParticles) return;

    uint valueLeft = SpatialIndices[indexLeft].z; // z = key
    uint valueRight = SpatialIndices[indexRight].z;

    // Swap entries if value is descending
    if (valueLeft > valueRight)
    {
        uint3 temp = SpatialIndices[indexLeft];
        SpatialIndices[indexLeft] = SpatialIndices[indexRight];
        SpatialIndices[indexRight] = temp;
    }
}

// Calculate offsets into the sorted Entries buffer (used for spatial hashing).
// For example, given an Entries buffer sorted by key like so: {2, 2, 2, 3, 6, 6, 9, 9, 9, 9}
// The resulting Offsets calculated here should be:            {-, -, 0, 3, -, -, 4, -, -, 6}
// (where '-' represents elements that won't be read/written)
//
// Usage example:
// Say we have a particular particle P, and we want to know which particles are in the same grid cell as it.
// First we would calculate the Key of P based on its position. Let's say in this example that Key = 9.
// Next we can look up Offsets[Key] to get: Offsets[9] = 6
// This tells us that SortedEntries[6] is the first particle that's in the same cell as P.
// We can then loop until we reach a particle with a different cell key in order to iterate over all the particles in the cell.
//
// NOTE: offsets buffer must filled with values equal to (or greater than) its length to ensure that this works correctly
[numthreads(64, 1, 1)]
void CalculateOffsets(uint3 id: SV_DispatchThreadID)
{
    if (id.x >= nbParticles) { return; }

    uint i = id.x;
    uint null = nbParticles;

    uint key = SpatialIndices[i].z; // .z = key
    uint keyPrev = i == 0 ? null : SpatialIndices[i - 1].z; // .z = key.

    if (key != keyPrev)
    {
        SpatialOffsets[key] = i;
    }
}

void ApplyForce(inout ParticleBody Particlebody, float3 force) {
    Particlebody.Force += force;
}

void SimulationBounds(inout ParticleBody Particlebody) 
{
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

    int3 originCell = GetCell3D(particlebody.Position, smoothingRadius);
    float density = 0;
    
    // Neighbour search.
    for (int i = 0; i < 27; i++)
    {
        uint hash = HashCell3D(originCell + offsets3D[i]);
        uint key = KeyFromHash(hash, nbParticles);
        uint currIndex = SpatialOffsets[key];

        while (currIndex < nbParticles)
        {
            uint3 indexData = SpatialIndices[currIndex];
            currIndex++;
            // Exit if no longer looking at correct bin
            if (indexData[2] != key) break;
            // Skip if hash does not match
            if (indexData[1] != hash) continue;

            uint neighbourIndex = indexData[0];

            ParticleBody otherParticleBody = bodies[neighbourIndex];

            float dx = particlebody.Position.x - otherParticleBody.Position.x;
            float dy = particlebody.Position.y - otherParticleBody.Position.y;
            float dz = particlebody.Position.z - otherParticleBody.Position.z;

            float distance = sqrt(dx * dx + dy * dy + dz * dz);

            float influence = SmoothingKernel(smoothingRadius, distance);
            density += influence;
        
            // float3 neighbourPos = bodies[neighbourIndex].Position;
            // float3 offsetToNeighbour = neighbourPos - particlebody.Position;
            // float sqrDstToNeighbour = dot(offsetToNeighbour, offsetToNeighbour);

            // // Skip if not within radius
            // if (sqrDstToNeighbour > sqrRadius) continue;

            // // Calculate density and near density
            // float dst = sqrt(sqrDstToNeighbour);
            // density += DensityKernel(dst, smoothingRadius);
            // nearDensity += NearDensityKernel(dst, smoothingRadius);
        }
    }

    // particlebody.Density = density;
    // bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer

    // float density = 0;
    // for (int i = 0; i < nbParticles; i++) {
        
    //     // if (i == threadID.x)
    //     //     continue;
        
    //     ParticleBody otherParticleBody = bodies[i];

    //     float dx = particlebody.Position.x - otherParticleBody.Position.x;
    //     float dy = particlebody.Position.y - otherParticleBody.Position.y;
    //     float dz = particlebody.Position.z - otherParticleBody.Position.z;

    //     float distance = sqrt(dx * dx + dy * dy + dz * dz);

    //     float influence = SmoothingKernel(smoothingRadius, distance);
    //     density += influence;
    // }

    particlebody.Density = density; //max(density, 0.01f);

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

    int3 originCell = GetCell3D(particlebody.Position, smoothingRadius);
    float3 pressureForce = float3(0.f, 0.f, 0.f);

    // Neighbour search
    for (int i = 0; i < 27; i++)
    {
        uint hash = HashCell3D(originCell + offsets3D[i]);
        uint key = KeyFromHash(hash, nbParticles);
        uint currIndex = SpatialOffsets[key];

        while (currIndex < nbParticles)
        {
            uint3 indexData = SpatialIndices[currIndex];
            currIndex++;
            // Exit if no longer looking at correct bin
            if (indexData[2] != key) break;
            // Skip if hash does not match
            if (indexData[1] != hash) continue;

            uint neighbourIndex = indexData[0];
            // Skip if looking at self
            if (neighbourIndex == threadID.x) continue;

            ParticleBody otherParticleBody = bodies[neighbourIndex];

            float3 offset = particlebody.Position - otherParticleBody.Position;
            float distance = length(offset);
            float3 direction = distance == 0 ? float3(0.f, 1.f, 0.f) : offset / distance;
            float slope = SmoothingKernelDerivative(smoothingRadius, distance);
            float otherDensity = otherParticleBody.Density;
            float sharedPressure = CalculateSharedPressure(otherDensity, particlebody.Density);
            pressureForce += sharedPressure * direction * slope * otherParticleBody.Mass / otherDensity;
        }
    }

    // for (int i = 0; i < nbParticles; i++)
    // {
    //     if (i == threadID.x)
    //         continue;

    //     ParticleBody otherParticleBody = bodies[i];
    //     float3 offset = particlebody.Position - otherParticleBody.Position;
    //     float distance = length(offset);
    //     float3 direction = distance == 0 ? float3(0.f, 1.f, 0.f) : offset / distance;
    //     float slope = SmoothingKernelDerivative(smoothingRadius, distance);
    //     float otherDensity = otherParticleBody.Density;
    //     float sharedPressure = CalculateSharedPressure(otherDensity, particlebody.Density);
    //     pressureForce += sharedPressure * direction * slope * otherParticleBody.Mass / otherDensity;
    // }

    ApplyForce(particlebody, pressureForce / particlebody.Density);
   
    bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer
}

[numthreads(64, 1, 1)]
void computeNeighborsViscosity(uint3 threadID: SV_DispatchThreadID)
{
    if (threadID.x >= nbParticles) return; // Ensure we do not go out of bounds

    ParticleBody particlebody = bodies[threadID.x]; // Access the Particlebody at the calculated index

    int3 originCell = GetCell3D(particlebody.Position, smoothingRadius);
    float3 viscosityForce = float3(0, 0, 0);

    // Neighbour search
    for (int i = 0; i < 27; i++)
    {
        uint hash = HashCell3D(originCell + offsets3D[i]);
        uint key = KeyFromHash(hash, nbParticles);
        uint currIndex = SpatialOffsets[key];

        while (currIndex < nbParticles)
        {
            uint3 indexData = SpatialIndices[currIndex];
            currIndex++;
            // Exit if no longer looking at correct bin
            if (indexData[2] != key) break;
            // Skip if hash does not match
            if (indexData[1] != hash) continue;

            uint neighbourIndex = indexData[0];
            // Skip if looking at self
            if (neighbourIndex == threadID.x) continue;

            ParticleBody otherParticleBody = bodies[neighbourIndex];

            float3 r = particlebody.Position - otherParticleBody.Position;
            float distance = length(r);

            if (distance > smoothingRadius) continue;

            float influence = SmoothingKernel(smoothingRadius, distance);
            viscosityForce += (otherParticleBody.Velocity - particlebody.Velocity) * influence * particlebody.Mass;
        }
    }

    // float3 viscosityForce = float3(0, 0, 0);

    // for (int i = 0; i < nbParticles; i++)
    // {
    //     if (i == threadID.x)
    //         continue;

    //     ParticleBody otherParticleBody = bodies[i];

    //     float3 r = particlebody.Position - otherParticleBody.Position;
    //     float distance = length(r);

    //     if (distance > smoothingRadius) continue;

    //     float influence = SmoothingKernel(smoothingRadius, distance);
    //     viscosityForce += (otherParticleBody.Velocity - particlebody.Velocity) * influence * particlebody.Mass;
    // }

    ApplyForce(particlebody, viscosityForce * viscosityStrength);

    bodies[threadID.x] = particlebody; // Write the updated Particlebody back to the buffer
}

float CalculateDensityAtPosition(float3 posiiton)
{
    int3 originCell = GetCell3D(posiiton, smoothingRadius);

    float density = 0;

    // Neighbour search
    for (int i = 0; i < 27; i++)
    {
        uint hash = HashCell3D(originCell + offsets3D[i]);
        uint key = KeyFromHash(hash, nbParticles);
        uint currIndex = SpatialOffsets[key];

        while (currIndex < nbParticles)
        {
            uint3 indexData = SpatialIndices[currIndex];
            currIndex++;
            // Exit if no longer looking at correct bin
            if (indexData[2] != key) break;
            // Skip if hash does not match
            if (indexData[1] != hash) continue;

            uint neighbourIndex = indexData[0];

            ParticleBody otherParticleBody = bodies[neighbourIndex];

            float dx = posiiton.x - otherParticleBody.Position.x;
            float dy = posiiton.y - otherParticleBody.Position.y;
            float dz = posiiton.z - otherParticleBody.Position.z;

            float distance = sqrt(dx * dx + dy * dy + dz * dz);

            //if (distance > smoothingRadius) continue;

            float influence = SmoothingKernel(smoothingRadius, distance);
            density += influence;
        }
    }

    // for (int i = 0; i < nbParticles; i++)
    // {
    //     // if (i == threadID.x)
    //     //     continue;

    //     ParticleBody otherParticleBody = bodies[i];

    //     float dx = posiiton.x - otherParticleBody.Position.x;
    //     float dy = posiiton.y - otherParticleBody.Position.y;
    //     float dz = posiiton.z - otherParticleBody.Position.z;

    //     float distance = sqrt(dx * dx + dy * dy + dz * dz);

    //     if (distance > smoothingRadius) continue;

    //     float influence = SmoothingKernel(smoothingRadius, distance);
    //     density += influence;
    // }

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