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
import MarchTables;

struct Vertex {
    float3 position = float3(-20000.f, 0.f, 0.f);
    // float3 normal;
    // int2 id;
};

struct Triangle {
    Vertex vertexA;
    Vertex vertexB;
    Vertex vertexC;
};

AppendStructuredBuffer<Triangle> triangles;
RWTexture3D<float> DensityTexture;

cbuffer PerFrameCB
{
    int numPointsPerAxis;
    int isoLevel;
    int textureSize;
    float boundSize;
    float SphereRadius;
}



[numthreads(8, 8, 8)]
void ComputeDensityTexture(uint3 threadID: SV_DispatchThreadID)
{
    float3 worldPos = (threadID / (textureSize - 1.f) - 0.5f) * boundSize;
    float dstFromWorldCenter = length(worldPos);
    float mapValue = dstFromWorldCenter - SphereRadius;

    DensityTexture[threadID] = mapValue;
}

float sampleDensity(int3 position)
{
    // Bounds check
    if (any(position < 0) || any(position >= textureSize))
        return 0;

    return DensityTexture.Load(position); // mip level 0
}

[numthreads(8, 8, 8)]
void ProcessCube(uint3 threadID: SV_DispatchThreadID)
{
    int numCubesPerAxis = numPointsPerAxis - 1;

    if (any(threadID) >= numCubesPerAxis) return;

    int3 cornerCoords[8];
    cornerCoords[0] = threadID + int3(0, 0, 0);
    cornerCoords[1] = threadID + int3(1, 0, 0);
    cornerCoords[2] = threadID + int3(1, 0, 1);
    cornerCoords[3] = threadID + int3(0, 0, 1);
    cornerCoords[4] = threadID + int3(0, 1, 0);
    cornerCoords[5] = threadID + int3(1, 1, 0);
    cornerCoords[6] = threadID + int3(1, 1, 1);
    cornerCoords[7] = threadID + int3(0, 1, 1);

    // Calculate unique index for each cube configuration.
    // There are 256 possible values (cube has 8 corners, so 2^8 possibilites).
    // A value of 0 means cube is entirely inside the surface; 255 entirely outside.
    // The value is used to look up the edge table, which indicates which edges of the cube the surface passes through.
    int cubeConfiguration = 0;
    for (int i = 0; i < 8; i++) {
        // Think of the configuration as an 8-bit binary number (each bit represents the state of a corner point).
        // The state of each corner point is either 0: above the surface, or 1: below the surface.
        // The code below sets the corresponding bit to 1, if the point is below the surface.
        if (sampleDensity(cornerCoords[i]) < isoLevel) {
            cubeConfiguration |= (1 << i);
        }
    }

    // Get array of the edges of the cube that the surface passes through.
    int edgeIndices[16] = triangulation[cubeConfiguration];

    // Create triangles for the current cube configuration
    for (int i = 0; i < 16; i += 3) {
        // If edge index is -1, then no further vertices exist in this configuration
        if (edgeIndices[i] == -1) { break; }

        // Get indices of the two corner points defining the edge that the surface passes through.
        // (Do this for each of the three edges we're currently looking at).
        int edgeIndexA = edgeIndices[i];
        int a0 = cornerIndexAFromEdge[edgeIndexA];
        int a1 = cornerIndexBFromEdge[edgeIndexA];

        int edgeIndexB = edgeIndices[i + 1];
        int b0 = cornerIndexAFromEdge[edgeIndexB];
        int b1 = cornerIndexBFromEdge[edgeIndexB];

        int edgeIndexC = edgeIndices[i + 2];
        int c0 = cornerIndexAFromEdge[edgeIndexC];
        int c1 = cornerIndexBFromEdge[edgeIndexC];

        // Calculate positions of each vertex.
        Vertex vertexA;
        Vertex vertexB;
        Vertex vertexC;
        vertexA.position = (cornerCoords[a0] + cornerCoords[a1]) * 0.5f;
        vertexB.position = (cornerCoords[b0] + cornerCoords[b1]) * 0.5f;
        vertexC.position = (cornerCoords[c0] + cornerCoords[c1]) * 0.5f;

        // Create triangle
        Triangle tri;
        tri.vertexA = vertexA;
        tri.vertexB = vertexB;
        tri.vertexC = vertexC;
        triangles.Append(tri);
    }
}