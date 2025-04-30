module SDF_Functions;

// --------------------------------------------//
//               Noise Functions
// --------------------------------------------//

// Taken from Inigo Quilez's Rainforest ShaderToy:
// https://www.shadertoy.com/view/4ttSWf
float hash1( float n )
{
    return fract( n*17.0*fract( n*0.3183099 ) );
}

// Taken from Inigo Quilez's Rainforest ShaderToy:
// https://www.shadertoy.com/view/4ttSWf
float noise( in float3 x )
{
    float3 p = floor(x);
    float3 w = fract(x);
    
    float3 u = w*w*w*(w*(w*6.0-15.0)+10.0);
    
    float n = p.x + 317.0*p.y + 157.0*p.z;
    
    float a = hash1(n+0.0);
    float b = hash1(n+1.0);
    float c = hash1(n+317.0);
    float d = hash1(n+318.0);
    float e = hash1(n+157.0);
	float f = hash1(n+158.0);
    float g = hash1(n+474.0);
    float h = hash1(n+475.0);

    float k0 =   a;
    float k1 =   b - a;
    float k2 =   c - a;
    float k3 =   e - a;
    float k4 =   a - b - c + d;
    float k5 =   a - c - e + g;
    float k6 =   a - b - e + f;
    float k7 = - a + b + c - d + e - f - g + h;

    return -1.0+2.0*(k0 + k1*u.x + k2*u.y + k3*u.z + k4*u.x*u.y + k5*u.y*u.z + k6*u.z*u.x + k7*u.x*u.y*u.z);
}

static const float3x3 m3 = float3x3
                      ( 0.00,  0.80,  0.60,
                      -0.80,  0.36, -0.48,
                      -0.60, -0.48,  0.64 );

// Taken from Inigo Quilez's Rainforest ShaderToy:
// https://www.shadertoy.com/view/4ttSWf
float fbm( in float3 x, int iterations, int iFrame)
{
    float f = 2.0;
    float s = 0.5;
    float a = 0.0;
    float b = 0.5;
    for( int i=min(0, iFrame); i<iterations; i++ )
    {
        float n = noise(x);
        a += b*n;
        b *= s;
        //x = f*m3*x;
        x = f * mul(m3, x);
    }
	return a;
}

// Taken from Inigo Quilez's Rainforest ShaderToy:
// https://www.shadertoy.com/view/4ttSWf
float fbm_4( in float3 x, int iFrame)
{
    return fbm(x, 4, iFrame);
}

float GetWaterWavesDisplacement(float3 position, float time)
{
    // const float scale = 1.0 / 5.0;

    // float height = 7.0 * sin(position.x / 15.0 + time * 1.3) +
    //                6.0 * cos(position.z / 15.0 + time / 1.1);

    // return height * (5.0 / 13.0); // scales from max ±13 to ±5

    return 7.0 * sin(position.x / 15.0 + time * 1.3) +
         6.0 * cos(position.z / 15.0 + time / 1.1);
}

float GetWaterNoise(float3 position, float waterTurbulence, float time, int iFrame)
{
    return waterTurbulence * fbm_4(position / 15.0 + time / 3.0, iFrame);
}

float QueryOceanDistanceField(float3 pos, float waterTurbulence, float time, int iFrame)
{
    return GetWaterWavesDisplacement(pos, time)
         + GetWaterNoise(pos, waterTurbulence, time, iFrame);
}

// Taken from https://iquilezles.org/articles/distfunctions
float sdPlane(float3 p)
{
    return p.y;
}

// Taken from https://iquilezles.org/articles/distfunctions
float sdSmoothSubtraction(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
    return mix(d2, -d1, h) + k * h * (1.0 - h);
}

float QueryVolumetricDistanceField(in float3 pos, float waterTurbulence, float time, int iFrame)
{
    float waveY = QueryOceanDistanceField(pos, waterTurbulence, time, iFrame); // gives you the surface height at this XZ point
    return pos.y - waveY;                               // signed distance from current point to the surface

    // float minDist = QueryOceanDistanceField(pos, waterTurbulence, time, iFrame);

    // minDist = sdSmoothSubtraction(sdSphere(pos, float3(0.0, 0.0, -3), 35.0) + 5.0 * 
    //     fbm_4(pos / float3(12, 20, 12) - time / 5.0, iFrame), minDist, 12.0);

    // float groundLevel = 5;

    // minDist = opSmoothUnion(minDist, sdPlane(pos - float3(0, groundLevel - 1.0, 0)), 13.0);

    // return minDist;
}

// --------------------------------------------//
//               SDF Functions
// --------------------------------------------//

public float sdSphere(float3 p, float3 sphereCenter, float sphereRadius)
{
    return length(p - sphereCenter) - sphereRadius;
}

public float sdSphere(float3 p, float sphereRadius)
{
    return length(p) - sphereRadius;
}

public float sdBox(float3 p, float3 b)
{
    float3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

// Define mix function in Slang (similar to GLSL mix or HLSL lerp)
public float mix(float x, float y, float a) {
    return x * (1.0 - a) + y * a;
}

public float opSmoothUnion(float d1, float d2, float k)
{
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

public float map(float3 p, float waterTurbulence, float time, int iFrame)
{
    float water = QueryVolumetricDistanceField(p, waterTurbulence, time, iFrame);
    //float water = GetWaterNoise(p, waterTurbulence, time, iFrame);

    // float3 c = float3(0, 10, 0);
    // float bound = sdBox(p - c, float3(10, 10, 10));

    // return max(water, bound); // Only visible where both inside

    float3 c = float3(water, 10, water);
    float sphere = sdSphere(p - c, 1);

    return water;

    // c = float3(0, 10, 0);
    // float sphere2 = sdSphere(p - c, .8);
    // // float box = sdBox(p - c, float3(0.75, 0.75, 0.75));

    // return opSmoothUnion(sphere, sphere2, 2);
}

// Taken from the Microsoft DirectX-Graphics-Samples git repo:
// https://github1s.com/microsoft/DirectX-Graphics-Samples/blob/master/Samples/Desktop/D3D12Raytracing/src/D3D12RaytracingProceduralGeometry/SignedDistancePrimitives.hlsli#L275-L283
public float3 sdCalculateNormal(in float3 pos, float waterTurbulence, float time, int iFrame)
{
    float2 e = float2(1.0, -1.0) * 0.5773 * 0.0001;
    return normalize(
        e.xyy * map(pos + e.xyy, waterTurbulence, time, iFrame) +
        e.yyx * map(pos + e.yyx, waterTurbulence, time, iFrame) +
        e.yxy * map(pos + e.yxy, waterTurbulence, time, iFrame) +
        e.xxx * map(pos + e.xxx, waterTurbulence, time, iFrame) );
}

public float mapDebug(float3 p, float time)
{
    float3 c = float3(sin(time * 3) * 3, 15, 0);
    float sphere = sdSphere(p - c, 1);

    c = float3(0, 15, 0);
    float box = sdBox(p - c, float3(0.75, 0.75, 0.75));

    return opSmoothUnion(sphere, box, 2); // Distance to a sphere of radius 1.
}