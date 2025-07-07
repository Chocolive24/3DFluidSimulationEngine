#pragma once
#include <array>

#include "../../Physics/include/SPH.h"

namespace Metrics // Meter in physical world != meter irl
{
    constexpr float MeterRatio = 100.f;

    [[nodiscard]] constexpr float PixelsToMeters(float pixels) noexcept
    {
        return pixels / MeterRatio;
    }

    [[nodiscard]] constexpr float MetersToPixels(float meters) noexcept
    {
        return meters * MeterRatio;
    }

    static float fixedDeltaTimeDiviser = 120.f;
    static float kFixedDeltaTime = 1.f / fixedDeltaTimeDiviser;
    

    constexpr int Width = MetersToPixels(19.20f);
    constexpr int Height = MetersToPixels(10.80f);

    inline int NbParticles = 10'000;
    
    static constexpr float WALLSIZE = Metrics::MetersToPixels(0.25f);
    static constexpr float WALLDIST = Metrics::MetersToPixels(0.3f);
    static float PARTICLESIZE = SPH::SmoothingRadius * 0.5f;
    //static constexpr float PARTICLESIZE = Metrics::MetersToPixels(0.005f);
    static constexpr float PARTICLESPACING = Metrics::MetersToPixels(0.05f);

  /*  inline int NbParticles = 40'000;

    static constexpr float WALLSIZE = Metrics::MetersToPixels(0.25f);
    static constexpr float WALLDIST = Metrics::MetersToPixels(0.3f);
    static float PARTICLESIZE = SPH::SmoothingRadius * 0.5f;*/
    // static constexpr float PARTICLESIZE = Metrics::MetersToPixels(0.005f);


    inline float densityGraphicsMultiplier = 100.f;
    inline int density_map_size = 150;
    inline float sim_bounds = (WALLDIST) * 2.f;
    inline float voxelSize = sim_bounds / float(density_map_size);

    constexpr std::array<int, 3> voxelGridResolution{8, 8, 8};
    constexpr int voxelGridTotalResolution =
        voxelGridResolution[0] * voxelGridResolution[1] * voxelGridResolution[2];
    //test
}
