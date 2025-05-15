#pragma once

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

    static constexpr float kFixedDeltaTime = 1.f / 60.f;

    constexpr int Width = MetersToPixels(19.20f);
    constexpr int Height = MetersToPixels(10.80f);

    inline int NbParticles = 10'000;
    
    static constexpr float WALLSIZE = Metrics::MetersToPixels(0.1f);
    static constexpr float WALLDIST = Metrics::MetersToPixels(2.f);
    static constexpr float PARTICLESIZE = Metrics::MetersToPixels(0.05f);
    static constexpr float PARTICLESPACING = Metrics::MetersToPixels(0.05f);

    inline float densityGraphicsMultiplier = 100.f;
    inline int density_map_size = 128;
    inline float sim_bounds = (WALLDIST) * 2.f;
    inline float voxelSize = sim_bounds / float(density_map_size);

    //test
}
