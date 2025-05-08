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

    constexpr int Width = MetersToPixels(19.20f);
    constexpr int Height = MetersToPixels(10.80f);

    inline int NbParticles = 5'000;
    
    static constexpr float WALLSIZE = Metrics::MetersToPixels(0.1f);
    static constexpr float WALLDIST = Metrics::MetersToPixels(2.0f);
    static constexpr float PARTICLESIZE = Metrics::MetersToPixels(0.05f);

    //test
}
