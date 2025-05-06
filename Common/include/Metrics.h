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

    int NbParticles = 1024;

    //test
}
