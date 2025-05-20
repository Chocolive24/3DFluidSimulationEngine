#pragma once
#include "Utility.h"
namespace SPH {

	inline float SmoothingRadius = 10;
	inline float TargetDensity = 50;
	inline float PressureMultiplier = 25;
	inline float nearPressureMultiplier = 0.5;
	inline float ViscosityStrength = 5000;

      /*  inline float SmoothingRadius = 1.5f;
        inline float TargetDensity = 120.f;
        inline float PressureMultiplier = 5.f;
        inline float nearPressureMultiplier = 2.16f;
        inline float ViscosityStrength = 500;*/

	//inline float SpikyPow2ScalingFactor = 6 / (PI * Pow(SmoothingRadius, 4));;
	//inline float SpikyPow3ScalingFactor = 10 / (PI * Pow(SmoothingRadius, 5));;
	//inline float SpikyPow2DerivativeScalingFactor = 12 / (PI * Pow(SmoothingRadius, 4));
	//inline float SpikyPow3DerivativeScalingFactor = 30 / (PI * Pow(SmoothingRadius, 5));
	//inline float Poly6ScalingFactor = 4 / (PI * Pow(SmoothingRadius, 8));
	inline float SpikyPow2ScalingFactor = 15.0f / (PI * Pow(SmoothingRadius, 5));
	inline float SpikyPow3ScalingFactor = 315.0f / (64.0f * PI * Pow(SmoothingRadius, 9));
	inline float SpikyPow2DerivativeScalingFactor = -30.0f / (PI * Pow(SmoothingRadius, 5));
	inline float SpikyPow3DerivativeScalingFactor = -945.0f / (32.0f * PI * Pow(SmoothingRadius, 9));
	inline float Poly6ScalingFactor = 315.0f / (64.0f * PI * Pow(SmoothingRadius, 9));


	float SmoothingKernelPoly6(float radius, float dst);

	float SpikyKernelPow3(float radius, float dst);
	float SpikyKernelPow2(float radius, float dst);

	float DerivativeSpikyPow3(float radius, float dst);

	float DerivativeSpikyPow2(float radius, float dst);

	float DensityKernel(float radius, float dst);

	float NearDensityKernel(float radius, float dst);

	float DensityDerivative(float radius, float dst);

	float NearDensityDerivative(float radius, float dst);

	float ViscosityKernel(float radius, float dst);
}
