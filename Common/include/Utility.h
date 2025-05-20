#pragma once

//constexpr double PI = 3.141592653589793238462643383279502884197;

static constexpr float PI = 3.1415926535f;

template<typename T>
[[nodiscard]] constexpr T Clamp(T x, T min, T max) noexcept
{
	if (min > max)
	{
		T tmp = min;
		min = max;
		max = tmp;
	}

	if (x < min)
	{
		x = min;
	}
	else if (x > max)
	{
		x = max;
	}

	return x;
}

template<typename T>
[[nodiscard]] constexpr T Abs(T nbr) noexcept
{
	return nbr < 0 ? -nbr : nbr;
}

template<typename T>
[[nodiscard]] constexpr T Min(const T x, const T y) noexcept
{
	return x < y ? x : y;
}

template<typename T>
[[nodiscard]] constexpr T Max(const T x, const T y) noexcept
{
	return x > y ? x : y;
}
template<typename T>
[[nodiscard]] constexpr T Fact(const T n) noexcept
{
	T result = 1;

	for (int i = 2; i <= Abs(n); i++)
	{
		result *= i;
	}

	if (n < 0)
	{
		result = -result;
	}

	return result;
}

template<typename T>
[[nodiscard]] constexpr T Pow(const T x, const unsigned int power)
{
	if (x == 0)
	{
		return 0;
	}

	T result = 1;

	for (unsigned int i = 1; i <= power; i++)
	{
		result *= x;
	}

	return result;
}
