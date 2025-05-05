#pragma once

#include <iostream>
#include "PhysicsSample.h"
#include "Random.h"

static constexpr float WALLSIZE = Metrics::MetersToPixels(0.1f);
static constexpr float WALLDIST = Metrics::MetersToPixels(1.0f);
static constexpr float PARTICLESIZE = Metrics::MetersToPixels(0.05f);

class WaterBathSample : public PhysicsSample, public ContactListener
{
private:
	std::vector<GraphicsData> _quadTreeGraphicsData;
        std::vector<float> _particle_densities{};
public:

	int NbParticles = 1024;

	std::string GetName() noexcept override;
	std::string GetDescription() noexcept override;
	void DrawImgui() noexcept override;

	void OnTriggerEnter(ColliderRef col1, ColliderRef col2) noexcept override {};

	void OnTriggerExit(ColliderRef col1, ColliderRef col2) noexcept override {};

	void OnCollisionEnter(ColliderRef col1, ColliderRef col2) noexcept override;

	void OnCollisionExit(ColliderRef col1, ColliderRef col2) noexcept override;

protected:
	void SampleSetUp() noexcept override;

	void SampleUpdate() noexcept override;

	void SampleTearDown() noexcept override;

private:
	void CreateBall(XMVECTOR position, float radius, BodyType type) noexcept;

	void CreateWall(XMVECTOR position, XMVECTOR minBound, XMVECTOR maxBound, bool isFilled) noexcept;

	void DrawQuadtree(const BVHNode& node) noexcept;

};
