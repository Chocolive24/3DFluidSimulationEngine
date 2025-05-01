#pragma once

#include "DefaultSample.h"
#include "StarSystemSample.h"
#include "TriggerSample.h"
#include "BouncingCollisionSample.h"
#include "GroundCollisionSample.h"
#include "WaterBathSample.h"

class SampleManager
{
private:
	std::vector<std::unique_ptr<PhysicsSample>> _samples;

	std::size_t _sampleIdx = 0;

public:
	[[nodiscard]] int GetSampleNbr() noexcept { return _samples.size(); }
	[[nodiscard]] int GetCurrentIndex() noexcept { return _sampleIdx; }

	[[nodiscard]] std::string GetSampleName(int idx) { return _samples[idx]->GetName(); }
	[[nodiscard]] std::string GetSampleDescription(int idx) { return _samples[idx]->GetDescription(); }
	[[nodiscard]] void DrawImgui(int idx) { _samples[idx]->DrawImgui(); }

	void SetUp() noexcept;

	void UpdateSample() const noexcept;

	void ChangeSample(int idx) noexcept;

	void NextSample() noexcept;

	void PreviousSample() noexcept;

	void RegenerateSample() const noexcept;

	[[nodiscard]] std::vector<GraphicsData>& GetSampleData() const noexcept;

        void GiveMousePositionToSample(XMVECTOR mousePosition) const noexcept;

	void GiveLeftMouseClickToSample() const noexcept;

	void GiveRightMouseClickToSample() const noexcept;
	void StopSample() const noexcept;

        [[nodiscard]] World& GetWorldRef() const noexcept { return _samples[_sampleIdx]->GetWorldRef(); }
};
