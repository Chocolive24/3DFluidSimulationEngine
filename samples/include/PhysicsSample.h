#pragma once

#include <string>

#include "Metrics.h"
#include "Timer.h"
#include "World.h"

#include "imgui.h"

struct Color {
  int r = 255, g = 255, b = 255, a = 255;
};
// for the graphics renderer to draw
struct GraphicsData {
  std::variant<SphereF, CuboidF> Shape{
      SphereF(XMVectorZero(), 1)};
  bool Filled = true;
  Color Color;
};

class PhysicsSample {
 public:
  std::vector<GraphicsData> AllGraphicsData;

 protected:
  World _world;

  std::vector<BodyRef> _bodyRefs;
  std::vector<ColliderRef> _colRefs;

  XMVECTOR _mousePos;

  bool _mouseLeftReleased = false;
  bool _mouseRightReleased = false;

 private:
  Timer _timer;
  bool _isActive = true;

 public:
  virtual std::string GetName() noexcept = 0;
  virtual std::string GetDescription() noexcept = 0;
  virtual void DrawImgui() noexcept = 0;

  void SetUp() noexcept;

  void TearDown() noexcept;

  void Update() noexcept;

  void GetMousePos(XMVECTOR mousePos) noexcept;

  void OnLeftClick() noexcept;

  void OnRightClick() noexcept;

  void OnStop() noexcept;

  virtual ~PhysicsSample() noexcept = default;

  [[nodiscard]] World& GetWorldRef() noexcept { return _world; }
  [[nodiscard]] std::vector<BodyRef>& GetBodyRefs() noexcept { return _bodyRefs; }

 protected:
  virtual void SampleSetUp() noexcept = 0;

  virtual void SampleTearDown() noexcept = 0;

  virtual void SampleUpdate() noexcept = 0;
};
