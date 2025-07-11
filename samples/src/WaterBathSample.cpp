#include "WaterBathSample.h"

std::string WaterBathSample::GetName() noexcept {
	return "Water Bath";
}

std::string WaterBathSample::GetDescription() noexcept {
	return "CONTROLS: WASD to move in space, CTRL to go down, SPACE to go up. "
		"\n\nSPH simulation ";
}

void WaterBathSample::DrawImgui() noexcept
{
    if (ImGui::SliderFloat("Gravity", &_world.Gravity, 0.0f, 500.0f))
    {
        _world.Gravity = _world.Gravity;
    }
    static int numParticles = NbParticles;
    if (ImGui::SliderInt("Number of Particles", &numParticles, 0, 10000))
    {
        //NbParticles = numParticles;
    }
    if (ImGui::SliderFloat("Smoothing radius", &SPH::SmoothingRadius, 0.1f, 50.0f))
    {
        SPH::SmoothingRadius = SPH::SmoothingRadius;
    }
    if (ImGui::SliderFloat("Target density", &SPH::TargetDensity, 0.0f, 1000.0f))
    {
        SPH::TargetDensity = SPH::TargetDensity;
    }
    if (ImGui::SliderFloat("Pressure multiplier", &SPH::PressureMultiplier, 0.0f, 1000.0f))
    {
        SPH::PressureMultiplier = SPH::PressureMultiplier;
    }
    if (ImGui::SliderFloat("Viscosity strength", &SPH::ViscosityStrength, 0.0f, 10000.0f))
    {
        SPH::ViscosityStrength = SPH::ViscosityStrength;
    }
    if (ImGui::SliderFloat("Fixed DeltaTime Diviser", &Metrics::fixedDeltaTimeDiviser, 0.0f, 200.f))
    {
        Metrics::fixedDeltaTimeDiviser = Metrics::fixedDeltaTimeDiviser;
        Metrics::kFixedDeltaTime = 1.f / Metrics::fixedDeltaTimeDiviser;
        // std::cout << kFixedDeltaTime << '\n' << std::endl;
    }
    if (ImGui::SliderFloat("Collision Damping", &SPH::collisionDamping, 0.0f, 1.f))
    {
        SPH::collisionDamping = SPH::collisionDamping;
    }
}

void WaterBathSample::OnCollisionEnter(ColliderRef col1,
									   ColliderRef col2) noexcept {
}

void WaterBathSample::OnCollisionExit(ColliderRef col1,
									  ColliderRef col2) noexcept {
}

void WaterBathSample::SampleSetUp() noexcept {
	_world.SetContactListener(this);
	GraphicsData gd;
        // Ground
        //CreateWall({0, -WALLDIST - WALLSIZE, 0}, {-WALLDIST, -WALLSIZE, -WALLDIST}, {WALLDIST *2, WALLSIZE, WALLDIST}, true);

        //// Wall 1
        //CreateWall({-WALLDIST - WALLSIZE, 0, 0}, {-WALLSIZE, -WALLDIST, -WALLDIST}, {WALLSIZE, WALLDIST, WALLDIST}, false);

        //// Wall 2
        //CreateWall({WALLDIST * 2 + WALLSIZE, 0, 0}, {-WALLSIZE, -WALLDIST, -WALLDIST}, {WALLSIZE, WALLDIST, WALLDIST}, false);

        //// Wall 3
        //CreateWall({0, 0, -WALLDIST - WALLSIZE}, {-WALLDIST, -WALLDIST, -WALLSIZE}, {WALLDIST* 2, WALLDIST, WALLSIZE}, false);

        //// Wall 4
        //CreateWall({0, 0, WALLDIST + WALLSIZE}, {-WALLDIST, -WALLDIST, -WALLSIZE}, {WALLDIST* 2, WALLDIST, WALLSIZE}, false);

        //// Roof
        //CreateWall({0, WALLDIST + WALLSIZE, 0}, {-WALLDIST, -WALLSIZE, -WALLDIST}, {WALLDIST* 2, WALLSIZE, WALLDIST}, false);

        //float cubeSize = 30; // Taille du cube sur chaque axe
        //int particlesPerAxis = static_cast<int>(std::ceil(std::cbrt(NbParticles)));
        //float spacing = cubeSize / static_cast<float>(particlesPerAxis);

        //int currentNbrParticles = 0;

        //for (int x = 0; x < particlesPerAxis; ++x)
        //    for (int y = 0; y < particlesPerAxis; ++y)
        //        for (int z = 0; z < particlesPerAxis; ++z)
        //        {
        //            if (currentNbrParticles >= NbParticles)
        //                break;

        //            float px = x * spacing;
        //            float py = y * spacing;
        //            float pz = z * spacing;

        //            XMVECTOR pos = XMVectorSet(px, py, pz, 0.0f);
        //            CreateBall(pos, PARTICLESIZE, BodyType::FLUID);

        //            currentNbrParticles++;
        //        }


        float spacing = SPH::SmoothingRadius * 0.8f; // ensures neighbor support
        int particlesPerAxis = static_cast<int>((2 * WALLDIST) / spacing);

        int currentNbrParticles = 0;

        for (int x = 0; x < particlesPerAxis; ++x)
            for (int y = 0; y < particlesPerAxis; ++y)
                for (int z = 0; z < particlesPerAxis; ++z)
                {
                    if (currentNbrParticles >= NbParticles)
                        break;

                    XMVECTOR pos = XMVectorSet(-WALLDIST + x * spacing, -WALLDIST + y * spacing, -WALLDIST + z * spacing, 0.0f);

                    CreateBall(pos, PARTICLESIZE, BodyType::FLUID);
                    currentNbrParticles++;
                }

        //for (size_t i = 0; i < NbParticles;)
        //{
        //    /*XMVECTOR pos = XMVectorSet(
        //        Random::Range(-WALLDIST / 4, WALLDIST),
        //        Random::Range(-WALLDIST / 4, WALLDIST),
        //        Random::Range(-WALLDIST / 4, WALLDIST),
        //        0.0f
        //    );*/

        //    XMVECTOR pos = XMVectorSet(
        //        Random::Range(-WALLDIST, WALLDIST),
        //        Random::Range(-WALLDIST, WALLDIST),
        //        Random::Range(-WALLDIST, WALLDIST),
        //        0.0f
        //    );

        //    bool overlaps = false;
        //    for (const auto& existing : particlePositions)
        //    {
        //        if (XMVectorGetX(XMVector3LengthSq(pos - existing)) < (PARTICLESIZE * PARTICLESIZE * 4))
        //        {
        //            overlaps = true;
        //            break;
        //        }
        //    }

        //    if (!overlaps)
        //    {
        //        particlePositions.push_back(pos);
        //        CreateBall(pos, PARTICLESIZE, BodyType::FLUID);
        //        ++i;
        //    }
        //}

       /* for (size_t i = 0; i < NbParticles; i++)
        {
            CreateBall({
                0.f, 0.f, 0.f, 0.f},
                PARTICLESIZE, BodyType::FLUID);
        }*/

        //int particlesPerAxis = static_cast<int>(cbrt(Metrics::NbParticles)) + 1;
        //float spacing = 5.f; //PARTICLESIZE * 2;

        //for (int i = 0; i < NbParticles; i++)
        //{
        //    int xIndex = i % particlesPerAxis;
        //    int yIndex = (i / particlesPerAxis) % particlesPerAxis;
        //    int zIndex = i / (particlesPerAxis * particlesPerAxis);

        //    float x = (xIndex - particlesPerAxis / 2.f + 0.5f) * spacing;
        //    float y = (yIndex - particlesPerAxis / 2.f + 0.5f) * spacing;
        //    float z = (zIndex - particlesPerAxis / 2.f + 0.5f) * spacing;

        //    x += Random::Range(0.f, 1.f);
        //    y += Random::Range(0.f, 1.f);
        //    z += Random::Range(0.f, 1.f);

        //    XMVECTOR pos{x, y, z};
        //    CreateBall(pos, PARTICLESIZE, BodyType::FLUID);
        //}


	/*for (size_t i = 0; i < NbParticles; i++) {
		CreateBall({ Random::Range(-WALLDIST * 0.8f, WALLDIST * 0.8f),
					 Random::Range(-WALLDIST * 0.8f, WALLDIST * 0.8f),
					 Random::Range(-WALLDIST * 0.8f, WALLDIST * 0.8f) }, PARTICLESIZE, BodyType::FLUID);
	}*/

}
void WaterBathSample::SampleUpdate() noexcept {
	//if (NbParticles + 5 < AllGraphicsData.size()) {
	//	AllGraphicsData.erase(AllGraphicsData.begin() + NbParticles + 5,
	//						  AllGraphicsData.end());
	//}

	//if (_mouseLeftReleased) {
	//	CreateBall({ 0,1000,0 }, PARTICLESIZE * 5, BodyType::DYNAMIC);
	//}
	/*else if (_mouseRightReleased) {
		CreateRect(_mousePos);
	}*/

	for (std::size_t i = 0; i < _colRefs.size(); ++i) {
		const auto& col = _world.GetCollider(_colRefs[i]);

		const auto& shape = _world.GetCollider(_colRefs[i]).Shape;

		auto& body = _world.GetBody(col.BodyRef);

		switch (shape.index()) {
		case static_cast<int>(ShapeType::Sphere):

		if (XMVectorGetY(col.BodyPosition) <= -WALLDIST * 2)//fix to reduce quadtree size
		{
			_world.GetBody(col.BodyRef).Position = XMVectorZero();
			_world.GetBody(col.BodyRef).Velocity = XMVectorZero();
		}

		if (XMVectorGetX(col.BodyPosition) <= -WALLDIST)
		{
			body.Velocity = XMVectorSetX(body.Velocity, Abs(XMVectorGetX(body.Velocity)));
		}
		else if (XMVectorGetX(col.BodyPosition) >= WALLDIST)
		{
			body.Velocity = XMVectorSetX(body.Velocity, -Abs(XMVectorGetX(body.Velocity)));
		}
		if (XMVectorGetY(col.BodyPosition) <= -WALLDIST)
		{
			body.Position = XMVectorSetY(body.Position, -WALLDIST);
			body.Velocity = XMVectorSetY(body.Velocity, Abs(XMVectorGetY(body.Velocity)));
		}
		else if (XMVectorGetY(col.BodyPosition) >= WALLDIST)
		{
			body.Velocity = XMVectorSetY(body.Velocity, -Abs(XMVectorGetY(body.Velocity)));
		}
		if (XMVectorGetZ(col.BodyPosition) <= -WALLDIST)
		{
			body.Velocity = XMVectorSetZ(body.Velocity, Abs(XMVectorGetZ(body.Velocity)));
		}
		else if (XMVectorGetZ(col.BodyPosition) >= WALLDIST)
		{
			body.Velocity = XMVectorSetZ(body.Velocity, -Abs(XMVectorGetZ(body.Velocity)));
		}

		AllGraphicsData[i].Shape = std::get<SphereF>(shape) + col.BodyPosition;

		break;
		case static_cast<int>(ShapeType::Cuboid):
		AllGraphicsData[i].Shape = std::get<CuboidF>(shape) + col.BodyPosition;
		break;
		default:
		break;
		}
	}

	//_quadTreeGraphicsData.clear();
	//DrawQuadtree(_world.OctTree.Nodes[0]);
	//AllGraphicsData.insert(AllGraphicsData.end(), _quadTreeGraphicsData.begin(),
	//					   _quadTreeGraphicsData.end());
}

void WaterBathSample::SampleTearDown() noexcept {}

void WaterBathSample::CreateBall(XMVECTOR position, float radius, BodyType type) noexcept {
	const auto sphereBodyRef = _world.CreateBody(type);
	_bodyRefs.push_back(sphereBodyRef);
	auto& sphereBody = _world.GetBody(sphereBodyRef);

	sphereBody.Mass = 1.f;//(type == BodyType::FLUID) ? 1 : 10;

	sphereBody.Position = position;

        
        float x = Random::Range(-1, 1);
        float y = Random::Range(-1, 1);
        float z = Random::Range(-1, 1);

        XMVECTOR vel = XMVectorSet(x, y, z, 0.0f);

        sphereBody.Velocity = vel;

	const auto sphereColRef = _world.CreateCollider(sphereBodyRef);
	_colRefs.push_back(sphereColRef);
	auto& sphereCol = _world.GetCollider(sphereColRef);
	sphereCol.Shape = Sphere(XMVectorZero(), radius);
	sphereCol.BodyPosition = sphereBody.Position;
	sphereCol.Restitution = 0.f;
	sphereCol.IsTrigger = false;

	GraphicsData gd;

	gd.Color = (type == BodyType::FLUID) ? Color{ 170, 213, 219 } : Color{ 100, 213, 100 };

	AllGraphicsData.emplace_back(gd);
}

void WaterBathSample::CreateWall(XMVECTOR position, XMVECTOR minBound, XMVECTOR maxBound, bool isFilled) noexcept
{
	const auto wallRef = _world.CreateBody();
	_bodyRefs.push_back(wallRef);
	auto& wallBody = _world.GetBody(wallRef);
	wallBody.Type = BodyType::STATIC;
	wallBody.Mass = 1;

	wallBody.Position = position;

	const auto wallColRef = _world.CreateCollider(wallRef);
	_colRefs.push_back(wallColRef);
	auto& wallCol = _world.GetCollider(wallColRef);
	wallCol.Shape = CuboidF(minBound, maxBound);
	wallCol.BodyPosition = position;
	wallCol.Restitution = 0.f;

	GraphicsData gd;
	gd.Color = { 160,160,160 };
	gd.Filled = isFilled;
        gd.Shape = wallCol.Shape;

	AllGraphicsData.emplace_back(gd);
}

void WaterBathSample::DrawQuadtree(const BVHNode& node) noexcept {
	if (node.Children[0] == nullptr) {
		_quadTreeGraphicsData.push_back({ CuboidF(node.Bounds), false });
	}
	else {
		for (int i = 0; i < 8; i++) {
			DrawQuadtree(*node.Children[i]);
		}
	}
}
