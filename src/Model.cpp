#include "Model.h"
#include "Integrator.h"
#include "Perturbation.h"

#include <algorithm>
#include <functional>

using namespace std;

Model::Model()
{
	m_Time = 0;
}

void Model::AddLargeParticle(const Particle &p)
{
	m_LargeParticles.push_back(p);
	UpdateParticleView();
}

void Model::AddSmallParticle(const Particle &p)
{
	m_SmallParticles.push_back(p);
	UpdateParticleView();
}

Particle *Model::GetParticle(int id)
{
	for (reference_wrapper<Particle> &i : m_AllParticles)
	{
		if (i.get().id == id)
			return &(i.get());
	}
	return nullptr;
}

void Model::RemoveParticle(int id)
{
	Particle p = Particle(id, 0, Vec(), Vec(), 0);
	
	m_LargeParticles.erase(std::remove(m_LargeParticles.begin(), m_LargeParticles.end(), p), m_LargeParticles.end());
	m_SmallParticles.erase(std::remove(m_SmallParticles.begin(), m_SmallParticles.end(), p), m_SmallParticles.end());
	
	UpdateParticleView();
}

size_t Model::GetNumLargeParticles()
{
	return m_LargeParticles.size();
}

size_t Model::GetNumSmallParticles()
{
	return m_SmallParticles.size();
}

size_t Model::GetNumParticles()
{
	return m_SmallParticles.size() + m_LargeParticles.size();
}

ParticleSet &Model::GetLargeParticles()
{
	return m_LargeParticles;
}

ParticleSet &Model::GetSmallParticles()
{
	return m_SmallParticles;
}

ParticleSetView &Model::GetAllParticles()
{
	return m_AllParticles;
}

Real Model::GetTime()
{
	return m_Time;
}

void Model::SetTime(Real time)
{
	m_Time = time;
}

Real Model::GetNewtonianPotentialEnergy()
{
	Real pot = 0;

	for (Particle &i : m_AllParticles)
	{
		for (Particle &j : m_AllParticles)
		{
			if (i != j)
			{
				Real r = sqrt((i.pos - j.pos).SquaredNorm());
				pot -= i.mass * j.mass / r;
			}
		}
	}
	return 0.5 * pot; // Counting double
}

Real Model::GetNewtonianKineticEnergy()
{
	Real kin = 0;

	for (Particle &p : m_AllParticles)
	{
		 kin += p.mass * p.vel.SquaredNorm();
	}
	return 0.5 * kin;
}

Real Model::GetNewtonianEnergy()
{
	return GetNewtonianKineticEnergy() + GetNewtonianPotentialEnergy();
}

Vec Model::GetNewtonianLinearMomentum()
{
	Vec p = Vec::Zero();
	
	for (Particle &i : m_AllParticles)
		p += i.mass * i.vel;
	
	return p;
}

void Model::UpdateParticleView()
{
	m_AllParticles = MakeParticleSetView(m_LargeParticles, m_SmallParticles);
}
