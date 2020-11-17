#ifndef MODEL_H
#define MODEL_H

#include "Common.h"
#include "Particle.h"

class Integrator;
class Interaction;

class Model
{
public:
	Model();
	
	void AddLargeParticle(const Particle &p);
	void AddSmallParticle(const Particle &p);
	
	Particle *GetParticle(int id);
	void RemoveParticle(int id);
	
	size_t GetNumLargeParticles();
	size_t GetNumSmallParticles();
	size_t GetNumParticles();
	
	// These functions provide access to particle data, but may not be
	// used for adding or removing particles.
	ParticleSet &GetLargeParticles();
	ParticleSet &GetSmallParticles();
	ParticleSetView &GetAllParticles();
	
	Real GetTime();
	void SetTime(Real time);
	
	Real GetNewtonianPotentialEnergy();
	Real GetNewtonianKineticEnergy();
	
	Real GetNewtonianEnergy();
	Vec GetNewtonianLinearMomentum();
	
private:
	void UpdateParticleView();
	
	Real m_Time;
	
	ParticleSet m_LargeParticles;
	ParticleSet m_SmallParticles;
	ParticleSetView m_AllParticles;
};

#endif // MODEL_H
