#ifndef PERTURBATION_H
#define PERTURBATION_H

#include "Common.h"
#include "Particle.h"

class Model;

class Perturbation
{
public:
	virtual ~Perturbation() {}
	virtual void CalculateAccelerations(Model *model, int thread_id, int num_threads);
	virtual void CalculateJerks(Model *model, Real dt, int thread_id, int num_threads);
	
	virtual void CalculateAccelerationsAndJerks(Model *model, Real dt, int thread_id, int num_threads);
	
	virtual Real GetEnergy(Model *model);
	virtual Vec GetLinearMomentum(Model *model);
	
	Real GetSpeedOfLight();
	void SetSpeedOfLight(Real speed_of_light);
	
protected:
	Real m_SpeedOfLight;
};

#endif // PERTURBATION_H
