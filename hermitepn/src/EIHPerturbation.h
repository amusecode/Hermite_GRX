#ifndef EIH_PERTURBATION
#define EIH_PERTURBATION

#include "Perturbation.h"

class EIHPerturbation
 : public Perturbation
{
	void CalculateAccelerations(Model *model, int thread_id, int num_threads);
	
	Real GetEnergy(Model *model);
	Vec GetLinearMomentum(Model *model);
};

#endif // EIH_PERTURBATION
