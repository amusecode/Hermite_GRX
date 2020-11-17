#ifndef PAIRWISE_1PN_PERTURBATION
#define PAIRWISE_1PN_PERTURBATION

#include "Perturbation.h"

class Pairwise1PNPerturbation
 : public Perturbation
{
	void CalculateAccelerations(Model *model, int thread_id, int num_threads);
	
	Real GetEnergy(Model *model);
	Vec GetLinearMomentum(Model *model);
};

#endif // PAIRWISE_1PN_PERTURBATION
