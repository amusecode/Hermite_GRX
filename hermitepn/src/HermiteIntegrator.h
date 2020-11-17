#include "Integrator.h"

#include "StoppingConditions.h"

#include <thread>
#include <vector>
#include <valarray>

class HermiteIntegrator : public Integrator
{
public:
	HermiteIntegrator(bool symmetrized, size_t num_threads = 1);
	virtual ~HermiteIntegrator();
	
private:
	void EvolveWorker(size_t thread_id);
	void NumThreadsChanged();
	
	void Predict(Real dt, size_t thread_id);
	void Evaluate(Real dt, size_t thread_id);
	void Correct(Real dt, size_t thread_id);
	
	void CalculateTimeStep(size_t thread_id);
	
	bool m_Symmetrized;
	
	std::valarray<Real> m_TimeSteps;
	Barrier m_InnerBarrier;
	
	StoppingConditions m_StoppingCond;
	bool m_ConditionsSet;
};
