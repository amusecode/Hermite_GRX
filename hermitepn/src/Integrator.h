#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "Common.h"
#include "Particle.h"
#include "Barrier.h"

#include <thread>

class Perturbation;
class Model;

class Integrator
{
public:
	Integrator(size_t num_threads);
	virtual ~Integrator();
	
	void Evolve(Model *model, Real t_end);
	
	Real GetTimeStepParameter();
	void SetTimeStepParameter(Real dtparam);
	
	size_t GetNumThreads();
	void SetNumThreads(size_t num_threads);
	
	Perturbation *GetPerturbation();
	void SetPerturbation(Perturbation *pert);
	
protected:
	virtual void EvolveWorker(size_t thread_id) = 0;
	virtual void NumThreadsChanged();
	
	Real m_DtParam;
	Perturbation *m_Perturbation;
	
	size_t m_NumThreads;
	std::vector<std::thread> m_Threads;
	
	Barrier m_BeginBarrier;
	Barrier m_EndBarrier;
	
	bool m_Stop;
	Model *m_Model;
	Real m_TimeEnd;
};

#endif // INTEGRATOR_H
