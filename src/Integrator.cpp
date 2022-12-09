#include "Integrator.h"

#include <iostream>

using namespace std;

Integrator::Integrator(size_t num_threads)
	: m_DtParam(0.04), m_Perturbation(nullptr),
	  m_BeginBarrier(num_threads + 1), m_EndBarrier(num_threads + 1),
	  m_Stop(false), m_Model(nullptr)
{
	SetNumThreads(num_threads);
}

Integrator::~Integrator()
{
	SetNumThreads(0);
}

void Integrator::Evolve(Model *model, Real t_end)
{
	m_Model = model;
	m_TimeEnd = t_end;

	if (m_NumThreads > 1)
	{
		m_BeginBarrier.Wait();
		m_EndBarrier.Wait();
	}
	else
	{
		EvolveWorker(0);
	}
}

Real Integrator::GetTimeStepParameter()
{
	return m_DtParam;
}

void Integrator::SetTimeStepParameter(Real dtparam)
{
	m_DtParam = dtparam;
}

size_t Integrator::GetNumThreads()
{
	return m_NumThreads;
}

void Integrator::SetNumThreads(size_t num_threads)
{
	if (m_NumThreads > 1)
	{
		// Stop the running threads
		m_Stop = true;
		m_BeginBarrier.Enforce(false);

		for (thread &th : m_Threads)
			th.join();
		m_Threads.clear();

		m_Stop = false;
		m_BeginBarrier.Enforce(true);
	}

	// Set the new number of threads
	m_NumThreads = num_threads;
	m_BeginBarrier.SetNumThreads(m_NumThreads + 1);
	m_EndBarrier.SetNumThreads(m_NumThreads + 1);
	NumThreadsChanged();

	if (m_NumThreads > 1)
	{
		// Start these threads with the worker program
		for (size_t i = 0; i < num_threads; ++i)
		{
			m_Threads.push_back(thread([this, i]()
									   {
				while (true)
				{
					m_BeginBarrier.Wait();
			
					if (m_Stop)
						break;
					else
						EvolveWorker(i);
					
					m_EndBarrier.Wait();
				} }));
		}
	}
}

void Integrator::NumThreadsChanged()
{
}

Perturbation *Integrator::GetPerturbation()
{
	return m_Perturbation;
}

void Integrator::SetPerturbation(Perturbation *pert)
{
	m_Perturbation = pert;
}
