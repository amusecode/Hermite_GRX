#include "HermiteIntegrator.h"
#include "Perturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

HermiteIntegrator::HermiteIntegrator(bool symmetrized, size_t num_threads)
	: Integrator(num_threads), m_Symmetrized(symmetrized), m_InnerBarrier(num_threads)
{
}

HermiteIntegrator::~HermiteIntegrator()
{
}

void HermiteIntegrator::EvolveWorker(size_t thread_id)
{
	if (thread_id == 0)
	{
		m_StoppingCond.Reset();
		m_ConditionsSet = false;
	}

	// Reevaluate the accelerations and jerks because they may have been
	// changed inbetween calls to Evolve.
	Evaluate(1e300, thread_id);

	m_InnerBarrier.Wait();
	CalculateTimeStep(thread_id);

	m_InnerBarrier.Wait();
	Real time_step_curr = m_TimeSteps.min();

	Real time = m_Model->GetTime();
	while (time < m_TimeEnd && !m_ConditionsSet)
	{
		for (Particle &p : Strided(m_Model->GetAllParticles(), thread_id, m_NumThreads))
			p.Save();

		Real dt = time_step_curr;
		std::cout << "Timestep " << dt;

		Predict(dt, thread_id);

		m_InnerBarrier.Wait();
		Evaluate(dt, thread_id);

		m_InnerBarrier.Wait();
		CalculateTimeStep(thread_id);

		m_InnerBarrier.Wait();
		Real time_step_next = m_TimeSteps.min();

		if (m_Symmetrized)
		{
			dt = 0.5 * (time_step_curr + time_step_next);

			m_InnerBarrier.Wait();
			Predict(dt, thread_id);

			m_InnerBarrier.Wait();
			Evaluate(dt, thread_id);

			m_InnerBarrier.Wait();
			CalculateTimeStep(thread_id);

			m_InnerBarrier.Wait();
			time_step_next = m_TimeSteps.min();
		}

		Correct(dt, thread_id);

		time += dt;
		time_step_curr = time_step_next;

		m_InnerBarrier.Wait();
		if (thread_id == 0)
			m_ConditionsSet = m_StoppingCond.Check(m_Model);
		m_InnerBarrier.Wait();
	}

	// Set the model time only once
	if (thread_id == 0)
		m_Model->SetTime(time);
}

void HermiteIntegrator::NumThreadsChanged()
{
	m_InnerBarrier.SetNumThreads(m_NumThreads);
	m_TimeSteps.resize(m_NumThreads);
}

void HermiteIntegrator::Predict(Real dt, size_t thread_id)
{
	ParticleSetView &all = m_Model->GetAllParticles();

	Real dt2_over_2 = dt * dt / 2;
	Real dt3_over_6 = dt * dt * dt / 6;

	for (Particle &p : Strided(all, thread_id, m_NumThreads))
	{
		p.pos = p.old_pos + p.old_vel * dt + (p.old_acc_newton + p.old_acc_pert) * dt2_over_2 + (p.old_jerk_newton + p.old_jerk_pert) * dt3_over_6;
		p.vel = p.old_vel + (p.old_acc_newton + p.old_acc_pert) * dt + (p.old_jerk_newton + p.old_jerk_pert) * dt2_over_2;
	}
}

void HermiteIntegrator::Evaluate(Real dt, size_t thread_id)
{
	ParticleSetView &all = m_Model->GetAllParticles();

	for (Particle &i : Strided(all, thread_id, m_NumThreads))
	{
		// Set initial accs and jerks at zero.
		i.acc_newton = i.jerk_newton = Vec::Zero();
		i.acc_pert = i.jerk_pert = Vec::Zero();

		// Calculate the Newtonian accs and jerks
		for (Particle &j : all)
		{
			if (i == j)
				continue;

			Vec r_ij = i.pos - j.pos;
			Vec v_ij = i.vel - j.vel;

			Real r2 = r_ij.SquaredNorm();
			Real r = sqrt(r2);
			Real r3 = r2 * r;

			Real m_over_r3 = j.mass / r3;

			i.acc_newton -= m_over_r3 * r_ij;
			i.jerk_newton -= m_over_r3 * (v_ij - (3 * r_ij.Dot(v_ij) / r2) * r_ij);
		}
	}

	// Calculate the perturbing accs and jerks
	m_Perturbation->CalculateAccelerationsAndJerks(m_Model, dt, thread_id, m_NumThreads);
}

void HermiteIntegrator::Correct(Real dt, size_t thread_id)
{
	ParticleSetView &all = m_Model->GetAllParticles();

	Real dt_over_2 = dt / 2;
	Real dt2_over_12 = dt * dt / 12;

	for (Particle &p : Strided(all, thread_id, m_NumThreads))
	{
		p.vel = p.old_vel + (p.old_acc_newton + p.old_acc_pert + p.acc_newton + p.acc_pert) * dt_over_2 + (p.old_jerk_newton + p.old_jerk_pert - p.jerk_newton - p.jerk_pert) * dt2_over_12;

		p.pos = p.old_pos + (p.old_vel + p.vel) * dt_over_2 + (p.old_acc_newton + p.old_acc_pert - p.acc_newton - p.acc_pert) * dt2_over_12;
	}
}

void HermiteIntegrator::CalculateTimeStep(size_t thread_id)
{
	ParticleSetView &all = m_Model->GetAllParticles();
	Real coll_time_q = 1e300;

	for (Particle &i : Strided(all, thread_id, m_NumThreads))
	{
		for (Particle &j : all)
		{
			// Ignore identical particles.
			if (i == j)
				continue;

			Real r2 = (j.pos - i.pos).SquaredNorm();
			Real v2 = (j.vel - i.vel).SquaredNorm();

			// First estimate: unaccelerated linear motion.
			Real coll_time_q_est = (r2 * r2) / (v2 * v2);
			coll_time_q = min(coll_time_q, coll_time_q_est);

			Real a2 = (j.acc_newton + j.acc_pert - i.acc_newton - i.acc_pert).SquaredNorm();
			Real m = i.mass + j.mass;

			// Second estimate: free fall time.
			coll_time_q_est = r2 / (a2 * (m * m));
			coll_time_q = min(coll_time_q, coll_time_q_est);
		}
	}

	m_TimeSteps[thread_id] = m_DtParam * sqrt(sqrt(coll_time_q));
}
