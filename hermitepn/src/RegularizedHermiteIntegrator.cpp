#include "RegularizedHermiteIntegrator.h"
#include "Perturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <utility>
#include <set>
#include <list>

using namespace std;

typedef pair<Real, pair<Particle *, Particle *>> RegularizationRank;

bool operator > (const RegularizationRank &l, const RegularizationRank &r)
{
	return l.first > r.first;
}

Real GetRegularizationRank(const Particle &p1, const Particle &p2)
{
	Vec r_21 = p2.pos - p1.pos;
	//Vec v_21 = p2.vel - p1.vel;

	//Real v2 = v_21.SquaredNorm();
	//Real r = r_21.Norm();
	Real mu = p1.mass + p2.mass;

	//Vec e = ((v2 - mu / r) * r_21 - r_21.Dot(v_21) * v_21) / mu;
	//Real ecc = e.Norm();
	
	//Real E = 0.5 * v2 - mu / r;
	//Real a = -mu / (2 * E);

	//Real p = a * (1 - ecc*ecc);

	return mu / r_21.SquaredNorm();
}

RegularizedPair::RegularizedPair(Particle &part1, Particle &part2)
{
	p1 = &part1;
	p2 = &part2;
	
	Quat x = Quat(p2->pos - p1->pos);
	
	// Swap the particles if the denominator becomes close to zero.
	Real xnorm = x.Norm();
	if (((xnorm + x[0]) / xnorm) < 1e-5)
	{
		p1 = &part2;
		p2 = &part1;
		x = -x;
	}
	
	Quat v = Quat(p2->vel - p1->vel);

	U = (x + xnorm) / (sqrt(2*(xnorm + x[0])));
	V = 0.5 * v * U.Conj().StarConj();
	
	Real mu = p1->mass + p2->mass;
	h = (mu - 2*V.SquaredNorm()) / xnorm;
	
	pos_cm = (p1->mass * p1->pos + p2->mass * p2->pos) / (p1->mass + p2->mass);
	vel_cm = (p1->mass * p1->vel + p2->mass * p2->vel) / (p1->mass + p2->mass);

	Evaluate();
	Save();
}

bool RegularizedPair::Contains(const Particle &part1, const Particle &part2) const
{
	if (part1 == *p1)
		return part2 == *p2;
	if (part1 == *p2)
		return part2 == *p1;
	return false;
}

bool RegularizedPair::Contains(const Particle &p) const
{
	return (p == *p1) || (p == *p2);
}

void RegularizedPair::Save()
{
	old_U = U;
	old_V = V;
	old_acc = acc;
	old_jerk = jerk;

	old_pos_cm = pos_cm;
	old_vel_cm = vel_cm;
	old_acc_cm = acc_cm;
	old_jerk_cm = jerk_cm;

	old_h = h;
	old_h_prime = h_prime;
	old_h_prime2 = h_prime2;
}

void RegularizedPair::Evaluate()
{
	Vec x_prime = (V * U.StarConj() + U * V.StarConj()).GetVecPart();

	Real r = U.SquaredNorm();
	Real r_prime = (V * U.Conj() + U * V.Conj())[0];

	Vec P = p2->acc_newton + p2->acc_pert - p1->acc_newton - p1->acc_pert;
	Vec P_prime = (p2->jerk_newton + p2->jerk_pert - p1->jerk_newton - p1->jerk_pert) * r;

	h_prime = -x_prime.Dot(P);

	acc = -h / 2 * U + r * Quat(P) * U.Conj().StarConj() / 2;
	jerk = 0.5 * (-h_prime * U - h * V + r_prime * Quat(P) * U.Conj().StarConj() + r * Quat(P_prime) * U.Conj().StarConj() + r * Quat(P) * V.Conj().StarConj());

	Vec x_prime2 = (acc * U.StarConj() + 2 * V * V.StarConj() + U * acc.StarConj()).GetVecPart();
	h_prime2 = -x_prime2.Dot(P) - x_prime.Dot(P_prime);

	acc_cm = (p1->mass * (p1->acc_newton + p1->acc_pert) + p2->mass * (p2->acc_newton + p2->acc_pert)) / (p1->mass + p2->mass);
	jerk_cm = (p1->mass * (p1->jerk_newton + p1->jerk_pert) + p2->mass * (p2->jerk_newton + p2->jerk_pert)) / (p1->mass + p2->mass);
}

void RegularizedPair::UpdateParticles()
{
	Vec x = (U * U.StarConj()).GetVecPart();
	Vec v = (2 * V * U.StarConj()).GetVecPart() / U.SquaredNorm();
	
	p1->pos = pos_cm - (p2->mass / (p1->mass + p2->mass)) * x;
	p2->pos = pos_cm + (p1->mass / (p1->mass + p2->mass)) * x;
	
	p1->vel = vel_cm - (p2->mass / (p1->mass + p2->mass)) * v;
	p2->vel = vel_cm + (p1->mass / (p1->mass + p2->mass)) * v;
}

// Using Funato et al (1995)
Real RegularizedPair::ToRealTime(Real ds)
{
	Quat old_snap = -6 * (old_acc - acc) / (ds*ds) - 2 * (2 * old_jerk + jerk) / ds;
	Quat old_crackle = 12 * (old_acc - acc) / (ds*ds*ds) + 6 * (old_jerk + jerk) / (ds*ds);
	//Quat snap = old_snap + ds * old_crackle;
	//Quat crackle = old_crackle;

	Quat half_U = old_U + old_V * ds / 2 + old_acc * (ds*ds) / 8 + old_jerk * (ds*ds*ds) / 48 + old_snap * (ds*ds*ds*ds) / 384 + old_crackle * (ds*ds*ds*ds*ds) / 3840;
	Quat half_V = old_V + old_acc * ds / 2 + old_jerk * (ds*ds) / 8 + old_snap * (ds*ds*ds) / 48 + old_crackle * (ds*ds*ds*ds) / 384;
	Quat half_acc = old_acc + old_jerk * ds / 2 + old_snap * (ds*ds) / 8 + old_crackle * (ds*ds*ds) / 48;
	Quat half_jerk = old_jerk + old_snap * ds / 2 + old_crackle * (ds*ds) / 8;
	Quat half_snap = old_snap + old_crackle * ds / 2;

	Real half_t1 = half_U.SquaredNorm();
	Real half_t3 = (half_acc * half_U.Conj() + 2 * half_V.SquaredNorm() + half_U * half_acc.Conj())[0];
	Real half_t5 = (half_snap * half_U.Conj() + 4 * half_jerk * half_V.Conj() + 6 * half_acc.SquaredNorm() + 4 * half_V * half_jerk.Conj() + half_U * half_snap.Conj())[0];

	Real dt = half_t1 * ds + half_t3 * (ds*ds*ds) / 24 + half_t5 * (ds*ds*ds*ds*ds) / 1920;
	return dt;
}

// Using Newton-Rapson iteration to invert the above.
Real RegularizedPair::ToRegularizedTime(Real dt)
{
	const int N_iter = 4;
	Real ds = dt / U.SquaredNorm();

	for (int i = 0; i < N_iter; ++i)
	{
		Real f = ToRealTime(ds) - dt;
		Quat half_U = old_U + old_V * ds / 2 + old_acc * (ds*ds) / 8 + old_jerk * (ds*ds*ds) / 48;
		Real fderiv = half_U.SquaredNorm();

		ds -= f / fderiv;
	}

	return ds;
}

RegularizedHermiteIntegrator::RegularizedHermiteIntegrator(bool symmetrized, size_t num_threads)
	: Integrator(num_threads), m_Symmetrized(symmetrized)
{
}

RegularizedHermiteIntegrator::~RegularizedHermiteIntegrator()
{
}

bool RegularizedHermiteIntegrator::IsRegularized(const Particle &p1, const Particle &p2) const
{
	for (const RegularizedPair &pair : m_RegularizedPairs)
	{
		if (pair.Contains(p1, p2))
			return true;
	}
	return false;
}

bool RegularizedHermiteIntegrator::IsRegularized(const Particle &p) const
{
	for (const RegularizedPair &pair : m_RegularizedPairs)
	{
		if (pair.Contains(p))
			return true;
	}
	return false;
}

bool RegularizedHermiteIntegrator::SelectRegularizedParticles()
{
	ParticleSetView &all = m_Model->GetAllParticles();

	typedef pair<Real, pair<Particle *, Particle *>> RegularizationRank;
	vector<RegularizationRank> ranks;

	Real min_rank = -1e300;
	const size_t MAX_REGULARIZED = 10;

	for (Particle &i : all)
	{
		for (Particle &j : all)
		{
			if (i.id <= j.id)
				continue;

			Real rank = GetRegularizationRank(i, j);

			if (rank > min_rank)
			{
				RegularizationRank rrank = make_pair(rank, make_pair(&i, &j));
				auto pos = lower_bound(ranks.begin(), ranks.end(), rrank, greater<RegularizationRank>());
				ranks.insert(pos, rrank);

				if (ranks.size() > MAX_REGULARIZED)
				{
					ranks.pop_back();
					min_rank = ranks.back().first;
				}
			}
		}
	}

	set<Particle *> regularized_particles;
	bool haschanged = false;

	for (auto &p : ranks)
	{
		// Check if the particles were already regularized
		if (IsRegularized(*(p.second.first), *(p.second.second)))
		{
			//cout << "Pair (" << p.second.first->id << ", " << p.second.second->id << ") already regularized." << endl;
			regularized_particles.emplace(p.second.first);
			regularized_particles.emplace(p.second.second);
			continue;
		}

		haschanged = true;
		
		// Check if the particles were used handled this time step
		if (regularized_particles.count(p.second.first) + regularized_particles.count(p.second.second))
		{
			//cout << "Pair (" << p.second.first->id << ", " << p.second.second->id << ") already handled this timestep." << endl;
			continue;
		}

		// We have a new pair that has to be regularized
		// First remove all currently regularized pairs that contain one of these particles
		m_RegularizedPairs.remove_if([&, this](RegularizedPair &pair) { return pair.Contains(*(p.second.first)) || pair.Contains(*(p.second.second)); });

		// Add a new regularized pair
		m_RegularizedPairs.push_back(RegularizedPair(*(p.second.first), *(p.second.second)));
		regularized_particles.emplace(p.second.first);
		regularized_particles.emplace(p.second.second);

		cout << "Regularized (" << p.second.first->id << ", " << p.second.second->id << ")" << endl;
	}

	return haschanged;
}

void RegularizedHermiteIntegrator::EvolveWorker(size_t thread_id)
{
	if (thread_id != 0)
		return;
	
	m_StoppingCond.Reset();
	m_ConditionsSet = false;
	
	//m_RegularizedPairs.clear();
	//m_RegularizedPairs.push_back(RegularizedPair(m_Model->GetLargeParticles()[0], m_Model->GetSmallParticles()[0]));
	
	// Reevaluate the accelerations and jerks because they may have been
	// changed inbetween calls to Evolve.
	
	Real dt = 1e300;

	Evaluate(dt);
	for (Particle &p : m_Model->GetAllParticles())
	{
		p.jerk_pert = Vec::Zero();
		p.Save();
	}
	
	Real time_step_curr = CalculateTimeStep();
	
	Real time = m_Model->GetTime();
	while (time < m_TimeEnd && !m_ConditionsSet)
	{
		bool hasregularizationchanged = SelectRegularizedParticles();
		if (hasregularizationchanged)
			Evaluate(dt);

		for (Particle &p : m_Model->GetAllParticles())
			p.Save();
		for (RegularizedPair &p : m_RegularizedPairs)
			p.Save();
		
		dt = time_step_curr;
		
		Predict(dt);
		Evaluate(dt);
		
		Real time_step_next = CalculateTimeStep();
		
		if (m_Symmetrized)
		{
			dt = 0.5 * (time_step_curr + time_step_next);
			
			Predict(dt);
			Evaluate(dt);

			time_step_next = CalculateTimeStep();
		}
		
		Correct(dt);
		
		time += dt;
		time_step_curr = time_step_next;
		
		m_ConditionsSet = m_StoppingCond.Check(m_Model);
	}
	
	m_Model->SetTime(time);
}

void RegularizedHermiteIntegrator::NumThreadsChanged()
{
}

void RegularizedHermiteIntegrator::Predict(Real dt)
{
	ParticleSetView &all = m_Model->GetAllParticles();
	
	Real dt2_over_2 = dt*dt / 2;
	Real dt3_over_6 = dt*dt*dt / 6;
	
	for (Particle &p : all)
	{
		if (IsRegularized(p))
			continue;
		
		p.pos = p.old_pos
			+ p.old_vel * dt 
			+ (p.old_acc_newton + p.old_acc_pert) * dt2_over_2 
			+ (p.old_jerk_newton + p.old_jerk_pert) * dt3_over_6;
		p.vel = p.old_vel
			+ (p.old_acc_newton + p.old_acc_pert) * dt 
			+ (p.old_jerk_newton + p.old_jerk_pert) * dt2_over_2;
	}
	
	for (RegularizedPair &p : m_RegularizedPairs)
	{
		Real ds = p.ToRegularizedTime(dt);
		
		p.U = p.old_U + ds * p.V + ds*ds/2 * p.old_acc + ds*ds*ds/6 * p.old_jerk;
		p.V = p.old_V + ds * p.old_acc + ds*ds/2 * p.old_jerk;
		
		p.pos_cm = p.old_pos_cm + dt * p.old_vel_cm + dt2_over_2 * p.old_acc_cm + dt3_over_6 * p.old_jerk_cm;
		p.vel_cm = p.old_vel_cm + dt * p.old_acc_cm + dt2_over_2 * p.old_jerk_cm;

		p.h = p.old_h + p.old_h_prime * ds + p.old_h_prime2 * ds*ds / 2;
		
		p.UpdateParticles();
	}
}

void RegularizedHermiteIntegrator::Evaluate(Real dt)
{
	ParticleSetView &all = m_Model->GetAllParticles();
	
	for (Particle &i : all)
	{
		// Set initial accs and jerks at zero.
		i.acc_newton = i.jerk_newton = Vec::Zero();
		i.acc_pert = i.jerk_pert = Vec::Zero();
		
		// Calculate the Newtonian accs and jerks
		for (Particle &j : all)
		{
			if (i == j)
				continue;
			
			if (IsRegularized(i, j))
				continue;

			Vec r_ij = i.pos - j.pos;
			Vec v_ij = i.vel - j.vel;

			Real r2 = r_ij.SquaredNorm();			
			Real m_over_r3 = j.mass / (r2*sqrt(r2));

			i.acc_newton -= m_over_r3 * r_ij;
			i.jerk_newton -= m_over_r3 * (v_ij - (3*r_ij.Dot(v_ij) / r2) * r_ij);
		}
	}
	
	// Calculate the perturbing accs and jerks
	m_Perturbation->CalculateAccelerationsAndJerks(m_Model, dt, 0, 1);
	
	// Evaluate Regularized Particles
	for (RegularizedPair &p : m_RegularizedPairs)
	{
		p.Evaluate();
	}
}

void RegularizedHermiteIntegrator::Correct(Real dt)
{
	ParticleSetView &all = m_Model->GetAllParticles();
	
	Real dt_over_2 = dt / 2;
	Real dt2_over_12 = dt*dt / 12;
	
	for (Particle &p : all)
	{
		if (IsRegularized(p))
			continue;
		
		p.vel = p.old_vel 
			+ (p.old_acc_newton + p.old_acc_pert + p.acc_newton + p.acc_pert) * dt_over_2
			+ (p.old_jerk_newton + p.old_jerk_pert - p.jerk_newton - p.jerk_pert) * dt2_over_12;

		p.pos = p.old_pos 
			+ (p.old_vel + p.vel) * dt_over_2
			+ (p.old_acc_newton + p.old_acc_pert - p.acc_newton - p.acc_pert) * dt2_over_12;
	}
	
	for (RegularizedPair &p : m_RegularizedPairs)
	{
		Real ds = p.ToRegularizedTime(dt);

		p.V = p.old_V + (p.old_acc + p.acc) * ds / 2 + (p.old_jerk - p.jerk) * ds * ds / 12;
		p.U = p.old_U + (p.old_V + p.V) * ds / 2 + (p.old_acc - p.acc) * ds * ds / 12;

		p.vel_cm = p.old_vel_cm + (p.old_acc_cm + p.acc_cm) * dt_over_2 + (p.old_jerk_cm - p.jerk_cm) * dt2_over_12;
		p.pos_cm = p.old_pos_cm + (p.old_vel_cm + p.vel_cm) * dt_over_2 + (p.old_acc_cm - p.acc_cm) * dt2_over_12;

		p.h = p.old_h + (p.old_h_prime + p.h_prime) * ds / 2 + (p.old_h_prime2 - p.h_prime2) * ds*ds / 12;
		
		p.UpdateParticles();
	}
}

Real RegularizedHermiteIntegrator::CalculateTimeStep()
{
	ParticleSetView &all = m_Model->GetAllParticles();
	Real coll_time_q = 1e300;

	for (Particle &i : all)
	{
		for (Particle &j : all)
		{
			// Ignore identical particles.
			if (i == j)
				continue;
			
			if (IsRegularized(i, j))
				continue;
			
			Real r2 = (j.pos - i.pos).SquaredNorm();
			Real v2 = (j.vel - i.vel).SquaredNorm();

			// First estimate: unaccelerated linear motion.
			Real coll_time_q_est = (r2*r2) / (v2*v2);
			coll_time_q = min(coll_time_q, coll_time_q_est);

			Real a2 = (j.acc_newton + j.acc_pert - i.acc_newton - i.acc_pert).SquaredNorm();
			Real m = i.mass + j.mass;

			// Second estimate: free fall time.
			coll_time_q_est = r2 / (a2 * (m * m));
			coll_time_q = min(coll_time_q, coll_time_q_est);
		}
	}
	Real dt_other = m_DtParam * sqrt(sqrt(coll_time_q));
	
	Real dt_regularized = 1e300;
	for (RegularizedPair &pair : m_RegularizedPairs)
	{
		//Real ds = m_DtParam * sqrt((pair.acc.Norm() * pair.U.Norm() + pair.V.SquaredNorm()) / (pair.jerk.Norm() * pair.V.Norm() + pair.acc.SquaredNorm()));
		Real temp = (pair.acc.Norm() * pair.U.Norm()) / (pair.jerk.Norm() * pair.V.Norm());
		temp = min(temp, pair.V.SquaredNorm() / pair.acc.SquaredNorm());
		
		Real ds = m_DtParam * sqrt(temp);

		Real dt = pair.ToRealTime(ds);
		
		dt_regularized = min(dt_regularized, dt);
	}
	
	return min(dt_other, dt_regularized);
}
