#include "Pairwise1PNPerturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>
#include <iostream>

using namespace std;

void Pairwise1PNPerturbation::CalculateAccelerations(Model *model, int thread_id, int num_threads)
{
	ParticleSetView &all = model->GetAllParticles();
	Real c2_recipr = 1 / (m_SpeedOfLight * m_SpeedOfLight);
	
	for (Particle &i : Strided(all, thread_id, num_threads))
	{
		Real v_i2 = i.vel.SquaredNorm();
		
		for (Particle &j : all)
		{
			if (i == j)
				continue;

			Vec r_ij = i.pos - j.pos;
			Vec v_ij = i.vel - j.vel;

			Real inv_r2 = 1 / r_ij.SquaredNorm();
			Real inv_r = sqrt(inv_r2);

			Real m_i_over_r = i.mass * inv_r;
			Real m_j_over_r = j.mass * inv_r;
			Real m_j_over_r2 = j.mass * inv_r2;

			Real v_j2 = j.vel.SquaredNorm();
			Real v_i_dot_v_j = i.vel.Dot(j.vel);

			Vec n_ij = r_ij * inv_r;
			Real n_ij_dot_v_i = n_ij.Dot(i.vel);
			Real n_ij_dot_v_j = n_ij.Dot(j.vel);

			i.acc_pert += c2_recipr * m_j_over_r2 * (-v_i2 - 2 * v_j2 + 4 * v_i_dot_v_j + 1.5 * n_ij_dot_v_j*n_ij_dot_v_j + 5 * m_i_over_r + 4 * m_j_over_r) * n_ij;
			i.acc_pert += c2_recipr * m_j_over_r2 * (4 * n_ij_dot_v_i - 3 * n_ij_dot_v_j) * v_ij;
		}
	}
}

Real Pairwise1PNPerturbation::GetEnergy(Model *model)
{
	ParticleSetView &all = model->GetAllParticles();
	Real energy = 0;
	
	for (Particle &a : all)
	{
		Real v_a2 = a.vel.SquaredNorm();
		
		Real energypart = 0.375 * v_a2*v_a2;
		
		for (Particle &b : all)
		{
			if (a == b)
				continue;
			
			Real r_ab = (a.pos - b.pos).Norm();
			Real m_a_over_r_ab = a.mass / r_ab;
			Real m_b_over_r_ab = b.mass / r_ab;
			Vec n_ab = (a.pos - b.pos) / r_ab;
			
			energypart += 0.5 * m_a_over_r_ab*m_b_over_r_ab;
			
			Real temp = 1.5 * v_a2;
			temp -= 1.75*a.vel.Dot(b.vel);
			temp -= 0.25*a.vel.Dot(n_ab) * b.vel.Dot(n_ab);
			energypart += m_b_over_r_ab * temp;
		}
		
		energy += energypart * a.mass;
	}
	
	energy /= m_SpeedOfLight * m_SpeedOfLight;

	return energy;
}

Vec Pairwise1PNPerturbation::GetLinearMomentum(Model *model)
{
	// TODO
	return Vec::Zero();
}
