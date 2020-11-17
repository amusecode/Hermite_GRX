#include "EIHPerturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>

using namespace std;

void EIHPerturbation::CalculateAccelerations(Model *model, int thread_id, int num_threads)
{
	ParticleSetView &all = model->GetAllParticles();
	Real c2_recipr = 1 / (m_SpeedOfLight * m_SpeedOfLight);
	
	for (Particle &a : Strided(all, thread_id, num_threads))
	{
		for (Particle &b : all)
		{
			if (a == b)
				continue;
			
			Vec x_ab = a.pos - b.pos;
			Vec v_ab = a.vel - b.vel;
			
			Real r_ab2 = x_ab.SquaredNorm();
			Real r_ab = sqrt(r_ab2);
			Real r_ab3 = r_ab * r_ab2;
			
			Vec n_ab = x_ab / r_ab;
			Real m_a_over_r_ab = a.mass / r_ab;
			Real m_b_over_r_ab = b.mass / r_ab;
			Real m_b_over_r_ab3 = b.mass / r_ab3;
			
			Real v_b_dot_n_ab = b.vel.Dot(n_ab);
			
			a.acc_pert += c2_recipr * v_ab * (m_b_over_r_ab3 * x_ab.Dot(4*a.vel - 3*b.vel));
			
			Real temp = 4 * m_b_over_r_ab;
			temp += 5 * m_a_over_r_ab;
			temp -= a.vel.SquaredNorm();
			temp += 4 * a.vel.Dot(b.vel);
			temp -= 2 * b.vel.SquaredNorm();
			temp += 1.5 * v_b_dot_n_ab * v_b_dot_n_ab;
			
			for (Particle &c : all)
			{
				if (a == c || b == c)
					continue;
				
				Vec x_bc = b.pos - c.pos;
				Vec x_ac = a.pos - c.pos;
				
				Real r_bc2 = x_bc.SquaredNorm();
				Real r_bc = sqrt(r_bc2);
				Real r_bc3 = r_bc2 * r_bc;
				Real r_ac = x_ac.Norm();
				
				Real m_c_over_r_bc3 = c.mass / r_bc3;
				
				temp += c.mass / r_bc;
				temp += 4 * c.mass / r_ac;
				temp -= 0.5 * m_c_over_r_bc3 * x_ab.Dot(x_bc);
				
				a.acc_pert -= c2_recipr * x_bc * (3.5 * m_b_over_r_ab * m_c_over_r_bc3);
			}
			
			a.acc_pert += c2_recipr * x_ab * (m_b_over_r_ab3 * temp);		
		}
	}
}

Real EIHPerturbation::GetEnergy(Model *model)
{
	ParticleSetView &all = model->GetAllParticles();
	Real energy = 0;
	
	for (Particle &a : all)
	{
		Vec v_a = a.vel;
		Real v_a2 = v_a.SquaredNorm();
		
		Real energypart = 0.375 * v_a2*v_a2;
		for (Particle &b : all)
		{
			if (a == b)
				continue;
			
			Real r_ab = sqrt((a.pos - b.pos).SquaredNorm());
			Real m_b_over_r_ab = b.mass / r_ab;
			Vec n_ab = (a.pos - b.pos) / r_ab;
			
			energypart += 1.5 * v_a2 * m_b_over_r_ab;
			
			Real temp = 7*a.vel.Dot(b.vel);
			temp += a.vel.Dot(n_ab) * b.vel.Dot(n_ab);
			energypart -= 0.25 * m_b_over_r_ab * temp;
			
			for (Particle &c : all)
			{
				// Note: it could be b == c, which is not clear in Will (2014a)
				if (a == c)
					continue;
				
				Real r_ac = (a.pos - c.pos).Norm();
				energypart += 0.5 * m_b_over_r_ab * c.mass / r_ac;
			}
		}
		
		energy += energypart * a.mass;
	}
	
	energy /= m_SpeedOfLight * m_SpeedOfLight;

	return energy;
}

Vec EIHPerturbation::GetLinearMomentum(Model *model)
{
	// TODO
	ParticleSetView &all = model->GetAllParticles();
	Vec momentum = Vec::Zero();
	
	for (Particle &p : all)
	{
		momentum += p.mass * p.vel;
	}
	
	return momentum;
}
