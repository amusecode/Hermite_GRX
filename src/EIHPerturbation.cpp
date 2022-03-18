#include "EIHPerturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>

using namespace std;

// This is the computation of the full EIH equations in 1PN limit. The equations are based
// on (Will 2014) section 3. See Arxiv:1312.1289v3
// The higher order terms that are added are found in REF

void EIHPerturbation::CalculateAccelerations(Model *model, int thread_id, int num_threads)
{
	ParticleSetView &all = model->GetAllParticles();
	Real c2_recipr = 1 / (m_SpeedOfLight * m_SpeedOfLight);
	Real c4_recipr = c2_recipr * c2_recipr;
	Real c5_recipr = c4_recipr / m_SpeedOfLight;
	// The primary double loop includes all direct interactions between the masses
	for (Particle &a : Strided(all, thread_id, num_threads))
	{
		for (Particle &b : all)
		{
			// No summation over the particle interacting with itself
			if (a == b)
				continue;

			//
			// Compute the differences in position and velocity
			Vec x_ab = a.pos - b.pos;
			Vec v_ab = a.vel - b.vel;
			// compute the distances and its powers
			Real r_ab2 = x_ab.SquaredNorm(); // PN0.0
			Real r_ab = sqrt(r_ab2);		 // PN1.0
			Real r_ab3 = r_ab * r_ab2;		 // PN1.0
			// Compute the normal vector, and mass to distance ratios
			Vec n_ab = x_ab / r_ab;				// PN0.0
			Real m_a_over_r_ab = a.mass / r_ab; // PN1.0
			Real m_b_over_r_ab = b.mass / r_ab; // PN1.0
			// Real m_b_over_r_ab3 = b.mass / r_ab3; // PN1.0
			Real m_b_over_r_ab2 = b.mass / r_ab2; // PN1.0

			// Compute the inner products used in the equations
			Real v_b_dot_n_ab = b.vel.Dot(n_ab);	// PN1.0
			Real v_a_dot_v_b = a.vel.Dot(b.vel);	// PN1.0
			Real v_a_dot_v_a = a.vel.SquaredNorm(); // PN1.0
			Real v_b_dot_v_b = b.vel.SquaredNorm(); // PN1.0
			// Real square_v_a_dot_v_a = v_a_dot_v_a * v_a_dot_v_a;				   // PN1.0
			Real square_v_b_dot_v_b = b.vel.SquaredNorm() * b.vel.SquaredNorm();   // PN1.0
			Real v_a_dot_n_ab = a.vel.Dot(n_ab);								   // PN2.0
			Real square_v_a_dot_v_b = a.vel.Dot(b.vel) * a.vel.Dot(b.vel);		   // PN2.0
			Real square_v_b_dot_n_ab = v_b_dot_n_ab * v_b_dot_n_ab;				   // PN2.0
			Real square_v_a_dot_n_ab = v_a_dot_n_ab * v_a_dot_n_ab;				   // PN2.0
			Real cubed_v_b_dot_n_ab = square_v_b_dot_n_ab * v_b_dot_n_ab;		   // PN2.0
			Real quarted_v_b_dot_n_ab = square_v_b_dot_n_ab * square_v_b_dot_n_ab; // PN2.0

			// Where does the PN 0.0 term apply?
			// PN 1.0
			Real temp = 5.0 * m_a_over_r_ab;
			temp += 4.0 * m_b_over_r_ab;
			temp += 1.5 * v_a_dot_n_ab;
			temp -= v_a_dot_v_a;
			temp += 4.0 * v_a_dot_v_b;
			temp -= 2.0 * v_b_dot_v_b;

			// Include the cross terms between all particle a and particles b and c.
			for (Particle &c : all)
			{
				// No summation over the particles interacting with themselves
				if (a == c || b == c)
					continue;

				Vec x_bc = b.pos - c.pos;
				Vec x_ac = a.pos - c.pos;

				Real r_bc2 = x_bc.SquaredNorm();
				Real r_bc = sqrt(r_bc2);
				Real r_bc3 = r_bc2 * r_bc;
				Real r_ac = x_ac.Norm();
				// Vec n_bc = x_bc / r_bc;

				Real m_c_over_r_bc3 = c.mass / r_bc3;
				// Add the 3 cross terms
				temp += c.mass / r_bc;
				temp += 4.0 * c.mass / r_ac;
				temp -= 0.5 * m_c_over_r_bc3 * x_ab.Dot(x_bc);
				// Seperately add the 7/2 cross term (since it does not scale with m_b_over_r_ab3)
				a.acc_pert -= c2_recipr * x_bc * (3.5 * m_b_over_r_ab * m_c_over_r_bc3);
			}
			// Add the 1st PN in the v_ab direction
			a.acc_pert += c2_recipr * v_ab * (m_b_over_r_ab2 * (4.0 * v_a_dot_n_ab - 3.0 * v_b_dot_n_ab));
			// Add the sum of the pair and cross terms in the n_ab direction
			a.acc_pert += c2_recipr * n_ab * (m_b_over_r_ab2 * temp);

			// PN 2.0
			Real inv_r_ab2 = 1.0 / r_ab2;					  // computing inverse and reusing it
			temp = -14.25 * a.mass * a.mass * inv_r_ab2;	  // 57/4
			temp -= 34.5 * a.mass * b.mass * inv_r_ab2;		  // 69/2
			temp -= 9.0 * b.mass * b.mass * inv_r_ab2;		  // 9
			temp -= 1.875 * quarted_v_b_dot_n_ab;			  // 16/8
			temp += 1.5 * square_v_b_dot_n_ab * v_a_dot_v_a;  // 3/2
			temp -= 6.0 * square_v_b_dot_v_b * v_a_dot_v_b;	  // 6
			temp -= 2.0 * square_v_a_dot_v_b;				  // 2
			temp += 4.5 * square_v_b_dot_n_ab * v_b_dot_n_ab; // 9/2
			temp += 4.0 * v_a_dot_v_b * v_b_dot_v_b;		  // 4
			temp -= 2.0 * v_b_dot_v_b;						  // 2
			// Compute longer terms with new dummy variable
			Real temp2 = 19.5 * square_v_a_dot_n_ab;	 // 39/2
			temp2 -= 39.0 * v_a_dot_n_ab * v_b_dot_n_ab; // 39
			temp2 += 8.5 * square_v_b_dot_n_ab;			 // 17/2
			temp2 -= 3.75 * v_a_dot_v_a;				 // 15/4
			temp2 += 2.5 * v_a_dot_v_b;					 // 5/2
			temp2 += 1.25 * v_b_dot_v_b;				 // 5/4
			temp += a.mass / r_ab * temp2;
			// Next long term
			temp2 = 2.0 * square_v_a_dot_n_ab;			// 2
			temp2 -= 4.0 * v_a_dot_n_ab * v_b_dot_n_ab; // 4
			temp2 -= 6.0 * square_v_b_dot_n_ab;			// 6
			temp2 -= 8.0 * v_a_dot_v_b;					// 8
			temp2 += 4.0 * v_b_dot_v_b;					// 4
			temp += b.mass / r_ab * temp2;
			// Add the 2.0 PN contributions in n_ab
			a.acc_pert += c4_recipr * n_ab * (m_b_over_r_ab2 * temp);
			// Compute the v_ab contributions
			temp = m_b_over_r_ab * (-2.0 * v_a_dot_n_ab - 2.0 * v_b_dot_n_ab);		// 2, 2
			temp += m_a_over_r_ab * (-15.75 * v_a_dot_n_ab + 13.75 * v_b_dot_n_ab); // -63/4, 55/4
			temp += 6.0 * v_a_dot_n_ab * square_v_b_dot_n_ab;						// 6
			temp += 4.5 * cubed_v_b_dot_n_ab;										// 9/2
			temp += v_b_dot_n_ab * v_a_dot_v_a;										// 1
			temp -= 4.0 * v_a_dot_n_ab * v_a_dot_v_b;								// -4
			temp += 4.0 * v_b_dot_n_ab * v_a_dot_v_b;								// +4
			temp += 4.0 * v_a_dot_n_ab * v_b_dot_v_b;								// +4
			temp -= 5.0 * v_b_dot_n_ab * v_b_dot_v_b;								// -5
			a.acc_pert += c4_recipr * v_ab * (m_b_over_r_ab2 * temp);

			// PN 2.5
			Real m_a_m_b_over_r_ab3 = a.mass * b.mass / r_ab3;
			Real v_ab_dot_n_ab = v_ab.Dot(n_ab);
			Real v_ab_dot_v_ab = v_ab.SquaredNorm();
			// Compute the n_ab terms for PN 2.5
			temp = 52.0 / 3.0 * b.mass / r_ab * v_ab_dot_n_ab; // 52/3
			temp -= 6.0 * a.mass / r_ab * v_ab_dot_n_ab;	   // 6
			temp += 3.0 * (v_ab_dot_n_ab)*v_ab_dot_v_ab;	   // 3
			a.acc_pert += c5_recipr * n_ab * (0.8 * m_a_m_b_over_r_ab3 * temp);
			temp = 2.0 * a.mass / r_ab;	 // 2
			temp -= 8.0 * b.mass / r_ab; // 8
			temp -= v_ab_dot_v_ab;		 // 1
			a.acc_pert += c5_recipr * v_ab * (0.8 * m_a_m_b_over_r_ab3 * temp);
		}
	}
}

Real EIHPerturbation::GetEnergy(Model *model)
{
	ParticleSetView &all = model->GetAllParticles();
	Real energy = 0;
	Real c2_recipr = 1 / (m_SpeedOfLight * m_SpeedOfLight);
	Real c4_recipr = c2_recipr * c2_recipr;

	for (Particle &a : all)
	{
		Vec v_a = a.vel;
		Real v_a_dot_v_a = v_a.SquaredNorm();
		Real square_v_a_dot_v_a = v_a_dot_v_a * v_a_dot_v_a; // PN1.0
		// 1PN kinetic term
		energy += c2_recipr * a.mass * 0.375 * v_a_dot_v_a * v_a_dot_v_a;
		// 2PN kinetic term
		energy += c4_recipr * a.mass * 0.3125 * v_a_dot_v_a * v_a_dot_v_a * v_a_dot_v_a;
		for (Particle &b : all)
		{
			if (a == b)
				continue;

			// Compute the differences in position and velocity
			Vec x_ab = a.pos - b.pos;
			// Vec v_ab = a.vel - b.vel;
			// compute the distances and its powers
			Real r_ab2 = x_ab.SquaredNorm(); // PN0.0
			Real r_ab = sqrt(r_ab2);		 // PN1.0
			// Real r_ab3 = r_ab * r_ab2;		 // PN1.0
			// Compute the normal vector, and mass to distance ratios
			Vec n_ab = x_ab / r_ab;				// PN0.0
			Real m_a_over_r_ab = a.mass / r_ab; // PN1.0
			// Real ma_am_b_over_r_ab = a.mass * b.mass / r_ab;  // PN1.0
			// Real ma_m_b_over_r_ab3 = a.mass * b.mass / r_ab3; // PN1.0
			// Real m_a_m_b_over_r_ab2 = a.mass * b.mass / r_ab2; // PN1.0
			Real m_b_over_r_ab = b.mass / r_ab; // PN1.0

			// Compute the inner products used in the equations
			Real v_b_dot_n_ab = b.vel.Dot(n_ab);	// PN1.0
			Real v_a_dot_v_b = a.vel.Dot(b.vel);	// PN1.0
			Real v_b_dot_v_b = b.vel.SquaredNorm(); // PN1.0
			// Real square_v_b_dot_v_b = b.vel.SquaredNorm() * b.vel.SquaredNorm();   // PN1.0
			Real v_a_dot_n_ab = a.vel.Dot(n_ab);						  // PN2.0
			Real square_v_b_dot_n_ab = v_b_dot_n_ab * v_b_dot_n_ab;		  // PN2.0
			Real square_v_a_dot_n_ab = v_a_dot_n_ab * v_a_dot_n_ab;		  // PN2.0
			Real cubed_v_a_dot_n_ab = square_v_a_dot_n_ab * v_a_dot_n_ab; // PN2.0
			// Real quarted_v_b_dot_n_ab = square_v_b_dot_n_ab * square_v_b_dot_n_ab; // PN2.0

			// 1PN terms
			Real energypart = 1.5 * v_a_dot_v_a * m_b_over_r_ab;
			Real energytemp = 7.0 * a.vel.Dot(b.vel);
			energytemp += a.vel.Dot(n_ab) * b.vel.Dot(n_ab);
			energypart -= 0.25 * m_b_over_r_ab * energytemp;
			// 1PN cross terms
			for (Particle &c : all)
			{
				// Note: it could be b == c, which is not clear in Will (2014a)
				if (a == c)
					continue;

				Real r_ac = (a.pos - c.pos).Norm();
				energypart += 0.5 * m_b_over_r_ab * c.mass / r_ac;
			}
			// ADD 1PN terms
			energy += c2_recipr * energypart * a.mass;
			// 2PN terms
			energypart = -0.5 * m_a_over_r_ab * m_a_over_r_ab;				  // 1/2
			energypart -= 2.375 * m_a_over_r_ab * m_b_over_r_ab;			  // 9/8
			energypart += 0.375 * cubed_v_a_dot_n_ab * v_b_dot_n_ab;		  // 3/8
			energypart += 0.1875 * square_v_a_dot_n_ab * square_v_b_dot_n_ab; // 3/16
			energypart -= 1.125 * v_a_dot_n_ab * v_b_dot_n_ab * v_a_dot_v_a;  // 9/8
			energypart -= 1.625 * square_v_b_dot_n_ab * v_a_dot_v_a;		  // 13/8
			energypart += 2.625 * square_v_a_dot_v_a;						  // 21/8
			energypart += 1.625 * square_v_a_dot_n_ab * v_a_dot_v_b;		  // 13/8
			energypart += 0.75 * v_a_dot_n_ab * v_b_dot_n_ab * v_a_dot_v_b;	  // 3/4
			energypart -= 6.875 * v_a_dot_v_a * v_a_dot_v_b;				  // 55/8
			energypart += 2.125 * v_a_dot_v_b * v_a_dot_v_b;				  // 17/8
			energypart += 1.9375 * v_a_dot_v_a * v_b_dot_v_b;				  // 31/16
			energytemp = 7.25 * square_v_a_dot_n_ab;						  // 29/4
			energytemp -= 3.25 * v_a_dot_n_ab * v_b_dot_n_ab;				  // 13/4
			energytemp += 0.5 * square_v_b_dot_n_ab;						  // 1/2
			energytemp -= 1.5 * v_a_dot_v_a;								  // 3/2
			energytemp += 1.75 * v_b_dot_v_b;								  // 7/4
			energypart += energytemp * m_a_over_r_ab;
			energy += a.mass * c4_recipr * energypart * m_b_over_r_ab;
		}
	}
	return energy;
}

Vec EIHPerturbation::GetLinearMomentum(Model *model)
{
	// This term is computed up to 1PN with cross terms
	ParticleSetView &all = model->GetAllParticles();
	Vec momentum = Vec::Zero();
	Real c2_recipr = 1 / (m_SpeedOfLight * m_SpeedOfLight);

	for (Particle &a : all)
	{
		// 0PN
		momentum += a.mass * a.vel;
		// 1PN
		Real v_a_dot_v_a = a.vel.SquaredNorm();
		Vec momentumpart = v_a_dot_v_a * a.vel;
		for (Particle &b : all)
		{
			if (b == a)
				continue;
			Vec x_ab = a.pos - b.pos;
			Real r_ab = x_ab.Norm();
			Vec n_ab = x_ab / r_ab;
			momentumpart -= b.mass / r_ab * a.vel;
			momentumpart -= b.mass / r_ab * (a.vel.Dot(n_ab)) * n_ab;
		}
		momentum += 0.5 * c2_recipr * a.mass * momentumpart;
	}
	return momentum;
}
