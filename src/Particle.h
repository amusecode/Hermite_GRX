#ifndef PARTICLE_H
#define PARTICLE_H

#include "Common.h"

#include <vector>
#include <functional>

class Particle
{
public:
	Particle()
	{
	}

	Particle(int i, Real m, const Vec &r, const Vec &v, Real rad)
	{
		id = i;
		mass = m;
		pos = r;
		vel = v;
		radius = rad;
		acc_newton = jerk_newton = acc_pert = jerk_pert = Vec::Zero();
		old_pos = old_vel = old_acc_newton = old_jerk_newton = old_acc_pert = old_jerk_pert = Vec::Zero();
	}

	int id;
	Real mass;
	Real radius;
	Vec pos, vel;
	
	Vec acc_newton, jerk_newton;
	Vec acc_pert, jerk_pert;
	
	Vec old_pos, old_vel;
	Vec old_acc_newton, old_jerk_newton;
	Vec old_acc_pert, old_jerk_pert;
	
	bool operator==(const Particle &b) const
	{
		return id == b.id;
	}

	bool operator!=(const Particle &b) const
	{
		return id != b.id;
	}
	
	void Save()
	{
		old_pos = pos;
		old_vel = vel;
		old_acc_newton = acc_newton;
		old_acc_pert = acc_pert;
		old_jerk_newton = jerk_newton;
		old_jerk_pert = jerk_pert;
	}
};

typedef std::vector<Particle> ParticleSet;
typedef std::vector<std::reference_wrapper<Particle>> ParticleSetView;

inline ParticleSetView MakeParticleSetView(ParticleSet &a, ParticleSet &b)
{
	ParticleSetView res;
	for (Particle &p : a)
		res.push_back(std::ref(p));
	for (Particle &p : b)
		res.push_back(std::ref(p));
	return res;
}

#endif // PARTICLE_H
