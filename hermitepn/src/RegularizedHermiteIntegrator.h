#include "Integrator.h"

#include "StoppingConditions.h"
#include "Quat.h"

#include <list>

struct RegularizedPair
{
	RegularizedPair(Particle &part1, Particle &part2);
	
	bool Contains(const Particle &part1, const Particle &part2) const;
	bool Contains(const Particle &p) const;
	void Save();
	
	void Evaluate();
	void UpdateParticles();

	Real ToRealTime(Real ds);
	Real ToRegularizedTime(Real dt);
	
	Particle *p1, *p2;

	Real h;
	Real h_prime, h_prime2;

	Real old_h;
	Real old_h_prime, old_h_prime2;

	Quat U, V;
	Quat acc, jerk;
	
	Vec pos_cm, vel_cm;
	Vec acc_cm, jerk_cm;
	
	Quat old_U, old_V;
	Quat old_acc, old_jerk;
	
	Vec old_pos_cm, old_vel_cm;
	Vec old_acc_cm, old_jerk_cm;
};

class RegularizedHermiteIntegrator : public Integrator
{
public:
	RegularizedHermiteIntegrator(bool symmetrized, size_t num_threads = 1);
	virtual ~RegularizedHermiteIntegrator();
	
private:
	bool IsRegularized(const Particle &p1, const Particle &p2) const;
	bool IsRegularized(const Particle &p) const;

	bool SelectRegularizedParticles();
	
	void EvolveWorker(size_t thread_id);
	void NumThreadsChanged();
	
	void Predict(Real dt);
	void Evaluate(Real dt);
	void Correct(Real dt);
	
	Real CalculateTimeStep();
	
	bool m_Symmetrized;
	
	StoppingConditions m_StoppingCond;
	bool m_ConditionsSet;
	
	std::list<RegularizedPair> m_RegularizedPairs;
};


