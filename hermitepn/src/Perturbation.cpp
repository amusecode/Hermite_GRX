#include "Perturbation.h"
#include "StridedContainer.h"
#include "Model.h"

#include <cmath>

using namespace std;

void Perturbation::CalculateAccelerations(Model *model, int thread_id, int num_threads)
{
}

void Perturbation::CalculateJerks(Model *model, Real dt, int thread_id, int num_threads)
{
}

void Perturbation::CalculateAccelerationsAndJerks(Model *model, Real dt, int thread_id, int num_threads)
{
	CalculateAccelerations(model, thread_id, num_threads);
	CalculateJerks(model, dt, thread_id, num_threads);
}

Real Perturbation::GetEnergy(Model *model)
{
	return 0;
}

Vec Perturbation::GetLinearMomentum(Model *model)
{
	return Vec::Zero();
}

Real Perturbation::GetSpeedOfLight()
{
	return m_SpeedOfLight;
}

void Perturbation::SetSpeedOfLight(Real speed_of_light)
{
	m_SpeedOfLight = speed_of_light;
}
