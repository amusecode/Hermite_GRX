#include "StoppingConditions.h"

#include "Model.h"

#include <stopcond.h>

using namespace std;
using namespace chrono;

StoppingConditions::StoppingConditions()
{
	set_support_for_condition(COLLISION_DETECTION);
	set_support_for_condition(TIMEOUT_DETECTION);
	set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
}

void StoppingConditions::Reset()
{
	reset_stopping_conditions();
	m_NumSteps = 0;
	m_TimeStart = steady_clock::now();
}

bool StoppingConditions::Check(Model *model)
{
	int timeout_detection, num_steps_detection, collision_detection;
	
	is_stopping_condition_enabled(TIMEOUT_DETECTION, &timeout_detection);
	is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, &num_steps_detection);
	is_stopping_condition_enabled(COLLISION_DETECTION, &collision_detection);

	if (timeout_detection)
	{
		double max_time;
		get_stopping_condition_timeout_parameter(&max_time);

		steady_clock::time_point now = steady_clock::now();
		auto duration = duration_cast<seconds>(now - m_TimeStart);
		if (duration.count() > max_time)
		{
			int stopping_index = next_index_for_stopping_condition();
			set_stopping_condition_info(stopping_index, TIMEOUT_DETECTION);
		}
	}

	if (num_steps_detection)
	{
		int max_steps;
		get_stopping_condition_number_of_steps_parameter(&max_steps);

		m_NumSteps++;
		if (m_NumSteps > max_steps)
		{
			int stopping_index = next_index_for_stopping_condition();
			set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION);
		}
	}

	if (collision_detection)
	{
		ParticleSetView &all = model->GetAllParticles();

		for (Particle &i : all)
		{
			for (Particle &j : all)
			{
				if (i.id >= j.id)
					continue;

				Real r2 = (i.pos - j.pos).SquaredNorm();
				Real sum_rad = i.radius + j.radius;

				if (r2 < sum_rad*sum_rad)
				{
					int stopping_index = next_index_for_stopping_condition();
					set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
					set_stopping_condition_particle_index(stopping_index, 0, i.id);
					set_stopping_condition_particle_index(stopping_index, 1, j.id);
				}
			}
		}
	}

	return set_conditions;
}
