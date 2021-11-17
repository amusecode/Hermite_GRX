#ifndef STOPPING_CONDITIONS_H
#define STOPPING_CONDITIONS_H

#include <chrono>

class Model;

class StoppingConditions
{
public:
	StoppingConditions();
	
	void Reset();
	bool Check(Model *model);

protected:
	int m_NumSteps;
	std::chrono::steady_clock::time_point m_TimeStart;
};

#endif // STOPPING_CONDITIONS_H
