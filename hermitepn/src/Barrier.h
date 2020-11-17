#include <mutex>
#include <condition_variable>
#include <iostream>

class Barrier
{
public:
	inline Barrier(size_t num_threads)
		: m_NumThreads(num_threads), m_Counter(num_threads), m_IsEnforced(true)
	{
	}
	
	inline ~Barrier()
	{
		Enforce(false);
	}
	
	inline void Wait()
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		
		if (!m_IsEnforced)
			return;
		
		if (--m_Counter)
		{
			m_Condition.wait(lock);
		}
		else
		{
			m_Condition.notify_all();
			m_Counter = m_NumThreads;
		}
	}
	
	inline void Release()
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		
		m_Condition.notify_all();
		m_Counter = m_NumThreads;
	}
	
	inline void Enforce(bool enforced)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);
		
		if (enforced)
		{
			m_IsEnforced = true;
		}
		else
		{
			m_Condition.notify_all();
			m_Counter = m_NumThreads;
			m_IsEnforced = false;
		}
	}

	inline void SetNumThreads(size_t num_threads)
	{
		std::unique_lock<std::mutex> lock(m_Mutex);

		if (m_NumThreads - m_Counter > num_threads)
		{
			m_Condition.notify_all();
			m_Counter = m_NumThreads;
		}

		m_Counter += num_threads - m_NumThreads;
		m_NumThreads = num_threads;
	}
	
private:
	size_t m_NumThreads;
	size_t m_Counter;
	bool m_IsEnforced;
	
	std::condition_variable m_Condition;
	std::mutex m_Mutex;
};
