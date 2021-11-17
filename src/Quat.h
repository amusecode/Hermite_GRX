#ifndef QUATERNION_H
#define QUATERNION_H

#include "Common.h"

#include <cmath>
#include <iostream>

class Quat
{
public:
	Quat()
	{
	}
	
	Quat(const Vec &vec)
	{
		m_Data[0] = vec[0];
		m_Data[1] = vec[1];
		m_Data[2] = vec[2];
		m_Data[3] = 0;
	}

	Quat(const Quat &quat)
	{
		m_Data[0] = quat.m_Data[0];
		m_Data[1] = quat.m_Data[1];
		m_Data[2] = quat.m_Data[2];
		m_Data[3] = quat.m_Data[3];
	}
	
	Quat(Real c)
	{
		m_Data[0] = c;
		m_Data[1] = m_Data[2] = m_Data[3] = 0;
	}

	Quat(Real q0, Real q1, Real q2, Real q3)
	{
		m_Data[0] = q0;
		m_Data[1] = q1;
		m_Data[2] = q2;
		m_Data[3] = q3;
	}
	
	Quat operator + (const Real &b) const
	{
		Quat res(*this);
		res += b;
		return res;
	}
	
	Quat &operator += (const Real &b)
	{
		m_Data[0] += b;
		return *this;
	}

	Quat operator + (const Quat &b) const
	{
		Quat res(*this);
		res += b;
		return res;
	}

	Quat &operator += (const Quat &b)
	{
		m_Data[0] += b.m_Data[0];
		m_Data[1] += b.m_Data[1];
		m_Data[2] += b.m_Data[2];
		m_Data[3] += b.m_Data[3];
		return *this;
	}
	
	Quat operator - (const Real &b) const
	{
		Quat res(*this);
		res -= b;
		return res;
	}
	
	Quat &operator -= (const Real &b)
	{
		m_Data[0] -= b;
		return *this;
	}

	Quat operator - (const Quat &b) const
	{
		Quat res(*this);
		res -= b;
		return res;
	}

	Quat &operator -= (const Quat &b)
	{
		m_Data[0] -= b.m_Data[0];
		m_Data[1] -= b.m_Data[1];
		m_Data[2] -= b.m_Data[2];
		m_Data[3] -= b.m_Data[3];
		return *this;
	}
	
	Quat operator - () const
	{
		return Quat(-m_Data[0], -m_Data[1], -m_Data[2], -m_Data[3]);
	}
	
	Quat operator * (const Real &b) const
	{
		Quat res(*this);
		res *= b;
		return res;
	}
	
	Quat &operator *= (const Real &b)
	{
		m_Data[0] *= b;
		m_Data[1] *= b;
		m_Data[2] *= b;
		m_Data[3] *= b;
		return *this;
	}
	
	Quat operator * (const Quat &b) const
	{
		Quat res(*this);
		res *= b;
		return res;
	}
	
	Quat &operator *= (const Quat &b)
	{
		Real t0, t1, t2;
		t0 = m_Data[0]*b.m_Data[0] - m_Data[1]*b.m_Data[1] - m_Data[2]*b.m_Data[2] - m_Data[3]*b.m_Data[3];
		t1 = m_Data[0]*b.m_Data[1] + m_Data[1]*b.m_Data[0] + m_Data[2]*b.m_Data[3] - m_Data[3]*b.m_Data[2];
		t2 = m_Data[0]*b.m_Data[2] - m_Data[1]*b.m_Data[3] + m_Data[2]*b.m_Data[0] + m_Data[3]*b.m_Data[1];
		
		m_Data[3] = m_Data[0]*b.m_Data[3] + m_Data[1]*b.m_Data[2] - m_Data[2]*b.m_Data[1] + m_Data[3]*b.m_Data[0];
		m_Data[0] = t0;
		m_Data[1] = t1;
		m_Data[2] = t2;
		return *this;
	}
	
	Quat operator / (const Real &b) const
	{
		Quat res(*this);
		res /= b;
		return res;
	}
	
	Quat &operator /= (const Real &b)
	{
		m_Data[0] /= b;
		m_Data[1] /= b;
		m_Data[2] /= b;
		m_Data[3] /= b;
		return *this;
	}
	
	Real &operator [] (int i)
	{
		return m_Data[i];
	}
	
	const Real &operator [] (int i) const
	{
		return m_Data[i];
	}
	
	Real SquaredNorm()
	{
		return m_Data[0] * m_Data[0] + m_Data[1] * m_Data[1] + m_Data[2] * m_Data[2] + m_Data[3] * m_Data[3];
	}
	
	Real Norm()
	{
		return std::sqrt(SquaredNorm());
	}
	
	Quat Conj()
	{
		return Quat(m_Data[0], -m_Data[1], -m_Data[2], -m_Data[3]);
	}
	
	Quat StarConj()
	{
		return Quat(m_Data[0], m_Data[1], m_Data[2], -m_Data[3]);
	}
	
	Vec GetVecPart()
	{
		return Vec(m_Data[0], m_Data[1], m_Data[2]);
	}
	
	static Quat Zero()
	{
		return Quat(0, 0, 0, 0);
	}
	
private:
	Real m_Data[4];
};

inline Quat operator + (Real a, const Quat &b)
{
	return b + a;
}

inline Quat operator - (Real a, const Quat &b)
{
	return -b + a;
}

inline Quat operator * (Real a, const Quat &b)
{
	return b * a;
}

inline std::ostream &operator << (std::ostream &out, const Quat &q)
{
	out << "[" << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << "]";
	return out;
}

#endif // QUATERNION_H
