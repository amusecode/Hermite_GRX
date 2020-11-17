#ifndef VEC_H
#define VEC_H

#include "Common.h"

#include <iostream>
#include <cmath>

class Vec
{
public:
	Vec()
	{
	}

	Vec(const Vec &vec)
	{
		m_Data[0] = vec.m_Data[0];
		m_Data[1] = vec.m_Data[1];
		m_Data[2] = vec.m_Data[2];
	}

	Vec(Real c)
	{
		m_Data[0] = m_Data[1] = m_Data[2] = c;
	}

	Vec(Real x, Real y, Real z)
	{
		m_Data[0] = x;
		m_Data[1] = y;
		m_Data[2] = z;
	}

	Vec operator + (const Vec &b) const
	{
		Vec res(*this);
		res += b;
		return res;
	}

	Vec &operator += (const Vec &b)
	{
		m_Data[0] += b.m_Data[0];
		m_Data[1] += b.m_Data[1];
		m_Data[2] += b.m_Data[2];
		return *this;
	}
	
	Vec operator - (const Vec &b) const
	{
		Vec res(*this);
		res -= b;
		return res;
	}

	Vec &operator -= (const Vec &b)
	{
		m_Data[0] -= b.m_Data[0];
		m_Data[1] -= b.m_Data[1];
		m_Data[2] -= b.m_Data[2];
		return *this;
	}

	Vec operator * (const Real &b) const
	{
		Vec res(*this);
		res *= b;
		return res;
	}

	Vec &operator *= (const Real &b)
	{
		m_Data[0] *= b;
		m_Data[1] *= b;
		m_Data[2] *= b;
		return *this;
	}

	Vec operator / (const Real &b) const
	{
		Vec res(*this);
		res /= b;
		return res;
	}

	Vec &operator /= (const Real &b)
	{
		m_Data[0] /= b;
		m_Data[1] /= b;
		m_Data[2] /= b;
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
	
	Real Dot(const Vec &b)
	{
		return m_Data[0] * b.m_Data[0] + m_Data[1] * b.m_Data[1] + m_Data[2] * b.m_Data[2];
	}

	Real SquaredNorm()
	{
		return Dot(*this);
	}

	Real Norm()
	{
		return sqrt(SquaredNorm());
	}

	static Vec Zero()
	{
		return Vec(0);
	}

private:
	Real m_Data[3];
};

inline Vec operator * (Real a, const Vec &b)
{
	return b * a;
}

inline std::ostream &operator << (std::ostream &out, const Vec &v)
{
	out << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
	return out;
}

#endif // VEC_H
