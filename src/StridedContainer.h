#include <iterator>

template <typename Iterator>
class StridedIterator
{
public:
	inline StridedIterator(Iterator iter, Iterator end, size_t stride)
		: m_Iterator(iter), m_End(end), m_Stride(stride)
	{
	}

	inline StridedIterator<Iterator> &operator++()
	{
		if (std::distance(m_Iterator, m_End) > (int)m_Stride)
			std::advance(m_Iterator, m_Stride);
		else
			m_Iterator = m_End;
		return *this;
	}

	inline typename Iterator::value_type &operator*()
	{
		return *m_Iterator;
	}

	inline bool operator!=(const StridedIterator<Iterator> &it)
	{
		return m_Iterator != it.m_Iterator;
	}

protected:
	Iterator m_Iterator;
	Iterator m_End;
	size_t m_Stride;
};

template <typename Container>
class StridedContainer
{
public:
	StridedContainer(Container &con, size_t start, size_t stride)
		: m_Container(con), m_Start(start), m_Stride(stride)
	{
	}

	inline StridedIterator<typename Container::iterator> begin()
	{
		return StridedIterator<typename Container::iterator>(m_Container.begin() + m_Start, m_Container.end(), m_Stride);
	}

	inline StridedIterator<typename Container::iterator> end()
	{
		return StridedIterator<typename Container::iterator>(m_Container.end(), m_Container.end(), m_Stride);
	}

protected:
	Container &m_Container;
	size_t m_Start;
	size_t m_Stride;
};

template <typename Container>
inline StridedContainer<Container> Strided(Container &con, size_t start, size_t stride)
{
	return StridedContainer<Container>(con, start, stride);
}
