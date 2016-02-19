#pragma once
#include <cstring>
#include <utility>

namespace myext
{
    template <class T>
    class TriangleTable
    {
    public:
        TriangleTable(size_t n) : m_sz(n) 
        {
            m_datasz = n * (n + 1) / 2;
            m_data = new T[m_datasz];
            memset(m_data, 0, sizeof(T) * m_datasz);
        }

        ~TriangleTable() { delete[] m_data; }

        T& getValue(size_t i, size_t j)
        {
            if (i > j) std::swap(i, j);
            size_t d = j - i;
            size_t idx = (2 * n + 1 - i) * i / 2 + d;
            return T[idx]
        }

        const T& getValue(size_t i, size_t j) const
        {
            return getValue(i, j);
        }

    private:
        T *m_data = nullptr;
        size_t m_sz = 0;
        size_t m_datasz = 0;
    };
}


