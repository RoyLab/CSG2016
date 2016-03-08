#pragma once
#include <cstring>
#include <utility>

namespace myext
{
    enum GraphType { AdjacencyList, Matrix };

    template <class T, GraphType GT = Matrix>
    class UndirectedGraph;

    template <class T>
    class UndirectedGraph<T, Matrix>
    {
    public:
        UndirectedGraph(size_t n) : m_sz(n)
        {
            m_datasz = n * (n + 1) / 2;
            m_data = new T[m_datasz];
            memset(m_data, 0, sizeof(T) * m_datasz);
        }

        ~UndirectedGraph() { delete[] m_data; }

        void setValue(size_t i, size_t j, T val)
        {
            getValue(i, j) = val;
        }

        T& getValue(size_t i, size_t j)
        {
            return m_data[idx(i, j)];
        }

        const T& getValue(size_t i, size_t j) const
        {
            return getValue(i, j);
        }

    private:

        size_t idx(size_t i, size_t j) const
        {
            if (i > j) std::swap(i, j);
            size_t d = j - i;
            return (2 * n + 1 - i) * i / 2 + d;
        }

        T *m_data = nullptr;
        size_t m_sz = 0;
        size_t m_datasz = 0;
    };
}


