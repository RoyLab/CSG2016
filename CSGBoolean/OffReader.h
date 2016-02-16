#pragma once
#include <fstream>
#include <string>

namespace CSG
{
    class COffReader
    {
        enum STATE { NA, READY, INVALID, END, VERTEX, FACE };
    public:
        COffReader(const char* fileName);
        COffReader(const wchar_t* fileName);
        ~COffReader(){}

        void release();
        int getVertexCount() const { return m_nv; }
        int getFaceCount() const { return m_nf; }

        bool getFirstVertex(float* v);
        bool getNextVertex(float* v);

        bool getFirstFacet(int* ids, int* count = 0, float* colors = 0);
        bool getNextFacet(int* ids, int* count = 0, float* colors = 0);

    private:

        bool readFromMyFile();
        bool getNextValidLine(std::string& line);

        std::ifstream   m_file;
        STATE           m_state = NA;
        int             m_nv = -1, m_nf = -1, m_ne = -1;
        int             m_vstart = -1, m_fstart = -1;
        int             m_vcount = 0, m_fcount = 0;
    };

}

