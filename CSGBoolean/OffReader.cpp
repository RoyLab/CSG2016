#include "OffReader.h"
#include "stringmanip.h"
#include "macro_util.h"
#include <cassert>
#include <sstream>


namespace CSG
{
    COffReader::COffReader(const char* fileName)
    {
        m_file.open(fileName);
        readFromMyFile();
    }

    COffReader::COffReader(const wchar_t* fileName)
    {
        m_file.open(fileName);
        readFromMyFile();
    }

    bool COffReader::readFromMyFile()
    {
        std::string line;

        if (m_file)
        {
            getNextValidLine(line);
            myext::trim(line);
            if (myext::strcmp_nocase(line, "off") != 0)
            {
                assert(0);
                m_state = INVALID;
                return false;
            }
            else
            {
                getNextValidLine(line);
                std::istringstream ss(line);
                ss >> m_nv;
                ss >> m_nf;
                ss >> m_ne;

                if (ss.good())
                {
                    m_state = READY;
                    m_vstart = m_file.tellg();
                    return true;
                }
                else
                {
                    m_state = INVALID;
                    return false;
                }
            }
        }
    }


    void COffReader::release()
    {
        m_file.close();
        m_state = INVALID;
        m_vstart = -1;
        m_fstart = -1;
    }

    bool COffReader::getFirstVertex(float* v)
    {
        bool hr;
        V_RETURN(m_vstart == READY);

        m_file.seekg(m_vstart);
        m_state = VERTEX;

        return getNextVertex(v);
    }

    bool COffReader::getNextVertex(float* v)
    {
        std::string line;
        bool hr;

        V_RETURN(getNextValidLine(line));

        float x;
        std::istringstream ss(line);

        for (size_t i = 0; i < 3; i++)
        {
            ss >> x;
            v[i] = x;
        }

        return true;
    }

    bool COffReader::getFirstFacet(int* ids, int* count, float* colors)
    {
    }

    bool COffReader::getNextFacet(int* ids, int* count, float* colors)
    {

    }

}