#include "ItstGraph.h"


namespace CSG
{
    ItstGraph::ItstGraph(FH fh)
    {
        m_bValid = false;

        if (fh->isSimple())
        {
            for (int i = 0; i < 3; i++)
            {
                if (fh->edges[i]->data && fh->edges[i]->data->vertices.size())
                    m_bValid = true;
            }

            if (!m_bValid)
            {
                ReportError("Cannot build a graph!");
                return;
            }
        }


    }


    ItstGraph::~ItstGraph()
    {
    }
}
