#pragma once
#include <vector>

#include "csgdefs.h"
#include "CGALext.h"
#include "MyMesh.h"


namespace CSG
{
    //typedef std::list<CGAL::Point_3<K>> PointList;

    /* check if it is setup */

    class SharedZone
    {
    public:
        SharedZone(std::vector<MyMesh*> &meshList){}
        ~SharedZone(){}

        void setupIsectFacet(MyMesh::Face_handle fh);

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void setupPonits(MyMesh::Face_handle fh0, MyMesh::Face_handle fh1, const myext::TriTriIsectResult<K>& result);

    private:
        const std::vector<MyMesh*>* meshList = nullptr;
        //PointList                   pointEntity;
    };

}
