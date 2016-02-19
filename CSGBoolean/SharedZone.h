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

        /*  �ϲ����ȼ���
        �Ѿ��Ǽ�Ϊ�����ģ���˳���Ŵ�С���ȵ���
        */
        void setupPonits(MyMesh::Face_handle fh0, MyMesh::Face_handle fh1, const myext::TriTriIsectResult<K>& result);

    private:
        const std::vector<MyMesh*>* meshList = nullptr;
        //PointList                   pointEntity;
    };

}
