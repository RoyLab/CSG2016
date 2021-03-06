/**  Forward declaration only */
#pragma once
#include <cstdint>
#include <memory>


namespace Boolean
{
    class Octree;
    class BSPTree;

	class XPlane;
	class PlaneLine;
	class PlanePoint;

    class MyVertex;
    class MyEdge;

	class IPolygon;
	class Triangle;
    class SubPolygon;
	class RegularMesh;

	class GlobalData;

    struct NeighborInfo;
    struct FacePbi;
    struct PbiRep;
    typedef PbiRep EdgePbi;

    template <class Pbi> class InsctData;
    template <class Mesh> class CSGTree;

    class EdgeInsctData;
    class FaceInsctData;

    template <class T> struct AutoPtr : std::shared_ptr<T> {};
}
