#pragma once
#include <CGAL\Polyhedron_3.h>

namespace CSG
{
    template <class Refs, class Plane>
    struct FaceWithId : public CGAL::HalfedgeDS_face_base<Refs, CGAL::Tag_true, Plane> {
        int idx;
        FaceWithId() {}
    };

    struct Items_FaceWithId : public CGAL::Polyhedron_items_3 {
        template <class Refs, class Traits>
        struct Face_wrapper {
            typedef typename Traits::Plane_3 Plane;
            typedef FaceWithId<Refs, Plane> Face;
        };
    };

    template <class Kernel>
    class HalfedgeMesh :
        public CGAL::Polyhedron_3 < Kernel, Items_FaceWithId >
    {};
}

