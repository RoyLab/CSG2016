#define _USE_MATH_DEFINES
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MPMesh.h"
#include <vector>
#include <string>
#include <ext\vld.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

using namespace CSG;
using namespace OpenMesh;

int main()
{
    MPMesh mesh;
    mesh.request_face_colors();

    OpenMesh::IO::Options ropt;
    ropt += IO::Options::FaceColor;

    IO::read_mesh(mesh, "a.off");
    IO::write_mesh(mesh, "b.off");

    system("pause");
    return 0;
}