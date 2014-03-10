#ifndef MESH_DEFINITIONS_H
#define MESH_DEFINITIONS_H

#include <cassert>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

const double EPSILON = 0.00001;

struct MyTraits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> Mesh;

#endif
