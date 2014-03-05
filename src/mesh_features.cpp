#include "mesh_features.h"
using namespace OpenMesh;

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
	// CHECK IF e IS A SILHOUETTE HERE -----------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------
    return false;
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
	// CHECK IF e IS SHARP HERE ------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------------------------------

    return false;
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
	return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

