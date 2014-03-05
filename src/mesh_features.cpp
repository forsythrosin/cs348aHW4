#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "mesh_features.h"
using namespace OpenMesh;
using namespace Eigen;

bool isSilhouette(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos)  {
    // CHECK IF e IS A SILHOUETTE HERE ------
    Mesh::VertexHandle vh = mesh.to_vertex_handle(mesh.halfedge_handle(e,0));
    Vec3f tc = cameraPos - mesh.point(vh);
    Vector3d toCamera(tc.values_[0],tc.values_[1],tc.values_[2]); 

    Mesh::FaceHandle fh0 = mesh.face_handle(mesh.halfedge_handle(e,0));
    Mesh::FaceHandle fh1 = mesh.face_handle(mesh.halfedge_handle(e,1));

    Vec3f mesh_n0 = mesh.normal(fh0);
    Vector3d normal0(mesh_n0.values_[0],mesh_n0.values_[1],mesh_n0.values_[2]); 
    Vec3f mesh_n1 = mesh.normal(fh1);
    Vector3d normal1(mesh_n1.values_[0],mesh_n1.values_[1],mesh_n1.values_[2]); 

    return((toCamera.transpose() * normal0) * (toCamera.transpose() * normal1) < 0);
    // --------------------------------------
}

bool isSharpEdge(Mesh &mesh, const Mesh::EdgeHandle &e) {
    // CHECK IF e IS SHARP HERE -------------

    Mesh::FaceHandle fh0 = mesh.face_handle(mesh.halfedge_handle(e,0));
    Mesh::FaceHandle fh1 = mesh.face_handle(mesh.halfedge_handle(e,1));

    Vec3f mesh_n0 = mesh.normal(fh0);
    Vector3d normal0(mesh_n0.values_[0],mesh_n0.values_[1],mesh_n0.values_[2]); 
    Vec3f mesh_n1 = mesh.normal(fh1);
    Vector3d normal1(mesh_n1.values_[0],mesh_n1.values_[1],mesh_n1.values_[2]); 

    return(normal0.transpose() * normal1 < 0.5);
    // --------------------------------------
}

bool isFeatureEdge(Mesh &mesh, const Mesh::EdgeHandle &e, Vec3f cameraPos) {
	return mesh.is_boundary(e) || isSilhouette(mesh,e, cameraPos) || isSharpEdge(mesh,e);
}

