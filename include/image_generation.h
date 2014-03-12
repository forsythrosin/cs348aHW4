#ifndef IMAGE_GENERATION_H
#define IMAGE_GENERATION_H

#include "mesh_definitions.h"
#include <string>
#include <list>

//A contour edge of interest for rendering for Part 3...
//Includes any feature edge or suggestive contour.
//Edges with the same endpoints are considered equal.
typedef struct ContourEdge {
  OpenMesh::Vec3f mp1;
  OpenMesh::Vec3f mp2;
  bool silhouette;

  ContourEdge (const OpenMesh::Vec3f& p1, const OpenMesh::Vec3f& p2, bool isSilhouette) {
    mp1 = p1;
    mp2 = p2;
    silhouette = isSilhouette;
  }

  inline const OpenMesh::Vec3f& source() const {
    return mp1;
  }
  inline const OpenMesh::Vec3f& target() const {
    return mp2;
  }
} ContourEdge;

bool isVisible(OpenMesh::Vec3f point, float margin);
void writeImage(Mesh &mesh, int width, int height, std::string filename, OpenMesh::Vec3f camPos, const std::list<ContourEdge>& contourEdges);

#endif
