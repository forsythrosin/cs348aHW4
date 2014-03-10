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

  ContourEdge (const OpenMesh::Vec3f& p1, const OpenMesh::Vec3f& p2) {
    mp1 = p1;
    mp2 = p2;
  }

  inline const OpenMesh::Vec3f& source() const {
    return mp1;
  }
  inline const OpenMesh::Vec3f& target() const {
    return mp2;
  }
} ContourEdge;

/*
//Less-than for contour edges
typedef struct ContourLess : public binary_function <ContourEdge, ContourEdge, bool> {
  bool operator() (const ContourEdge &c1, const ContourEdge &c2) const {
    if (c1[0] < c2[0]) {
      return true;
    }
    else if (c1[0] > c2[0]) {
      return false;
    }

    if (c1[1] < c2[1]) {
      return true;
    }
    else if (c1[1] > c2[1]) {
      return false;
    }

    if (c1[2] < c2[2]) {
      return true;
    }
    else if (c1[2] > c2[2]) {
      return false;
    }

    return false; // equal
  }
}
*/

bool isVisible(OpenMesh::Vec3f point);
void writeImage(Mesh &mesh, int width, int height, std::string filename, OpenMesh::Vec3f camPos, const std::list<ContourEdge>& contourEdges);

#endif
