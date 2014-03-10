#include "image_generation.h"
#include "mesh_features.h"
#include <GL/glut.h>
#include <fstream>
#include <map>
#include <set>
#include <list>
#include <sstream>
using namespace OpenMesh;
using namespace std;

struct VertexCompare : public binary_function <Vec3f, Vec3f, bool> {
  bool operator() (const Vec3f& v1, const Vec3f &v2) const {
    if (v1[0] < v2[0]) return true;
    else if (v1[0] > v2[0]) return false;

    if (v1[1] < v2[1]) return true;
    else if (v1[1] > v2[1]) return false;

    if (v1[2] < v2[2]) return true;
    else if (v1[2] > v2[2]) return false;

    return false;
  }
};

//Multimap: Vertex -> all of its neighbors
typedef map<Vec3f, const Vec3f*, VertexCompare> VertexMap;
typedef map<const Vec3f*, list<const Vec3f*> > EdgeMap;
struct ContourGraph {
  VertexMap vertices; //to find unique copy of a vertex
  EdgeMap edges; //vertex -> list of neighbors

  void insertVertex(const Vec3f& v) {
    if (vertices.find(v) == vertices.end()) {
      const Vec3f* newVertex = new Vec3f(v);
      vertices[v] = newVertex;
      edges[newVertex].clear();
    }
  }

  const Vec3f* getVertex(const Vec3f& v) const {
    if (vertices.find(v) != vertices.end()) {
      return vertices.find(v)->second;
    }
    return NULL;
  }
};

// Split vertex V into 2 vertices. The new vertex is
// connected to one of V's neigbhors, which is cut off
// from V.
void splitVertex(ContourGraph &graph, const Vec3f* v) {
  //Add slightly perturbed vertex
  Vec3f newVert = *v;
  while (graph.getVertex(newVert) != NULL) {
    newVert[0] += EPSILON;
    newVert[1] += EPSILON;
    newVert[2] += EPSILON;
  }
  graph.insertVertex(newVert);
  const Vec3f* newVertPtr = graph.getVertex(newVert);

  //Remove edge from neighbor -> *bvit, add edge from neighbor -> newVert
  const Vec3f *neighbor = graph.edges[v].front();
  graph.edges[v].remove(neighbor);
  graph.edges[neighbor].remove(v);
  graph.edges[newVertPtr].push_back(neighbor);
  graph.edges[neighbor].push_back(newVertPtr);
}

//Build graph of contour edges
void buildGraph(ContourGraph &graph, const list<ContourEdge>& contourEdges) {
  graph.vertices.clear();
  graph.edges.clear();

  //Collect all vertices and their neighbors
  for (list<ContourEdge>::const_iterator cit = contourEdges.begin();
       cit != contourEdges.end();
       ++cit) {
    const ContourEdge &c = *cit;
    const Vec3f &v1 = c.source();
    const Vec3f &v2 = c.target();

    if (v1 == v2) continue; //no self loops allowed
    graph.insertVertex(v1);
    graph.insertVertex(v2);
    graph.edges[graph.getVertex(v1)].push_back(graph.getVertex(v2));
    graph.edges[graph.getVertex(v2)].push_back(graph.getVertex(v1));
  }

  //Remove duplicate vertices from neighbor lists, and find vertices with 3+ neighbors.
  set<const Vec3f*> breakVerts;
  for (EdgeMap::iterator eit = graph.edges.begin(); eit != graph.edges.end(); ++eit) {
    eit->second.unique();
    if (eit->second.size() > 2) {
      breakVerts.insert(eit->first);
    }
  }

  //Some vertices (silhouette edges?) have 3 or more neighbors.
  //Break those apart, because our algorithm to find chains assumes no vertex
  //has 3 or more neighbors.
  for (set<const Vec3f*>::const_iterator bvit = breakVerts.begin(); bvit != breakVerts.end(); ++bvit) {
    const Vec3f* v = *bvit;
    list<const Vec3f *> &vnbrs = graph.edges[v];
    while (vnbrs.size() > 2) {
      splitVertex(graph, v);
    }
  }
}

//List of chains of vertices
typedef list< list <const Vec3f*> > ChainList;

//Explore a chain starting at currVertex, and traveling away from lastVertex.
void exploreChain(ContourGraph &graph, const Vec3f *lastVertex,
                  const Vec3f *currVertex, list <const Vec3f*>& chain) {
  assert (graph.vertices.find(*lastVertex) != graph.vertices.end());
  assert (graph.vertices.find(*currVertex) != graph.vertices.end());
  while (currVertex) {
    chain.push_back (currVertex);

    assert(graph.edges.find(currVertex) != graph.edges.end());
    //if (graph.edges.find(currVertex) == graph.edges.end())
    //  return;
    
    list <const Vec3f *>* currNeighbors = &graph.edges[currVertex];
    assert(currNeighbors->size() <= 2);
    currNeighbors->remove (lastVertex);
    assert(currNeighbors->size() <= 1);
    lastVertex = currVertex;
    if (!currNeighbors->empty()) {
      currVertex = currNeighbors->front();
      currNeighbors->pop_front();
    }
    else {
      currVertex = NULL; //reached the end
    }
    graph.vertices.erase(*lastVertex);
    graph.edges.erase(lastVertex);
  }
}

//Eat ContourGraph and build list of chains
void buildChains(ContourGraph &graph, ChainList &chains) {
  //Collect vertices with one neighbor first-- easy starting points.
  set<const Vec3f*> endPoints;
  set<const Vec3f*> isolatedPoints; //no connecting edge
  for (VertexMap::const_iterator vit = graph.vertices.begin(); vit != graph.vertices.end(); ++vit) {
    const Vec3f *v = vit->second;
    if (graph.edges[v].size() == 1) {
      endPoints.insert(v);
    }
  }
  
  //Chomp chains that end in a single vertex
  while (!endPoints.empty()) {
    const Vec3f* v = *endPoints.begin();
    list <const Vec3f*> chain;
    assert (graph.edges.find(v) != graph.edges.end());
    list <const Vec3f *>& vnbrs = graph.edges[v];
    assert (vnbrs.size() == 1);
    exploreChain(graph, v, vnbrs.front(), chain);
    chain.push_front(v);
    chains.push_back(chain);
    graph.vertices.erase(*v);
    graph.edges.erase(v);

    //remove all points on the chain
    for (list<const Vec3f*>::const_iterator cit = chain.begin(); cit != chain.end(); ++cit) {
      endPoints.erase(*cit);
    }
  }

  //Chomp what's left-- these are chains that loop back on themselves.
  //Cut them and eat them just like the normal chains
  cout << graph.vertices.size() << " vertices belong to loops" << endl;
  while (!graph.vertices.empty()) {
    const Vec3f *v = graph.vertices.begin()->second;
    assert (graph.edges.find(v) != graph.edges.end());
    assert (graph.edges[v].size() == 2);
    splitVertex(graph, v); //cut the loop

    //now chomp it just like the well-behaved chains
    list <const Vec3f*> chain;
    assert (graph.edges.find(v) != graph.edges.end());
    list <const Vec3f *>& vnbrs = graph.edges[v];
    assert (vnbrs.size() == 1);
    exploreChain(graph, v, vnbrs.front(), chain);
    chain.push_front(v);
    chains.push_back(chain);
    graph.vertices.erase(*v);
    graph.edges.erase(v);
  }
}

Vec3f toImagePlane(Vec3f point) {
	GLdouble point3DX = point[0], point3DY = point[1], point3DZ = point[2];

	GLdouble modelMatrix[16], projMatrix[16];
	GLint viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	GLdouble point2DX, point2DY, point2DZ;
	gluProject(point3DX, point3DY, point3DZ, modelMatrix, projMatrix, viewport, &point2DX, &point2DY, &point2DZ);
	
	return Vec3f(point2DX,point2DY,point2DZ);
}

// Adapted from
// http://stackoverflow.com/questions/1311869/opengl-how-to-determine-if-a-3d-rendered-point-is-occluded-by-other-3d-rende
bool isVisible(Vec3f point) {
	Vec3f projected = toImagePlane(point);

	GLfloat bufDepth = 0.0;
	glReadPixels(static_cast<GLint>( projected[0] ), static_cast<GLint>( projected[1] ),
		     1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &bufDepth);

	GLdouble EPSILON = 0.01;
	return (bufDepth - projected[2]) > -EPSILON; // check sign!
}

void writeImage(Mesh &mesh, int width, int height, string filename, Vec3f camPos,
                const list<ContourEdge>& contourEdges) {
  cout << "building graph" << endl;
  ContourGraph graph;
  buildGraph(graph, contourEdges);

  cout << "building chains" << endl;
  ChainList chains;
  buildChains(graph, chains);
  cout << "done building stuff" << endl;

  cout << "writing image: " << filename << endl;
  ofstream outfile(filename.c_str());
  outfile << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
  outfile << "<svg width=\"5in\" height=\"5in\" viewBox=\"0 0 " << width << ' ' << height << "\">\n";
  outfile << "<g stroke=\"black\" fill=\"black\">\n";

  // Sample code for generating image of the entire triangle mesh:
  /*for (Mesh::ConstEdgeIter it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
    Mesh::HalfedgeHandle h0 = mesh.halfedge_handle(it,0);
    Mesh::HalfedgeHandle h1 = mesh.halfedge_handle(it,1);
    Vec3f source(mesh.point(mesh.from_vertex_handle(h0)));
    Vec3f target(mesh.point(mesh.from_vertex_handle(h1)));
	
    if (!isVisible(source) || !isVisible(target)) continue;
		
    Vec3f p1 = toImagePlane(source);
    Vec3f p2 = toImagePlane(target);
    outfile << "<line ";
    outfile << "x1=\"" << p1[0] << "\" ";
    outfile << "y1=\"" << height-p1[1] << "\" ";
    outfile << "x2=\"" << p2[0] << "\" ";
    outfile << "y2=\"" << height-p2[1] << "\" stroke-width=\"1\" />\n";
    }*/
  
  // WRITE CODE HERE TO GENERATE A .SVG OF THE MESH --------------------------------------------------------------
  
  // Dumb code to render all contour edges
  /*
    for (list<ContourEdge>::const_iterator cit = contourEdges.begin(); cit != contourEdges.end(); ++cit) {
    const ContourEdge& c = *cit;
    
    if (!isVisible(c.source()) || !isVisible(c.target())) continue;
    
    Vec3f p1 = toImagePlane(c.source());
    Vec3f p2 = toImagePlane(c.target());
    outfile << "<line ";
    outfile << "x1=\"" << p1[0] << "\" ";
    outfile << "y1=\"" << height-p1[1] << "\" ";
    outfile << "x2=\"" << p2[0] << "\" ";
    outfile << "y2=\"" << height-p2[1] << "\" stroke-width=\"1\" />\n";
    }
  */
  
  // Render chains of edges
  //static const char* COLORS[] = {"red", "black", "blue", "green", "orange", "magenta", "cyan"};
  //static const int N_COLORS = sizeof(COLORS) / sizeof(const char*);
  for (ChainList::const_iterator lcit = chains.begin(); lcit != chains.end(); ++lcit) {
    const list<const Vec3f*>& chain = *lcit;
    const Vec3f* lastVertex = NULL;
    //const char* color = COLORS[rand() % N_COLORS];
    stringstream sscolor;
    sscolor << "rgb(" << rand() % 256 << "," << rand() % 256 << "," << rand() % 256 << ")";
    for (list<const Vec3f*>::const_iterator cit = chain.begin(); cit != chain.end(); ++cit) {
      const Vec3f* currVertex = *cit;
      if (lastVertex) {
        //TODO: why doesn't visibility test work?
        //if (isVisible(*lastVertex) && isVisible(*currVertex)) {
          Vec3f p1 = toImagePlane(*lastVertex);
          Vec3f p2 = toImagePlane(*currVertex);
          outfile << "<line ";
          outfile << "x1=\"" << p1[0] << "\" ";
          outfile << "y1=\"" << height-p1[1] << "\" ";
          outfile << "x2=\"" << p2[0] << "\" ";
          outfile << "y2=\"" << height-p2[1] << "\" stroke=\"" << sscolor.str() << "\" stroke-width=\"1\" />\n";
          //}
      }
      lastVertex = currVertex;
    }
  }
  
  // -------------------------------------------------------------------------------------------------------------
  
  outfile << "</g>\n";
  outfile << "</svg>\n";
  cout << "done" << endl;
}
