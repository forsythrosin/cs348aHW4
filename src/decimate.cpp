#include "decimate.h"
#include <iostream>
#include <set>
#include <float.h>
using namespace OpenMesh;
using namespace std;

const double EPSILON = 0.00001;

VPropHandleT<Quadricd> vquadric;
VPropHandleT<float> vprio;
VPropHandleT<Mesh::HalfedgeHandle> vtarget;

void initDecimation(Mesh & mesh);
bool is_collapse_legal(Mesh &mesh, Mesh::HalfedgeHandle _hh);
float priority(Mesh &mesh, Mesh::HalfedgeHandle _heh);
void decimate(Mesh &mesh, unsigned int _n_vertices);
void enqueue_vertex(Mesh &mesh, Mesh::VertexHandle vh); 


// access quadric of vertex _vh
Quadricd& quadric(Mesh& mesh, Mesh::VertexHandle _vh) {
	return mesh.property(vquadric, _vh);
}

// access priority of vertex _vh
float& priority(Mesh& mesh, Mesh::VertexHandle _vh) {
	return mesh.property(vprio, _vh);
}

// access target halfedge of vertex _vh
Mesh::HalfedgeHandle& target(Mesh& mesh, Mesh::VertexHandle _vh) {
	return mesh.property(vtarget, _vh);
}


// NOTE:  We're making a global pointer to the mesh object here for notational convenience.
//    This will NOT work if you make your code multithreaded and is horrible programming practice.
//    But, it works and we're not asking you to implement this part.

Mesh *meshPtr;

// compare functor for priority queue
struct VertexCmp {
	bool operator()(Mesh::VertexHandle _v0, Mesh::VertexHandle _v1) const {
		Mesh &mesh = *meshPtr;
		// std::set needs UNIQUE keys -> handle equal priorities
		return ((priority(mesh,_v0) == priority(mesh,_v1)) ?
				(_v0.idx() < _v1.idx()) : (priority(mesh,_v0) < priority(mesh,_v1)));
	}
};

typedef std::set<Mesh::VertexHandle, VertexCmp> PQueue;
PQueue queue;

void simplify(Mesh &mesh, float percentage) {
	meshPtr = &mesh; // NEVER EVER DO THIS IN REAL LIFE

	// add required properties
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();
	mesh.request_face_normals();
	mesh.add_property(vquadric);
	mesh.add_property(vprio);
	mesh.add_property(vtarget);

	// compute normals & quadrics
	initDecimation(mesh);

	// decimate
	decimate(mesh, (int) (percentage * mesh.n_vertices()));
        std::cout << "Simplifying to #vertices: " << mesh.n_vertices() << std::endl;
}

void computeQuadric (Mesh &mesh, Mesh::VertexHandle vi) {
  priority(mesh, vi) = -1.0;
  quadric(mesh, vi).clear();

  // INSERT CODE HERE FOR PART 1
  // calc vertex quadrics from incident triangles
  Mesh::Point n;
  double a, b, c, d, length;
  Mesh::VertexFaceIter vf_it;      // To iterate through incident faces
  for(vf_it = mesh.vf_iter(vi); vf_it; ++vf_it) {
    Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle());
    n = mesh.normal(vf_it);
    a = n[0];
    b = n[1];
    c = n[2];
    length = sqrt(a*a + b*b + c*c);
    a /= length;
    b /= length;
    c /= length;
    assert (a*a + b*b + c*c - 1.0 < EPSILON);
    d = -1;
    quadric(mesh, vi) += Quadricd(a, b, c, d);
  }
}

void initDecimation(Mesh &mesh) {
	// compute face normals
	mesh.update_face_normals();

	Mesh::VertexIter v_it, v_end = mesh.vertices_end();
	Mesh::Scalar sum;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
          priority(mesh, v_it) = -1.0;
          sum = 0;                            // Reset for each iteration
          computeQuadric (mesh, v_it.handle());
	}
    std::cout << "Finished init" << std::endl;
}

bool is_collapse_legal(Mesh &mesh, Mesh::HalfedgeHandle _hh) {
	// collect vertices
	Mesh::VertexHandle v0, v1;
	v0 = mesh.from_vertex_handle(_hh);
	v1 = mesh.to_vertex_handle(_hh);

	// collect faces
	Mesh::FaceHandle fl = mesh.face_handle(_hh);
	Mesh::FaceHandle fr = mesh.face_handle(
			mesh.opposite_halfedge_handle(_hh));

	// backup point positions
	Mesh::Point p0 = mesh.point(v0);
	Mesh::Point p1 = mesh.point(v1);

	// topological test
	if (!mesh.is_collapse_ok(_hh))
		return false;

	// test boundary stuff
	if (mesh.is_boundary(v0) && !mesh.is_boundary(v1))
		return false;

    // Test for normal flipping
    for(Mesh::VertexFaceIter vf_i = mesh.vf_iter(v0); vf_i; ++vf_i){
        if(fl == mesh.face_handle(vf_i) || fr == mesh.face_handle(vf_i))
            continue;
        Mesh::Normal orig = mesh.normal(vf_i);
        
        Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_i);
        Vec3f p1 = mesh.point(fv_it);
        Vec3f p2 = mesh.point(++fv_it);
        Vec3f p3 = mesh.point(++fv_it);
        Mesh::Normal res;

        if(p1 == mesh.point(v0))
            p1 = mesh.point(v1);
        else if(p2 == mesh.point(v0))
            p2 = mesh.point(v0);
        else
            p3 = mesh.point(v0);

        res = (p2-p1) % (p3-p1);
        if((orig | res) / (orig.length()*res.length()) < 1/sqrt(2.0))
            return false;                           // Flipped triangle
    }

	// collapse passed all tests -> ok
	return true;
}

float priority(Mesh &mesh, Mesh::HalfedgeHandle _heh) {
  // INSERT CODE HERE FOR PART 2---------------------------------------------------------------------------------
  // return priority: the smaller the better
  // use quadrics to estimate approximation error
  // -------------------------------------------------------------------------------------------------------------
  Mesh::VertexHandle vi_handle = mesh.from_vertex_handle(_heh);
  Mesh::VertexHandle vj_handle = mesh.to_vertex_handle(_heh);
  Mesh::Point vj = mesh.point(vj_handle);
  Quadricd Qi = quadric(mesh, vi_handle);
  Quadricd Qj = quadric(mesh, vj_handle);
  Quadricd Q = Qi;
  Q += Qj;

  return Q(OpenMesh::VectorT<double, 3>(vj[0], vj[1], vj[2]));
}

void enqueue_vertex(Mesh &mesh, Mesh::VertexHandle _vh) {
	float prio, min_prio(FLT_MAX);
	Mesh::HalfedgeHandle min_hh;

	// find best out-going halfedge
	for (Mesh::VOHIter vh_it(mesh, _vh); vh_it; ++vh_it) {
		if (is_collapse_legal(mesh,vh_it)) {
			prio = priority(mesh, vh_it);
			if (prio != -1.0 && prio < min_prio) {
				min_prio = prio;
				min_hh = vh_it.handle();
			}
		}
	}

	// update queue
	if (priority(mesh, _vh) != -1.0) {
		queue.erase(_vh);
		priority(mesh, _vh) = -1.0;
	}

	if (min_hh.is_valid()) {
		priority(mesh, _vh) = min_prio;
		target(mesh, _vh) = min_hh;
		queue.insert(_vh);
	}
}

void decimate(Mesh &mesh, unsigned int _n_vertices) {
	unsigned int nv(mesh.n_vertices());
        std::cout << "Got to decimate" << std::endl;

	Mesh::HalfedgeHandle hh;
	Mesh::VertexHandle to, from;
	Mesh::VVIter vv_it;

	std::vector < Mesh::VertexHandle > one_ring;
	std::vector<Mesh::VertexHandle>::iterator or_it, or_end;

	// build priority queue
	Mesh::VertexIter v_it = mesh.vertices_begin(), v_end =
			mesh.vertices_end();

	queue.clear();
	for (; v_it != v_end; ++v_it)
          enqueue_vertex(mesh, v_it.handle());


	// INSERT CODE HERE FOR PART 3-------
	// DECIMATE THE MESH
	// ----------------------------------
        while (nv >= _n_vertices) {
          if (queue.empty()) {
            std::cout << "decimate: ran out of vertices" << std::endl;
            break;
          }
          // RIP OUT THE FIRST ELEMENT OF THE PRIORITY QUEUE
          PQueue::iterator pit = queue.begin();
          Mesh::VertexHandle vi_handle = *pit;

          // COLLAPSE THE EDGE
          Mesh::HalfedgeHandle vhe_handle = target(mesh, vi_handle);
          if (is_collapse_legal(mesh, vhe_handle)) {
            Mesh::VertexHandle vj_handle = mesh.to_vertex_handle(vhe_handle);
            //gather affected vertices and update their quadrics.
            set<Mesh::VertexHandle> affectedVerts;
            for (Mesh::VertexVertexIter vv_it = mesh.vv_iter(vi_handle); vv_it; ++vv_it) {
              affectedVerts.insert(vv_it.handle());
            }
            Quadricd Qi = quadric(mesh, vi_handle);
            quadric(mesh, vj_handle) += Qi;
            Qi.renormalize();

            //collapse the half-edge
            mesh.collapse(vhe_handle);
            --nv;

            //update neighboring vertices and their quadrics
            for (set<Mesh::VertexHandle>::iterator vv_it = affectedVerts.begin(); vv_it != affectedVerts.end(); ++vv_it) {
              //computeQuadric(mesh, *vv_it);
              //quadric(mesh, *vv_it) += Qi;
              //quadric(mesh, *vv_it).renormalize();
            }
            for (set<Mesh::VertexHandle>::iterator vv_it = affectedVerts.begin(); vv_it != affectedVerts.end(); ++vv_it) {
              enqueue_vertex(mesh, *vv_it);
            }
          }
          //--nv;
          queue.erase(pit);
        }

	// clean up after decimation
	queue.clear();

	// now, delete the items marked to be deleted
	mesh.garbage_collection();
        std::cout << "Out of decimate" << std::endl;
}

