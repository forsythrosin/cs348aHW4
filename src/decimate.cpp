#include <Eigen/Core>
#include <Eigen/LU>
#include "decimate.h"
#include <iostream>
#include <set>
#include <float.h>
using namespace OpenMesh;
using namespace Eigen;
using namespace std;

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

void initDecimation(Mesh &mesh) {
	// compute face normals
	mesh.update_face_normals();

	Mesh::VertexIter v_it, v_end = mesh.vertices_end();
	Mesh::Point n;
	Mesh::VertexFaceIter vf_it;      // To iterate through incident faces
	double a, b, c, d, length;
	Mesh::Scalar sum;

	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
          priority(mesh, v_it) = -1.0;
          quadric(mesh, v_it).clear();
          sum = 0;                            // Reset for each iteration
          
          // INSERT CODE HERE FOR PART 1
          // calc vertex quadrics from incident triangles
          // TODO: do we need a^2 + b^2 + c^2 = 1?
          Mesh::Point v = mesh.point(v_it.handle());
          for(vf_it = mesh.vf_iter(v_it.handle()); vf_it; ++vf_it) {
            Mesh::FaceVertexIter fv_it = mesh.fv_iter(vf_it.handle());
            Mesh::Point fp = mesh.point(fv_it); //a vertex on face
            n = mesh.normal(vf_it);
            a = n[0];
            b = n[1];
            c = n[2];
            d = -a*fp[0] - b*fp[1] - c*fp[2];
            quadric(mesh, v_it) += Quadricd(a, b, c, d);
          }
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
  Mesh::Point vi = mesh.point(vi_handle);
  Mesh::Point vj = mesh.point(vj_handle);
  Quadricd Qi = quadric(mesh, vi_handle);
  Quadricd Qj = quadric(mesh, vj_handle);
  Quadricd Q = Qi;
  Q += Qj;

  //solve for vbar: the vertex minimizing error on this edge
  Matrix4d mQ = Q.toMatrix();
  mQ(3, 0) = 0;
  mQ(3, 1) = 0;
  mQ(3, 2) = 0;
  mQ(3, 3) = 1;

  //TODO do the inverse?
  //TODO check if mQ is singular first?
  //Vector4d vbar = mQ.inverse() * Vector4d(0, 0, 0, 1);
  
  Vector4d vbar = Vector4d(vj[0], vj[1], vj[2], 1);

  return vbar.transpose() * Q.toMatrix() * vbar;
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
        while(nv >= _n_vertices) {
          if (queue.empty()) {
            std::cout << "decimate: ran out of vertices" << std::endl;
            break;
          }
          // RIP OUT THE FIRST ELEMENT OF THE PRIORITY QUEUE
          PQueue::iterator pit = queue.begin();
          Mesh::VertexHandle vi_handle = *pit;
          queue.erase(pit);

          // CRUSH THIS HALFEDGE LIKE A CHICKEN BONE
          Mesh::HalfedgeHandle vhe_handle = target(mesh, vi_handle);
          if (is_collapse_legal(mesh, vhe_handle)) {
            Mesh::VertexHandle vj_handle = mesh.to_vertex_handle(vhe_handle);
            //Remove target vertex from queue
            //PQueue::iterator qt_it = queue.find(mesh.to_vertex_handle(vhe_handle));
            //if (qt_it != queue.end())
            //  queue.erase(qt_it);
            //collapse the 'from' vertex vi to the 'to' vertex vj
            mesh.collapse(vhe_handle);
            --nv;

            //TODO: move the vertex to vbar or v2?
            //TODO: ?

            //TODO: put it back in or leave it?
            //enqueue_vertex(mesh, vi_handle);

            // SWEEP STORM OF DESTRUCTION OVER NEIGHBORING VERTICES
            for(Mesh::VertexVertexIter vv_it = mesh.vv_iter(vj_handle); vv_it; ++vv_it) {
              enqueue_vertex(mesh, vv_it.handle());
            }
          }
        }

	// clean up after decimation
	queue.clear();

	// now, delete the items marked to be deleted
	mesh.garbage_collection();
        std::cout << "Out of decimate" << std::endl;
}

