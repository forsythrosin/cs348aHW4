#define DECIMATE_HH

//== INCLUDES =================================================================

#include "mesh_definitions.h"

void simplify(Mesh &mesh, float percentage);

//== CLASS DEFINITION =========================================================

/** /class QuadricT

 Stores a quadric as a 4x4 symmetrix matrix. Used by the
 error quadric based mesh decimation algorithms.
 **/

template<class Scalar>
class QuadricT {
public:

	/// construct with upper triangle of symmetric 4x4 matrix
	QuadricT(Scalar _a, Scalar _b, Scalar _c, Scalar _d, Scalar _e, Scalar _f,
			Scalar _g, Scalar _h, Scalar _i, Scalar _j) :
			a(_a), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h), i(_i), j(_j) {
	}

	/// constructor from given plane equation: ax+by+cz+d=0
	QuadricT(Scalar _a = 0.0, Scalar _b = 0.0, Scalar _c = 0.0, Scalar _d = 0.0) :
			a(_a * _a), b(_a * _b), c(_a * _c), d(_a * _d), e(_b * _b), f(
					_b * _c), g(_b * _d), h(_c * _c), i(_c * _d), j(_d * _d) {
	}

	/// set all entries to zero
	void clear() {
		a = b = c = d = e = f = g = h = i = j = 0.0;
	}

	/// add quadrics
	QuadricT<Scalar>& operator+=(const QuadricT<Scalar>& _q) {
		a += _q.a;
		b += _q.b;
		c += _q.c;
		d += _q.d;
		e += _q.e;
		f += _q.f;
		g += _q.g;
		h += _q.h;
		i += _q.i;
		j += _q.j;
		return *this;
	}

	/// multiply by scalar
	QuadricT<Scalar>& operator*=(Scalar _s) {
		a *= _s;
		b *= _s;
		c *= _s;
		d *= _s;
		e *= _s;
		f *= _s;
		g *= _s;
		h *= _s;
		i *= _s;
		j *= _s;
		return *this;
	}

	/// evaluate quadric Q at vector v: v*Q*v
	template<typename T>
	Scalar operator()(const OpenMesh::VectorT<T, 3> _v) const {
		Scalar x(_v[0]), y(_v[1]), z(_v[2]);
		return a * x * x + 2.0 * b * x * y + 2.0 * c * x * z + 2.0 * d * x
				+ e * y * y + 2.0 * f * y * z + 2.0 * g * y + h * z * z
				+ 2.0 * i * z + j;
	}

        /// build matrix from this quadric
        Eigen::Matrix4d toMatrix() const {
          Eigen::Matrix4d mQ;
          mQ(0, 0) = a;
          mQ(0, 1) = b;
          mQ(0, 2) = c;
          mQ(0, 3) = d;
          mQ(1, 0) = b;
          mQ(1, 1) = e;
          mQ(1, 2) = f;
          mQ(1, 3) = g;
          mQ(2, 0) = c;
          mQ(2, 1) = f;
          mQ(2, 2) = h;
          mQ(2, 3) = i;
          mQ(3, 0) = d;
          mQ(3, 1) = g;
          mQ(3, 2) = i;
          mQ(3, 3) = j;
          return mQ;
        }

private:

	Scalar a, b, c, d, e, f, g, h, i, j;
};

/// Quadric using floats
typedef QuadricT<float> Quadricf;

/// Quadric using double
typedef QuadricT<double> Quadricd;
