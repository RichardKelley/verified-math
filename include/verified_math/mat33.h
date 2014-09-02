#ifndef MAT33_H
#define MAT33_H

#include "verified_math/vec3.h"

namespace verified_math {
  
  /*
    A 3x3 matrix.
   */
  template<typename Scalar>
  class Mat33 {
  public:
    Scalar x11 {0.0}; Scalar x12 {0.0}; Scalar x13 {0.0};
    Scalar x21 {0.0}; Scalar x22 {0.0}; Scalar x23 {0.0};
    Scalar x31 {0.0}; Scalar x32 {0.0}; Scalar x33 {0.0};

    Mat33<Scalar> (Scalar _x11, Scalar _x12, Scalar _x13,
		   Scalar _x21, Scalar _x22, Scalar _x23,
		   Scalar _x31, Scalar _x32, Scalar _x33) 
    : x11{_x11}, x12{_x12}, x13{_x13},
      x21{_x21}, x22{_x22}, x23{_x32},
      x31{_x31}, x32{_x32}, x33{_x33} { }
  };


  template<typename Scalar>
  Mat33<Scalar> operator+(const Mat33<Scalar>& m1, const Mat33<Scalar>& m2) {
    return Mat33<Scalar>{
      m1.x11 + m2.x11, m1.x12 + m2.x12, m1.x13 + m2.x13,
	m1.x21 + m2.x21, m1.x22 + m2.x22, m1.x23 + m2.x23,
	m1.x31 + m2.x31, m1.x32 + m2.x32, m1.x33 + m2.x33
    };
  }

  template<typename Scalar>
  Mat33<Scalar> operator-(const Mat33<Scalar>& m1, const Mat33<Scalar>& m2) {
    return Mat33<Scalar> {
      m1.x11 - m2.x11, m1.x12 - m2.x12, m1.x13 - m2.x13,
	m1.x21 - m2.x21, m1.x22 - m2.x22, m1.x23 - m2.x23,
	m1.x31 - m2.x31, m1.x32 - m2.x32, m1.x33 - m2.x33
    };
  }

  template<typename Scalar>
  Mat33<Scalar> operator*(const Scalar c, const Mat33<Scalar> m) {
    return Mat33<Scalar>{
      c * m.x11, c * m.x12, c * m.x13,
	c * m.x21, c * m.x22, c * m.x23,
	c * m.x31, c * m.x32, c * m.x33
    };
  }

  template<typename Scalar>
  Mat33<Scalar> operator*(const Mat33<Scalar> m, const Scalar c) {
    return Mat33<Scalar>{
      m.x11 * c, m.x12 * c, m.x13 * c,
	m.x21 * c, m.x22 * c, m.x23 * c,
	m.x31 * c, m.x32 * c, m.x33 * c
    };
  }

  template<typename Scalar>
  Vec3<Scalar> operator*(const Mat33<Scalar>& m, const Vec3<Scalar>& v) { 
    return Vec3<Scalar> {
      m.x11 * v.x1 + m.x12 * v.x2 + m.x13 * v.x3,
	m.x21 * v.x1 + m.x22 * v.x2 + m.x23 * v.x3,
	m.x31 * v.x1 + m.x32 * v.x2 + m.x33 * v.x3
    };
  }

  template<typename Scalar>
  Mat33<Scalar> operator*(const Mat33<Scalar>& m1, const Mat33<Scalar>& m2) {
    return Mat33<Scalar> { m1.x11 * m2.x11 + m1.x12 * m2.x21 + m1.x13 * m2.x31,
			  m1.x11 * m2.x12 + m1.x12 * m2.x22 + m1.x13 * m2.x32,
			  m1.x11 * m2.x13 + m1.x12 * m2.x23 + m1.x13 * m2.x33,
			  
			  m1.x21 * m2.x11 + m1.x22 * m2.x21 + m1.x23 * m2.x31,
			  m1.x21 * m2.x12 + m1.x22 * m2.x22 + m1.x23 * m2.x32,
			  m1.x21 * m2.x13 + m1.x22 * m2.x23 + m1.x23 * m2.x33,
			  
			  m1.x31 * m2.x11 + m1.x32 * m2.x21 + m1.x33 * m2.x31,
			  m1.x31 * m2.x12 + m1.x32 * m2.x22 + m1.x33 * m2.x32,
			  m1.x31 * m2.x13 + m1.x32 * m2.x23 + m1.x33 * m2.x33
    };
  }

  template<typename Scalar>
  Mat33<Scalar> transpose(const Mat33<Scalar>& m) {
    return Mat33<Scalar> {
      m.x11, m.x21, m.x31,
	m.x12, m.x22, m.x32,
	m.x13, m.x23, m.x33
    };
  }

  template<typename Scalar>
  Scalar det(const Mat33<Scalar>& m) {
    return m.x11 * (m.x22 * m.x33 - m.x23 * m.x32) +
      m.x12 * (m.x23 * m.x31 - m.x21 * m.x33) +
      m.x13 * (m.x21 * m.x32 - m.x31 * m.x22);
  }

  template<typename Scalar>
  Mat33<Scalar> inverse(const Mat33<Scalar>& m) {
    return (1.0 / det(m)) * Mat33<Scalar>{
      (m.x22 * m.x33 - m.x23 * m.x32), -(m.x12 * m.x33 - m.x13 * m.x32), (m.x12 * m.x23 - m.x13 * m.x22),
	-(m.x21 * m.x33 - m.x23 * m.x31), (m.x11 * m.x33 - m.x13 * m.x31), -(m.x11 * m.x23 - m.x13 * m.x21),
	(m.x21 * m.x32 - m.x22 * m.x31), -(m.x11 * m.x32 - m.x12 * m.x31), (m.x11 * m.x22 - m.x12 * 21)
    };
  }

  template<typename Scalar>
  Scalar trace(const Mat33<Scalar>& m) {
    return m.x11 + m.x22 + m.x33;
  }

}

#endif // MAT33_H
