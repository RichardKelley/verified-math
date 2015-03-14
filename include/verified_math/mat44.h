#ifndef MAT44_H
#define MAT44_H

namespace verified_math {

  /*
    A 4x4 matrix.
   */
  template<typename Scalar>
  class Mat44 {
  public:
    Scalar x11{0}; Scalar x12{0}; Scalar x13{0}; Scalar x14{0};
    Scalar x21{0}; Scalar x22{0}; Scalar x23{0}; Scalar x24{0};
    Scalar x31{0}; Scalar x32{0}; Scalar x33{0}; Scalar x34{0};
    Scalar x41{0}; Scalar x42{0}; Scalar x43{0}; Scalar x44{0};

    Mat44<Scalar> (Scalar _x11, Scalar _x12, Scalar _x13, Scalar _x14,
		   Scalar _x21, Scalar _x22, Scalar _x23, Scalar _x24,
		   Scalar _x31, Scalar _x32, Scalar _x33, Scalar _x34,
		   Scalar _x41, Scalar _x42, Scalar _x43, Scalar _x44)
      : x11{_x11}, x12{_x12}, x13{_x13}, x14{_x14},
        x21{_x21}, x22{_x22}, x23{_x23}, x24{_x24},
        x31{_x31}, x32{_x32}, x33{_x33}, x34{_x34},
        x41{_x41}, x42{_x42}, x43{_x43}, x44{_x44} { }

	Scalar l2_norm() const {
	  return x11 * x11 + x12 * x12 + x13 * x13 + x14 * x14 +
	    x21 * x21 + x22 * x22 + x23 * x23 + x24 * x24 +
	    x31 * x31 + x32 * x32 + x33 * x33 + x34 * x34 +
	    x41 * x41 + x42 * x42 + x43 * x43 + x44 * x44;
	}
  };

  template<typename Scalar>
  Mat44<Scalar> operator+(const Mat44<Scalar>& m1, const Mat44<Scalar>& m2) {
    return Mat44<Scalar>{
      m1.x11 + m2.x11, m1.x12 + m2.x12, m1.x13 + m2.x13, m1.x14 + m2.x14,
	m1.x21 + m2.x21, m1.x22 + m2.x22, m1.x23 + m2.x23, m1.x24 + m2.x24,
	m1.x31 + m2.x31, m1.x32 + m2.x32, m1.x33 + m2.x33, m1.x34 + m2.x34,
	m1.x41 + m2.x41, m1.x42 + m2.x42, m1.x43 + m2.x43, m1.x44 + m2.x44
    };
  }

  template<typename Scalar>
  Mat44<Scalar> operator-(const Mat44<Scalar>& m1, const Mat44<Scalar>& m2) {
    return Mat44<Scalar> {
      m1.x11 - m2.x11, m1.x12 - m2.x12, m1.x13 - m2.x13, m1.x14- m2.x14,
	m1.x21 - m2.x21, m1.x22 - m2.x22, m1.x23 - m2.x23, m1.x24 - m2.x24,
	m1.x31 - m2.x31, m1.x32 - m2.x32, m1.x33 - m2.x33, m1.x34 - m2.x34,
	m1.x41 - m2.x41, m1.x42 - m2.x42, m1.x43 - m2.x43, m1.x44 - m2.x44
    };
  }

  template<typename Scalar>
  Mat44<Scalar> operator*(const Scalar c, const Mat44<Scalar> m) {
    return Mat44<Scalar>{
      c * m.x11, c * m.x12, c * m.x13, c * m.x14,
	c * m.x21, c * m.x22, c * m.x23, c * m.x24,
	c * m.x31, c * m.x32, c * m.x33, c * m.x34,
	c * m.x41, c * m.x42, c * m.x43, c * m.x44
     };
  }

  template<typename Scalar>
  Mat44<Scalar> operator*(const Mat44<Scalar> m, const Scalar c) {
    return Mat44<Scalar>{
      m.x11 * c, m.x12 * c, m.x13 * c, m.x14 * c,
	m.x21 * c, m.x22 * c, m.x23 * c, m.x24 * c,
	m.x31 * c, m.x32 * c, m.x33 * c, m.x34 * c,
	m.x41 * c, m.x42 * c, m.x43 * c, m.x44 * c
    };
  }

  template<typename Scalar>
  Vec4<Scalar> operator*(const Mat44<Scalar>& m, const Vec4<Scalar>& v) { 
    return Vec3<Scalar> {
      m.x11 * v.x1 + m.x12 * v.x2 + m.x13 * v.x3 + m.x14 * v.x4,
	m.x21 * v.x1 + m.x22 * v.x2 + m.x23 * v.x3 + m.x24 * v.x4,
	m.x31 * v.x1 + m.x32 * v.x2 + m.x33 * v.x3 + m.x34 * v.x4,
	m.x41 * v.x1 + m.x42 * v.x2 + m.x43 * v.x3 + m.x44 * v.x4
    };
  }

    template<typename Scalar>
  Mat44<Scalar> operator*(const Mat44<Scalar>& m1, const Mat44<Scalar>& m2) {
    return Mat44<Scalar> {m1.x11 * m2.x11 + m1.x12 * m2.x21 + m1.x13 * m2.x31 + m1.x14 * m2.x41,
			  m1.x11 * m2.x12 + m1.x12 * m2.x22 + m1.x13 * m2.x32 + m1.x14 * m2.x42,
			  m1.x11 * m2.x13 + m1.x12 * m2.x23 + m1.x13 * m2.x33 + m1.x14 * m2.x43,
                          m1.x11 * m2.x14 + m1.x12 * m2.x24 + m1.x13 * m2.x34 + m1.x14 * m2.x44,
			  
			  m1.x21 * m2.x11 + m1.x22 * m2.x21 + m1.x23 * m2.x31 + m1.x24 * m2.x41,
			  m1.x21 * m2.x12 + m1.x22 * m2.x22 + m1.x23 * m2.x32 + m1.x24 * m2.x42,
			  m1.x21 * m2.x13 + m1.x22 * m2.x23 + m1.x23 * m2.x33 + m1.x24 * m2.x43,
			  m1.x21 * m2.x14 + m1.x22 * m2.x24 + m1.x23 * m2.x34 + m1.x24 * m2.x44,
			  
			  m1.x31 * m2.x11 + m1.x32 * m2.x21 + m1.x33 * m2.x31 + m1.x34 * m2.x41,
			  m1.x31 * m2.x12 + m1.x32 * m2.x22 + m1.x33 * m2.x32 + m1.x34 * m2.x42,
                          m1.x31 * m2.x13 + m1.x32 * m2.x23 + m1.x33 * m2.x33 + m1.x34 * m2.x43,
	                  m1.x31 * m2.x14 + m1.x32 * m2.x24 + m1.x33 * m2.x34 + m1.x34 * m2.x44,

			  m1.x41 * m2.x11 + m1.x42 * m2.x21 + m1.x43 * m2.x31 + m1.x44 * m2.x41,
			  m1.x41 * m2.x12 + m1.x42 * m2.x22 + m1.x43 * m2.x32 + m1.x44 * m2.x42,
                          m1.x41 * m2.x13 + m1.x42 * m2.x23 + m1.x43 * m2.x33 + m1.x44 * m2.x43,
	                  m1.x41 * m2.x14 + m1.x42 * m2.x24 + m1.x43 * m2.x34 + m1.x44 * m2.x44
    };
  }

  template<typename Scalar>
  Mat44<Scalar> transpose(const Mat44<Scalar>& m) {
    return Mat44<Scalar> {
      m.x11, m.x21, m.x31, m.x41,
	m.x12, m.x22, m.x32, m.x42,
	m.x13, m.x23, m.x33, m.x43,
	m.x14, m.x24, m.x34, m.x44
    };
  }


  template<typename Scalar>
  Scalar det(const Mat44<Scalar>& m) {
    // TODO fixme
    return 0.0;
  }

  template<typename Scalar>
  Mat44<Scalar> inverse(const Mat44<Scalar>& m) {
    // TODO fixme
    return (1.0 / det(m)) * Mat44<Scalar>{
        0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0
    };
  }

  template<typename Scalar>
  Scalar condition_number(const Mat44<Scalar>& m) {
    auto norm = m.l2_norm();
    auto inv = inverse(m);
    auto norm_inv = inv.l2_norm();
    return norm * norm_inv;
  }

  template<typename Scalar>
  Scalar trace(const Mat44<Scalar>& m) {
    return m.x11 + m.x22 + m.x33 + m.x44;
  }

}


#endif // MAT44_H
