#ifndef VEC4_H
#define VEC4_H

namespace verified_math {

  template <typename Scalar>
  class Vec4 {
  public:
    Scalar x1 {0.0};
    Scalar x2 {0.0};
    Scalar x3 {0.0};
    Scalar x4 {0.0};

    Vec4<Scalar>(Scalar _x1, Scalar _x2, Scalar _x3, Scalar _x4)
      : x1{_x1}, x2{_x2}, x3{_x3}, x4{_x4} { }
  };

  /*
    Basic vector algebra.
  */
  template<typename Scalar>
  Vec4<Scalar> operator+(const Vec4<Scalar>& x, const Vec4<Scalar>& y) {
    return Vec4<Scalar>(x.x1 + y.x1, x.x2 + y.x2, 
			x.x3 + y.x3, x.x4 + y.x4);
  }

  template<typename Scalar>
  Vec4<Scalar> operator-(const Vec4<Scalar>& x, const Vec4<Scalar>& y) {
    return Vec4<Scalar>(x.x1 - y.x1, x.x2 - y.x2, 
			x.x3 - y.x3, x.x4 - y.x4);
  }

  template<typename Scalar>
  Vec4<Scalar> operator*(const Scalar c, const Vec4<Scalar>& x) {
    return Vec4<Scalar>(c * x.x1, c * x.x2, 
			c * x.x3, c * x.x4);
  }

  template<typename Scalar>
  Vec4<Scalar> operator*(const Vec4<Scalar>& x, Scalar c) {
    return Vec4<Scalar>(x.x1 * c, x.x2 * c, 
			x.x3 * c, x.x4 * c);
  }

  /*
    Vector multiplications
  */
  // dot product
  template<typename Scalar>
  Scalar dot(const Vec4<Scalar>& x, const Vec4<Scalar>& y) {
    return (x.x1 * y.x1 + x.x2 * y.x2 + 
	    x.x3 * y.x3 + x.x4 * y.x4);
  }

}

#endif // VEC3_H
