#ifndef VEC3_H
#define VEC3_H

namespace verified_math {

  template <typename Scalar>
  class Vec3 {
  public:
    Scalar x1 {0.0};
    Scalar x2 {0.0};
    Scalar x3 {0.0};

    Vec3<Scalar>(Scalar _x1, Scalar _x2, Scalar _x3)
      : x1{_x1}, x2{_x2}, x3{_x3} { }
  };

  /*
    Basic vector algebra.
  */
  template<typename Scalar>
  Vec3<Scalar> operator+(Vec3<Scalar> x, Vec3<Scalar> y) {
    return Vec3<Scalar>(x.x1 + y.x1, x.x2 + y.x2, x.x3 + y.x3);
  }

  template<typename Scalar>
  Vec3<Scalar> operator-(Vec3<Scalar> x, Vec3<Scalar> y) {
    return Vec3<Scalar>(x.x1 - y.x1, x.x2 - y.x2, x.x3 - y.x3);
  }

  template<typename Scalar>
  Vec3<Scalar> operator*(Scalar c, Vec3<Scalar> x) {
    return Vec3<Scalar>(c * x.x1, c * x.x2, c * x.x3);
  }

  template<typename Scalar>
  Vec3<Scalar> operator*(Vec3<Scalar> x, Scalar c) {
    return Vec3<Scalar>(x.x1 * c, x.x2 * c, x.x3 * c);
  }

  /*
    Vector multiplications
  */
  // dot product
  template<typename Scalar>
  Scalar dot(Vec3<Scalar> x, Vec3<Scalar> y) {
    return (x.x1 * y.x1 + x.x2 * y.x2 + x.x3 * y.x3);
  }

  // cross product
  template<typename Scalar>
  Vec3<Scalar> cross(Vec3<Scalar> x, Vec3<Scalar> y) {
    return Vec3<Scalar>(x.x2 * y.x3 - x.x3 * y.x2,
			x.x3 * y.x1 - x.x1 * y.x3,
			x.x1 * y.x2 - x.x2 * y.x1);
  }


}

#endif // VEC3_H
