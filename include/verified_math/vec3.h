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

  template<typename Scalar>
  Vec3<Scalar> operator+(Vec3<Scalar> x, Vec3<Scalar> y) {
    return Vec3<Scalar>(x.x1 + y.x1, x.x2 + y.x2, x.x3 + y.x3);
  }
}

#endif // VEC3_H
