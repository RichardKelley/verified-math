#include "verified_math/vec3.h"
#include "gtest/gtest.h"
#include "checkpp/checkpp.h"

/*
TEST(TestVec3, TestIdentityVecCreation) {
  verified_math::Vec3<double> v{};
  ASSERT_FLOAT_EQ(v.x1, 0.0);
  ASSERT_FLOAT_EQ(v.x2, 0.0);
  ASSERT_FLOAT_EQ(v.x3, 0.0);
}
*/

TEST(TestVec3, TestAdditionCommutativity) {
  auto commutative_property = [](float x1, float x2, float x3,
				 float y1, float y2, float y3) {
    auto v1 = verified_math::Vec3<float>{x1, x2, x3};
    auto v2 = verified_math::Vec3<float>{y1, y2, y3};

    verified_math::Vec3<float> sum1 = v1 + v2;
    verified_math::Vec3<float> sum2 = v2 + v1;

    return (abs(sum1.x1 - sum2.x1) < 0.001 &&
	    abs(sum1.x2 - sum2.x2) < 0.001 &&
	    abs(sum1.x3 - sum2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<float, float, float, float, float, float> {
	commutative_property
	  }, 10000
      )
    );
}

TEST(TestVec3, TestAdditionAssociativity) {
  auto associative_property = [](float x1, float x2, float x3,
				 float y1, float y2, float y3,
				 float z1, float z2, float z3) {
    auto v1 = verified_math::Vec3<float>{x1, x2, x3};
    auto v2 = verified_math::Vec3<float>{y1, y2, y3};
    auto v3 = verified_math::Vec3<float>{z1, z2, z3};

    auto sum1 = (v1 + v2) + v3;
    auto sum2 = v1 + (v2 + v3);

    return (abs(sum1.x1 - sum2.x1) < 0.001 &&
	    abs(sum1.x2 - sum2.x2) < 0.001 &&
	    abs(sum1.x3 - sum2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(
			     checkpp::Property<float, float, float,
			     float, float, float,
			     float, float, float> {
			       associative_property
				 }, 10000
			     )
	      );
}
