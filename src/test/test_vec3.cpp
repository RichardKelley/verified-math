#include "verified_math/vec3.h"
#include "gtest/gtest.h"
#include "checkpp/checkpp.h"

TEST(TestVec3, TestAdditionDoubleCommutativity) {
  auto commutative_property = [](double x1, double x2, double x3,
				 double y1, double y2, double y3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    auto v2 = verified_math::Vec3<double>{y1, y2, y3};

    verified_math::Vec3<double> sum1 = v1 + v2;
    verified_math::Vec3<double> sum2 = v2 + v1;

    return (abs(sum1.x1 - sum2.x1) < 0.001 &&
	    abs(sum1.x2 - sum2.x2) < 0.001 &&
	    abs(sum1.x3 - sum2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double> {
	commutative_property
	  }, 10000
      )
    );
}

TEST(TestVec3, TestAdditionDoubleAssociativity) {
  auto associative_property = [](double x1, double x2, double x3,
				 double y1, double y2, double y3,
				 double z1, double z2, double z3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    auto v2 = verified_math::Vec3<double>{y1, y2, y3};
    auto v3 = verified_math::Vec3<double>{z1, z2, z3};

    auto sum1 = (v1 + v2) + v3;
    auto sum2 = v1 + (v2 + v3);

    return (abs(sum1.x1 - sum2.x1) < 0.001 &&
	    abs(sum1.x2 - sum2.x2) < 0.001 &&
	    abs(sum1.x3 - sum2.x3) < 0.001);
  };

  EXPECT_FALSE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double> { associative_property }, 10000 ));
  
}

TEST(TestVec3, TestAdditionDoubleIdentity) {
  auto identity_addition_property = [](double x1, double x2, double x3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    auto zero = verified_math::Vec3<double>{0.0, 0.0, 0.0};

    auto sum1 = v1 + zero;

    return (abs(sum1.x1 - v1.x1) < 0.001 &&
	    abs(sum1.x2 - v1.x2) < 0.001 &&
	    abs(sum1.x3 - v1.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double> {
	identity_addition_property}, 10000)
    );
}

TEST(TestVec3, TestAdditionDoubleInverse) {
  auto inverse_addition_property = [](double x1, double x2, double x3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    auto neg_v1 = -1.0 * v1;
    auto sum1 = v1 + neg_v1;
    auto zero = verified_math::Vec3<double>{0.0f, 0.0f, 0.0f};

    return (abs(sum1.x1 - zero.x1) < 0.001 &&
	    abs(sum1.x2 - zero.x2) < 0.001 &&
	    abs(sum1.x3 - zero.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double> {
	inverse_addition_property} , 10000)
    );
      
}

TEST(TestVec3, TestScalarMultCommutative) {
  auto scalar_commutative_prop = [](double x1, double x2, double x3, double c) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    
    auto prod1 = c * v1;
    auto prod2 = v1 * c;

    return (abs(prod1.x1 - prod2.x1) < 0.001 &&
	    abs(prod1.x2 - prod2.x2) < 0.001 &&
	    abs(prod1.x3 - prod2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double> {
	scalar_commutative_prop
	  }, 10000)
    );
}

TEST(TestVec3, TestFieldMultCompatibility) {
  auto field_mult_compatibility = [](double a, double b, double x1, double x2, double x3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};

    auto prod1 = a * (b * v1);
    auto prod2 = (a * b) * v1;

    return (abs(prod1.x1 - prod2.x1) < 0.001 &&
	    abs(prod1.x2 - prod2.x2) < 0.001 &&
	    abs(prod1.x3 - prod2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double> {
	field_mult_compatibility
      }, 10000)
    );
}

TEST(TestVec3, TestMultIdentity) {
  auto test_mult_identity = [](double x1, double x2, double x3) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};

    auto prod1 = 1.0 * v1;

    return (abs(prod1.x1 - v1.x1) < 0.001 &&
	    abs(prod1.x2 - v1.x2) < 0.001 &&
	    abs(prod1.x3 - v1.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double> {
	test_mult_identity
      }, 10000)
    );
}

TEST(TestVec3, TestDistributivityWRTVectorAddition) {
  auto test_distributivity_wrt_vec_addition = [](double x1, double x2, double x3, 
						 double y1, double y2, double y3,
						 double a) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};
    auto v2 = verified_math::Vec3<double>{y1, y2, y3};

    auto prod1 = a * (v1 + v2);
    auto prod2 = (a * v1) + (a * v2);

    return (abs(prod1.x1 - prod2.x1) < 0.001 &&
	    abs(prod1.x2 - prod2.x2) < 0.001 &&
	    abs(prod1.x3 - prod2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double> {
	test_distributivity_wrt_vec_addition
      }, 10000)
    );

}

TEST(TestVec3, TestDistributivityWRTFieldAddition) {
  auto test_distributivity_wrt_field_addition = [](double x1, double x2, double x3, 
						 double a, double b) {
    auto v1 = verified_math::Vec3<double>{x1, x2, x3};


    auto prod1 = (a + b) * v1;
    auto prod2 = (a * v1) + (b * v1);

    return (abs(prod1.x1 - prod2.x1) < 0.001 &&
	    abs(prod1.x2 - prod2.x2) < 0.001 &&
	    abs(prod1.x3 - prod2.x3) < 0.001);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double> {
	test_distributivity_wrt_field_addition
      }, 10000)
    );
}
