#include "verified_math/mat33.h"
#include "gtest/gtest.h"
#include "checkpp/checkpp.h"

#include <iostream>
#include <cmath>

#define epsilon 0.1

TEST(TestMat33, TestEyeCommutative) {
  auto eye_commutative = [](double x11, double x12, double x13,
			    double x21, double x22, double x23,
			    double x31, double x32, double x33) {
    auto m = verified_math::Mat33<double> {
      x11, x12, x13,
      x21, x22, x23,
      x31, x32, x33      
    };

    auto eye = verified_math::Mat33<double> {
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
    };

    auto prod1 = eye * m;
    auto prod2 = m * eye;

    return (fabs(prod1.x11 - prod2.x11) < epsilon &&
	    fabs(prod1.x12 - prod2.x12) < epsilon &&
	    fabs(prod1.x13 - prod2.x13) < epsilon &&
	    fabs(prod1.x21 - prod2.x21) < epsilon &&
	    fabs(prod1.x22 - prod2.x22) < epsilon &&
	    fabs(prod1.x23 - prod2.x23) < epsilon &&
	    fabs(prod1.x31 - prod2.x31) < epsilon &&
	    fabs(prod1.x32 - prod2.x32) < epsilon &&
	    fabs(prod1.x33 - prod2.x33) < epsilon);
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, 
			     double, double, double, 
			     double, double, double> {
			       eye_commutative }, 10000
			     )
	      );
}

TEST(TestMat33, TestMatInv) {
  auto mat_inv = [](double x11, double x12, double x13,
		    double x21, double x22, double x23,
		    double x31, double x32, double x33) {
    auto m = verified_math::Mat33<double> {
	x11, x12, x13,
	x21, x22, x23,
	x31, x32, x33      
    };

    auto kappa = condition_number(m);
    
    auto inv_m = verified_math::inverse(m);

    auto eye = verified_math::Mat33<double> {
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
    };
      
    auto prod = m * inv_m;

    auto truth_val = (fabs(prod.x11 - eye.x11) < epsilon &&
		      fabs(prod.x12 - eye.x12) < epsilon &&
		      fabs(prod.x13 - eye.x13) < epsilon &&
		      fabs(prod.x21 - eye.x21) < epsilon &&
		      fabs(prod.x22 - eye.x22) < epsilon &&
		      fabs(prod.x23 - eye.x23) < epsilon &&
		      fabs(prod.x31 - eye.x31) < epsilon &&
		      fabs(prod.x32 - eye.x32) < epsilon &&
		      fabs(prod.x33 - eye.x33) < epsilon);

      return truth_val || kappa > 1.1;
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double> {
	mat_inv }, 1000)
  );      
}

/*
  Tests related to the properties of the determinant.
*/
TEST(TestMat33, TestDetEyeIsOne) {
  auto m = verified_math::Mat33<double> {
    1.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
  };

  auto determinant = det(m);
  EXPECT_TRUE(fabs(determinant - 1.0) < epsilon);
}

TEST(TestMat33, TestDetTransposeInvariant) {
  auto det_transpose_invariant = [](double x11, double x12, double x13,
	      double x21, double x22, double x23,
	      double x31, double x32, double x33) {
    auto m = verified_math::Mat33<double> {
      x11, x12, x13,
      x21, x22, x23,
      x31, x32, x33
    };

    auto kappa = condition_number(m);
        
    auto m_transpose = verified_math::transpose(m);
    
    auto det1 = verified_math::det(m);
    auto det2 = verified_math::det(m_transpose);

    return fabs(det1 - det2) < epsilon || kappa > 1.1;
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double,
			     double, double, double,
			     double, double, double> {
			       det_transpose_invariant
				 }, 10000)
	      );
}

    TEST(TestMat33, TestMatInvInv) {
      auto inv_inv = [](double x11, double x12, double x13,
			double x21, double x22, double x23,
			double x31, double x32, double x33) {
	auto m = verified_math::Mat33<double> {
	    x11, x12, x13,
	    x21, x22, x23,
	    x31, x32, x33
	  };

	auto kappa = condition_number(m);

	  auto inv_m = verified_math::inverse(m);
	  auto inv_inv = verified_math::inverse(inv_m);

	  auto truth_val = (fabs(inv_inv.x11 - m.x11) < epsilon &&
			    fabs(inv_inv.x12 - m.x12) < epsilon &&
			    fabs(inv_inv.x13 - m.x13) < epsilon &&
			    fabs(inv_inv.x21 - m.x21) < epsilon &&
			    fabs(inv_inv.x22 - m.x22) < epsilon &&
			    fabs(inv_inv.x23 - m.x23) < epsilon &&
			    fabs(inv_inv.x31 - m.x31) < epsilon &&
			    fabs(inv_inv.x32 - m.x32) < epsilon &&
			    fabs(inv_inv.x33 - m.x33) < epsilon);
	  return truth_val || kappa > 1.1;
	  };

	EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double> {
	      inv_inv
	    }, 10000)
	  );
    }

TEST(TestMat33, TestDetInverse) {
  auto det_inverse = [](double x11, double x12, double x13,
	      double x21, double x22, double x23,
	      double x31, double x32, double x33) {
    auto m = verified_math::Mat33<double> {
      x11, x12, x13,
      x21, x22, x23,
      x31, x32, x33
    };

    auto kappa = condition_number(m);

    auto det1 = verified_math::det(m);
    auto det2 = verified_math::det(verified_math::inverse(m));

    return fabs(det1 - det2) < epsilon || kappa > 1.1;
  };

  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double> {
	det_inverse}, 10000
      )
    );
}

TEST(TestMat33, TestDetIsHomomorphic) {
  auto det_homomorphic = [](double x11, double x12, double x13, double x21, double x22, double x23, double x31, double x32, double x33,
			    double y11, double y12, double y13, double y21, double y22, double y23, double y31, double y32, double y33) {
    
    verified_math::Mat33<double> m1 = verified_math::Mat33<double>{x11, x12, x13, x21, x22, x23, x31, x32, x33 };
    verified_math::Mat33<double> m2 = verified_math::Mat33<double> {y11, y12, y13, y21, y22, y23, y31, y32, y33 };

    auto kappa1 = condition_number(m1);
    auto kappa2 = condition_number(m2);
    
    auto det1 = verified_math::det(m1 * m2);
    auto det2 = verified_math::det(m1);
    auto det3 = verified_math::det(m2);
	
    return fabs(det1 - (det2 * det3)) < epsilon || kappa1 > 1.1 || kappa2 > 1.1;
	
  };
	
  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double, 
			     double, double, double, double, double, double, double, double, double> {
		det_homomorphic }, 10000)
  );      
}
      
TEST(TestMat33, TestScalarPowerInDet) {
  auto scalar_power_in_det = [](double x11, double x12, double x13,
				double x21, double x22, double x23,
				double x31, double x32, double x33, double c) {
    auto m = verified_math::Mat33<double> {
      x11, x12, x13,
      x21, x22, x23,
      x31, x32, x33
    };

    auto kappa = condition_number(m);

    auto det1 = verified_math::det(c * m);
    auto det2 = (c * c * c) * verified_math::det(m);

    return fabs(det1 - det2) < epsilon || kappa > 1.1;
  };
  
  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double, double>{
	scalar_power_in_det
	  }, 10000)
    );
}

/*
  Tests related to the trace.
 */
TEST(TestMat33, TestMatMatTraceCommutative) {
  auto mat_trace_commutative = 
    [](double x11, double x12, double x13,
       double x21, double x22, double x23,
       double x31, double x32, double x33,
       double y11, double y12, double y13,
       double y21, double y22, double y23,
       double y31, double y32, double y33) {
    
    verified_math::Mat33<double> m1 = verified_math::Mat33<double>{
      x11, x12, x13,
      x21, x22, x23,
      x31, x32, x33
    };

    verified_math::Mat33<double> m2 = verified_math::Mat33<double> {
      y11, y12, y13,
      y21, y22, y23,
      y31, y32, y33
    };

    auto kappa1 = condition_number(m1);
    auto kappa2 = condition_number(m2);

    auto truth_val = fabs(verified_math::trace(m1 * m2) - verified_math::trace(m2 * m1)) < epsilon;
    
    return truth_val || kappa1 > 1.1 || kappa2 > 1.1;
    
  };
  
  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double> {
	mat_trace_commutative
	  }, 10000
      )
    );
}

