#include "verified_math/mat33.h"
#include "gtest/gtest.h"
#include "checkpp/checkpp.h"

#include <iostream>
#include <cmath>

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

    return (fabs(prod1.x11 - prod2.x11) < 0.001 &&
	    fabs(prod1.x12 - prod2.x12) < 0.001 &&
	    fabs(prod1.x13 - prod2.x13) < 0.001 &&
	    fabs(prod1.x21 - prod2.x21) < 0.001 &&
	    fabs(prod1.x22 - prod2.x22) < 0.001 &&
	    fabs(prod1.x23 - prod2.x23) < 0.001 &&
	    fabs(prod1.x31 - prod2.x31) < 0.001 &&
	    fabs(prod1.x32 - prod2.x32) < 0.001 &&
	    fabs(prod1.x33 - prod2.x33) < 0.001);
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

      if (verified_math::det(m) < 0.0001) {
	return true;
      }
      
      auto inv_m = verified_math::inverse(m);

    auto eye = verified_math::Mat33<double> {
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
    };
      

      auto prod = m * inv_m;

    auto truth_val = (fabs(prod.x11 - eye.x11) < 0.01 &&
		      fabs(prod.x12 - eye.x12) < 0.01 &&
		      fabs(prod.x13 - eye.x13) < 0.01 &&
		      fabs(prod.x21 - eye.x21) < 0.01 &&
		      fabs(prod.x22 - eye.x22) < 0.01 &&
		      fabs(prod.x23 - eye.x23) < 0.01 &&
		      fabs(prod.x31 - eye.x31) < 0.01 &&
		      fabs(prod.x32 - eye.x32) < 0.01 &&
		      fabs(prod.x33 - eye.x33) < 0.01);

      return truth_val;
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
  EXPECT_TRUE(fabs(determinant - 1.0) < 0.0001);
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

    auto m_transpose = verified_math::transpose(m);
    
    auto det1 = verified_math::det(m);
    auto det2 = verified_math::det(m_transpose);

    return fabs(det1 - det2) < 0.001;
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

	  if (fabs(verified_math::det(m)) < 0.001) {
	    return true;
	  }

	  auto inv_m = verified_math::inverse(m);
	  auto inv_inv = verified_math::inverse(inv_m);

	  std::cout << inv_inv.x11 << " " << m.x11 << ", " 
	            << inv_inv.x12 << " " << m.x12 << ", "
	            << inv_inv.x13 << " " << m.x13 << ", "
	            << inv_inv.x21 << " " << m.x21 << ", "
	            << inv_inv.x22 << " " << m.x22 << ", "
	            << inv_inv.x23 << " " << m.x23 << ", "
	            << inv_inv.x31 << " " << m.x31 << ", "
	            << inv_inv.x32 << " " << m.x32 << ", "
	  << inv_inv.x33 << " " << m.x33 << ", " << std::endl << std::endl;


	  auto truth_val = (fabs(inv_inv.x11 - m.x11) < 0.01 &&
			    fabs(inv_inv.x12 - m.x12) < 0.01 &&
			    fabs(inv_inv.x13 - m.x13) < 0.01 &&
			    fabs(inv_inv.x21 - m.x21) < 0.01 &&
			    fabs(inv_inv.x22 - m.x22) < 0.01 &&
			    fabs(inv_inv.x23 - m.x23) < 0.01 &&
			    fabs(inv_inv.x31 - m.x31) < 0.01 &&
			    fabs(inv_inv.x32 - m.x32) < 0.01 &&
			    fabs(inv_inv.x33 - m.x33) < 0.01);
	  return truth_val;
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

    auto det1 = verified_math::det(m);
    auto det2 = verified_math::det(verified_math::inverse(m));

    return fabs(det1 - det2) < 0.0001;
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
	
    auto det1 = verified_math::det(m1 * m2);
    auto det2 = verified_math::det(m1);
    auto det3 = verified_math::det(m2);
	
    return fabs(det1 - (det2 * det3)) < 0.001;
	
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

    auto det1 = verified_math::det(c * m);
    auto det2 = (c * c * c) * verified_math::det(m);
    
    return fabs(det1 - det2) < 0.001;
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

    return fabs(verified_math::trace(m1 * m2) - verified_math::trace(m2 * m1)) < 0.000001;
  };
  
  EXPECT_TRUE(checkpp::check(checkpp::Property<double, double, double, double, double, double, double, double, double, 
			     double, double, double, double, double, double, double, double, double> {
			       mat_trace_commutative
				 }, 10000
			     )
	      );
}

