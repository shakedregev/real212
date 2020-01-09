/**
 * @file Main file provided as a part of starter kit 
 * for 2019 CME 212 final exam.
 *
 * @warning You are not allowed to modify this file except to set
 * parameter `tolerance`.
 *
 */

#include <iostream>
#include <cmath>

#include "Legendre.hpp" // <~ Written by you
#include "PolynomialTest.hpp" // <~ Provided in starter kit

int main()
{
  using real_type     = double;
  using integer_type  = int;

  real_type tolerance = 1.0e-15;/* Specify tolerance here */
  
  Legendre<real_type, integer_type> P(10);

  PolynomialTest<real_type, integer_type> test(tolerance);


  test("P0' ",  P.der(0, 0.0));
  test("P1' ",  P.der(1, 0.0) - 1.0);
  test("P2  ",  P(2, 0.577350269189626));
  test("P5  ",  P(5, 0.906179845938664));
  test("P5' ",  P.der(5, 1.0) - 15.0);
  test("P6  ",  P(6,  0.661209386466265));
  test("P7  ",  P(7, -0.405845151377397));
  test("P8  ",  P(8, -0.183434642495650));
  test("P14 ",  P(14, 0.108054948707344));
  test("P16'",  P.der(16, 1.0) - 136.0);
  
  std::cout << test.getNumPasses() << " out of "
            << test.getNumTests()  << " tests passed\n";

  return 0;
}
