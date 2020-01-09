/**
 * @file Mini unit test framework for verifying functionality
 * of the Legendre polynomial class. This file is 
 * provided as a part of starter kit for 2019 CME 212 final exam.
 *
 * @warning You are not allowed to modify this file.
 *
 */

#ifndef _POLYNOMIAL_TEST_HPP_
#define _POLYNOMIAL_TEST_HPP_

#include <iostream>

template <class RealT, class IdxT>
class PolynomialTest
{
public:
  PolynomialTest(const RealT tol)
    : score_(0),
      num_tests_(0),
      tol_(tol)
  {}

  ~PolynomialTest() = default;

  /**
   * @brief Checks if absolute value of `zero` is less than `tol_`
   */
  void operator()(const std::string testName, const RealT& zero)
  {
    std::cout << "Testing " << testName << " ... ";
    if(abs(zero) < tol_)
    {
      std::cout << "  PASS\n";
      ++score_;
      ++num_tests_;
    }
    else
    {
      std::cout << "  FAIL\n";
      ++num_tests_;
    }
  }

  int getNumTests()
  {
    return num_tests_;
  }

  int getNumPasses()
  {
    return score_;
  }

private:
  int score_;
  int num_tests_;
  const RealT tol_;
};

#endif // _POLYNOMIAL_TEST_HPP_
