/**
 * @file The file containing orthogonal polynomial interface is 
 * provided as a part of starter kit for 2019 CME 212 final exam.
 *
 * @warning You are not allowed to modify this file.
 *
 */

#ifndef _ORTHOGONAL_POLYNOMIAL_HPP_
#define _ORTHOGONAL_POLYNOMIAL_HPP_

/**
 * @brief Pure virtual class defining orthogonal polynomial interface.
 *
 * @tparam RealT - real number representation
 * @tparam IdxT  - integer number representation (signed or unsigned)
 *
 */
template <typename RealT, typename IdxT>
class OrthogonalPolynomial
{
public:
  OrthogonalPolynomial(){}
  virtual ~OrthogonalPolynomial(){}

  /**
   * @brief Prototype of an operator returning value of an orthogonal
   * polynomial.
   *
   * @param[in] n polynomial order
   * @param[in] x polynomial argument
   * @return polynomial value evaluated at x
   *
   */
  virtual RealT operator()(IdxT n, RealT x) = 0;

  /**
   * @brief Prototype of a method returning first derivative of
   * an orthogonal polynomial.
   *
   * @param[in] n polynomial order
   * @param[in] x polynomial argument
   * @return polynomial derivative evaluated at x
   *
   */
  virtual RealT der(IdxT n, RealT x) = 0;
};


#endif // _ORTHOGONAL_POLYNOMIAL_HPP_
