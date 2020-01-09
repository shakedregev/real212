
#ifndef _LEGENDRE_HPP_
#define _LEGENDRE_HPP_

#include <vector>
#include <cmath>
/**
 * @brief A class defining Legendre polynomial interface.
 *  Includes: a constructor for the class
 *  comp_coef - which computes the polynomial coefficients up to morder
 * "overloaded" () - which evaluates a polynomial at a location
 * der() - evaluates a derivative of a polynomial at a location
 * @tparam RealT - real number representation
 * @tparam IdxT  - integer number representation (signed or unsigned)
 *
 */
template <typename RealT, typename IdxT>
class Legendre
{
public:
  Legendre(){}
  std::vector<std::vector<RealT> > coefs;

    /**
   * @brief Constructor for a Legendre polynomial
   *Precomputes the coefficients of the polynomials
   *up to max order.
   * @param[in] maxorder polynomial order
   */
  Legendre(IdxT  maxorder)
  {
    comp_coef(maxorder);
  }
  /**
   * @brief Precomputes the coefficients of the polynomials
   *up to max order using the possibly initialized
   *coefs vector of vectors as a starting point
   * @param[in] maxorder polynomial order
   */
  void comp_coef(IdxT morder)
  {
    for(unsigned int i= (unsigned int) coefs.size(); i<=(unsigned int) morder; ++i)
    {
      std::vector<RealT> blankvec;
      coefs.push_back(blankvec);
      if (i==0)
      {
        coefs[i].push_back(1.0);
      }
      if (i==1)
      {
        coefs[i].push_back(0.0);
        coefs[i].push_back(1.0);
      }
      if (i>1) //use recurrence
      {
        for (unsigned int j = 0; j <= i-2; ++j) // from two places back
        {
          coefs[i].push_back((1.0-i)*coefs[i-2][j]/i);
        }
        coefs[i].push_back(0.0);
        coefs[i].push_back(0.0);
        for (unsigned int j = 0; j <= i-1; ++j) //from one place back
        {
          coefs[i][j+1]+=(2*i-1.0)*coefs[i-1][j]/i;
        }
      }

    }
  }
  /**
   * @brief Operator returning value of a Legendre
   * polynomial.
   *
   * @param[in] n polynomial order
   * @param[in] x polynomial argument
   * @return polynomial value evaluated at x
   *
   */
  RealT operator()(IdxT n, RealT x)
  {
    RealT y=0;
    comp_coef(n); //if coefs.size()>n the for loop doesn't occur
    for (unsigned int i = 0; i <= (unsigned int) n; ++i)
    {
      y+=coefs[n][i]*pow(x,i);
    }
    return y;
  }

  /**
   * @brief Method returning first derivative of
   * a Legendre polynomial.
   *
   * @param[in] n polynomial order
   * @param[in] x polynomial argument
   * @return polynomial derivative evaluated at x
   *
   */
  RealT der(IdxT n, RealT x)
  {
    RealT y=0;
    comp_coef(n); //if coefs.size()>n the for loop doesn't occur
    for (unsigned int i = 1; i <= (unsigned int)n; ++i)
    {
      y+=coefs[n][i]*pow(x,i-1)*i;
    }
    return y;
  }
};


#endif // _LEGENDRE_HPP_
