This is a writeup for the final exam of CME 212 by Shaked Regev

I implemented the Legendre class with the following methods:

Legendre- a constructor for the class

comp_coef - which computes the polynomial coefficients up to "morder"

"overloaded" () - which evaluates a polynomial of order n at a location x

der() - evaluates a derivative of a polynomial of order n at a location x

The code performs the desired operations. I created a vector of vectors "coefs"
to store the coefficients with size maxorder+1, and each internal size for a 
row i is i+1. I realize that since Legendre polynomials of even order have odd
coefficients equal to 0, and vice-versa, I could have saved up to half the
storage/computations. However, that would require dealing separately with even
and odd cases, so I decided to make the code simpler and more readable since
the compuation time is still asymptotically the same. 

I decided not to use the constructor directly, because all of the methods can 
potentially require computing coefficients. So I made a separate method for
that. Had I not done this, the same lines would have been repeated three times.