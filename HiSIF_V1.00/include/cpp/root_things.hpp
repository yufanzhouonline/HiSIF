#include <cmath>
#include <limits>

#define PI       3.14159265358979323846264338328      /* pi */
#define kMAXLGM 2.556348e305
#define kMACHEP  1.11022302462515654042363166809e-16
#define kMAXLOG  709.782712893383973096206318587

static double kBiginv =  2.22044604925031308085e-16;
static double kBig = 4.503599627370496e15;
static double LS2PI  =  0.91893853320467274178;

/* Logarithm of gamma function */
/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */

static double A[] = {
   8.11614167470508450300E-4,
   -5.95061904284301438324E-4,
   7.93650340457716943945E-4,
   -2.77777777730099687205E-3,
   8.33333333333331927722E-2
};

static double B[] = {
   -1.37825152569120859100E3,
   -3.88016315134637840924E4,
   -3.31612992738871184744E5,
   -1.16237097492762307383E6,
   -1.72173700820839662146E6,
   -8.53555664245765465627E5
};
   
static double C[] = {
/* 1.00000000000000000000E0, */
   -3.51815701436523470549E2,
   -1.70642106651881159223E4,
   -2.20528590553854454839E5,
   -1.13933444367982507207E6,
   -2.53252307177582951285E6,
   -2.01889141433532773231E6
};

/*
 * calculates a value of a polynomial of the form:
 * a[0]x^N+a[1]x^(N-1) + ... + a[N]
*/
double Polynomialeval(double x, double* a, unsigned int N)
{
   if (N==0) return a[0];
   else
   {
      double pom = a[0];
      for (unsigned int i=1; i <= N; i++)
         pom = pom *x + a[i];
      return pom;
   }
}

/*
 * calculates a value of a polynomial of the form:
 * x^N+a[0]x^(N-1) + ... + a[N-1]
*/
double Polynomial1eval(double x, double* a, unsigned int N)
{
   if (N==0) return a[0];
   else
   {
      double pom = x + a[0];
      for (unsigned int i=1; i < N; i++)
         pom = pom *x + a[i];
      return pom;
   }
}


double lgam( double x )
{
   double p, q, u, w, z;
   int i;

   int sgngam = 1;

   if (x >= std::numeric_limits<double>::infinity())
      return(std::numeric_limits<double>::infinity());

   if( x < -34.0 )
   {
      q = -x;
      w = lgam(q); 
      p = std::floor(q);
      if( p==q )//_unur_FP_same(p,q)
         return (std::numeric_limits<double>::infinity());
      i = (int) p;
      if( (i & 1) == 0 )
         sgngam = -1;
      else
         sgngam = 1;
      z = q - p;
      if( z > 0.5 )
      {
         p += 1.0;
         z = p - q;
      }
      z = q * std::sin( PI * z );
      if( z == 0 )
         return (std::numeric_limits<double>::infinity());
      z = std::log(PI) - std::log( z ) - w;
      return( z );
   }

   if( x < 13.0 )
   {
      z = 1.0;
      p = 0.0;
      u = x;
      while( u >= 3.0 )
      {
         p -= 1.0;
         u = x + p;
         z *= u;
      }
      while( u < 2.0 )
      {
         if( u == 0 )
            return (std::numeric_limits<double>::infinity());
         z /= u;
         p += 1.0;
         u = x + p;
      }
      if( z < 0.0 )
      {
         sgngam = -1;
         z = -z;
      }
      else
         sgngam = 1;
      if( u == 2.0 )
         return( std::log(z) );
      p -= 2.0;
      x = x + p;
      p = x * Polynomialeval(x, B, 5 ) / Polynomial1eval( x, C, 6);
      return( std::log(z) + p );
   }

   if( x > kMAXLGM )
      return( sgngam * std::numeric_limits<double>::infinity() );

   q = ( x - 0.5 ) * std::log(x) - x + LS2PI;
   if( x > 1.0e8 )
      return( q );

   p = 1.0/(x*x);
   if( x >= 1000.0 )
      q += ((   7.9365079365079365079365e-4 * p
                - 2.7777777777777777777778e-3) *p
            + 0.0833333333333333333333) / x;
   else
      q += Polynomialeval( p, A, 4 ) / x;
   return( q );
}


// prototype used here
double igamc( double a, double x );

double igam( double a, double x )
{
   double ans, ax, c, r;

   // LM: for negative values returns 1.0 instead of zero
   // This is correct if a is a negative integer since Gamma(-n) = +/- inf
   if (a <= 0)  return 1.0;  

   if (x <= 0)  return 0.0;

   if( (x > 1.0) && (x > a ) )
      return( 1.0 - igamc(a,x) );

/* Compute  x**a * exp(-x) / gamma(a)  */
   ax = a * std::log(x) - x - lgam(a);
   if( ax < -kMAXLOG )
      return( 0.0 );

   ax = std::exp(ax);

/* power series */
   r = a;
   c = 1.0;
   ans = 1.0;

   do
   {
      r += 1.0;
      c *= x/r;
      ans += c;
   }
   while( c/ans > kMACHEP );

   return( ans * ax/a );
}


double igamc( double a, double x ) 
{

   double ans, ax, c, yc, r, t, y, z;
   double pk, pkm1, pkm2, qk, qkm1, qkm2;

   // LM: for negative values returns 0.0
   // This is correct if a is a negative integer since Gamma(-n) = +/- inf
   if (a <= 0)  return 0.0;  

   if (x <= 0) return 1.0;

   if( (x < 1.0) || (x < a) )
      return( 1.0 - igam(a,x) );

   ax = a * std::log(x) - x - lgam(a);
   if( ax < -kMAXLOG )
      return( 0.0 );

   ax = std::exp(ax);

/* continued fraction */
   y = 1.0 - a;
   z = x + y + 1.0;
   c = 0.0;
   pkm2 = 1.0;
   qkm2 = x;
   pkm1 = x + 1.0;
   qkm1 = z * x;
   ans = pkm1/qkm1;

   do
   { 
      c += 1.0;
      y += 1.0;
      z += 2.0;
      yc = y * c;
      pk = pkm1 * z  -  pkm2 * yc;
      qk = qkm1 * z  -  qkm2 * yc;
      if(qk)
      {
         r = pk/qk;
         t = std::abs( (ans - r)/r );
         ans = r;
      }
      else
         t = 1.0;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;
      if( std::abs(pk) > kBig )
      {
         pkm2 *= kBiginv;
         pkm1 *= kBiginv;
         qkm2 *= kBiginv;
         qkm1 *= kBiginv;
      }
   }
   while( t > kMACHEP );

   return( ans * ax );
}


double inc_gamma_c( double a, double x) { 
   return igamc(a,x); 
} 


double gamma_cdf_c(double x, double alpha, double theta, double x0) {

      return inc_gamma_c(alpha, (x-x0)/theta);

   }

// double gamma_cdf_c(double x, double alpha, double theta, double x0 = 0);

double poisson_cdf(unsigned int n, double mu) {
      // mu must be >= 0  . Use poisson - gamma relation
      //  Pr ( x <= n) = Pr( y >= a)   where x is poisson and y is gamma distributed ( a = n+1)
      double a = (double) n + 1.0; 
      return gamma_cdf_c(mu, a, 1.0, 0);
   }


// root function
double Poisson(double x, double par)
{
  // compute the Poisson distribution function for (x,par)
  // The Poisson PDF is implemented by means of Euler's Gamma-function
  // (for the factorial), so for any x integer argument it is correct.
  // BUT for non-integer x values, it IS NOT equal to the Poisson distribution.
  // see TMath::PoissonI to get a non-smooth function.
  // Note that for large values of par, it is better to call
  //     TMath::Gaus(x,par,sqrt(par),kTRUE)
//Begin_Html
/*
<img src="gif/Poisson.gif">
*/
//End_Html

   if (x<0)
      return 0;
   else if (x == 0.0)
      return 1./exp(par);
   else {
      // double lnpoisson = x*log(par)-par-LnGamma(x+1.);
      double lnpoisson = x*log(par)-par-lgamma(x+1.);
      return exp(lnpoisson);
   }
   // An alternative strategy is to transition to a Gaussian approximation for
   // large values of par ...
   //   else {
   //     return Gaus(x,par,Sqrt(par),kTRUE);
   //   }
}

/*
Double_t TMath::LnGamma(Double_t z)
{
   // Computation of ln[gamma(z)] for all z.
   //
   // C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
   //
   // The accuracy of the result is better than 2e-10.
   //
   //--- Nve 14-nov-1998 UU-SAP Utrecht

   return ::ROOT::Math::lgamma(z);
}*/
