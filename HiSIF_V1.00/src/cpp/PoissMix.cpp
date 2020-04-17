#ifndef POISSMIX_CPP
#define POISSMIX_CPP

#include "CPP_Utilities.hpp"
#include "root_things.hpp"

//-------------------------------------------------------------------------
PoissMix::PoissMix(int N){
  	for(int k=0; k<N; k++){
  		Wt.push_back(0.); 
  		Mut.push_back(0.);
  	}
  }
//-------------------------------------------------------------------------
// Returns the value log(Gamma(y)) for y > 0. See http://www.nr.com/.
double PoissMix::Gammaln(double y)
{
    double x, z, Tmp, Ser;
    int   j;

    static double Cof[6] = { (double)76.18009172947146, -(double)86.50532032941677,
        (double)24.01409824083091, -(double)1.231739572450155,
        (double)0.1208650973866179E-2, -(double)0.5395239384953E-5 };

    static double Stp = (double)2.5066282746310005;

    z = x = y; Tmp = x + (double)5.5; Tmp -= (x + (double)0.5) * (double)log(Tmp);

    Ser = (double)1.000000000190015;

    for (j = 0; j < 6; j++) Ser += Cof[j] / ++z;

    return (-Tmp + (double)log(Stp * Ser / x));
} // Gammaln
//-------------------------------------------------------------------------
double PoissMix::PoissPDF(double x,double lambda){
	return exp(x * log(lambda) - Gammaln(x + 1.0) - lambda);
}
//-------------------------------------------------------------------------
double PoissMix::pSum(int j){
  double res = 0.;
  for(int k = 0; k < N; k++){
      res = res + PoissPDF(x[j],Mu[k])*W[k];
    }
  return res+10e-20;
}
//-------------------------------------------------------------------------
double PoissMix::pIJ(int i, int j){
  double res;
  res = PoissPDF(x[j],Mu[i])*W[i]/pSum(j);
  return res;
}
//-------------------------------------------------------------------------
void PoissMix::W_update(){
  double res;
  for(int i = 0; i < N; i++){
    res = 0.;
    for(int j = 0; j < M; j++){
      res = res + pIJ(i,j);
    }
    Wt[i] = res/double(M);
  }
}
//-------------------------------------------------------------------------
void PoissMix::Mu_update(){
  double resu, resd;
  
  for(int i = 0; i < N; i++){
    resu = 0.;
    resd = 0.;
    for(int j = 0; j < M; j++){
      resu = resu + x[j]*pIJ(i,j);
      resd = resd + pIJ(i,j);
    }
    Mut[i] = resu/resd;
  }
 }
//-------------------------------------------------------------------------
void PoissMix::update(){
  for(int i=0; i<N; i++)
    {
      W[i] = Wt[i];
      Mu[i] = Mut[i];
    }
}
//-------------------------------------------------------------------------
long double PoissMix::lhood(){
  long double res2, pdf;
  for(int j=0;j<M;j++){
      res2=0.;
      for(int i=0; i<N; i++){
		  pdf = W[i]*PoissPDF(x[j],Mu[i]);
	     // printf("PoissMix::lhood: pdf == %.17g\n",  pdf);

        res2 = res2 + pdf;
	     // printf("PoissMix::lhood: res2 == %.17g\n",  res2);
      }
      res2 += log(res2);
    }
	  // printf("PoissMix::lhood: final res2 == %.17g\n",  res2);
  return(res2);    
}
//-------------------------------------------------------------------------
double PoissMix::BIC(){
  double res;
  res=-2.0*log(exp(lhood()))+double(3.*N-1.)*log(double(M));
  return(res);
}
//-------------------------------------------------------------------------
double PoissMix::AIC(){
  double res;
  res=-2.0*log(exp(lhood()))+double(3.*N-1.)*2.;
  return(res);
}
//-------------------------------------------------------------------------
double PoissMix::PoissCDF(double x,double lambda){
	double value = poisson_cdf(x,lambda);
	return value;
}

double PoissMix::PValue(double t,double lambda1,double lambda2,double alpha1,double alpha2){
	double pt = PoissPDF(t, lambda1)*alpha1/(alpha1 + alpha2) + PoissPDF(t, lambda2)*alpha2 / (alpha1 + alpha2);
        double pv = 0;
	for (int j = 0; j < 200; j++) {
		double pj = PoissPDF(j, lambda1)*alpha1 / (alpha1 + alpha2) + PoissPDF(j, lambda2)*alpha2 / (alpha1 + alpha2);
                if (pj <= pt) {
			pv = pv + pj;
		}
	}
	return pv;
}
//-------------------------------------------------------------------------
double PoissMix::FDRThreshold(double t,double lambda1,double lambda2,double alpha1,double alpha2){
        double numrator = (1-PoissCDF(t-1,lambda1))*alpha1;
	double denrator = (1-PoissCDF(t-1,lambda1))*alpha1+(1-PoissCDF(t-1,lambda2))*alpha2;
	return numrator/denrator;
        
        /*
        double arr[200];
	double arr_p[200];
        int t = int(t1);
	for (int i = 0; i < 200; i++) {
		arr_p[i] = PoissPDF(i, lambda1)*alpha1 / (alpha1 + alpha2) + PoissPDF(i, lambda2)*alpha2 / (alpha1 + alpha2);
	}
	for (int i = 0; i < 200; i++) {
		double pv = 0;
		for (int j = 0; j < 200; j++) {
			if (arr_p[j] <= arr_p[i]) {
				pv = pv + arr_p[j];
			}
		}
		arr[i] = pv;
	}
	int rk = 1;
	for (int i = 0; i < 200; i++) {
		if (arr[i] < arr[t]) {
			rk++;
		}
	}
	double fdr = arr[t] * rk / 200;
	return fdr;
        */
}
//-------------------------------------------------------------------------
#endif
