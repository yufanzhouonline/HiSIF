#ifndef POISSMIX_H
#define POISSMIX_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TRandom.h"
#include <TApplication.h>
#include "Math/ProbFunc.h"
#include "TMath.h"

using namespace std;

class PoissMix{
 public: 
  PoissMix(int N);
  double Gammaln(double y);
  double PoissPDF(double x,double lambda);
  double pSum(int j);
  double pIJ(int i, int j);
  void W_update();
  void Mu_update();
  void update();  //replace W and Mu with the new one
  double lhood(); // Calculate the likelihood
  double BIC();   // Calculate the BIC of the maximum likelihood
  double AIC();   // Calculate the AIC of the maximum likelihood
  double PoissCDF(double x,double lambda);
  double FDRThreshold(double t,double lambda1,double lambda2,double alpha1,double alpha2);
  double PValue(double t,double lambda1,double lambda2,double alpha1,double alpha2);

  vector<double> x, W, Mu;
  int N;  // number of mixture
  int M;  // number of data points
  double ep; // precision goal
  map<double, map<double, vector<int> > > mixpara ;

  private: 
  vector<double> Wt, Mut; // store the update
   
};
#endif
