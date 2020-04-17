// #ifndef POISSMIX_H
// #define POISSMIX_H
// #ifndef FRAGINFOWITHCUTTINGSITE_H
// #define FRAGINFOWITHCUTTINGSITE_H
// #ifndef PAIR_FRAG_2D_PEAK_INFO
// #define PAIR_FRAG_2D_PEAK_INFO

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

// #include "TTree.h"
// #include "TROOT.h"
// #include "TCanvas.h"
// #include "TGaxis.h"
// #include "TH1F.h"
// #include "TRandom.h"
// #include <TApplication.h>
// #include "Math/ProbFunc.h"
// #include "TMath.h"

// C things
#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

using namespace std;


// primary c++ functions */
void chrSort(char **filelist, int count);
int extractChrNum(char *str);

int Get_Chr_Len(const string&);
int Map_Chr_Len(const string&, map<int, int>&, const string&);
int Parser_Enzyme_Cutting_Site(const string&, map< int, vector<int> >&, const string&);



double FindQT(vector<double>*);


/****************************CLASS DEFINITIONS*****************************************/

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
  long double lhood(); // Calculate the likelihood
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

// Frag_Info_With_Cutting_Site class definition
class Frag_Info_With_Cutting_Site
{
	public:
		int 	end1_chr;
		int 	end1_pos;
		int 	end1_cuttingSite;
		bool 	end1_strand;

		int 	end2_chr;
		int 	end2_pos;
		int 	end2_cuttingSite;
		bool 	end2_strand;

		Frag_Info_With_Cutting_Site():
		end1_chr(0),
		end1_pos(0),
		end1_cuttingSite(0),
		end1_strand(false),
		end2_chr(0),
		end2_pos(0),
		end2_cuttingSite(0),
		end2_strand(false){}
};

// used when combining the chrfiles
struct cuttingSite{
	unsigned int chr; unsigned int length;
	
	Frag_Info_With_Cutting_Site frag;
};
typedef struct cuttingSite cBundle;

// when all hybrid frags cannot fit in memory, and are in a file
// PoissMix& Get_Prox_Ligation_Frags_File(int, unsigned int, vector<double>&, int, int, ofstream&, int);

PoissMix& Get_Prox_Ligation_Frags_File(int , vector<double>& , float , int , ofstream&, ofstream&);

//int Pair_Frags_Interaction_Freq_File(map<int, < map <int, int> > &, 
//			   map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > &);
int Pair_Frags_Interaction_Freq_File(map<int, map< int, string > >&, 
	                                      map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >&, map<int, int>&);

int Gen_Ran_Dis_File(map<int, map< int, string > >&, 
								 map< int, vector<int> >& , 
	                      map<int, int>&, 
	                      int);


//concept of fragPair (a pair of enzyme digested frags) is different then frag (reads pair)
class Pair_Frag_2D_Peak_Info
{
	public:
		int chrNo_domain1;
		int chrNo_domain2;
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > region_info;
		int fragPair_num;
		int total_frag;
		vector<Frag_Info_With_Cutting_Site> frag_info;
		pair< pair<int, int>, pair<int, int> > peakRegion;

		bool end1_strand;
		bool end2_strand;

		int end1_regionStart;
		int end1_regionEnd;

		int end1_max;
		int end1_min;
		int end2_max;
		int end2_min;

		double sumLogProb;


		Pair_Frag_2D_Peak_Info(): 
		chrNo_domain1(0),
		chrNo_domain2(0),
		fragPair_num(0),
		total_frag(0),
		end1_max(0),
		end1_min(0),
		end2_max(0),
		end2_min(0),
		sumLogProb(0.0){}

	// used during standard use of this program
	int getSumLogProb(const map<int, double>& ranDisMap_local, 
	map<int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragMap_local, 
	const double minProb_local, const double aveSiteFrag_local, const int binSize_local)
	{
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = frag_info.begin();
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt_end = frag_info.end();
		for ( ; fragIt != fragIt_end; ++fragIt )
		{
			map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator end1Site = cuttingSiteFragMap_local[fragIt->end1_chr].find(fragIt->end1_cuttingSite);
			double correctionCoeff = (end1Site->second.size()) * 1.0 / aveSiteFrag_local;
			double minProb_corrected = minProb_local * correctionCoeff;
			double logMinProb = log(minProb_corrected);
			if (fragIt->end1_chr == fragIt->end2_chr)
			{
				int dis = abs(fragIt->end2_cuttingSite - fragIt->end1_cuttingSite);
				int binIndex = dis / binSize_local;

				map<int, double>::const_iterator disProbIt = ranDisMap_local.find(binIndex);
				if (disProbIt != ranDisMap_local.end())
				{
					sumLogProb += log((disProbIt->second)*correctionCoeff);
					return 0;
				}
				else
				{
					sumLogProb += logMinProb;
				}
			}
			else
			{
				sumLogProb += logMinProb;
			}
		}

	return 0;
	}

	// this is used during memory saving
	int getSumLogProbFile(const map<int, double>& ranDisMap_local, 
	map< int, map < int, string> >& files,
	const double minProb_local, const double aveSiteFrag_local, const int binSize_local)
	{
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = frag_info.begin();
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt_end = frag_info.end();
		for ( ; fragIt != fragIt_end; ++fragIt )
		{
			map< int, string >::const_iterator end1Site = files[fragIt->end1_chr].find(fragIt->end1_cuttingSite);

			
			// double correctionCoeff = (end1Site->second.size()) * 1.0 / aveSiteFrag_local;
			
			int fd;
			if ((fd = open(end1Site->second.c_str(), O_RDONLY) < 0)){
				printf("getSumLogProb...opening file error.\n");
				printf("Could not open file %s\n", end1Site->second.c_str());
				exit(-1);
			}

			unsigned int filesize = lseek(fd, 0, SEEK_END);
			// get # of elements
			filesize /= sizeof(cBundle);
			double correctionCoeff = filesize * 1.0 / aveSiteFrag_local;
			double minProb_corrected = minProb_local * correctionCoeff;
			double logMinProb = log(minProb_corrected);
			if (fragIt->end1_chr == fragIt->end2_chr)
			{
				int dis = abs(fragIt->end2_cuttingSite - fragIt->end1_cuttingSite);
				int binIndex = dis / binSize_local;

				map<int, double>::const_iterator disProbIt = ranDisMap_local.find(binIndex);
				if (disProbIt != ranDisMap_local.end())
				{
					sumLogProb += log((disProbIt->second)*correctionCoeff);
					return 0;
				}
				else
				{
					sumLogProb += logMinProb;
				}
			}
			else
			{
				sumLogProb += logMinProb;
			}
		}

	return 0;
	}


	/*
	int getSumLogProb(const map<int, double>&, 
		map<int, map< int, vector<Frag_Info_With_Cutting_Site> > >&,
		const double, const double, const int);
	*/
};

void fragCopy(cBundle&, const Frag_Info_With_Cutting_Site&);

void printcBundle(cBundle&);


void cbundleCopy(Frag_Info_With_Cutting_Site&, cBundle&);
