#include <string>
#include <vector>
#include <map>

#include "Frag_Info_With_Cutting_Site.hpp"

#ifndef PAIR_FRAG_2D_PEAK_INFO
#define PAIR_FRAG_2D_PEAK_INFO

using std::pair;
using std::string;
using std::vector;
using std::map;


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

	int getSumLogProb(const map<int, double>&, map<int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, const double, const double, const int);
};

int Pair_Frag_2D_Peak_Info::getSumLogProb(const map<int, double>& ranDisMap_local, map<int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragMap_local, const double minProb_local, const double aveSiteFrag_local, const int binSize_local)
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

#endif
