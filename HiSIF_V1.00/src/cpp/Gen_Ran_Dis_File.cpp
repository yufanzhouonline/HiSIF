#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cstring>
#include <string>

#include <cstdlib>
#include <fcntl.h>
#include <unistd.h>

#include <CPP_Utilities.hpp>
#include <C_Utilities.hpp>


//construct distribution of random interaction genomic distance
int Gen_Ran_Dis_File(map<int, map< int, string > >&files, 
								 map< int, vector<int> >& cuttingSitesMap_local, 
	                      map<int, int>& ranDisMap_local, 
	                      int binSize_local)
{
	int counter_frag = 0;

	//sumup genomic length
	int totalCuttingSite = 0;

	// map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin();
	// map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end = cuttingSiteFragsMap_local.end();
	
	map< int, map <int, string > >::const_iterator chrIt = files.begin(), chrIt_end = files.end();
	
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		// how many cuttingSites, size of the map
		totalCuttingSite += files[chrIt->first].size();
	}	

	chrIt = files.begin();
	for ( ; chrIt != chrIt_end; ++chrIt ){
		int sitesMinusChr = totalCuttingSite - files[chrIt->first].size();
		//site on each chromosome should have different chance to have a random inter-chromosomal ligation, the smaller the chromosome is,
		// larger the rest of the genome is, thus for each particular site on a smaller chromosome will have a better chance to form a inter chromosomal ligation
		int interChrFrag = 0; 
		// map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		map< int, string>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();

		
		int bytesread;
		cBundle data;
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			// this is the file we're looking through
			int fd = open(siteIt->second.c_str(), O_RDONLY);

			// vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();

			// while there are entries in this file
			while ((bytesread = r_read(fd, (char*)&data, sizeof(cBundle))) > 0){

				if (data.frag.end1_chr != data.frag.end2_chr)
				{
					++interChrFrag;
				}
				else
				{
					int dis = abs(data.frag.end2_cuttingSite - data.frag.end1_cuttingSite);
					int binIndex = dis / binSize_local;

					ranDisMap_local[binIndex]++;
					++counter_frag;
				}
			}

			// close after each file is opened
			close(fd);

			/*
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				if (fragIt->end1_chr != fragIt->end2_chr)
				{
					++interChrFrag;
				}
				else
				{
					int dis = abs(fragIt->end2_cuttingSite - fragIt->end1_cuttingSite);
					int binIndex = dis / binSize_local;

					ranDisMap_local[binIndex]++;
					++counter_frag;
				}
			}*/
		}

		double interChrFragProb = interChrFrag * 1.0 / sitesMinusChr;
		//cout<<chrIt->first<<"\tinter chromosomal frag prob:\t"<<interChrFragProb;
	}

	return counter_frag;
}

