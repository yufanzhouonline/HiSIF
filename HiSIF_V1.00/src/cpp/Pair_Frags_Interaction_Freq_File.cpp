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
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// a map of open files (that correspong to that nested map structure), will be used to access the 'vectors' there.
// note that an element of the vector before will be an entry into a file here
//
// note the structure of the data type we're reading:
// cBundle:
// 	unsigned int chr;
// 	unsigned int length;
//
// 	Frag_Info_With_Cutting_Site frag;
// 	
// reads from a file containing these structures, and does things


void cbundleCopy(Frag_Info_With_Cutting_Site&, cBundle&);

using namespace std;
int Pair_Frags_Interaction_Freq_File(map<int, map< int, string > >&files, 
	                                      map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >& freqMaptoWrite,
													  map<int, int> &Chr_Len_Map)
{
	freqMaptoWrite.clear();
	// map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(), chrIt_end = cuttingSiteFragsMap_local.end();
	map< int, map <int, string > >::const_iterator chrIt = files.begin(), chrIt_end = files.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		// map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt_begin = chrIt->second.begin(), siteIt_frag1End = chrIt->second.begin(), siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		map<int, string>::const_iterator siteIt_begin = chrIt->second.begin(), siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end(), siteIt_frag1End = chrIt->second.begin();

		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			// vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();

			// this is the file we're looking through
			int fd = open(siteIt->second.c_str(), O_RDONLY);

			int bytesread;
			cBundle data;
			
						
			while ((bytesread = r_read(fd, (char*)&data, sizeof(cBundle))) > 0){
				int frag1_start = 0;
				int frag1_end = 0;

				if (data.frag.end1_strand){
					frag1_end = data.chr;
					if (siteIt != siteIt_begin)
					{
						siteIt_frag1End--;
						frag1_start = siteIt_frag1End->first;
					}
					else
					{
						frag1_start = 1;
					}

				}
				else{
					frag1_start = siteIt->first;

					siteIt_frag1End++;
					if (siteIt_frag1End != siteIt_end)
					{
						frag1_end = siteIt_frag1End->first;
					}
					else
					{
						frag1_end = Chr_Len_Map[data.frag.end1_chr];
					}
				}

				// end2
				map< int, map< int, string > >::const_iterator chrIt_end2 = files.find(data.frag.end2_chr);

				// if this is in here
				if (chrIt_end2 != chrIt_end)
				{
					map< int, string >::const_iterator siteIt_end2_begin = chrIt_end2->second.begin(), siteIt_frag2End = chrIt_end2->second.begin(), siteIt_end2 = chrIt_end2->second.find(data.frag.end2_cuttingSite), siteIt_end2_end = chrIt_end2->second.end();
					int frag2_start = 0;
					int frag2_end = 0;
					siteIt_frag2End = siteIt_end2;
					if (data.frag.end2_strand)
					{
						frag2_end = siteIt_end2->first;
						if (siteIt_end2 != siteIt_end2_begin)
						{
							siteIt_frag2End--;
							frag2_start = siteIt_frag2End->first;
						}
						else
						{
							frag2_start = 1;
						}
					}
					else
					{
						frag2_start = siteIt_end2->first;

						siteIt_frag2End++;
						if (siteIt_frag2End != siteIt_end2_end)
						{
							frag2_end = siteIt_frag2End->first;
						}
						else
						{
							frag2_end = Chr_Len_Map[data.frag.end2_chr];
						}
					}
					// convert this element to a Frag_Info_With_Cutting_Site object
					Frag_Info_With_Cutting_Site frag;
					cbundleCopy(frag, data);
					freqMaptoWrite[make_pair(make_pair(chrIt->first, make_pair(frag1_start, frag1_end)), make_pair(chrIt_end2->first, make_pair(frag2_start, frag2_end)))].push_back(frag);
				}


			}

			// error checking
			if (bytesread == -1){
				fprintf(stderr, "Pair_Frags...read error\n");
				exit(-1);
			}

			/*
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				int frag1_start = 0; //a pair of interacting frag, frag1 is end1, frag2 is end2
				int frag1_end = 0;
				siteIt_frag1End = siteIt;

				// positive strand
				if (fragIt->end1_strand)
				{
					frag1_end = siteIt->first;
					if (siteIt != siteIt_begin)
					{
						siteIt_frag1End--;
						frag1_start = siteIt_frag1End->first;
					}
					else
					{
						frag1_start = 1;
					}
				}// negative strand
				else
				{
					frag1_start = siteIt->first;

					siteIt_frag1End++;
					if (siteIt_frag1End != siteIt_end)
					{
						frag1_end = siteIt_frag1End->first;
					}
					else
					{
						frag1_end = Chr_Len_Map[fragIt->end1_chr];
					}
				}

				// second end
				map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end2 = cuttingSiteFragsMap_local.find(fragIt->end2_chr);
				if (chrIt_end2 != chrIt_end)
				{
					map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt_end2_begin = chrIt_end2->second.begin(), siteIt_frag2End = chrIt_end2->second.begin(), siteIt_end2 = chrIt_end2->second.find(fragIt->end2_cuttingSite), siteIt_end2_end = chrIt_end2->second.end();
					int frag2_start = 0;
					int frag2_end = 0;
					siteIt_frag2End = siteIt_end2;
					if (fragIt->end2_strand)
					{
						frag2_end = siteIt_end2->first;
						if (siteIt_end2 != siteIt_end2_begin)
						{
							siteIt_frag2End--;
							frag2_start = siteIt_frag2End->first;
						}
						else
						{
							frag2_start = 1;
						}
					}
					else
					{
						frag2_start = siteIt_end2->first;

						siteIt_frag2End++;
						if (siteIt_frag2End != siteIt_end2_end)
						{
							frag2_end = siteIt_frag2End->first;
						}
						else
						{
							frag2_end = Chr_Len_Map[fragIt->end2_chr];
						}
					}
					freqMaptoWrite[make_pair(make_pair(chrIt->first, make_pair(frag1_start, frag1_end)), make_pair(chrIt_end2->first, make_pair(frag2_start, frag2_end)))].push_back(*fragIt);
				}
			}*/
		}
	}
	return 0;
}

