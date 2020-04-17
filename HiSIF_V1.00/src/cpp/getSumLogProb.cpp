// function prototypes
	int getSumLogProb(const map<int, double>& ranDisMap_local, 
	map<int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragMap_local, 
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
				printf("Could not open file %s\n", siteIt->second.c_str());
				exit(-1);
			}

			unsigned int filesize = lseek(fd, 0, SEEK_END);
			// get # of elements
			filesize =/ sizeof(cBundle);
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

