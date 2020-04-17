//#ifndef HICPEAK_CPP
//#define HICPEAK_CPP

#include <algorithm>										// allows std::sort
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <cstring>
#include <string>
#include <vector>


// #include <TDirectory.h>
// #include <TSystem.h>
// #include <TAxis.h>
// use of Poisson function
// #include <TMath.h>
// #include <TF1.h>
// #include <TGraph.h>
// #include <TApplication.h>

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/sysinfo.h>

// used to offset the large structure
#include <FragDumper.hpp>
// the functions we use in root here
// #include <root_things.hpp>


double Poisson(double x, double par);

using namespace std;
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// comparison function for Frag_Info_With_Cutting_Site data object, starting with checking End1 */
// TODO possible improvement of boolean logic
bool Comp_BothEndsMappedFrag_End1(const Frag_Info_With_Cutting_Site& lhs, const Frag_Info_With_Cutting_Site& rhs)
{
	if (lhs.end1_chr!=rhs.end1_chr)
	{
		return lhs.end1_chr<rhs.end1_chr;
	}
	else if (lhs.end1_pos!=rhs.end1_pos)
	{
		return lhs.end1_pos<rhs.end1_pos;
	}
	else if (lhs.end2_chr!=rhs.end2_chr)
	{
		return lhs.end2_chr<rhs.end2_chr;
	}
	else
	{
		return lhs.end2_pos<rhs.end2_pos;
	}
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// comparison function for Frag_Info_With_Cutting_Site data object, starting with checking End2 */
// TODO possible improvement of boolean logic
bool Comp_BothEndsMappedFrag_End2(const Frag_Info_With_Cutting_Site& lhs, const Frag_Info_With_Cutting_Site& rhs)
{
	if (lhs.end2_chr!=rhs.end2_chr)
	{
		return lhs.end2_chr<rhs.end2_chr;
	}
	else if (lhs.end2_pos!=rhs.end2_pos)
	{
		return lhs.end2_pos<rhs.end2_pos;
	}
	else if (lhs.end1_chr!=rhs.end1_chr)
	{
		return lhs.end1_chr<rhs.end1_chr;
	}
	else
	{
		return lhs.end1_pos<rhs.end1_pos;
	}
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void swap_num(Frag_Info_With_Cutting_Site &a, Frag_Info_With_Cutting_Site &b)  // NEW
{

    if (a.end1_pos > b.end2_pos){
    	a.end1_pos = a.end1_pos + b.end2_pos;
    	b.end2_pos = a.end1_pos - b.end2_pos;
    	a.end1_pos = a.end1_pos - b.end2_pos;
    }    
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class HiC_Peak{
private:
	string dir_genomeFasta;
	string file_source;
	string outputfile;
	string file_cuttingSites;
    
    //size and number of reads will be determined after reading file(readSize:parser_eland,readsNum:reader_gff)
	int readSize;
	int extent_cuting_size;
	bool flag_thres;
	int modelBinSize;

	vector<int> binsizeVec;
	vector<double> poisVec;
	vector<double> percVec;

	// the length of each chromosome in the form: <1, length>, <2, length>, ..., <24, length>
	map<int, int> Chr_Len_Map;
	
	PoissMix* two;            // Create Poisson mixture model
	// TH1F* hist;
	// TF1* pois;

	inline void usage() const
	{
		// cout <<"Usage: HiC_Peak -w read_length cuttingSiteExtent fragmentextent -T peakthreshold -c cuttingSitesMap -o outputfile inputDirectory"<<endl;
		
		fprintf(stderr, "\nProgram: HiSIF - HiC Significant Interaction Fragments\n");
		fprintf(stderr, "Version: 1.0.0\nHiSIF [options] <inputDirectory>\n\n");
		
		// options
		fprintf(stderr, "\t-g <DIR>\t\treference genome directory\n");
		fprintf(stderr, "\t-c <FILE>\t\tcutting sites map .bed file\n");
		fprintf(stderr, "\t-p <INT> <INT>\t\tpoisson mixture model parameters\n");
		fprintf(stderr, "\t-w <INT> <INT> <INT>\treadLength, cuttingSiteExtent, binSize\n");
		fprintf(stderr, "\t-t <INT>\t\tpeakThreshold value, 1, 1.5, 2 and so on\n");
		fprintf(stderr, "\t-s <0.0-1.0>\t\tpercentage of dataset for bootstrapping, default is 1\n");
		fprintf(stderr, "\t-i <INT>\t\tbootstrapping iterations, default is 50\n");
		fprintf(stderr, "\t-f <INT>\t\toutput fragment size, default is the same as binSize\n");
		//fprintf(stderr, "\t-o <FILE>\t\toutputfile\n");
		fprintf(stderr, "\t-x\t\t\tlimit number of child processes used, default is 0 (no limit)\n\n");
		fprintf(stderr, "\tFor example:\n\n");
		fprintf(stderr, "\tRun the following for HindIII digested Hi-C experiments\n\n");		
		fprintf(stderr, "\tbin/HiSIF -g <hg19genome> -c <resources/hg19.HindIII.bed> -p 1 29 -w 50 500 20000 -t 1 -i 2 chrfiles\n\n");
		fprintf(stderr, "\tRun the following for MboI digested Hi-C experiments\n\n");		
		fprintf(stderr, "\tbin/HiSIF -g <hg19genome> -c <resources/hg19.MboI.bed> -p 1 29 -w 50 500 3000 -t 1 -i 2 chrfiles\n\n");
		fprintf(stderr, "\tThe significant interaction fragments will be output to the file chrfiles_t1_peak.txt with columns:\n\n");
		fprintf(stderr, "\tchr1\tstart1\tend1\tchr2\tstart2\tend2\tcounts\tFDR\n\n");
	}

public:
	int Extend_Cutting_Site(map< int, vector<int> >&, map< int, map< bool, vector< pair<int, int> > > >&, int&);
	void Make_Restriction_Vector(const string&, vector<Frag_Info_With_Cutting_Site>&);
	pair <int, int> Map_Hybrid_Frag2_CuttingSite_End1(vector<Frag_Info_With_Cutting_Site>& , map< int, map< bool, vector< pair<int, int> > > >&, const int);
	pair <int, int> Map_Hybrid_Frag2_CuttingSite_End2(vector<Frag_Info_With_Cutting_Site>& , map< int, map< bool, vector< pair<int, int> > > >&, const int);
	int Get_Cutting_Site_Frags(const map< int, vector<int> >&, const vector<Frag_Info_With_Cutting_Site>&, map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, int, ofstream&);
	int Pair_Frags_Interaction_Freq(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >&);
	int Pair_Frags_Interaction_Freq_With_Both_Side(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, map< int, map<int, int> >&, map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >&);
	int Site_Vec2_Frag_Map(const map< int, vector<int> >&, map< int, map<int, int> >&);
	
	//2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, p value, score etc> >
	void Get_Interacting_Sites(const map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >, map< int, map<int, int> >&, const double, vector<Pair_Frag_2D_Peak_Info>&); 

	// modified function to use file i/o class
	void Get_Interacting_Sites_File(const map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >, map< int, map<int, int> >&, const double, vector<Pair_Frag_2D_Peak_Info>&);
        void Get_Interacting_Sites_random(const map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >, vector<Pair_Frag_2D_Peak_Info>&); 
	int Gen_Ran_Dis(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, map< int, vector<int> >&, map<int, int>&, const int);
	int Freq_To_Prob(const map<int, int>&, const int, map<int, double>&);
	
	//frag here is enzyme cutted frag, map<pair<chr1,pair<frag1Start, frag1End> >, pair<chr2,pair<frag2STart, frag2End> > >, pair<number of interaction, vector<interaction frag> > >
	int Search_Neighbour_Frags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >&, map< int, map<int, int> >&, int, int, pair<int, int>, pair<int, int>, int, int, vector< pair< pair< pair<int, int>, pair<int, int> >, int > >&, double, int&, int&, vector<Frag_Info_With_Cutting_Site>&);
	// modified function to use file i/o class
	int Search_Neighbour_Frags(FragDumper*, map< int, map<int, int> >&, int, int, pair<int, int>, pair<int, int>, int, int, vector< pair< pair< pair<int, int>, pair<int, int> >, int > >&, double, int&, int&, vector<Frag_Info_With_Cutting_Site>&);
	int Find_Peak_Summit(vector<Pair_Frag_2D_Peak_Info>&);

	// new function
	double FindQT(vector<double> *para_values);

	int Writer_Frag_Interaction_Map(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >, const string);
	int Writer_Frag_Freq_Map(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local, const string fragFreqFiletoWrite);
	int Writer_Digested_Frags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local,const string fragInteractionFiletoWrite);
	void Writer_PoissMix(map<unsigned int,unsigned int>&, string&);
	void BinWriter_PoissMix(map<unsigned int,unsigned int>& , string&);
	void Writer_FragSize(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, map<unsigned int,unsigned int>&, map<unsigned int, unsigned int>&, string&);
	void Writer_LigateProb(map<int, int>&, string&);
	int Writer_Dis_Distribution(const map<int, int>&, const string&);
	int Writer_Dis_Distribution_Prob(const map<int, double>&, const string&);
	int Writer_Cut_Site_Frags(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, const string&);
	int Writer_2D_Peak(vector<Pair_Frag_2D_Peak_Info>&, const string&);
	int Writer_2D_Peak_with_FDR(vector<Pair_Frag_2D_Peak_Info>&,PoissMix*mix ,const string&);
    int Writer_2D_Peak_with_FDR_Random(vector<Pair_Frag_2D_Peak_Info>&,vector<Pair_Frag_2D_Peak_Info>&,const double,const string&);
	int Writer_2D_Peak_with_FDR_Random_And_Length_Limit(vector<Pair_Frag_2D_Peak_Info>&, vector<Pair_Frag_2D_Peak_Info>&,	const double, const string&, int);

	// remove large thing from memory to prevent crashes (dirty thing)
	void fragMemDump(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, string filename);
	void fragMemRecovery(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, string filename);


	// void fragCopy(cBundle& , const Frag_Info_With_Cutting_Site&);
	// void cbundleCopy(Frag_Info_With_Cutting_Site&, cBundle&);
	// void printcBundle(cBundle& bundle);
	double poissonf(double*x,double*par);


	
	HiC_Peak(int argc, char* argv[]): 
	// basically like calling a constructor
	readSize(0),
	extent_cuting_size(0),
	flag_thres(false),
	modelBinSize(0),
	binsizeVec(0),
	poisVec(0),
	percVec(0),
	two(0){
		if (argc < 2){
			usage();
			exit(1);
		}
	// hist(0),
	/*pois(0){
		}*/

		//number of option arguments (arguments that not source file name)
		int optionArg 		= 0;
		optionArg 			= argc - 2;
		int binSize 		= 0;
		int poisspa			= 0;
		float pSampleSize = 1;					// used to determine % of total dataset size
		int saveMemory = 0;
		int maxProc = 0;
		int iterations = 50;							// how many iterations?
		int fragmentSize = 0;
		double percentNo 	= 0.0;

		// Map Structures */



		// Get cutting sites according to chromosome to a map where map< chr, vector<pos> >
		map< int, vector<int> > cuttingSitesMap;

		//map< chr, map< strand, vector< pair<start, end> > > >
		map< int, map< bool, vector< pair<int, int> > > > cuttingSiteExtMap;


		// options parsing!!! */
		//options -w(bin size)[20,1000], -p(threshold percentage)[0.5,1.0), -s(input file is sorted), -f(calculate FDR);
		int argIt 			= 1;
		while ( argIt <= optionArg ){
			if (argv[argIt][0] == '-'){
				switch (argv[argIt][1]){

					case 'g':
						if ( strlen(argv[argIt]) != 2)
						{
							cout <<"Error: -g option can not be combined with other options."<<endl;
							exit(1);
						}

						++argIt;
						if ( argv[argIt][0] == '-' || argIt > optionArg )
						{
							cout <<"Specify the directory contains genome sequence information following -g option."<<endl;
							exit(1);
						}
						else
						{
							dir_genomeFasta = argv[argIt];
							++argIt;
							if ( argIt <= optionArg && argv[argIt][0] != '-' )
							{
								cout <<"Only one genome could be chosen for one data."<<endl;
								exit(1);
							}
						}
						break;
					
					case 'w':
						if ( strlen(argv[argIt]) != 2)
						{
							cout <<"Error: -w option can not be combined with other options."<<endl;
							exit(1);
						}

						++argIt;
						if ( argv[argIt][0] == '-' || argIt > optionArg )
						{
							cout <<"Error: no bin size specified following -w option."<<endl;
							exit(1);
						}
						while ( argv[argIt][0] != '-' && argIt <= optionArg ){
							binSize = atoi(argv[argIt]);
							binsizeVec.push_back(binSize);
							++argIt;
						}
						break;

					case 'p':
						if ( strlen(argv[argIt]) != 2)
						{
							cout <<"Error: -p option can not be combined with other options."<<endl;
							exit(1);
						}

						++argIt;
						if ( argv[argIt][0] == '-' || argIt > optionArg )
						{
							cout <<"Error: no mixture parameters specified after -p option."<<endl;
							exit(1);
						}
						while ( argv[argIt][0] != '-' && argIt <= optionArg ){
							poisspa = atof(argv[argIt]);
							poisVec.push_back(poisspa);
							++argIt;
						}
						break;

					//HiC_Peak -t 4.5 6.4 -w 50.....
					case 't':
						flag_thres = true;
						if ( strlen(argv[argIt]) != 2){
							cout <<"Error: -t option can not be combined with other options."<<endl;
							exit(1);
						}
						++argIt;
						if ( argv[argIt][0] == '-' || argIt > optionArg ){
							cout <<"Error: no threshold specified following -p option."<<endl;
							exit(1);
						}
						while ( argv[argIt][0] != '-' && argIt <= optionArg ){
							percentNo = atof(argv[argIt]);
							percVec.push_back(percentNo);
							++argIt;
						}
						break;

					// save memory flag
					case 'm':
						++argIt;
						saveMemory = atoi(argv[argIt]);

						if (saveMemory < 1)
							saveMemory = 1;
						if (saveMemory > 2)
							saveMemory = 2;
                                                ++argIt;
						break;
					// less processes flag
					case 'x':
						++argIt;
						// -1, for the main process
						maxProc = atoi(argv[argIt])-1;
						++argIt;
						break;

					case 'o':
						++argIt;
						outputfile = argv[argIt];
						++argIt;
						break;
				
					case 's':
						++argIt;
						pSampleSize == atof(argv[argIt]);
						++argIt;
						break;

					case 'i':
						++argIt;
						iterations = atoi(argv[argIt]);
						++argIt;
						break;

					case 'f':
						++argIt;
						fragmentSize = atoi(argv[argIt]);
						++argIt;
						break;

					case 'c':
						if ( strlen(argv[argIt]) != 2){
							cout <<"Error: -c option can not be combined with other options."<<endl;
							exit(1);
						}
						++argIt;
						if ( argv[argIt][0] == '-' || argIt > optionArg ){
							cout <<"Specify the name of the enzyme cutting site map file following -c option."<<endl;
							exit(1);
						}
						else{
							file_cuttingSites = argv[argIt];
							++argIt;
							if ( argIt <= optionArg && argv[argIt][0] != '-' ){
								cout <<"Only one enzyme cutting site map file could be used for one data."<<endl;
								exit(1);
							}
						}
						break;

										
					default:
					cout <<"Error: Undefined option: "<<endl;
					usage();
					exit(1);		
				}
			}
			else{
				cout <<"Error: Undefined argument: "<<argv[argIt]<<endl;
				usage();
				exit(1);
			}
		}

		// end of parsing */

		two= new PoissMix(2);


		// last of the arguments
		file_source = argv[argc-1];
		// remove trailing '/'
		if (file_source.find_last_of('/') == file_source.length()-1)
			file_source.erase(file_source.find_last_of('/'), 1);

		int start = file_source.find_last_of('/');
		int cuttoff = file_source.find('_');
		string prefix = file_source.substr(start+1, cuttoff-(start+1));
		if(dir_genomeFasta.size() > 0){
			// determine the length of each chromosome
			Map_Chr_Len(dir_genomeFasta, Chr_Len_Map, "build/chrlens.log");
		}

		// these are required to proceed
		if (!flag_thres){
			cout << "The -T switch is required." << endl;
			exit(-1);
		}
			
	//if (flag_thres){

		readSize 				= binsizeVec[0];
		extent_cuting_size 	    = binsizeVec[1]; //we use 500 bp according to Yaffe and Tany  
		modelBinSize 			= binsizeVec[2];
		if (fragmentSize == 0)
			fragmentSize = modelBinSize; 

		// Start data processing 
		cout <<"(=:...........Start processing files...........:=)"<<endl;

		// Get the current working directory
		char thisDirectory[1024];
		getcwd(thisDirectory, 1024);
		string cwd = thisDirectory;
		// gather the cutting site location (middle of 6 base pairs), store into map
		int cuttingSiteTotal = Parser_Enzyme_Cutting_Site(file_cuttingSites, cuttingSitesMap, "");
		printf("cuttingSiteTotal == %d\n", cuttingSiteTotal);
		cout <<"<-----Parsed enzyme cutting site map----->"<<endl;
		
		ostringstream thresBuffer;
		thresBuffer.precision(6);
		thresBuffer<<percVec[0];

		string SitesPerChr = cwd + "/" + prefix + "_t" + thresBuffer.str() + "_PerChr.txt";
		// printf("sites per chr %s\n", SitesPerChr.c_str());
		// cout << SitesPerChr << endl;
		ofstream EnzymeFragSizesFile(SitesPerChr.c_str());

		for (int i=1;i<=24;i++)
			EnzymeFragSizesFile << i << "\t" << cuttingSitesMap[i].size() << endl; 

		EnzymeFragSizesFile.close();

		Extend_Cutting_Site(cuttingSitesMap, cuttingSiteExtMap, extent_cuting_size);
		cout <<"<-----Extended cutting site region----->"<<endl;

		// variables used for each chromosome
		vector<Frag_Info_With_Cutting_Site> interactingFrags;
		vector<Frag_Info_With_Cutting_Site> bothEndMappedVec;
		pair<int, int> end1MapInfo, end2MapInfo;
		map< int, map< int, vector<Frag_Info_With_Cutting_Site> > > cuttingSiteFragsMap;

		double aveSiteFrag;
		int totalHybridFrag;

		
		// fprintf(stderr, "<-----Found %d files!----->\n", filecount);
		// sort the list! (just in case)
		// chrSort(filelist, filecount);

		string tmp = "";

		string poissFile;

		/*****************Multi-Processing Setup*********************/
		// pipes that allow communication to a parent process
		int sumpipe[2], cbundle[2], bootstrap_pipe[2];

		// keep track of children being created
		pid_t pid;

		// create the pipes
		pipe(sumpipe);
		pipe(cbundle);
		pipe(bootstrap_pipe);


		// used for limiting number of processes
		// 1) make n processes to start
		// 2) should any of them end, and we have
		// 	more files, make more processes
		// 3) wait on all children
		//
		// create the children to read the files
	
		char *dir = (char *)file_source.c_str();
		int i, filecount;
		char buf[NAME_MAX];

		int nextFile = 0, procNum = 0;
		char **filelist;
		// make 2 main procesess
		pid = fork();

		// create one process to control spawning processes
		// children fall out and begin their work
		if (pid == 0){
			// gather the files in the directory
			filelist = getfilelist(dir, &filecount);
			cout << "<-----Found " << filecount - 2 << " files----->" << endl;

			// no one is reading from these pipes
			close(bootstrap_pipe[0]);
			close(cbundle[0]);
			close(sumpipe[0]);

			// no process limit
			if (maxProc == 0){
				for (i = 0; i < filecount; i++){
					// used for child numbering
					nextFile = i;
					pid = fork();

					// children leave the forking loop
					if (pid == 0)
						break;
				}
				
				// main process waits on these
				if (pid != 0){
					// close write ends
					close(bootstrap_pipe[1]);
					close(cbundle[1]);
					close(sumpipe[1]);
					freeFileList(filelist, filecount);
					while (waitpid(-1, NULL, 0)){
						if (errno == ECHILD)
							break;
					}
					// cout << "Child spawner finished" << endl;
		
					exit(0);
				}
			// limited processes
			}else{
				// while we have files left
				while (nextFile < filecount){
					pid = fork();
					if (pid == 0)
						break;
					// used for files, and child numbering
					nextFile++;
					procNum++;

					// hit max
					if (procNum == maxProc){
						// wait on 1 to finish
						waitpid(-1, NULL, 0);
						procNum--;
					}
				}

				if (pid != 0){
					// close pipes
					close(bootstrap_pipe[1]);
					close(cbundle[1]);
					close(sumpipe[1]);
					freeFileList(filelist, filecount);

					// finish waiting
					while (waitpid(-1, NULL, 0)){
						if (errno == ECHILD)
							break;
					}
					// cout << "Child spawner finished" << endl;
					exit(0);
				}
			}
		}	// end of the child spawning parent
		// start of other main processes code
		else{
			/*****************************************************************************************************************/
			/* Now begins the code for the parent, to collect the information delivered by the children */


			//parent code, wait read from it's pipe and then
			// wait on it's children to finish
			// then finally, write it's contents to the new file

			// how many elements in the final combined file, counted as
			// this process reads from it's combining pipe

			// freeFileList(filelist, filecount);
			unsigned int totalElements = 0;
			map <unsigned int, unsigned int> final_res_count;

			// ofstream tempOut("tempOut");

			// close all the write ends
			close(sumpipe[1]);
			close(cbundle[1]);
			close(bootstrap_pipe[1]);

			unsigned int values[2];
			int bytesread;

			/*****************************************************************************************************************/
			/*****************************************************************************************************************/
			/* main process will become 2 processes, one for reading summation of __PoisMix.txt and writing it's own combined
				file, the other will be reading in the data structures written from the children, and proceeding to finish
				the algorithm. This is done to accomplish parallelism. */
			/*****************************************************************************************************************/
			/*****************************************************************************************************************/
			pid = fork();
				
			/*****************************************************************************************************************/
			/* 1 of 2 main processes, reading from the __PoisMix pipe */
			if (pid == 0){
				// not using these
				close(cbundle[0]);
				close(bootstrap_pipe[0]);

				printf("<-----Reading sum pipe----->\n");
				// first child reads in the sum pipe
				// read in everything from the sum pipe
				while ((bytesread = r_read(sumpipe[0], (char*)values, sizeof(values))) > 0){
					// copy over these read values, add them to the final map
					// is there anything here yet?
					// tempOut << values[0] << " " << values[1] << endl;

					if (final_res_count.count(values[0]) == 0)
						final_res_count[values[0]] = values[1];
					else	// there already exists a value for this, add it to the value
						final_res_count[values[0]] += values[1];
				}

				if (bytesread == -1){
					fprintf(stderr, "Error: there was an error for the main process reading from the poisMix sum pipe.\n");
					exit(-1);
				}

				close(sumpipe[0]);

				// filtration of removing the counts that are: 33 33, 234 234, as in there was only 1
				printf("<-----Performing filtration----->\n");

				// perform filtration
				map<unsigned int, unsigned int> ::const_iterator mapit_s = final_res_count.begin();
				map<unsigned int, unsigned int> ::const_iterator mapit_e = final_res_count.end();

				if (!final_res_count.empty()) mapit_s++;
				unsigned int temp;
				while (mapit_s != mapit_e){
					if (mapit_s->first == mapit_s->second){
						temp = mapit_s->first;
						mapit_s++;
						final_res_count.erase(temp);
					}else
					mapit_s++;
				}

				// write the final distribution to file
				poissFile = cwd + "/" + prefix + "_t" + thresBuffer.str() + "_PoisMix.txt";
				printf("<-----Writing to %s_PoisMix.txt----->\n", prefix.c_str());
				Writer_PoissMix(final_res_count, poissFile);

				final_res_count.clear();
				cout << "<-----Main Process 1 finished writing distrubitions----->" << endl;
				exit(0);
			}

			
			// not reading from this pipe
			close(sumpipe[0]);

			// are we writing to a file for saving memory?
			int largeFileFd;
			unsigned int size;

			// used for bootstrapping
			largeFileFd = open("largeTempFile", O_RDWR | O_TRUNC | O_CREAT, 0666);

			// create last 2 parents
			pid = fork();

			/*****************************************************************************************************************/
			/* process 2 of 3, collecting information from children into 1 file. Need to use pipes for atomic writes.*/
			if (pid == 0){

				// read the sizes of the vectors here for bootstrapping
				cout << "<-----Reading Vector Sizes for Bootstrapping----->" << endl;
				while ((bytesread = r_read(bootstrap_pipe[0], (char*)&size, sizeof(unsigned int))) > 0){
					// cout << "Writing to large File!" << endl;
					r_write(largeFileFd, (char*)&size, sizeof(unsigned int));
				}

				cout << "<-----Finished Vector Sizes for Bootstrapping----->" << endl;
				close(bootstrap_pipe[0]);
				close(largeFileFd);
				
				exit(0);
			}

			
			/*****************************************************************************************************************/
			/* process 3 of 3, will be re-constructing the data structures created by the children here, to proceed with
				finishing the algorithm. Note that the cBundle structure is used only for reading data from the child processes */

			
			// this will hold a file descriptor for each chr/length entry
			map< int, map < int, string > > FIWCS_files;
			int filePostfix = 0;

			// perhaps this can be set by the user
			char *fragprefix = "fiwcs";

			// read to reconstruct data structure
			// fprintf(stderr, "<-----Combining Data from Child Processes----->\n");
			
			
			// this is where we get the largest amount of information
			cout << "<-----Combining Data from Child Processes----->" << endl;
			cBundle indata;
			unsigned int combinedHybridFrags = 0;
			char file[1024];

			// holds the currently open file
			int openfd;
			while ((bytesread = r_read(cbundle[0], (char *)&indata, sizeof(cBundle))) > 0){
				combinedHybridFrags++;
				
				// add the F_I_W_C_S object to the main processes vector list
				if (!saveMemory){
					Frag_Info_With_Cutting_Site temp;
				
					// copy contents from cBundle data structure to F_I_W_C_S object
					cbundleCopy(temp, indata);

					cuttingSiteFragsMap[indata.chr][indata.length].push_back(temp);
				}
				else{		// perform file i/o to save memory

					// create a file for each entry that would be in the structure.
					// to not overflow the open file descriptor table, when writing
					// to a file, open, then close

					// cout << "chr: " << indata.chr << " length: " << indata.length << endl;
					// need to open this new file
					if (FIWCS_files.find(indata.chr) == FIWCS_files.end()
					|| FIWCS_files[indata.chr].find(indata.length) == FIWCS_files[indata.chr].end()){
						memset(file, 0, sizeof(file));
						sprintf(file, "fiwcs/%s%d", fragprefix, filePostfix);
						// create c++ string object
						string tmpobj = file;
						filePostfix++;
						FIWCS_files[indata.chr][indata.length] = tmpobj;
						// printf("%s\n", tmpobj.c_str());

						// create/truncate a tmp file
						openfd = open(file, O_WRONLY | O_TRUNC | O_CREAT, 0666);
						// cout << "Created file " << file << endl;
						// cout << "FD == " << FIWCS_files[indata.chr][indata.length] << endl;
					}else{
						openfd = open(file, O_WRONLY);
						lseek(openfd, 0, SEEK_END);
					}

					// cout << "Writing to " << FIWCS_files[indata.chr][indata.length] << endl;
					// write this content to the respective file
					r_write(openfd, (char*)&indata, sizeof(cBundle));

					// close the file
					close(openfd);
				}
			}
		
			// done reading from this pipe
			close(cbundle[0]);
			// fprintf(stderr, "Total Hybrid Frags: %u\n", combinedHybridFrags);

			// calculate the averageFragSite
			aveSiteFrag = combinedHybridFrags * 1.0 / cuttingSiteTotal;


			// need to wait on the other process, to ensure using the entire file
			while (waitpid(-1, NULL, 0)){
				if (errno == ECHILD)
					break;
			}

			cout << "<-----Child Processes Finished----->" << endl;
			// now the file is complete, proceed
			

			// uses too much memory here!
			// perhaps save off cuttingSiteFragsMap to a file?
			/**************DUMP**********************************/
			/*
			FragDumper *f = new FragDumper();
 
			cout<<"<-----Dumping large memory chunk----->"<<endl;
			f->dump(cuttingSiteFragsMap);
			cout<<"<-----Dumped large memory chunk----->"<<endl;
			cuttingSiteFragsMap.clear();*/

			// output file for bootstrapping
			string bootstrapFile= cwd + "/" + prefix + "_t" + thresBuffer.str() + "_BootStrapping.txt";
							ofstream oboot(bootstrapFile.c_str());

			bootstrapFile= cwd + "/" + prefix + "_t" + thresBuffer.str() + "_maxIteration.txt";
							ofstream iter(bootstrapFile.c_str());
			// int poissfd = open(poissFile.c_str(), O_RDONLY);
			// results fed in
			
			cout << "<-----Beginning Bootstrapping----->" << endl;
			// PoissMix poissOut = Get_Prox_Ligation_Frags(cuttingSiteFragsMap, poisVec, pSampleSize, iterations, oboot, iter);
			PoissMix poissOut = Get_Prox_Ligation_Frags_File(largeFileFd, poisVec, pSampleSize, iterations, oboot, iter);
			oboot.close();
			iter.close();

			cout << "<-----Proximate ligation events extracted from mixture model----->" << endl;
			
			/*
			cout<<"<-----Restoring large memory chunk----->"<<endl;
			f->restore(cuttingSiteFragsMap);
			cout<<"<-----Restored large memory chunk----->"<<endl;
			*/
			
			// string cuttFreqFile        = cwd + "/"+ file_source.c_str() + "_CuttedFragMap";
		        /*	
			string cuttFreqFile        = cwd + "/"+ "__CuttedFragMap";
			string fragFreqFile        = cwd + "/"+ "__FreqMap";
			string enzymeDigsetFile    = cwd + "/"+ "__DegestedMap";

			Writer_Frag_Interaction_Map(proximate_ligations, cuttFreqFile);
			Writer_Frag_Freq_Map(proximate_ligations, fragFreqFile);
			Writer_Digested_Frags(proximate_ligations,enzymeDigsetFile);
			cout<<"Frag interaction frequency map done."<<endl;
			*/
			
			map< int, map<int, int> > digestedFragMap;
			Site_Vec2_Frag_Map(cuttingSitesMap, digestedFragMap);
			cout<<"<-----Transform site vector to digested fragments map----->"<<endl;

			// map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionFreqMap;
			map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > proximate_ligations;

			cout << "<-----Starting Frequency Generation----->" << endl;
			// uses a file that was written to earlier
			// note the file one, passes Chr_Len_Map, since function definition is in another file
			if (!saveMemory)
				//Pair_Frags_Interaction_Freq(cuttingSiteFragsMap, proximate_ligations);
				Pair_Frags_Interaction_Freq_With_Both_Side(cuttingSiteFragsMap, digestedFragMap, proximate_ligations);
			else
				Pair_Frags_Interaction_Freq_File(FIWCS_files, proximate_ligations, Chr_Len_Map);

			cout << "<-----Finished Frequency Generation----->" << endl;

			//construct background model
			map<int, int> ranDisMap;
			int intraChrFrag;
			if (!saveMemory)
				intraChrFrag = Gen_Ran_Dis(cuttingSiteFragsMap, cuttingSitesMap, ranDisMap, modelBinSize);
			else
				intraChrFrag = Gen_Ran_Dis_File(FIWCS_files, cuttingSitesMap, ranDisMap, modelBinSize);

			cuttingSitesMap.clear();

			string rand_dist = cwd + "/" + prefix + "_t" + thresBuffer.str() + "_randomDis.txt";
			Writer_Dis_Distribution(ranDisMap, rand_dist);
			cout<<"<-----Constructed random distribution----->"<<endl;

			map<int, double> ranDisMap_prob;
			Freq_To_Prob(ranDisMap, intraChrFrag, ranDisMap_prob);
                        string rand_dist_prob = cwd + "/" + prefix + "_t" + thresBuffer.str() + "_randomDisProb.txt";
			Writer_Dis_Distribution_Prob(ranDisMap_prob, rand_dist_prob);
			cout<<"<-----Converted to distribution prob----->"<<endl;

			// uses too much memory here!
			// perhaps save off cuttingSiteFragsMap to a file?
			/**************DUMP**********************************/
			/*
			cout<<"<-----Dumping large memory chunk----->"<<endl;
			f->dump(cuttingSiteFragsMap);
			cout<<"<-----Dumped large memory chunk----->"<<endl;
			cuttingSiteFragsMap.clear();
			*/

			cout<<"<-----Begining Get_Interacting_Sites----->"<<endl;
			vector<Pair_Frag_2D_Peak_Info> peakMap_2D;
                        vector<Pair_Frag_2D_Peak_Info> peakMap_2D_random;
			Get_Interacting_Sites(proximate_ligations, digestedFragMap, percVec[0], peakMap_2D);
			Get_Interacting_Sites_random(proximate_ligations, peakMap_2D_random);
			proximate_ligations.clear();
			digestedFragMap.clear();

			cout<<"<-----Created significant interaction peaks----->"<<endl;

			/*
			cout<<"<-----Restoring large memory chunk----->"<<endl;
			f->restore(cuttingSiteFragsMap);
			cout<<"<-----Restored large memory chunk----->"<<endl;
			*/

			// restore here?
			//construct background model
			// map<int, int> ranDisMap;
			// int intraChrFrag = Gen_Ran_Dis(cuttingSiteFragsMap, cuttingSitesMap, ranDisMap, modelBinSize);

			// string rand_dist = cwd + "/" + prefix + "_randomDis";
			//Writer_Dis_Distribution(ranDisMap, rand_dist);
			//cout<<"<-----Constructed random distribution----->"<<endl;
			

			//map<int, double> ranDisMap_prob;
			//Freq_To_Prob(ranDisMap, intraChrFrag, ranDisMap_prob);
			//cout<<"<-----Converted to distribution prob----->"<<endl;

			double minProb = (ranDisMap_prob.rbegin())->second;

			vector<Pair_Frag_2D_Peak_Info>::iterator peakIt_2D = peakMap_2D.begin();
			vector<Pair_Frag_2D_Peak_Info>::iterator peakIt_2D_end = peakMap_2D.end();
			// call this function for every element
			// printf("AverageSiteFrag == %d\n", aveSiteFrag);
			for ( ; peakIt_2D != peakIt_2D_end; ++peakIt_2D ){
				if (!saveMemory)
					peakIt_2D->getSumLogProb(ranDisMap_prob, cuttingSiteFragsMap, minProb, aveSiteFrag, modelBinSize);
				else
					peakIt_2D->getSumLogProbFile(ranDisMap_prob, FIWCS_files, minProb, aveSiteFrag, modelBinSize);
			}
			cout<<"<-----Calculated sum log probability----->"<<endl;

			int num_peak = Find_Peak_Summit(peakMap_2D);
			cout<<"<-----Acquired peak summit----->"<<endl;

			string file_2Dpeak = cwd + "/" + prefix + "_t" + thresBuffer.str() + "_peak.txt";
			//Writer_2D_Peak_with_FDR_Random(peakMap_2D, peakMap_2D_random, percVec[0], file_2Dpeak);
			Writer_2D_Peak_with_FDR_Random_And_Length_Limit(peakMap_2D, peakMap_2D_random, percVec[0], file_2Dpeak, fragmentSize);
			peakMap_2D.clear();
                        peakMap_2D_random.clear();

			cuttingSiteExtMap.clear();
			if (!saveMemory)
				cuttingSiteFragsMap.clear();
			//else
			// delete the directory containing all the temporary files
			char cmd[1024];
			memset(cmd, 0, sizeof(cmd));
			sprintf(cmd, "rm -r %s", fragprefix);
			//printf("%s\n", cmd);
			// system(cmd);

			cout << "<-----Algorithm Finished Successfully----->" << endl;
			remove("largeTempFile");
			exit(0);

		} // end of the last parent code

		// child code starts here
		
		// child process reading the files, performing filtration, and sending
		// information/data out to the main process to collect
		/* maxProc child processes will be reading the 25 chr files here */
		// save an id for this process
		int child_id = nextFile+1;

		// close read ends of pipes, don't need
		ofstream outfile(outputfile.c_str());

		memset(buf, 0, sizeof(buf));
		strcpy(buf, dir);

		// ensure '/' at the end of the directory string
		addDirEnd(buf);

		// add filename to the directory string
		strcat(buf, filelist[i]);

		// free the rest of the files, we only use 1
		freeFileList(filelist, filecount);

		// make a string object instead of a c_str
		tmp = string(buf);

		// read the file, add to the vector
		// used to be Read_Interacting_Region!
		Make_Restriction_Vector(tmp, bothEndMappedVec);


		// use algo sort, with pre-defined function Comp_BothEndsMappedFrag_End1 */
		sort(bothEndMappedVec.begin(), bothEndMappedVec.end(), Comp_BothEndsMappedFrag_End1);
		// cout<<"Child " << child_id << ": sorted end 1"<<endl;

		
		// perform bin size and duplication filtering
		Frag_Info_With_Cutting_Site prevFrag = bothEndMappedVec[0];
		try{
			if((prevFrag.end1_chr != prevFrag.end2_chr) || (prevFrag.end1_pos - prevFrag.end2_pos > modelBinSize) || (prevFrag.end2_pos - prevFrag.end1_pos > modelBinSize))
				interactingFrags.push_back(prevFrag);
		} catch (bad_alloc& ba){
			// cout << "Child " << child_id << ": bad_alloc caught: " << ba.what() << endl;
			exit(-1);
		}

		vector<Frag_Info_With_Cutting_Site>::iterator bothEndMappedVecIt = bothEndMappedVec.begin(); 
		vector<Frag_Info_With_Cutting_Site>::iterator bothEndMappedVecIt_end= bothEndMappedVec.end();

		for ( ++bothEndMappedVecIt; bothEndMappedVecIt != bothEndMappedVecIt_end; ++bothEndMappedVecIt ){
			if ( ( (bothEndMappedVecIt->end1_chr != bothEndMappedVecIt->end2_chr) || (bothEndMappedVecIt->end1_pos - bothEndMappedVecIt->end2_pos > modelBinSize) || (bothEndMappedVecIt->end2_pos - bothEndMappedVecIt->end1_pos > modelBinSize) ) 
					&& 
				(bothEndMappedVecIt->end1_chr != prevFrag.end1_chr || bothEndMappedVecIt->end1_pos != prevFrag.end1_pos || bothEndMappedVecIt->end1_strand != prevFrag.end1_strand || bothEndMappedVecIt->end2_chr != prevFrag.end2_chr || bothEndMappedVecIt->end2_pos != prevFrag.end2_pos || bothEndMappedVecIt->end2_strand != prevFrag.end2_strand) ){
				try{
					// note that here we could write this to a file should we run low on memory
					interactingFrags.push_back(*bothEndMappedVecIt);
				} catch (bad_alloc& ba){
					// cout << "Child " << child_id << ": bad_alloc caught: " << ba.what() << endl;
					struct sysinfo info;
					sysinfo(&info);
					cout << "Memory available: " << info.freeram << endl;
					exit(-1);
				}
			}
			prevFrag = *bothEndMappedVecIt;
		}

		bothEndMappedVec.clear();
		// cout << "Child " << child_id << ": Filtered for bin size and dupicates" << endl;

		// this string variable is re-used for all ofstream files
		string file;
		end1MapInfo = Map_Hybrid_Frag2_CuttingSite_End1(interactingFrags, cuttingSiteExtMap, readSize/2.0);
		// cout<<"Child " << child_id <<": mapped end 1 to cutting site"<<endl;

		// use algo sort, with pre-defined function Comp_BothEndsMappedFrag_End2 */
		sort(interactingFrags.begin(), interactingFrags.end(), Comp_BothEndsMappedFrag_End2);
		// cout<<"Child " << child_id <<": sort end 2"<<endl;

		end2MapInfo = Map_Hybrid_Frag2_CuttingSite_End2(interactingFrags, cuttingSiteExtMap, readSize/2.0);
		// cout<<"Child " << child_id <<": mapped end 2 to cutting site"<<endl;

		// don't need this anymore, free memory
		cuttingSiteExtMap.clear();
		

		// count total hybrid fragments, as well as printing to outfile ofstream
		totalHybridFrag = Get_Cutting_Site_Frags(cuttingSitesMap, interactingFrags, cuttingSiteFragsMap, i+1, outfile);
		
		// cout<<"Child " << child_id <<": construct cutting site fragment map"<<endl;
		
		// clear out the information of this chromosome
		interactingFrags.clear();

		// fprintf(stderr, "Finished chr%d\n", i+1);

		// NEW Add a function to _PoisMix.txt _FragSize.txt _LigateProb.txt
		// string forpoismix = cwd + "/" + file_source.substr(0, file_source.size()-4) + "_PoisMix.txt";
		// string forpoismix = cwd + "/" + "__PoisMix.txt";
		map<unsigned int, unsigned int> count_res_frags, count_prob_ligation;
		// file = cwd + "/" + "fragsizes/" + prefix + to_string((long long int)i+1) + "__FragSize.txt";
		file = "";
		// this HAS to be called before Writer_PoissMix and Writer_LigateProb, for NOW
		// also needed to calculate the number of restriction fragments needed for summing up PoissMix

		Writer_FragSize(cuttingSiteFragsMap, count_res_frags, count_prob_ligation, file);
		// Writer_PoissMix(count_res_frags, file);
		// Writer_LigateProb(count_prob_ligation, cwd + "/" + "__LigateProb.txt");
	
		/**************************************************************************************************/
		/* Here the children will write their values for __PoisMix.txt to the main process to be added up */
		map<unsigned int, unsigned int> ::const_iterator mapit_s = count_res_frags.begin();
		map<unsigned int, unsigned int> ::const_iterator mapit_e = count_res_frags.end();
	
		// write to the PoissMix sum pipe (for __PoisMix.txt)
		unsigned int temp[2];

		for(; mapit_s != mapit_e; ++mapit_s){
			temp[0] = mapit_s->first;
			temp[1] = mapit_s->second;

			// write to the pipe
			ssize_t byteswritten;
			byteswritten = r_write(sumpipe[1], temp, sizeof(temp));
		}
		// printf("Child %d: finished sending distributions\n", child_id);

		// we're done writing to sumpipe, close the line
		close(sumpipe[1]);
		/*****************************************************************************************************************/
		/* children will send their fragment data structure information to the main process to be combined */

					
		// used for re-constructing this datastructure in main process
		unsigned int hybridFrags = 0;
		cBundle newdata;

		map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap.begin(); 
		map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end = cuttingSiteFragsMap.end();



		// send the vector sizes
		for(;chrIt != chrIt_end;++chrIt){
			map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It 		= chrIt->second.begin(); 
			map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It_end 	= chrIt->second.end();

			
			for(;CS_It != CS_It_end; ++CS_It){
				unsigned int size = CS_It->second.size();
				if (size != 0){
					// send the size of this vector, needed for bootstrapping later
					r_write(bootstrap_pipe[1], &size, sizeof(unsigned int));
				}
			}
		}
		//cout << "Child " << child_id << ": finished sending the vector sizes" << endl;

		// now the data structures
		chrIt = cuttingSiteFragsMap.begin();
		for(;chrIt != chrIt_end;++chrIt){
			newdata.chr = chrIt->first;
			map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It 		= chrIt->second.begin(); 
			map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It_end 	= chrIt->second.end();

			for(;CS_It != CS_It_end; ++CS_It){
				newdata.length = CS_It->first;
				if (CS_It->second.size() != 0){
						
					// send this vector element, plus two indices
					for (i = 0; i < CS_It->second.size(); i++){
						fragCopy(newdata, CS_It->second[i]);
						hybridFrags++;
						// write to the pipe
						ssize_t byteswritten;
						byteswritten = r_write(cbundle[1], &newdata, sizeof(cBundle));
					}
				}
			}
		}
		close(cbundle[1]);
		close(bootstrap_pipe[1]);
		// free memory used here!
		cuttingSiteFragsMap.clear();

		//printf("Child %d: sent %u hybrid frags\n", child_id, hybridFrags);

		outfile.close();
		// close the write end, we're done
		// printf("Child %d: finished\n", child_id);
		exit(0);
		// end of children code
	}
	//}
};
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//leave the read size alone, it will be taken cared of during mapping all the reads to the cutting site, in which reads should be shitfed half the read size
int HiC_Peak::Extend_Cutting_Site(map< int, vector<int> >& cuttingSiteMap_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, int& cuttingSiteExtent_local) 
{
	map< int, vector<int> >::const_iterator chrIt = cuttingSiteMap_local.begin(), chrIt_end = cuttingSiteMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		vector<int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();

		//if true strand, site should be second, if false then should be first
		//first site
		if (*siteIt <= cuttingSiteExtent_local)
		{
			cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(1, *siteIt));
		}
		else
		{
			cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair((*siteIt)-cuttingSiteExtent_local, *siteIt));
		}
		int preSite = *siteIt;
		++siteIt;

		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			if (*siteIt-preSite > cuttingSiteExtent_local)
			{
				cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(preSite, preSite+cuttingSiteExtent_local));
				cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(*siteIt-cuttingSiteExtent_local, *siteIt));
			}
			else
			{
				cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(preSite, *siteIt-1));
				cuttingSiteExtMap_local[chrIt->first][true].push_back(make_pair(preSite+1, *siteIt));
			}

			preSite = *siteIt;
		}

		//last site
		--siteIt;
		cuttingSiteExtMap_local[chrIt->first][false].push_back(make_pair(*siteIt, *siteIt+cuttingSiteExtent_local)); //it seems that it would not be a problem for the last region to exceed the chromosomal len
	}

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// reads the 6-column format, after filtering
void HiC_Peak::Make_Restriction_Vector(const string& fileToread, vector<Frag_Info_With_Cutting_Site>& interactingRegionMap_local)
{
	vector<Frag_Info_With_Cutting_Site>* pinteractingRegionMap_local;
	pinteractingRegionMap_local= &interactingRegionMap_local;

	// input data format chr1 pos1 strand1(1,0) chr2 pos2 strand2(1,0)
	ifstream inputFile(fileToread.c_str());
	if (!inputFile)
	{
		cout <<"\n"<< "Error opening " << fileToread << "." << endl;
		exit(1);
	}

	string line = "";
	int chr1    = 0;
	int pos1    = 0;
	int strand1 = 0;
	int chr2    = 0;
	int pos2    = 0;
	int strand2 = 0;
	Frag_Info_With_Cutting_Site fragInfo;

	while(getline(inputFile,line)){

		istringstream ss(line);
		ss >> chr1 >> pos1 >> strand1 >> chr2 >> pos2 >> strand2;

		// Get read 1 info
		fragInfo.end1_chr = chr1;
		fragInfo.end1_pos = pos1;

		if(strand1 == 1)
		fragInfo.end1_strand = true;
		else fragInfo.end1_strand = false;

		// Get read 2 info
		fragInfo.end2_chr = chr2;
		fragInfo.end2_pos = pos2;

		if(strand2 == 1)
		fragInfo.end2_strand = true;
		else fragInfo.end2_strand = false;

		pinteractingRegionMap_local->push_back(fragInfo);
	}
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//map hybrid reads to cutting sites
pair <int, int> HiC_Peak::Map_Hybrid_Frag2_CuttingSite_End1(vector<Frag_Info_With_Cutting_Site>& hybridFragVec_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, const int halfReadSize_local)
{
	int regionStart = 0;
	int regionEnd = 0;
	unsigned int totalFragsInRegions = 0;
	int totalRegionNum = 0;
	double totalAveDensity = 0.0;

	map< int, map< bool, vector< pair<int, int> > > >::iterator chrIt_regionMap = cuttingSiteExtMap_local.begin();
	map< int, map< bool, vector< pair<int, int> > > >::iterator chrIt_regionMap_end = cuttingSiteExtMap_local.end();

	vector<Frag_Info_With_Cutting_Site>::iterator fragIt_begin = hybridFragVec_local.begin();
	vector<Frag_Info_With_Cutting_Site>::iterator fragIt = hybridFragVec_local.begin();
	vector<Frag_Info_With_Cutting_Site>::iterator fragIt_chrStart = hybridFragVec_local.begin();
	vector<Frag_Info_With_Cutting_Site>::iterator fragIt_end = hybridFragVec_local.end();

	//true strand first, then false strand
	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end1_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end1_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end1_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[true].begin();
			vector< pair<int, int> >::const_iterator regionIt_end = chrIt_regionMap->second[true].end();
			while (regionIt != regionIt_end && fragIt->end1_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) < regionStart )
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end1_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end1_strand == true) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end1_cuttingSite = regionEnd; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}

	//false strand
	chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end1_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end1_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end1_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[false].begin(), regionIt_end = chrIt_regionMap->second[false].end();
			while (regionIt != regionIt_end && fragIt->end1_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end1_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end1_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end1_chr == chrIt_regionMap->first && (fragIt->end1_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end1_strand == false) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end1_cuttingSite = regionStart; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}

	return make_pair(totalRegionNum, totalFragsInRegions);
	//return totalFragsInPeaks;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//map hybrid reads to cutting sites
pair <int, int> HiC_Peak::Map_Hybrid_Frag2_CuttingSite_End2(vector<Frag_Info_With_Cutting_Site>& hybridFragVec_local, map< int, map< bool, vector< pair<int, int> > > >& cuttingSiteExtMap_local, const int halfReadSize_local)
{
	int regionStart = 0;
	int regionEnd = 0;
	int totalFragsInRegions = 0;
	int totalRegionNum = 0;
	double totalAveDensity = 0.0;

	map< int, map< bool, vector< pair<int, int> > > >::iterator chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	vector<Frag_Info_With_Cutting_Site>::iterator fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	//true strand first, then false strand
	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end2_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end2_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end2_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[true].begin(), regionIt_end = chrIt_regionMap->second[true].end();
			while (regionIt != regionIt_end && fragIt->end2_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end2_chr != chrIt_regionMap->first )
				{
					break;
				}


				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end2_strand == true) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end2_cuttingSite = regionEnd; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}



	//false strand
	chrIt_regionMap = cuttingSiteExtMap_local.begin(), chrIt_regionMap_end = cuttingSiteExtMap_local.end();
	fragIt_begin = hybridFragVec_local.begin(), fragIt = hybridFragVec_local.begin(), fragIt_chrStart = hybridFragVec_local.begin(), fragIt_end = hybridFragVec_local.end();

	for ( ; chrIt_regionMap != chrIt_regionMap_end; ++chrIt_regionMap )
	{
		//cout <<"\nentered"<<i<<endl;
		while (fragIt->end2_chr < chrIt_regionMap->first && fragIt != fragIt_end)
		{
			++fragIt;
		}

		if (fragIt == fragIt_end)
		{
			break;
		}
		else if (fragIt->end2_chr > chrIt_regionMap->first)
		{
			continue;
		}
		else if (fragIt->end2_chr == chrIt_regionMap->first)
		{
			fragIt_chrStart = fragIt;

			vector< pair<int, int> >::const_iterator regionIt = chrIt_regionMap->second[false].begin(), regionIt_end = chrIt_regionMap->second[false].end();
			while (regionIt != regionIt_end && fragIt->end2_chr == chrIt_regionMap->first)
			{
				//cout <<"entered while"<<endl;
				regionStart = regionIt->first;
				regionEnd = regionIt->second;

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) < regionStart ) //--fragIt is to counter the effect of the following ++fragIt, when it would not go into the following while loop
				{
					++fragIt;
					//cout <<"ET1W1\t";
				}

				bool moveBackFlag = false;
				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
					moveBackFlag = true;
				}

				while ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && fragIt != fragIt_chrStart && (fragIt->end2_pos + halfReadSize_local) >= regionStart )
				{
					--fragIt;
				}
				//cout <<"OT1W2\t";

				if ( fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && moveBackFlag )
				{
					++fragIt;
				}
				else if ( fragIt == fragIt_end || fragIt->end2_chr != chrIt_regionMap->first )
				{
					break;
				}

				//cout <<"OT1W2\t";
				int fragsCounter = 0;
				//int startingFrag = 0;

				for ( ; fragIt != fragIt_end && fragIt->end2_chr == chrIt_regionMap->first && (fragIt->end2_pos + halfReadSize_local) <= regionEnd ; ++fragIt )
				//for (startingRead_forward = *fragIt_forward; *fragIt_forward < fragRangeEnd_forward; ++fragIt_forward)
				{
					//cout <<"ET1F\t";
					if (fragIt->end2_strand == false) //this made the totalRegionNum meaningless. Even the for loop is entered, there might be no frags in this region since the frag might only on the other strand
					{
						fragIt->end2_cuttingSite = regionStart; //true strand = regionEnd and false strand = regionStart
						fragsCounter += 1;
					}
				}

				++totalRegionNum;
				totalFragsInRegions += fragsCounter;

				++regionIt;
			}
		}
	}

	return make_pair(totalRegionNum, totalFragsInRegions);
	//return totalFragsInPeaks;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//1) filter out 0 or 1 end mapped to cutting site frags, and keep both end mapped. 
//2) link frag to the cutting site map 
//3) remove two ends mapped to the same cutting site or adjacent cutting site if the two ends following the self loop characteristic (end mapped to the precedent cutting site is -, and the other end is +)
int HiC_Peak::Get_Cutting_Site_Frags(const map< int, vector<int> >& cuttingSitesMap_local, const vector<Frag_Info_With_Cutting_Site>& interactingFrags_local, map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap_local, int chrNum, ofstream& outfile)
{
	vector<Frag_Info_With_Cutting_Site> emptyFragVec;
	//initialize the cuttingSiteFragMap, add all cutting site to the map, with empty vector as value
	map< int, vector<int> >::const_iterator chrIt = cuttingSitesMap_local.begin(), chrIt_end = cuttingSitesMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		vector<int>::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			cuttingSiteFragsMap_local[chrIt->first][*siteIt] = emptyFragVec;
		}
	}

	unsigned int totalFrag = 0;
	int bothEndMapped = 0;
	unsigned int noBothEndMapped = 0;
	int interChr = 0;
	int intraChr = 0;
	int bothEndMapToOneSite = 0;
	int bothEndMapToDiffSite = 0;
	int intraHybridFrags = 0;
	int selfLoop = 0;

	vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = interactingFrags_local.begin(), fragIt_end = interactingFrags_local.end();
	for ( ; fragIt != fragIt_end; ++fragIt )
	{
		++totalFrag;
		if ( fragIt->end1_cuttingSite != 0 && fragIt->end2_cuttingSite != 0 )
		{
			++bothEndMapped;
			if ( fragIt->end1_chr != fragIt->end2_chr )
			{
				++interChr;
				cuttingSiteFragsMap_local[fragIt->end1_chr][fragIt->end1_cuttingSite].push_back(*fragIt);
			}
			else
			{
				++intraChr;
				//change the parsing frag step, allow each frag has two record in the fragmap, instead of only allow one for two ends dis < 2000
				if ( fragIt->end1_cuttingSite == fragIt->end2_cuttingSite )
				{
					++bothEndMapToOneSite;
				}
				else
				{
					++bothEndMapToDiffSite;
					if ( fragIt->end1_cuttingSite > fragIt->end2_cuttingSite && fragIt->end1_strand == true && fragIt->end2_strand == false ) //self loop smaller coordinates should be on negative strand and larger should be on positive strand
					{
						map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator end2SiteIt = cuttingSiteFragsMap_local[fragIt->end2_chr].find(fragIt->end2_cuttingSite);
						if ( end2SiteIt != cuttingSiteFragsMap_local[fragIt->end2_chr].end() )
						{
							++end2SiteIt;
							if ( end2SiteIt != cuttingSiteFragsMap_local[fragIt->end2_chr].end() )
							{
								if ( end2SiteIt->first != fragIt->end1_cuttingSite )
								{
									if( abs(fragIt->end1_pos - fragIt->end2_pos) < modelBinSize )
										continue;
									++intraHybridFrags;
									cuttingSiteFragsMap_local[fragIt->end1_chr][fragIt->end1_cuttingSite].push_back(*fragIt);
								}
								else
								{
									++selfLoop;
									// cout <<"Self loop found!\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
								}
							}
							else
							{
								cout <<"Error! End2 site is the end, end1 site not found!\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
							}
						}
						else
						{
							cout <<"Error! End2 site not found!\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
						}
					}
					else if ( fragIt->end1_cuttingSite < fragIt->end2_cuttingSite && fragIt->end1_strand == false && fragIt->end2_strand == true )
					{
						map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator end1SiteIt = cuttingSiteFragsMap_local[fragIt->end1_chr].find(fragIt->end1_cuttingSite);
						if ( end1SiteIt != cuttingSiteFragsMap_local[fragIt->end1_chr].end() )
						{
							++end1SiteIt;
							if ( end1SiteIt != cuttingSiteFragsMap_local[fragIt->end1_chr].end() )
							{
								if ( end1SiteIt->first != fragIt->end2_cuttingSite )
								{
									if( abs(fragIt->end1_pos - fragIt->end2_pos) < modelBinSize )
										continue;
									++intraHybridFrags;
									cuttingSiteFragsMap_local[fragIt->end1_chr][fragIt->end1_cuttingSite].push_back(*fragIt);
								}
								else
								{
									++selfLoop;
									// cout <<"Self loop found!\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
								}
							}
							else
							{
								cout <<"Error! End1 site is the end, end2 site not found\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
							}
						}
						else
						{
							cout <<"Error! End1 site not found!\t"<<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<endl;
						}
					}
					else
					{
						++intraHybridFrags;
						cuttingSiteFragsMap_local[fragIt->end1_chr][fragIt->end1_cuttingSite].push_back(*fragIt);
					}
				}
			}
		}
		else
		{
			++noBothEndMapped;
		}
	}

	int totalHybridFrag = interChr + intraHybridFrags;

	outfile<<"Chr"<<chrNum<<"\n\tTotal fragments:\t"<<totalFrag<<"\n\tBoth ends mapped to cutting sites:\t"<<bothEndMapped<<"\n\tNo end or 1 end mapped to cutting site:\t"<<noBothEndMapped<<"\n\tInter chromosomal frags:\t"<<interChr<<"\n\tIntra chromosomal frags:\t"<<intraChr<<"\n\tBoth ends mapped to one cutting site:\t"<<bothEndMapToOneSite<<"\n\tTwo ends mapped to different cutting siteds:\t"<<bothEndMapToDiffSite<<"\n\tIntra chromosomal hybrid frags:\t"<<intraHybridFrags<<"\n\tSelf loop:\t"<<selfLoop<<endl<<endl;
	return totalHybridFrag;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//find interaction frequency between each pair of ENZYME DIGESTED FRAGMENTS, then search for freq that above the threshold (including the neighbouring frag)
//map<pair<chr1,pair<frag1Start, frag1End> >, pair<chr2,pair<frag2STart, frag2End> > >, pair<number of interaction, vector<interaction frag> > >
int HiC_Peak::Pair_Frags_Interaction_Freq(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap_local, 
	                                      map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >& freqMaptoWrite)
{
	freqMaptoWrite.clear();
	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(), chrIt_end = cuttingSiteFragsMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt_begin = chrIt->second.begin(), siteIt_frag1End = chrIt->second.begin(), siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				int frag1_start = 0; //a pair of interacting frag, frag1 is end1, frag2 is end2
				int frag1_end = 0;
				siteIt_frag1End = siteIt;
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
				}
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
			}
		}
	}

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//construct distribution of random interaction genomic distance
int HiC_Peak::Gen_Ran_Dis(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap_local, 
	                      map< int, vector<int> >& cuttingSitesMap_local, 
	                      map<int, int>& ranDisMap_local, 
	                      int binSize_local)
{
	int counter_frag = 0;

	//sumup genomic length
	int totalCuttingSite = 0;

	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin();
	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end = cuttingSiteFragsMap_local.end();
	
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		totalCuttingSite += cuttingSitesMap_local[chrIt->first].size();
	}	

	chrIt = cuttingSiteFragsMap_local.begin();
	for ( ; chrIt != chrIt_end; ++chrIt ){
		int sitesMinusChr = totalCuttingSite - cuttingSitesMap_local[chrIt->first].size();
		//site on each chromosome should have different chance to have a random inter-chromosomal ligation, the smaller the chromosome is,
		// larger the rest of the genome is, thus for each particular site on a smaller chromosome will have a better chance to form a inter chromosomal ligation
		int interChrFrag = 0; 
		map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();
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
			}
		}

		double interChrFragProb = interChrFrag * 1.0 / sitesMinusChr;
		//cout<<chrIt->first<<"\tinter chromosomal frag prob:\t"<<interChrFragProb;
	}

	return counter_frag;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Freq_To_Prob(const map<int, int>& ranDisMap_local, 
	                       const int intraChrFrag_local, 
	                       map<int, double>& ranDisMap_prob_local)
{
	map<int, int>::const_iterator disIt = ranDisMap_local.begin(), disIt_end = ranDisMap_local.end();
	for ( ; disIt != disIt_end; ++disIt )
	{
		ranDisMap_prob_local[disIt->first] = disIt->second * 1.0 / intraChrFrag_local;
	}

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Site_Vec2_Frag_Map(const map< int, vector<int> >& cuttingSitesMap_local, 
	                             map< int, map<int, int> >& digestedFragMap_local)
{
	map< int, vector<int> >::const_iterator chrIt = cuttingSitesMap_local.begin();
	map< int, vector<int> >::const_iterator chrIt_end = cuttingSitesMap_local.end();

	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		int prevSite = 1;
		vector<int>::const_iterator siteIt = chrIt->second.begin();
		vector<int>::const_iterator siteIt_end = chrIt->second.end();

		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			digestedFragMap_local[chrIt->first][prevSite] = *siteIt;
			prevSite = *siteIt;
		}
		digestedFragMap_local[chrIt->first][prevSite] = Chr_Len_Map[chrIt->first];
	}

	return 0;
}

typedef unsigned long long int UINT64;
UINT64 getRandom(UINT64 const& min = 0, UINT64 const& max = 0)
{
    return (((UINT64)(unsigned int)rand() << 32) + (UINT64)(unsigned int)rand()) % (max - min) + min;
}

void HiC_Peak::Get_Interacting_Sites_random(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local, vector<Pair_Frag_2D_Peak_Info>& peakMaptoWrite_2D)
{
	vector<Pair_Frag_2D_Peak_Info> peak_selected_info;
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator fragPairIt = fragInteractionMap_local.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator fragPairIt_end = fragInteractionMap_local.end();
	
	long int it = 0; //number of interatcing sites
	while (fragPairIt != fragPairIt_end) 
        {
	        it++;
		int chr_domain1 = fragPairIt->first.first.first;
		int chr_domain2 = fragPairIt->first.second.first;
	
		pair<int, int> position_frag1 = fragPairIt->first.first.second;
		pair<int, int> position_frag2 = fragPairIt->first.second.second;
	
		int fragPairScore = fragPairIt->second.size();
		int totalScore = fragPairIt->second.size();
	
	        int fragPairNum = 1;//my $binNo = 0;
	
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > peakRegionsCoordinates;
		Pair_Frag_2D_Peak_Info peakInfo;
		vector<Frag_Info_With_Cutting_Site> peakFragInfo;
		peakFragInfo = fragPairIt->second;
	
		peakRegionsCoordinates.push_back(make_pair(make_pair(position_frag1, position_frag2), fragPairScore));
		fragInteractionMap_local.erase(fragPairIt);
	
		peakInfo.chrNo_domain1 = chr_domain1;
		peakInfo.chrNo_domain2 = chr_domain2;
		peakInfo.fragPair_num = fragPairNum;
		peakInfo.total_frag = totalScore;
		peakInfo.frag_info = peakFragInfo;
		peakInfo.region_info = peakRegionsCoordinates;
		peak_selected_info.push_back(peakInfo);
		fragPairIt = fragInteractionMap_local.begin();
	}
	
	long int num_random = it / 10; // Number of random
	vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st = peak_selected_info.begin(), peak_sp = peak_selected_info.end();
        for (long int i = 0; i < num_random; i++)
	{       
		long int num = getRandom(0, it);
                vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st = peak_selected_info.begin();
                //for (long int j = 0; j < num; j++)
                //{
                //        if (peak_st < peak_sp)
                //        {
                //                ++peak_st;
                //        }
                //        else
                //        {
                //                break;
                //        }
                //}
                peak_st += num;
                Pair_Frag_2D_Peak_Info peakInfo;
                peakInfo.chrNo_domain1 	= peak_st->chrNo_domain1;
		peakInfo.chrNo_domain2 	= peak_st->chrNo_domain2;
		peakInfo.fragPair_num 	= peak_st->fragPair_num;
		peakInfo.total_frag 	= peak_st->total_frag;
		peakInfo.frag_info 	= peak_st->frag_info;
		peakInfo.region_info 	= peak_st->region_info;
		peakMaptoWrite_2D.push_back(peakInfo);
																								}
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//2Dpeakmap: Map< pair<pair<domain1 chr, domain2 chr>, vector(coordinates and score of each region in this peak)< pair< pair<domain1 position, domain2 position>,score> >, vector<peak property, p value, score etc> >
// map<
// 	pair<
// 		pair<int, pair<int (domain), int (position)>
// 		>
// 	pair<
// 		pair<int, pair<int(domain), int (position)>
// 		>
// >
void HiC_Peak::Get_Interacting_Sites(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local,
									map< int, map<int, int> >& digestedFragMap_local, 
	                                const double thres, 
	                                vector<Pair_Frag_2D_Peak_Info>& peakMaptoWrite_2D)
{
	// ensure this vector is empty?
	// peakMaptoWrite_2D.clear();

	// count how many total peaks
	// unsigned int countPeaknum 			= 0;
	const double doubleThres 	= 2 * thres;

	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator fragPairIt = fragInteractionMap_local.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator fragPairIt_end = fragInteractionMap_local.end();

	unsigned int it = 0;
	//iteratively search for neighbour region that has a higher score than the threshold
	while(fragPairIt != fragPairIt_end){
		/********************SETUP**************************/
		// here the large structure will be dumped to a file,
		// and read from to save memory


		/*
		fprintf(stderr, "Dumping the structure!\n");
		FragDumper *fragDumper = new FragDumper();

		// dumps this structure to files
		fragDumper->saveThing(fragInteractionMap_local);
		fprintf(stderr, "Dumped the structure!\n");


		// fragInteractionMap_local.clear();

		exit(0);*/
                it++;
		int chr_domain1 = fragPairIt->first.first.first;
		int chr_domain2 = fragPairIt->first.second.first;

		pair<int, int> position_frag1 = fragPairIt->first.first.second;
		pair<int, int> position_frag2 = fragPairIt->first.second.second;

		int fragPairScore = fragPairIt->second.size();
		int totalScore    = fragPairIt->second.size();

		int fragPairNum = 1;//my $binNo = 0;

		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > peakRegionsCoordinates;
		Pair_Frag_2D_Peak_Info peakInfo;
		vector<Frag_Info_With_Cutting_Site> peakFragInfo;
		peakFragInfo = fragPairIt->second;

		peakRegionsCoordinates.push_back( make_pair(make_pair(position_frag1,position_frag2), fragPairScore) );
		fragInteractionMap_local.erase(fragPairIt);

		if (fragPairScore >= thres)
		{
			// cout << "Get_Interacting_Sites: running search_neighbor_frags" << endl;
			Search_Neighbour_Frags(fragInteractionMap_local, 
								   digestedFragMap_local, 
								   chr_domain1, 
								   chr_domain2, 
								   position_frag1, 
								   position_frag2, 
								   0, 
								   0, 
								   peakRegionsCoordinates, 
								   thres, 
								   fragPairNum, 
								   totalScore, 
								   peakFragInfo);
		}
        //if regionIt itself or the sum of all scores in the peak is larger than double threshold -> peak
		//if (totalScore >= doubleThres)
		{
			//cout << "Get_Interacting_Sites: Found a peak!" << endl;
			peakInfo.chrNo_domain1 	= chr_domain1;
			peakInfo.chrNo_domain2 	= chr_domain2;
			peakInfo.fragPair_num 	= fragPairNum;
			peakInfo.total_frag 	= totalScore;
			peakInfo.frag_info 		= peakFragInfo;
			peakInfo.region_info 	= peakRegionsCoordinates;

			peakMaptoWrite_2D.push_back(peakInfo);

			// ++countPeaknum;
		}

		fragPairIt = fragInteractionMap_local.begin();
	}
	cout << "Get_Interacting_Sites: Finished!" << endl;
	// return countPeaknum;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Search_Neighbour_Frags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >& freqMaptoSearch_local, 
	                                 map< int, map<int, int> >& digestedFragMap_search, 
	                                 int chr_domain1_local, 
	                                 int chr_domain2_local, 
	                                 pair<int, int> position_domain1_local, 
	                                 pair<int, int> position_domain2_local, 
	                                 int parentSetoff_x, 
	                                 int parentSetoff_y, 
	                                 vector< pair< pair< pair<int, int>, pair<int, int> >, int > >& peakRegionsCoordinates_local, 
	                                 double thres_local, 
	                                 int& pairFragNum_local, 
	                                 int& totalScore_local, 
	                                 vector<Frag_Info_With_Cutting_Site>& peakFragInfo_local)
{
	//cout <<chr_domain1_local<<"\t"<<chr_domain2_local<<"\t"<<position_domain1_local.first<<"\t"<<position_domain2_local.first<<"\t"<<parentSetoff_x<<"\t"<<parentSetoff_y<<"\t"<<thres_local<<"\t"<<pairFragNum_local<<"\t"<<totalScore_local<<endl;
	
	// cout << "Entered search_neighbor!" << endl;
	for (int x=-1; x<=1; ++x)
	{
		//cout<<"enter x."<<endl;
		map<int, int>::const_iterator targetFrag_x = digestedFragMap_search[chr_domain1_local].find(position_domain1_local.first);

		if ( x==-1 && targetFrag_x != digestedFragMap_search[chr_domain1_local].begin() && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() ) //search for the previous one
		{
			//cout<<x<<endl;
			targetFrag_x--;
		}
		else if ( x==1 && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() )
		{
			//cout<<x<<endl;
			targetFrag_x++;
			if (targetFrag_x == digestedFragMap_search[chr_domain1_local].end())
			{
				//cout<<"if continue"<<endl;
				continue;
			}
		}
		//else if ( x==0 && targetFrag_x != digestedFragMap_search[chr_domain1_local].end() ){}
		else
		{
			//cout<<"else continue"<<endl;
			continue;
		}

		for (int y=-1; y<=1; ++y)
		{
			//cout<<"enter y."<<endl;
			map<int, int>::const_iterator targetFrag_y = digestedFragMap_search[chr_domain2_local].find(position_domain2_local.first);

			if ( y==-1 && targetFrag_y != digestedFragMap_search[chr_domain2_local].begin() && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() ) //search for the previous one
			{
				targetFrag_y--;
			}
			else if ( y==1 && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() )
			{
				targetFrag_y++;
				if (targetFrag_y == digestedFragMap_search[chr_domain2_local].end())
				{
					continue;
				}
			}
			//else if ( y==0 && targetFrag_y != digestedFragMap_search[chr_domain2_local].end() ) {}
			else
			{
				continue;
			}

			if ( !( (x==0 && y==0) || ( x==(0-parentSetoff_x) && y==(0-parentSetoff_y) ) ) )
			{
				//cout<<"enter xy."<<endl;
				map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::iterator neighbourPairFrag = freqMaptoSearch_local.find( make_pair( make_pair(chr_domain1_local, *targetFrag_x), make_pair(chr_domain2_local, *targetFrag_y) ) );

				if (neighbourPairFrag != freqMaptoSearch_local.end())
				{
					//cout<<"enter recurrsive."<<endl;
					int neighbourScore = neighbourPairFrag->second.size();
					pair<int, int> neighbourPosition_domain1 = neighbourPairFrag->first.first.second, neighbourPosition_domain2 = neighbourPairFrag->first.second.second;
					vector<Frag_Info_With_Cutting_Site> neighbourPairFragInfo = neighbourPairFrag->second;
					freqMaptoSearch_local.erase(neighbourPairFrag);//delete the element before pass down to recurrsive checking, if not, for example, in 3*3 array, [0,0][0,1][1,1],has value > thres, [0,0]->[0,1]->[1,1]->[0,0]->.....

					if (neighbourScore >= thres_local)
					{
						totalScore_local += neighbourScore;
						++pairFragNum_local;
						peakRegionsCoordinates_local.push_back( make_pair( make_pair( neighbourPosition_domain1, neighbourPosition_domain2 ), neighbourScore ) );
						vector<Frag_Info_With_Cutting_Site>::iterator fragIt = neighbourPairFragInfo.begin(), fragIt_end = neighbourPairFragInfo.end();
						for ( ; fragIt != fragIt_end; ++fragIt)
						{
							peakFragInfo_local.push_back(*fragIt);
						}
						Search_Neighbour_Frags(freqMaptoSearch_local, 
											   digestedFragMap_search, 
											   chr_domain1_local, 
											   chr_domain2_local, 
											   neighbourPosition_domain1, 
											   neighbourPosition_domain2, 
											   x, 
											   y, 
											   peakRegionsCoordinates_local, 
											   thres_local, 
											   pairFragNum_local, 
											   totalScore_local, 
											   peakFragInfo_local);
					}
				}
			}
		}
	}
	// cout << "Finished search_neighbor_frags" << endl;
	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//find the summit of a peak, if there are multiple paired enzyme cut fragments, choose the one with highest score, 
//if there multiple with highest score then: connect them if continuous, make two peak if separated
int HiC_Peak::Find_Peak_Summit(vector<Pair_Frag_2D_Peak_Info>& peakMap_2D_local)
{
	int addPeakNum = 0;
	vector<Pair_Frag_2D_Peak_Info> tempPeakVec;

	vector<Pair_Frag_2D_Peak_Info>::iterator peakIt 	= peakMap_2D_local.begin();
	vector<Pair_Frag_2D_Peak_Info>::iterator peakIt_end = peakMap_2D_local.end();

	for ( ; peakIt != peakIt_end; ++peakIt ){
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > >::const_iterator regionIt 		= peakIt->region_info.begin();
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > >::const_iterator regionIt_end 	= peakIt->region_info.end();

		int maxScore =0;
		for ( ; regionIt != regionIt_end; ++regionIt ){
			if (regionIt->second > maxScore){
				maxScore = regionIt->second;
			}
		}

		vector< pair< pair< pair<int, int>, pair<int, int> >, int > > maxScorePairFrags;
		regionIt = peakIt->region_info.begin();
		for ( ; regionIt != regionIt_end; ++regionIt ){
			if ( regionIt->second == maxScore ){
				maxScorePairFrags.push_back(*regionIt);
			}
		}

		regionIt = maxScorePairFrags.begin();
		regionIt_end = maxScorePairFrags.end();

		peakIt->peakRegion.first.first 		= regionIt->first.first.first;
		peakIt->peakRegion.first.second 	= regionIt->first.first.second;
		peakIt->peakRegion.second.first 	= regionIt->first.second.first;
		peakIt->peakRegion.second.second 	= regionIt->first.second.second;

		++regionIt;
		for ( ; regionIt != regionIt_end; ++regionIt ){
			if ( peakIt->peakRegion.first.first > regionIt->first.first.second || peakIt->peakRegion.first.second < regionIt->first.first.first || peakIt->peakRegion.second.first > regionIt->first.second.second || peakIt->peakRegion.second.second < regionIt->first.second.first )
			{
				if (tempPeakVec.empty())
				{
					Pair_Frag_2D_Peak_Info tempPeakInfo;
					tempPeakInfo.chrNo_domain1 	= peakIt->chrNo_domain1;
					tempPeakInfo.chrNo_domain2 	= peakIt->chrNo_domain2;
					tempPeakInfo.fragPair_num 	= 1;
					tempPeakInfo.total_frag 	= regionIt->second;
					tempPeakInfo.frag_info 		= peakIt->frag_info; //added peak has no frag information
					tempPeakInfo.region_info.push_back(*regionIt);

					tempPeakInfo.peakRegion.first.first 	= regionIt->first.first.first;
					tempPeakInfo.peakRegion.first.second 	= regionIt->first.first.second;
					tempPeakInfo.peakRegion.second.first 	= regionIt->first.second.first;
					tempPeakInfo.peakRegion.second.second 	= regionIt->first.second.second;

					tempPeakInfo.sumLogProb = peakIt->sumLogProb;

					tempPeakVec.push_back(tempPeakInfo);

					++addPeakNum;
				}
				else
				{
					vector<Pair_Frag_2D_Peak_Info>::reverse_iterator lastAddedPeak = tempPeakVec.rbegin();
					if ( lastAddedPeak->peakRegion.first.first > regionIt->first.first.second || lastAddedPeak->peakRegion.first.second < regionIt->first.first.first || lastAddedPeak->peakRegion.second.first > regionIt->first.second.second || lastAddedPeak->peakRegion.second.second < regionIt->first.second.first )
					{
						Pair_Frag_2D_Peak_Info tempPeakInfo;
						tempPeakInfo.chrNo_domain1 	= peakIt->chrNo_domain1;
						tempPeakInfo.chrNo_domain2 	= peakIt->chrNo_domain2;
						tempPeakInfo.fragPair_num 	= 1;
						tempPeakInfo.total_frag 	= regionIt->second;
						tempPeakInfo.frag_info 		= peakIt->frag_info; //added peak has no frag information
						tempPeakInfo.region_info.push_back(*regionIt);

						tempPeakInfo.peakRegion.first.first 	= regionIt->first.first.first;
						tempPeakInfo.peakRegion.first.second 	= regionIt->first.first.second;
						tempPeakInfo.peakRegion.second.first 	= regionIt->first.second.first;
						tempPeakInfo.peakRegion.second.second 	= regionIt->first.second.second;

						tempPeakInfo.sumLogProb = peakIt->sumLogProb;
						//peakInfo.Find_Peak_Summit();

						tempPeakVec.push_back(tempPeakInfo);

						++addPeakNum;
					}
					else
					{
						if (regionIt->first.first.first < lastAddedPeak->peakRegion.first.first)
						{
							lastAddedPeak->peakRegion.first.first = regionIt->first.first.first;
						}
						if (regionIt->first.first.second > lastAddedPeak->peakRegion.first.second)
						{
							lastAddedPeak->peakRegion.first.second = regionIt->first.first.second;
						}
						if (regionIt->first.second.first < lastAddedPeak->peakRegion.second.first)
						{
							lastAddedPeak->peakRegion.second.first = regionIt->first.second.first;
						}
						if (regionIt->first.second.second > lastAddedPeak->peakRegion.second.second)
						{
							lastAddedPeak->peakRegion.second.second = regionIt->first.second.second;
						}
					}

				}

			}
			else
			{
				if (regionIt->first.first.first < peakIt->peakRegion.first.first)
				{
					peakIt->peakRegion.first.first = regionIt->first.first.first;
				}
				if (regionIt->first.first.second > peakIt->peakRegion.first.second)
				{
					peakIt->peakRegion.first.second = regionIt->first.first.second;
				}
				if (regionIt->first.second.first < peakIt->peakRegion.second.first)
				{
					peakIt->peakRegion.second.first = regionIt->first.second.first;
				}
				if (regionIt->first.second.second > peakIt->peakRegion.second.second)
				{
					peakIt->peakRegion.second.second = regionIt->first.second.second;
				}
			}
		}
	}
	peakMap_2D_local.insert(peakMap_2D_local.end(), tempPeakVec.begin(), tempPeakVec.end());

	return addPeakNum;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_2D_Peak(vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D, const string& peakFiletoWrite_2D)
{
	ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
	vector<Pair_Frag_2D_Peak_Info>::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
	for ( ; peakIt != peakIt_end; ++peakIt )
	{
		output2DpeakFile <<">\t"<<peakIt->chrNo_domain1<<"\t"<<peakIt->chrNo_domain2<<"\t"<<peakIt->fragPair_num<<"\t"<<peakIt->total_frag<<"\t"<<(peakIt->peakRegion.first.first)<<"\t"<<(peakIt->peakRegion.first.second)<<"\t"<<peakIt->peakRegion.second.first<<"\t"<<peakIt->peakRegion.second.second<<"\t"<<peakIt->sumLogProb;
		vector< pair< pair< pair<int, int>, pair<int, int> >, int > >::const_iterator regionIt = peakIt->region_info.begin(), regionIt_end = peakIt->region_info.end();
		for ( ; regionIt != regionIt_end; ++regionIt )
		{
			output2DpeakFile <<"\t"<<(regionIt->first.first.first)<<"\t"<<(regionIt->first.first.second)<<"\t"<<(regionIt->first.second.first)<<"\t"<<(regionIt->first.second.second)<<"\t"<<regionIt->second;
		}
		output2DpeakFile<<endl;

		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = peakIt->frag_info.begin(), fragIt_end = peakIt->frag_info.end();
		for ( ; fragIt != fragIt_end; ++fragIt)
		{
			output2DpeakFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<endl;
		}
	}
	output2DpeakFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_2D_Peak_with_FDR(vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D,PoissMix*mix,const string& peakFiletoWrite_2D)
{
        ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
        vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st = peakMaptoRead_2D.begin(), peak_sp = peakMaptoRead_2D.end();
	double peak_sum = 0;
	int peak_score_max = 0;
	for (; peak_st != peak_sp; ++peak_st) // find max score
	{
		peak_sum = peak_sum + 1;
		if (peak_score_max <= peak_st->total_frag)
		{
			peak_score_max = peak_st->total_frag;
		}
	}
	if (peak_score_max >= 1)
	{
		double rzlt[peak_score_max];
		double rzlt_Pvalue[peak_score_max];
		double rzlt_FDR[peak_score_max];
		for (int i = 0; i < peak_score_max; i++)
		{
			rzlt_Pvalue[i] = 0;
			rzlt_FDR[i] = 0;
			rzlt[i] = 0;
		}
		for (peak_st = peakMaptoRead_2D.begin(); peak_st != peak_sp; ++peak_st)
		{
			rzlt[peak_st->total_frag - 1] = rzlt[peak_st->total_frag - 1] + 1.0 / peak_sum;
		}
		for (int i = 0; i < peak_score_max; i++) // p-value
		{
			double tmp = rzlt[i];
			for (int j = 0; j < peak_score_max; j++)
			{
				if (tmp >= rzlt[j])
				{
					rzlt_Pvalue[i] = rzlt_Pvalue[i] + rzlt[j];
				}
			}
		}
		for (int i = 0; i < peak_score_max; i++) // q-value
		{
			int rk = 1;
			for (int k = 0; k < peak_score_max; k++) // get rank value 
			{
				if (rzlt_Pvalue[k] < rzlt_Pvalue[i]) {
					rk++;
				}
			}
			rzlt_FDR[i] = rzlt_Pvalue[i] * peak_score_max / rk;
		}
		
	
                vector<Pair_Frag_2D_Peak_Info>::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
                for ( ; peakIt != peakIt_end; ++peakIt )
                {
                        //output2DpeakFile <<peakIt->chrNo_domain1<<"\t"<< peakIt->peakRegion.first.first << "\t" << peakIt->peakRegion.first.second << "\t"<<peakIt->chrNo_domain2 <<"\t"<<peakIt->peakRegion.second.first<<"\t" << peakIt->peakRegion.second.second << "\t"<<peakIt->total_frag << "\t" <<mix->FDRThreshold(peakIt->total_frag,mix->Mu[0],mix->Mu[1],mix->W[0],mix->W[1])*100 << "\t"<<peakIt->sumLogProb << endl;
                        output2DpeakFile << peakIt->chrNo_domain1 << "\t" << peakIt->peakRegion.first.first << "\t" << peakIt->peakRegion.first.second << "\t" << peakIt->chrNo_domain2 << "\t" << peakIt->peakRegion.second.first << "\t" << peakIt->peakRegion.second.second << "\t" << peakIt->total_frag << "\t" << rzlt_Pvalue[peakIt->total_frag-1] << "\t" << rzlt_FDR[peakIt->total_frag-1] << endl;
                }
        }
        output2DpeakFile.close();
        return 0;
}

int HiC_Peak::Writer_2D_Peak_with_FDR_Random(vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D,vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D_Random,const double thres,const string& peakFiletoWrite_2D)
{
        ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
        vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st = peakMaptoRead_2D.begin(), peak_sp = peakMaptoRead_2D.end();
        double peak_sum = 0;
        int peak_score_max = 0;
        for (; peak_st != peak_sp; ++peak_st) // find max score
        {
                peak_sum = peak_sum + 1;
                if (peak_score_max <= peak_st->total_frag)
                {
                        peak_score_max = peak_st->total_frag;
                }
        }
        vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st_r = peakMaptoRead_2D_Random.begin(), peak_sp_r = peakMaptoRead_2D_Random.end();
        double peak_sum_r = 0;
        int peak_score_max_r = 0;
        for (; peak_st_r != peak_sp_r; ++peak_st_r) // find max score
        {
                peak_sum_r = peak_sum_r + 1;
                if (peak_score_max_r <= peak_st_r->total_frag)
                {
                        peak_score_max_r = peak_st_r->total_frag;
                }
        }
        int score_max = max(peak_score_max_r, peak_score_max);
        double rzlt[score_max];
        double rzlt_r[score_max];
        double rzlt_FDR[score_max];
        for (int i = 0; i < peak_score_max; i++)
        {
                rzlt_r[i] = 0;
                rzlt_FDR[i] = 0;
                rzlt[i] = 0;
        }
        for (peak_st = peakMaptoRead_2D.begin(); peak_st != peak_sp; ++peak_st)
        {
                rzlt[peak_st->total_frag - 1] = rzlt[peak_st->total_frag - 1] + 1.0 / peak_sum;
        }
        for (peak_st_r = peakMaptoRead_2D_Random.begin(); peak_st_r != peak_sp_r; ++peak_st_r)
        {
                rzlt_r[peak_st_r->total_frag - 1] = rzlt_r[peak_st_r->total_frag - 1] + 1.0 / peak_sum_r;
        }
        for (int i = 0; i < peak_score_max; i++)
        {
                if ((rzlt[i] + rzlt_r[i]) > 0)
                {
                        rzlt_FDR[i] = rzlt_r[i] / (rzlt[i] + rzlt_r[i]);
                }
        }

        if (peak_score_max >= 1)
        {
                vector<Pair_Frag_2D_Peak_Info>::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
                for ( ; peakIt != peakIt_end; ++peakIt )
                {
                        if (peakIt->total_frag > 2*thres)
                        {
                                output2DpeakFile << peakIt->chrNo_domain1 << "\t" << peakIt->peakRegion.first.first << "\t" << peakIt->peakRegion.first.second << "\t" << peakIt->chrNo_domain2 << "\t" << peakIt->peakRegion.second.first << "\t" << peakIt->peakRegion.second.second << "\t" << peakIt->total_frag << "\t" << rzlt_FDR[peakIt->total_frag-1] << endl;
                        }
                }
        }
        output2DpeakFile.close();
        return 0;
}


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Cut_Site_Frags(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap_local, 
	                                const string& siteFragFiletoWrite)
{
	ofstream outputFragFile(siteFragFiletoWrite.c_str());
	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(), chrIt_end = cuttingSiteFragsMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt = chrIt->second.begin(), siteIt_end = chrIt->second.end();
		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			outputFragFile <<">\t"<<chrIt->first<<"\t"<<siteIt->first<<endl;
			vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(), fragIt_end = siteIt->second.end();
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				outputFragFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_strand<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_strand<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<endl;
			}
		}
	}
	outputFragFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Frag_Interaction_Map(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local, 
	                                      const string fragInteractionFiletoWrite)
{
	ofstream outputFragInteractionFile(fragInteractionFiletoWrite.c_str());
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt = fragInteractionMap_local.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt_end = fragInteractionMap_local.end();

	for ( ; fragPairIt != fragPairIt_end; ++fragPairIt )
	{
		outputFragInteractionFile <<">\t"<<fragPairIt->first.first.first<<"\t"<<fragPairIt->first.first.second.first<<"\t"<<fragPairIt->first.first.second.second<<"\t"<<fragPairIt->first.second.first<<"\t"<<fragPairIt->first.second.second.first<<"\t"<<fragPairIt->first.second.second.second<<"\t"<<fragPairIt->second.size()<<endl;
		
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = fragPairIt->second.begin();
		vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt_end = fragPairIt->second.end();
		for ( ; fragIt != fragIt_end; ++fragIt )
		{
			outputFragInteractionFile <<fragIt->end1_chr<<"\t"<<fragIt->end1_pos<<"\t"<<fragIt->end1_cuttingSite<<"\t"<<fragIt->end2_chr<<"\t"<<fragIt->end2_pos<<"\t"<<fragIt->end2_cuttingSite<<"\t"<<endl;
		}
	}
	outputFragInteractionFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Digested_Frags(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local,
	                                const string fragInteractionFiletoWrite)
{
	ofstream outputFragInteractionFile(fragInteractionFiletoWrite.c_str());
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt = fragInteractionMap_local.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt_end = fragInteractionMap_local.end();
    
    map<int, int > digested_frag_len1;
    map<int, int > digested_frag_len2;

	for ( ; fragPairIt != fragPairIt_end; ++fragPairIt )
	{
		if(fragPairIt->first.first.first == fragPairIt->first.second.first){
			digested_frag_len1[fragPairIt->first.first.second.second-fragPairIt->first.first.second.first]++;
			digested_frag_len2[fragPairIt->first.second.second.second-fragPairIt->first.second.second.first]++;
		}	
	}

	digested_frag_len1.insert(digested_frag_len2.begin(), digested_frag_len2.end());

	map<int, int> ::const_iterator mapit_s = digested_frag_len1.begin();
	map<int, int> ::const_iterator mapit_e = digested_frag_len1.end();

	for(; mapit_s != mapit_e; ++mapit_s){
		outputFragInteractionFile << mapit_s->first << "\t" << mapit_s->second << endl;
	}

	outputFragInteractionFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Frag_Freq_Map(map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> > fragInteractionMap_local, 
	                                      const string fragFreqFiletoWrite)
{
	ofstream outputFragFreqFile(fragFreqFiletoWrite.c_str());
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt = fragInteractionMap_local.begin();
	map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >::const_iterator fragPairIt_end = fragInteractionMap_local.end();

	map<int, int> fragFreqmap;
	for ( ; fragPairIt != fragPairIt_end; ++fragPairIt )
	{
		fragFreqmap[fragPairIt->second.size()]++;
	}

	map<int, int> ::const_iterator mapit_s = fragFreqmap.begin();
	map<int, int> ::const_iterator mapit_e = fragFreqmap.end();

	for(; mapit_s != mapit_e; ++mapit_s){
		outputFragFreqFile << mapit_s->first << "\t" << mapit_s->second << endl;
	}

	outputFragFreqFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Dis_Distribution(const map<int, int>& distributionMap_local, const string& distributionFiletoWrite)
{
	ofstream outputFile(distributionFiletoWrite.c_str());

	map<int, int>::const_iterator binIt = distributionMap_local.begin(), binIt_end = distributionMap_local.end();
	for ( ; binIt != binIt_end; ++binIt )
	{
		outputFile<<binIt->first<<"\t"<<binIt->second<<endl;
	}

	outputFile.close();

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int HiC_Peak::Writer_Dis_Distribution_Prob(const map<int, double>& distributionMap_local, const string& distributionFiletoWrite)
{
	ofstream outputFile(distributionFiletoWrite.c_str());

	map<int, double>::const_iterator binIt = distributionMap_local.begin(), binIt_end = distributionMap_local.end();
	for ( ; binIt != binIt_end; ++binIt )
	{
		outputFile<<binIt->first<<"\t"<<binIt->second<<endl;
	}

	outputFile.close();

	return 0;
}
//--------------------------------------------------------------------------------------------
// should the string be NULL, just perform the counting
void HiC_Peak::Writer_FragSize(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap, map<unsigned int,unsigned int>& count_res_frags, map<unsigned int, unsigned int>& count_prob_ligation, string& fragFile){

				ofstream FragSizeFile;
				if (!fragFile.empty())
					FragSizeFile.open(fragFile.c_str());

				map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap.begin(); 
				map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end = cuttingSiteFragsMap.end();

				// iterate throught each chomosome
				for(;chrIt != chrIt_end;++chrIt){
					
					map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It 		= chrIt->second.begin(); 
					map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It_end 	= chrIt->second.end();

					for(;CS_It != CS_It_end; ++CS_It){

						// count restriction fragments
						if (CS_It->second.size() != 0){
							vector<Frag_Info_With_Cutting_Site> ::const_iterator frag_It = CS_It->second.begin();
							vector<Frag_Info_With_Cutting_Site> ::const_iterator frag_It_end = CS_It->second.end();
							
							
							count_res_frags[CS_It->second.size()] += CS_It->second.size();

							// count_res_frags[CS_It->second.size()] ++;

							vector<Frag_Info_With_Cutting_Site> ::const_iterator fragIt = CS_It->second.begin();
							vector<Frag_Info_With_Cutting_Site> ::const_iterator fragIt_end = CS_It->second.end();
							for(;fragIt != fragIt_end;++fragIt){
								if(fragIt->end1_chr==fragIt->end2_chr){
									if (!fragFile.empty()){
										FragSizeFile << fragIt->end1_chr << "\t" << fragIt->end1_pos << "\t" << fragIt->end2_chr << "\t" << fragIt->end2_pos <<"\t" << CS_It->second.size() << endl;
									}
									count_prob_ligation[abs(fragIt->end1_cuttingSite-fragIt->end2_cuttingSite)]++;
								}
							}	
						}
					}
				}
	FragSizeFile.close();
}

// write the binary of these values to the file
void HiC_Peak::BinWriter_PoissMix(map<unsigned int,unsigned int>& count_res_frags, string& pstring){
	int fd = open(pstring.c_str(), O_WRONLY | O_CREAT, 0666);

	map<unsigned int, unsigned int> ::const_iterator mapit_s = count_res_frags.begin();
	map<unsigned int, unsigned int> ::const_iterator mapit_e = count_res_frags.end();

	// output to the chr_file
	for(; mapit_s != mapit_e; ++mapit_s){
		// add to the final PoissMix
		unsigned int tmp = mapit_s->second;
		// fprintf(stderr, "%u\n", mapit_s->second);
		// unsigned int tmp[2];
		// tmp[0] = mapit_s->first;
		// tmp[1] = mapit_s->second;
		r_write(fd, &tmp, sizeof(unsigned int));
	}
	
	close(fd);
}

void HiC_Peak::Writer_PoissMix(map<unsigned int,unsigned int>& count_res_frags, string& pstring){
	ofstream PoissMixFile(pstring.c_str());
	map<unsigned int, unsigned int> ::const_iterator mapit_s = count_res_frags.begin();
	map<unsigned int, unsigned int> ::const_iterator mapit_e = count_res_frags.end();

				
	// output to the chr_file
	for(; mapit_s != mapit_e; ++mapit_s){
		// add to the final PoissMix
		PoissMixFile << mapit_s->first << "\t" << mapit_s->second << endl;

	}



	PoissMixFile.close();
}


void HiC_Peak::Writer_LigateProb(map<int, int>& count_prob_ligation, string& ProbLigationFile){
		ofstream ProbLigateFile(ProbLigationFile.c_str());

		map<int,int> ::const_iterator probit_s  = count_prob_ligation.begin();
		map<int,int> ::const_iterator probit_e	= count_prob_ligation.end();

		for(; probit_s != probit_e; ++probit_s){
			ProbLigateFile << probit_s->first << "\t" << probit_s->second << endl;
		}
		ProbLigateFile.close();
}


double HiC_Peak::poissonf(double*x,double*par)                                         
{                                                                              
  return par[0]*Poisson(x[0],par[1]) + par[2]*Poisson(x[0],par[3]);
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	HiC_Peak dataFile(argc, argv);

	return 0;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// utility to copy values
/*
void HiC_Peak::fragCopy(cBundle& bundle, const Frag_Info_With_Cutting_Site& frag){
	bundle.frag.end1_chr = frag.end1_chr;
	bundle.frag.end1_pos = frag.end1_pos;
	bundle.frag.end1_cuttingSite = frag.end1_cuttingSite;
	bundle.frag.end1_strand = frag.end1_strand;

	bundle.frag.end2_chr = frag.end2_chr;
	bundle.frag.end2_pos = frag.end2_pos;
	bundle.frag.end2_cuttingSite = frag.end2_cuttingSite;
	bundle.frag.end2_strand = frag.end2_strand;
}

void HiC_Peak::printcBundle(cBundle& bundle){
	fprintf(stderr, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", bundle.chr, bundle.length,
		bundle.frag.end1_chr, bundle.frag.end1_pos, bundle.frag.end1_cuttingSite,
		bundle.frag.end1_strand, bundle.frag.end2_chr, bundle.frag.end2_pos,
		bundle.frag.end2_cuttingSite, bundle.frag.end2_strand);
}

void HiC_Peak::cbundleCopy(Frag_Info_With_Cutting_Site& frag, cBundle& bundle){
	// need to make a copy, this is immutable
	frag.end1_chr = bundle.frag.end1_chr;
	frag.end1_pos = bundle.frag.end1_pos;
	frag.end1_cuttingSite = bundle.frag.end1_cuttingSite;
	frag.end1_strand = bundle.frag.end1_strand;

	frag.end2_chr = bundle.frag.end2_chr;
	frag.end2_pos = bundle.frag.end2_pos;
	frag.end2_cuttingSite = bundle.frag.end2_cuttingSite;
	frag.end2_strand = bundle.frag.end2_strand;

}*/

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double HiC_Peak::FindQT(vector<double> *para_values)
{
	typedef vector<double>::size_type vecSize;
    	vecSize N = para_values->size();
	
	// declare new variables
	vecSize NMod4 = (N % 4);  // identification of 1 of the 4 known datum distribution profiles
	string datumDistr = "";   // datum distribution profile
	vecSize M, ML, MU;        // core vector indices for quartile computation
	double m, ml, mu;         // quartile values are store here
	   
	sort(para_values->begin(),para_values->end());

        // compute quartiles for the 4 known patterns	
	if ( NMod4 == 0 ){
        	// Q1-Q3 datum distribution: [0 0 0]
        	datumDistr = "[0 0 0]";
        	M = N / 2;
        	ML = M / 2;
        	MU = M + ML;
                                          
        	// grab quartile values
        	ml= (para_values->at(ML) + para_values->at(ML-1)) / 2;     // datum: 0
        	m = (para_values->at(M) + para_values->at(M-1)) / 2;       // datum: 0
        	mu = (para_values->at(MU) + para_values->at(MU-1)) / 2;    // datum: 0
        }

	else if ( NMod4 == 1 ){
        	// Q1-Q3 datum distribution: [0 1 0]
                datumDistr = "[0 1 0]";
                M = N / 2;
                ML = M / 2;
                MU = M + ML + 1;
                                          
                // grab quartile values
                datumDistr = "[0 0 0]";
                ml= (para_values->at(ML) + para_values->at(ML-1)) / 2;      // datum: 0
                m = para_values->at(M);                       // datum: 1
                mu = (para_values->at(MU) + para_values->at(MU-1)) / 2;     // datum: 0
       }
	
	else if ( NMod4 == 2 ){
        	datumDistr = "[1 0 1]";
        	M = N / 2;
        	ML = M / 2;
        	MU = M + ML;
 
        	// grab quartile values
                ml= para_values->at(ML);                    // datum: 1
                m = (para_values->at(M) + para_values->at(M-1)) / 2;     // datum: 0
                mu = para_values->at(MU);                   // datum: 1
       }

	else if ( NMod4 == 3 ){
        
        	datumDistr = "[1 1 1]";
        	M = N / 2;
        	ML = M / 2;
        	MU = M + ML + 1;
 
       
        	ml= para_values->at(ML);                    // datum: 1
        	m = para_values->at(M);                     // datum: 0
        	mu = para_values->at(MU);                   // datum: 1
    	}

	else{
        	cout << "Unknown pattern discovered - new algorithm may be required.";
		return 0;
	}

	double IQR 	= mu-ml;
	double start 	= ml-1.5*IQR;
	double end  	= mu+1.5*IQR;
	int elements    = 0;
	double sum      = 0;
	for (vecSize i = 0; i<N; i++){
        	if(para_values->at(i) >  start && para_values->at(i) < end){
			sum+=para_values->at(i);
			elements++;
		}
		else continue;
    	}
	double average = sum/elements;

	return m;
}


// gather information from this data structure, write it to a file
// print this informatin to the given file
/*
void fragMemDump(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& map, string filename){
	ofstream out(filename.c_str());

	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = fragInteractionMap_local.begin(); 
	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt_end = fragInteractionMap_local.end();

	// iterate through printing each value for later re-construction
	for(;chrIt != chrIt_end;++chrIt){
		map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It 		= chrIt->second.begin(); 
		map< int, vector<Frag_Info_With_Cutting_Site> > ::const_iterator CS_It_end 	= chrIt->second.end();

		for(;CS_It != CS_It_end; ++CS_It){
			
			vector<Frag_Info_With_Cutting_Site>::const_iterator vecIt = CS_It->begin();
			vector<Frag_Info_With_Cutting_Site>::const_iterator vecIt_end = CS_It->end();
			for (;vecIt != vecIt_end; vecIt++){
					out << vecIt->end1_chr << " " << vecIt->end1_pos << " " << vecIt->end1_cuttingSite << " " << end1_strand;
					out << vecIt->end2_chr << " " << vecIt->end2_pos << " " << vecIt->end2_cuttingSite << " " << end2_strand;
				}
			}
		}
	}


	out.close();

	// remove from memory
	map.clear();
}


// recover data dumped to the file
void fragMemRecovery(map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >&, string filename){


}*/


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int HiC_Peak::Pair_Frags_Interaction_Freq_With_Both_Side(const map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >& cuttingSiteFragsMap_local,
										map< int, map<int, int> >& digestedFragMap_local,
	                                    map< pair< pair< int, pair<int, int> >, pair< int, pair<int, int> > >, vector<Frag_Info_With_Cutting_Site> >& freqMaptoWrite)
{
	freqMaptoWrite.clear();
	map< int, map< int, vector<Frag_Info_With_Cutting_Site> > >::const_iterator chrIt = cuttingSiteFragsMap_local.begin(),
		chrIt_end = cuttingSiteFragsMap_local.end();
	for ( ; chrIt != chrIt_end; ++chrIt )
	{
		map< int, vector<Frag_Info_With_Cutting_Site> >::const_iterator siteIt = chrIt->second.begin(),
			siteIt_end = chrIt->second.end();

		for ( ; siteIt != siteIt_end; ++siteIt )
		{
			vector<Frag_Info_With_Cutting_Site>::const_iterator fragIt = siteIt->second.begin(),
				fragIt_end = siteIt->second.end();
				
			for ( ; fragIt != fragIt_end; ++fragIt )
			{
				// find the cutting site of end1 and save to frag1_start, frag1_end
				int frag1_chr = 0;
				int frag1_start = 0; //a pair of interacting frag, frag1 is end1, frag2 is end2
				int frag1_end = 0;
				if (fragIt->end1_strand)
				{
					//find the cutting site of end1 in positive strand
					map< int, map<int, int> >::const_iterator cuttingchr = digestedFragMap_local.find(fragIt->end1_chr);
					frag1_chr = cuttingchr->first;
					map<int, int>::const_iterator cuttingsite = cuttingchr->second.find(fragIt->end1_cuttingSite);
					--cuttingsite;
					frag1_start = cuttingsite->first;
					frag1_end = cuttingsite->second;

				}
				else
				{
					//find the cutting site of end1 in negative strand
					map< int, map<int, int> >::const_iterator cuttingchr = digestedFragMap_local.find(fragIt->end1_chr);
					frag1_chr = cuttingchr->first;
					map<int, int>::const_iterator cuttingsite = cuttingchr->second.find(fragIt->end1_cuttingSite);
					frag1_start = cuttingsite->first;
					frag1_end = cuttingsite->second;
				}
				
				// find the cutting site of end2 and save to frag2_start, frag2_end
				int frag2_chr = 0;
				int frag2_start = 0;
				int frag2_end = 0;
				if (fragIt->end2_strand)
				{
					//find the cutting site of end1 in positive strand
					map< int, map<int, int> >::const_iterator cuttingchr = digestedFragMap_local.find(fragIt->end2_chr);
					frag2_chr = cuttingchr->first;
					map<int, int>::const_iterator cuttingsite = cuttingchr->second.find(fragIt->end2_cuttingSite);
					--cuttingsite;
					frag2_start = cuttingsite->first;
					frag2_end = cuttingsite->second;

				}
				else
				{
					//find the cutting site of end1 in negative strand
					map< int, map<int, int> >::const_iterator cuttingchr = digestedFragMap_local.find(fragIt->end2_chr);
					frag2_chr = cuttingchr->first;
					map<int, int>::const_iterator cuttingsite = cuttingchr->second.find(fragIt->end2_cuttingSite);
					frag2_start = cuttingsite->first;
					frag2_end = cuttingsite->second;
				}
				
				freqMaptoWrite[make_pair(make_pair(frag1_chr, make_pair(frag1_start, frag1_end)), make_pair(frag2_chr, make_pair(frag2_start, frag2_end)))].push_back(*fragIt);
			}
		}
	}

	return 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//function: Calculate FDR, Filtering the size of fragment , Remove the duplicate and Save to file
int HiC_Peak::Writer_2D_Peak_with_FDR_Random_And_Length_Limit(
	vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D,
	vector<Pair_Frag_2D_Peak_Info>& peakMaptoRead_2D_Random,
	const double thres,
	const string& peakFiletoWrite_2D,
	int fragSize_local)
{
    ofstream output2DpeakFile(peakFiletoWrite_2D.c_str());
    vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st = peakMaptoRead_2D.begin(), peak_sp = peakMaptoRead_2D.end();
    double peak_sum = 0;
    int peak_score_max = 0;
    for (; peak_st != peak_sp; ++peak_st) // find max score
    {
        peak_sum = peak_sum + 1;
        if (peak_score_max <= peak_st->total_frag)
        {
            peak_score_max = peak_st->total_frag;
        }
    }
    vector<Pair_Frag_2D_Peak_Info>::const_iterator peak_st_r = peakMaptoRead_2D_Random.begin(), peak_sp_r = peakMaptoRead_2D_Random.end();
    double peak_sum_r = 0;
    int peak_score_max_r = 0;
    for (; peak_st_r != peak_sp_r; ++peak_st_r) // find max score
    {
        peak_sum_r = peak_sum_r + 1;
        if (peak_score_max_r <= peak_st_r->total_frag)
        {
            peak_score_max_r = peak_st_r->total_frag;
        }
    }
    int score_max = max(peak_score_max_r, peak_score_max);
    double rzlt[score_max];
    double rzlt_r[score_max];
    double rzlt_FDR[score_max];
    for (int i = 0; i < peak_score_max; i++)
    {
        rzlt_r[i] = 0;
        rzlt_FDR[i] = 0;
        rzlt[i] = 0;
    }
    for (peak_st = peakMaptoRead_2D.begin(); peak_st != peak_sp; ++peak_st)
    {
        rzlt[peak_st->total_frag - 1] = rzlt[peak_st->total_frag - 1] + 1.0 / peak_sum;
    }
    for (peak_st_r = peakMaptoRead_2D_Random.begin(); peak_st_r != peak_sp_r; ++peak_st_r)
    {
        rzlt_r[peak_st_r->total_frag - 1] = rzlt_r[peak_st_r->total_frag - 1] + 1.0 / peak_sum_r;
    }
    for (int i = 0; i < peak_score_max; i++)
    {
        if ((rzlt[i] + rzlt_r[i]) > 0)
        {
            rzlt_FDR[i] = rzlt_r[i] / (rzlt[i] + rzlt_r[i]);
        }
    }
	//give the output frags to the vector outputfrags
	vector<pair<pair<pair<int, pair<int, int> >, pair<int, pair<int, int> > >, pair<int, double> > > outputfrags;
    if (peak_score_max >= 1)
    {
        vector<Pair_Frag_2D_Peak_Info>::const_iterator peakIt = peakMaptoRead_2D.begin(), peakIt_end = peakMaptoRead_2D.end();
        for ( ; peakIt != peakIt_end; ++peakIt )
        {
            if (peakIt->total_frag > 2*thres)
            {
				outputfrags.push_back(make_pair(make_pair(
				make_pair(peakIt->chrNo_domain1, make_pair(peakIt->peakRegion.first.first, peakIt->peakRegion.first.second)), 
				make_pair(peakIt->chrNo_domain2, make_pair(peakIt->peakRegion.second.first, peakIt->peakRegion.second.second))), 
				make_pair(peakIt->total_frag, rzlt_FDR[peakIt->total_frag-1])));
            }
        }
    }
	//transfer from outputfrags to dedupfrag for filtering the size to < 200K and give the mid of fragments
	vector<pair<pair<pair<int, int>, pair<int, int> >, pair<int, double> > > dedupfrag;
	vector<pair<pair<pair<int, pair<int, int> >, pair<int, pair<int, int> > >, pair<int, double> > >::const_iterator
		outputIt = outputfrags.begin(),
		outputIt_end = outputfrags.end();
	for (; outputIt != outputIt_end; ++outputIt)
	{
		int frag1_start = outputIt->first.first.second.first;
		int frag1_end = outputIt->first.first.second.second;
		int frag2_start = outputIt->first.second.second.first;
		int frag2_end = outputIt->first.second.second.second;
		if ((frag1_end - frag1_start < 200000) && (frag2_end - frag2_start < 200000))
		{
			dedupfrag.push_back(make_pair(make_pair(
			make_pair(outputIt->first.first.first, (frag1_start + frag1_end) / 2),
			make_pair(outputIt->first.second.first, (frag2_start + frag2_end) / 2)),
			make_pair(outputIt->second.first, outputIt->second.second)));
		}
	}
	// remove the duplicates and transfer to output2DpeakFile
	for (int i = 0; i < dedupfrag.size() - 1; ++i)
	{
		int overlaptag = 0;
		for (int j = i+1; j < dedupfrag.size(); ++j)
		{
			//chr1[i] == chr2[j] && mid1[i] == mid2[j] && chr2[i] == chr1 [j] && mid2[i] == mid1[j]
			if((dedupfrag[i].first.first.first == dedupfrag[j].first.second.first) && 
			   (dedupfrag[i].first.first.second == dedupfrag[j].first.second.second) &&
			   (dedupfrag[i].first.second.first == dedupfrag[j].first.first.first) &&
			   (dedupfrag[i].first.second.second == dedupfrag[j].first.first.second))
			{
				overlaptag = 1;
			}
		}
		if (overlaptag == 0)
		{
			output2DpeakFile << dedupfrag[i].first.first.first << "\t" 
							 << dedupfrag[i].first.first.second -  fragSize_local / 2 << "\t" 
							 << dedupfrag[i].first.first.second +  fragSize_local / 2 << "\t" 
							 << dedupfrag[i].first.second.first << "\t" 
							 << dedupfrag[i].first.second.second -  fragSize_local / 2 << "\t" 
							 << dedupfrag[i].first.second.second +  fragSize_local / 2 << "\t" 
							 << dedupfrag[i].second.first << "\t" 
							 << dedupfrag[i].second.second
							 << endl;
		}
	}
	// save the last peaks
	output2DpeakFile << dedupfrag[dedupfrag.size() - 1].first.first.first << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].first.first.second -  fragSize_local / 2 << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].first.first.second +  fragSize_local / 2 << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].first.second.first << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].first.second.second -  fragSize_local / 2 << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].first.second.second +  fragSize_local / 2 << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].second.first << "\t" 
					 << dedupfrag[dedupfrag.size() - 1].second.second
					 << endl;

    output2DpeakFile.close();
    return 0;
}
