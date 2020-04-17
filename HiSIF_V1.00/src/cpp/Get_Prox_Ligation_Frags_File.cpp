#include <cstdlib>
#include <ctime>
#include <stack>
#include <cstring>

// PoissMix definition
#include <CPP_Utilities.hpp>
#include <C_Utilities.hpp>

using namespace std;

struct iteration{
	int count;

	// how many lhoods to hold
	double lhoods[100];
};

typedef struct iteration iterTrack;


// an unbiased version of rand()
int unbiasedRand(int n){
	int x;

		do {
		 x = rand();
	} while (x >= (RAND_MAX - RAND_MAX % (n+1)));

	return x %= (n+1);
}


/****************************************************************************************
 * Purpose:
 * 	Given original dataset, get sampleSize number of random elements.
 * 	Push them into the vector of doubles x, for later use.
 *
 * Parameters:
 * 	PoissMix& ts					training set PM
 * 	vector<double>& dataset		entire dataset to use
 * 	int sampleSize					how many elements to get
 * 	int dataSize					size of entire dataset
 ***************************************************************************************/
/*void generateRandomDataset(PoissMix& ts, vector<double>& dataset, int sampleSize, int dataSize){
	int j;
	
	for (j = 0; j < sampleSize; j++){
		double val = dataset[unbiasedRand(dataSize)];
		ts.x.push_back(val);
	}
}*/


/****************************************************************************************
 * Purpose:
 * 	Given a file of entire dataset, get sampleSize of elements from it.
 *
 * Parameters:
 * 	PoissMix& ts					training set PM
 * 	int fd							file descriptor for data
 * 	int sampleSize					how many elements
 * 	int dataSize					size of the data (# elements)
 ***************************************************************************************/
void generateRandomDatasetFile(PoissMix& ts, int fd, int sampleSize, int dataSize){
	int j, bytesread;

		for (j = 0; j < sampleSize; j++){

		unsigned int val;
		off_t offset = sizeof(unsigned int) * unbiasedRand(dataSize);
		// fprintf(stderr, "val location: %d\n", offset / sizeof(unsigned int));
		// move file pointer to the location
		if (lseek(fd, offset, SEEK_SET) < 0){
			fprintf(stderr, "Error: lseek failed\n");
			return;
		}

		// read an unsigned int
		bytesread = r_read(fd, (char *)&val, sizeof(unsigned int));
		// fprintf(stderr, "bytesread == %d\n", bytesread);
		// fprintf(stderr, "random value: %u\n", val);
		ts.x.push_back(val);
	}
}


/****************************************************************************************
 * Purpose:
 * 	Modify values in the Poisson Mixture object. This sets up for updating with
 * 	the new values.
 *
 * Parameters:
 * 	PoissMix& ts					trainign set we're using
 * 	double W1						first value for W
 * 	double W2						second value for W
 * 	double Mu1						first value for Mu
 * 	double Mu2						second value for Mu
 * 	int M								size of trainig set (number of datapoints)
 * 	int N								number of mixtures
 ***************************************************************************************/
void modPM(PoissMix& ts, double W1, double W2, double Mu1, double Mu2, int M, int N){
	ts.W.push_back(W1);
	ts.W.push_back(W2);
	ts.Mu.push_back(Mu1);
	ts.Mu.push_back(Mu2);
	ts.M = M;
	ts.N = N;
}

/******************************************************************************
 * Purpose:
 * 	Iteratively call update methods for W, Mu, and entire PoissMix model, as 
 * 	well as checking whether the lhood values are a difference of under the 
 * 	expected precision (ep). Also ensures both alpha values are under ep.
 *
 * Parameters:
 * 	PoissMix& ts					reference to the training set we're using
 * 	int maxIter						maximum iterations to perform for this set
 * 	double ep						expected precision value
 * 	struct iterTrack it			holds max iteration count and lhood values
 *
 *
 * Notes:
 * 	we will keep the maximum iteration for the entire bootstrapping,
 * 	this will be checked here.
 *****************************************************************************/


void updatePM(PoissMix& ts, int maxIter, double ep, ofstream& out, iterTrack *it){
	double lam1, lam2, alp1, alp2;
	long double lhood;


	// create thing to track doubles
	iterTrack tempTrack;

	// zero out the values
	// memset(tempTrack.lhoods, 0, sizeof(tempTrack.lhoods));
	tempTrack.count = 0;

	for(int k=0; k < maxIter; k++){
    lhood = ts.lhood();
    ts.W_update();
    ts.Mu_update();
    ts.update();

	// here for RAO, we're getting nan
   cout << "\tIteration: "<< k+1<< "--> LLH: " << ts.lhood() << endl;
 
	// store the lhood value
	tempTrack.count++;
	tempTrack.lhoods[k] = lhood;

	if (trunc(10.*ts.lhood())-trunc(10.*lhood) <= ep && k > 1){
		lam1 = ts.Mu[0];
		lam2 = ts.Mu[1];
		alp1 = ts.W[0];
		alp2 = 1 - alp1;

      out << (int)lam1 << "\t\t\t" << (int)lam2 << "\t\t\t" << alp1 << "\t" << alp2 << endl;
      break;
    }
  }

	
	// is this bigger than what we have?
	if (tempTrack.count > it->count){
		memcpy(it, &tempTrack, sizeof(iterTrack));
	}
}

/************************************************************************************************************************
 * Purpose:
 * 	Get proximate ligation events using bootstrapping of n datasets.
 *
 * Parameters:
 * 	int fd																										file containing all vector sizes
 * 	vector<double> poisVec																					vector with start values
 * 	float pSampleSize																							% of total size dataset for training
 * 	int n																											number of total iterations
 * 	ofstream& out																								file for bootstrapping
 * 	ofstream& iterationFile																					file for iteration output
 *	
 *	Results:
 *		Creates n training datasets, and uses these to stabilize the initial parameters
 *
 * Notes:
 * 	Used to store entire dataset in mydata_set, now only stored in mix PM
 ***********************************************************************************************************************/
PoissMix& Get_Prox_Ligation_Frags_File(int fd, vector<double>& poisVec, float pSampleSize, int n, ofstream& out, ofstream& iterationFile)
{
	srand(time(NULL));						// seed random generator

	long double ep     		= 10e-20;
  	long maxIter  		= 100;

	int i;
	int sampleSize;							// # elements from total size

	iterTrack it;								// used to find the maximum iterations for this dataset

	// zero out
	it.count = 0;

	PoissMix ts(2);							// training set used
	PoissMix *mix = new PoissMix(2);		// final returned parameters
		
	// track each iteration values
	vector<double> vecAl1;
	vector<double> vecAl2;
	vector<double> vecLam1;
	vector<double> vecLam2;

	// how many elements in this file?
	unsigned int elements = lseek(fd, 0, SEEK_END);
	elements /= sizeof(unsigned int);

	cout << "<-----Using samplesize of " << elements << " elements----->" << endl;
	lseek(fd, 0, SEEK_SET);
	// use % of dataset size
	sampleSize = elements * pSampleSize;

	out << "Lamda1\tLamda2\tAlpha1\tAlpha2" << endl;

	// generate n random datasets, perform bootstrapping, save values
	for (i = 0; i < n; i++){
		printf("<-----Random Dataset %d----->\n", i+1);
		ts.x.clear();
		generateRandomDatasetFile(ts, fd, sampleSize, elements);
		modPM(ts, 0.1, 0.9, poisVec[0], poisVec[1], ts.x.size(), 2);
		updatePM(ts, maxIter, ep, out, &it);

		// store these values
		vecAl1.push_back(ts.W[0]);
		vecAl2.push_back(ts.W[1]);
		vecLam1.push_back((int)ts.Mu[0]);
		vecLam2.push_back((int)ts.Mu[1]);

		// reset W and Mu
		ts.W.clear();
		ts.Mu.clear();
	}

	// find quartiles, and use these
	double alpha1 = FindQT(&vecAl1);
	mix->W.push_back(alpha1);
	mix->W.push_back(1-alpha1);
	mix->Mu.push_back(FindQT(&vecLam1));
	mix->Mu.push_back(FindQT(&vecLam2));
	
	out << "Final Quartile Values\n" << endl;
	out << (int)mix->Mu[0] << "\t\t\t" << (int)mix->Mu[1] << "\t\t\t" << mix->W[0] << "\t" << mix->W[1] << endl;

	// print out the maximum iteration values
	// out << "\n\n\nMaximum Iteration Count:\n" << endl;
	for (i = 0; i < it.count; i++){
		iterationFile << i+1 << " " << it.lhoods[i] << endl;
	}


	// remove the values here in mix
	mix->x.clear();
	return (*mix);
}
