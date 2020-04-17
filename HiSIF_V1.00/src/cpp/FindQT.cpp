#include <vector>
#include <string>
#include <algorithm>
#include <iostream>


using namespace std;
// used to find quartiles

double FindQT(vector<double> *para_values)
{
	typedef vector<double>::size_type vecSize;
    	vecSize N = para_values->size();
	
	// declare new variables
	vecSize NMod4 = (N % 4);  // identification of 1 of the 4 known datum distribution profiles
	string datumDistr = "";   // datum distribution profile
	vecSize M, ML, MU;        // core vector indices for quartile computation
	double m, ml, mu;         // quartile values are store here
	   
	sort(para_values->begin(),para_values->end());

	// printf("FindQT Sorted!\n");
        // compute quartiles for the 4 known patterns	
	if ( NMod4 == 0 ){
			// printf("NMod4 == 0!\n");
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
			// printf("NMod4 == 1!\n");
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
			// printf("NMod4 == 2!\n");
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
			// printf("NMod4 == 3!\n");
        
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
	for (int i = 0; i<N; i++){
		// fprintf(stderr, "para_values->at(%d) == %f\n", para_values->at(i));
		if(para_values->at(i) >  start && para_values->at(i) < end){
			sum+=para_values->at(i);
			elements++;
		}
		else continue;
	}

	// not returning this??? What do we need this for?
	// double average = sum/elements;

	return m;
}
