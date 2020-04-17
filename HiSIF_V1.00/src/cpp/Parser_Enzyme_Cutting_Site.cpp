/*************************************************************************
 * Purpose:
 * 	Read through a .bed file given, read each line and gather up all
 * 	of the cutting sites, per chromosome, and store into a map.
 *
 * 	Makes use of the extractChrNum() function
 * 	Optional: print the results to a tab-delimited file
 *
 * Parameters:
 * 	const string&					path to the .bed file
 * 	map<int, vector<int>>&		reference to the pre-made map
 *		const string&					path to the output file
 *
 * Return Values:
 * 	-1									error occured
 * 	0									success
 *
 *	Notes:
 *
 ************************************************************************/
// C things
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


// C++ things
#include <cstdio>
#include <string>
#include <map>
#include <vector>

#include <CPP_Utilities.hpp>
#include <C_Utilities.hpp>

using namespace std;

void printMap(map< int, vector<int> >&, char *filename);

int Parser_Enzyme_Cutting_Site(const string& filestr, map< int, vector<int> >& map, const string& outstr){

																				// convert to C strings */
	char *filepath = (char *)filestr.c_str();
	char *out;
	if (!outstr.empty())
		out = (char *)outstr.c_str();
	else
		out = NULL;

																// file descriptor for the bed file */
	int bedfd;

	int nbytes;
	
	char buf[1024];

	int total = 0;
																							// param check	*/
	if (filepath == NULL){
		fprintf(stderr, "Error: passed a null filepath\n");
		return -1;
	}
	
	// attempt to open the file
	if ((bedfd = open(filepath, O_RDONLY)) == -1){
		fprintf(stderr, "Error: could not open file %s\n", filepath);
		perror(NULL);
		return -1;
	}
																								// vectors */
	vector<int> vectors[24];

	// fprintf(stderr, "Starting the vectors\n");
																	// read in each line, and parse */
	char *tmp, *chr;
	memset(buf, 0, sizeof(buf));
	while ((nbytes = readline(bedfd, buf, sizeof(buf))) > 0){
		total++;
																				// chr number field */
		chr = strtok(buf, " 	");

											  // read in the first position, we need the next */
		tmp = strtok(NULL, " 	");
		tmp = strtok(NULL, " 	");

														  // add value to it's respective vector */
		if ((extractChrNum(chr)-1) >= 0)
		{
			vectors[extractChrNum(chr)-1].push_back(atoi(tmp));
			// fprintf(stderr, "Chr%d\tpos == %d\n", extractChrNum(chr), atoi(tmp));
			memset(buf, 0, sizeof(buf));
		}
	}
																							// read error */
	if (nbytes == -1){
		perror("Error: there was a read error.\n");
		return -1;
	}

																		 // add vectors into the map */
	unsigned int i;
	for (i = 1; i <= 24; i++)
		map[i] = vectors[i-1];

	if (out != NULL){
		printMap(map, out);
	}
	return total;
}

// print a map, in tabular formatting to be imported to excel
void printMap(map< int, vector<int> >& map, char *filename){
	ofstream out;

	out.open(filename);

	// file check
	if (out.is_open() == false){
		cerr << "Error: unable to open file given.\n";
		return;
	}

	unsigned int i,j;
														  // content of the map */
	for (i = 1; i <= 22; i++){
		// fprintf(stderr, "%d\n", (*map)[i]);
		out << "chr" << i;

													// iterate through vector */
		// fprintf(stderr, "map size == %d\n", map[
		for (j = 0; j < map[i].size(); j++)
			out << "\t" << map[i][j];
		out << "\n";
	}

														// odd ball chromosomes */
	out << "chrX";
	for (j = 0; j < map[i].size(); j++)
		out << "\t" << map[i][j];
	out << "\n";

	out << "chrY";

	for (j = 0; j < map[i+1].size(); j++)
		out << "\t" << map[i+1][j];

	out.close();
}
