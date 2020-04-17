/********************************************************************
 * Purpose:
 *		Given a directory of .fa files, create individual processes for
 *		each file, generating their lengths and storing them into their
 *		respective maps.
 *
 *	Parameters:
 *		char *dirpath					path to directory containing .fa files
 *		std::map<int, int> *map		pointer to pre-defined std::map
 *		char *out						output filename, NULL if not using
 *
 * Return Value:
 * 	-1			there was an error
 * 	0			no error occured, the completed map is stored in map
 *
 * 	*map -> will have the chr (number only) as keys, with their
 * 	respective lengths
 *	Note:
 *		ChrX and ChrY should be stored as Chr23 and Chr24
 *		Improve performance by adding parrallelization!
 *
 *******************************************************************/
#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

// #include <string.h>

// C++ things
#include <string>
#include <cstring>

#include <CPP_Utilities.hpp>
#include <C_Utilities.hpp>

using namespace std;

// prototype functions
static void printMap(map<int, int>& map, char *filename);

int Map_Chr_Len(const string& dirstr, map<int, int>& map, const string& outstr){

	char *filepath = (char *)dirstr.c_str();
	char *out = (char *)outstr.c_str();
										 // list of files from the directory */
	char **filelist;
	// char buf[1024];
								 // number of files in the given directory */
	int numfiles;

	if (filepath == NULL){
		fprintf(stderr, "Error: passed a null filepath\n");
		return -1;
	}

								// make sure filepath has a '/' at the end */
	
	char *newdirpath = (char*)malloc(NAME_MAX * sizeof(char));
	strcpy(newdirpath, filepath);
	if (filepath[strlen(filepath) - 1] != '/'){
		// add on the new '/', as well as the null-terminator */

		newdirpath[strlen(filepath)] = '/';
		newdirpath[strlen(filepath)+1] = '\0';
	}
												 // string array of file names */

	// fprintf(stderr, "newdirpath = %s\n", newdirpath);
	filelist = getfilelist(newdirpath, &numfiles);


													 // don't continue if error */
	if (filelist == NULL){
		fprintf(stderr, "Error: getfilelist() returned a null array\n");
		exit(-1);
	}


	// fprintf(stderr, "Got file list\n");
	

	int i;
	// for (i = 0; i < numfiles; i++)
		// fprintf(stderr, "filelist[%d] == %s\n", i, filelist[i]);
	
	char *backuppath = (char *)malloc(NAME_MAX * sizeof(char));
	strcpy(backuppath, newdirpath);
	for (i = 0; i < numfiles; i++){
		// generate the string for the filename
		// memset(buf, sizeof(buf), 0);

											  // prepare the file for reading */
		// fprintf(stderr, "newdirpath == %s\nfilelist[i] == %s\nstrlen(filelist[i]) == %d\n",
			// newdirpath, filelist[i], strlen(filelist[i]));
		//strncat(newdirpath, filelist[i], strlen(filelist[i]));
		
		// fprintf(stderr, "newdirpath == %s\n", newfilepath);
		// fprintf(stderr, "file == %s\n", filelist[i]);

		strcat(newdirpath, filelist[i]);

		// fprintf(stderr, "Got the string %s\n", newdirpath);

		// test extractChrNum
		// fprintf(stderr, "ChrNum == %d\n", extractChrNum(filelist[i]));
		// fprintf(stderr, "ChrLen == %d\n", Get_Chr_Len(newdirpath));
		map[extractChrNum(filelist[i])] = Get_Chr_Len(newdirpath);
		// fprintf(stderr, "map test == %d\n", (*map)[extractChrNum(filelist[i])]);
	
		memset(newdirpath, 0, sizeof(char) * NAME_MAX);	
		strcpy(newdirpath, backuppath);
		
	}

	// fprintf(stderr, "output file == %s\n", out);
	// are we printing to a file?
	/*
	if (out != NULL)
		printMap(map, out);
	*/

																		  // cleanup */
	free(newdirpath);
	free(backuppath);
	freeFileList(filelist, numfiles);
	return 0;
}

// print a map, in tabular formatting to be imported to excel
static void printMap(map<int, int>& map, char *filename){
	ofstream out;

	out.open(filename);

	// file check
	if (out.is_open() == false){
		cerr << "printMap: Error: unable to open file given.\n";
		return;
	}
																	// file header */
	int i;
	for (i = 1; i <= 22; i++){
		out << "chr" << i << "\t";
	}
	out << "chrX\tchrY\n";

														  // content of the map */
	for (i = 1; i <= 24; i++){
		// fprintf(stderr, "%d\n", (*map)[i]);
		out << map[i] << "\t";
	}

	out.close();
}
