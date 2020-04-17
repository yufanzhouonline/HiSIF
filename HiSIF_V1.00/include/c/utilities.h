/******************************************************************************
 * Contains all the c subroutines created.
 * ***************************************************************************/
#include "regionStructs.h"

#if defined (__cplusplus)
extern "C"{
#endif


// utility to print the same formatted error messages
void funcErr(char *func, char *msg, int cntl);

// make sure string has '/'
void addDirEnd(char *dirpath);

// read the file, determine the chr length (.fa files)
int getchrlen(char *filepath);

// given a string chr3*, return 3. X is 23, Y is 24
int getChrNum(char *string);

// read the binary of interRegionPair from a file
int displayInteractingRegions(int fd);

// write an array to a file
int regionDump(struct regionIndex *index, int outfd, int nbytes);

// thread function for reading interacting regions
void *readInteractingRegionsThread(void *);

// read fragments from files using threads
int gen_frag_thread(void *indirpath, char *outdirpath);

// normal function for reading interacting regions
int readInteractingRegions(int infd, int pipes[25][2], int dist, char option);

// generate processes to read all of the files given
int frag_proc(char *indirpath, char *outdirpath, char *option);

// return a char pointer to an array of the filenames
char **getfilelist(char *path, int *count);