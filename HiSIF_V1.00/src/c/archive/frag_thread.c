/*****************************************************************************
 * Purpose:
 * 	Given a directory, gather the list of files. Setup global array structures,
 * 	mutexes, write file descriptors, and size variables for each file. Assign 
 * 	a thread to each file, and each thread will run readInteractingRegionsThread, 
 * 	to write each read to a shared array, and upon the array being close to
 * 	a block size, write this array to the file, and flush.
 *
 * Parameters:
 * 	char *indirpath						path to directory containing frag files
 * 	char *outdirpath						path to create(or re-use) directory
 *
 * Returns:
 * 	0 on success
 * 	-1 if any errors occured
 *
 *	Notes:
 *		One thread will read one file, and attempt to write to it's respective
 *		file based upon it's chr1 value.
 *
 *		Eventually add custom file extension functionality, currently does *.tmp
 *****************************************************************************/

#if defined (__cplusplus)
extern "C"{
#endif
#include <limits.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "regionStructs.h"


// globals for threads */

															// a mutex for each output file */
pthread_mutex_t mutexes[25];
					 // how many interRegionPairs can we fit into a system block? */
int nstructs;
														// the # of bytes for a full array */
int full_bytes;
												// file descriptors for each output file */
int writefds[25];

														 // array structures for each file */
struct regionIndex *regions[25];

												 // function to get files in a directory */
char **getfilelist(char *path, int *count);


// if we want to run this as a separate thread itself?
void *frag_thread(void *indirpath, char *outdirpath){
	// gen_frag_threads as a separate thread
}



/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int gen_frag_threads(char *indirpath, char *outdirpath){
														  // used for generating filenames */
	char buf[1024];


											 // used for /directory/filename generation */
	char *fullpath = malloc(NAME_MAX * sizeof(char));

	int i, filecount, error = 0;

												 // string list of filenames, without \n */
	char **filelist;

	// ensure these are both directories */
	if (!isdirectory(indirpath)){
		error = 1;
		funcErr("gen_frag_threads", "first argument is not a directory or does not exist", 1);
	}
	if (!isdirectory(outdirpath)){
		error = 1;
		funcErr("gen_frag_threads", "second argument is not a directory or does not exist", 1);
	}

	if (error)
		return -1;

	// read the directory, open each file
	filelist = getfilelist(indirpath, &filecount);

	fprintf(stderr, "---gen_frag_threads:: Found %d files in %s---\n", filecount, indirpath);
	struct fragArgs args[filecount];


	pthread_t tids[filecount];
	pthread_attr_t tattr;
	if (pthread_attr_init(&tattr)){
		funcErr("frag_thread", "failed to create thread attribute obj", 0);
		return -1;
	}

	if (pthread_attr_setscope(&tattr, PTHREAD_SCOPE_SYSTEM)){
		funcErr("frag_thread", "failed to set scope to system");
		return -1;
	}

	
											  // file descriptors for each file to read */
	int readfds[filecount];

	// open all the files
	char *newindirpath = (char*)malloc(NAME_MAX * sizeof(char));
	char *newoutdirpath = (char*)malloc(NAME_MAX * sizeof(char));
	for (i = 0; i < filecount; i++){
		
		// make sure filepath has a '/' at the end */
		memset(newindirpath, 0, NAME_MAX * sizeof(char));
		strcpy(newindirpath, indirpath);
		if (newindirpath[strlen(indirpath) - 1] != '/'){
			// add on the new '/', as well as the null-terminator */

			newindirpath[strlen(indirpath)] = '/';
			newindirpath[strlen(indirpath)+1] = '\0';
		}
		
		// concatenate name/dir
		strcat(newindirpath, filelist[i]);

		// ensure that it has a trailing '/'
		if ((readfds[i] = open(newindirpath, O_RDONLY)) == -1){
			fprintf(stderr, "For file %s\n", filelist[i]);
			perror("Error gen_frag_threads: read error for file.\n");
			return -1;
		}
	}
	
	// gather system block size
	struct stat file;
	stat("/", &file);

										// determine block size, for efficient writing */
	nstructs = (file.st_blksize / sizeof(struct interRegionPair));
	// fprintf(stderr, "System block size == %d\n", file.st_blksize);
	
	full_bytes = nstructs * sizeof(struct interRegionPair);
	fprintf(stderr, "---gen_frag_threads:: number of structs per block == %d---\n", nstructs);
	

	// TODO create a function for this
	// allocate space for the arrays
	for (i = 0; i <=24; i++){
		regions[i] = malloc(sizeof(struct regionIndex));
		regions[i]->start = malloc(nstructs * sizeof(struct interRegionPair));
		// fprintf(stderr, "regions[%d] address == %X\n", i, regions[i]);
		memset(regions[i]->start, 0, nstructs * sizeof(struct interRegionPair));
		regions[i]->cur = regions[i]->start;
		regions[i]->end = regions[i]->start + nstructs;

	}


	// ensure there is a '/'
	memset(newoutdirpath, 0, NAME_MAX * sizeof(char));
	strcpy(newoutdirpath, outdirpath);
	if (newoutdirpath[strlen(outdirpath) - 1] != '/'){
		// add on the new '/', as well as the null-terminator */

		newoutdirpath[strlen(outdirpath)] = '/';
		newoutdirpath[strlen(outdirpath)+1] = '\0';
	}

	
	// setup output fd's, along with mutexes
	for (i = 0; i < 24; i++){
				
		memset(buf, 0, sizeof(buf));
		
		sprintf(buf, "%schr%d.tmp", newoutdirpath, i+1);

		if ((writefds[i] = open(buf, O_CREAT|O_WRONLY|O_TRUNC, 0666)) == -1){
			perror("Error gen_frag_threads: could not open file for write\n");
			return -1;
		}
		
		// mutex setup
		pthread_mutex_init(&(mutexes[i]), NULL);
	}

	sprintf(buf, "%sinstra.tmp", newoutdirpath, i+1);
	// last file setup (intra-chromosomal)
	if ((writefds[i] = open(buf, O_CREAT|O_WRONLY|O_TRUNC, 0666)) == -1){
		perror("Error gen_frag_threads: could not open the intra file\n");
		return -1;
	}

	pthread_mutex_init(&(mutexes[i]), NULL);
	
	// setup args
	for (i = 0; i < filecount; i++){
		args[i].id = i+1;
		args[i].readfd = readfds[i];
	}


	// generate threads
	for (i = 0; i < filecount; i++){
		if (pthread_create(&tids[i], &tattr, readInteractingRegionsThread, &(args[i])) < 0){
			perror("Error gen_frag_threads: thread couldn't be created\n"); 
			return -1;
		}
	}

	// wait for all threads
	for (i = 0; i < filecount; i++){
		if (pthread_join(tids[i], NULL) < 0){
				perror("Error gen_frag_threads: could not join thread\n");
				return -1;
		}
	}

	// fprintf(stderr, "Write remaining from main thread\n");
	// write remaining array to file, and free the starts
	for (i = 0; i < 25; i++){
		regionDump(regions[i], writefds[i], 0);
		// fprintf(stderr, "Start: %X\tEnd: %X\n", regions[i]->start, regions[i]->end);
		free(regions[i]->start);
		free(regions[i]);
		// fprintf(stderr, "Freed region %d\n", i);
	}

	// fprintf(stderr, "Freeing memory in main thread\n");
	// free(args);
	freeFileList(filelist, filecount);
	free(fullpath);
	free(newindirpath);
	free(newoutdirpath);

	pthread_attr_destroy(&tattr);
	return 0;

}

#if defined (__cplusplus)
}
#endif
