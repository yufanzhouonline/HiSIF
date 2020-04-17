#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "lowstructs.h"

// global mutexes from main
extern pthread_mutex_t mutexes[25];
extern int nstructs;
extern int writefds[25];
extern struct regionIndex *regions[25];
extern int full_bytes;

void parseRAO(struct interRegionPair *pair, char *buf, char *saveptr);
void parseOWN(struct interRegionPair *pair, char *buf, char *saveptr);


/******************************************************************************
 * Purpose:
 * 	A thread that reads from 1 file, and writes to many. There will be other
 * 	threads doign the same, so thread synchronization is a must.
 *
 * Globals:
 * 	pthread_mutex_t mutexes[25]				a mutex for each array
 * 	int nstructs									how many structs in the array
 * 	int writefds[25]								fd's for each output file
 * 	struct regionIndex *regions[25]			the arrays themselves
 *		int full_bytes									# bytes to write for full array
 *
 * Parameters:
 * 	void *args										contains id and readfd
 *
 *	Note:
 *		**There needs to be one last array dump, once all of these threads
 *		are finished. This is handles by the calling thread, who knows that there
 *		could be information left in the arrays at thread exit.
 *
 *		strtok is NOT thread safe! Must use strtok_r
 *
 *****************************************************************************/

void *readInteractingRegionsThread(void *args)
{
	struct fragArgs *fargs = (struct fragArgs*)args;

	char *saveptr;


	fprintf(stderr, "Thread %d reading fd %d\n", fargs->id, fargs->readfd);
	struct interRegionPair pair;
	memset(&pair, 0, sizeof(struct interRegionPair));

	// written bytes, read bytes
	int wbytes = -1, rbytes, index;

	char buf[1024];

	memset(buf, 0, sizeof(buf));
		
	// main read loop
	while((rbytes = readline(fargs->readfd, buf, sizeof(buf))) > 0){
		// parse the string (traditional way)
		// parseOWN(&pair, buf, saveptr);
		parseRAO(&pair, buf, saveptr);
				
		// setup the index
		if (pair.end1.chr != pair.end2.chr)
			index = 24;
		else
			index = (pair.end1.chr) - 1;


		// attempt to gain control of array
		pthread_mutex_lock(&(mutexes[index]));

		// add to array
		memcpy(regions[index]->cur, &pair, sizeof(struct interRegionPair));

		// increment the current pointer
		(regions[index]->cur)++;

		// write entire chunk to the file, then erase the array, reset the pointer
		if (regions[index]->cur == regions[index]->end){
			// fprintf(stderr, "Thread %d is dumping array %d\n", fargs->id, index);
			regionDump(regions[index], writefds[index], full_bytes);

			regions[index]->cur = regions[index]->start;

			memset(regions[index]->cur, 0, full_bytes);
		}

		// return control
		pthread_mutex_unlock(&(mutexes[index]));
		

		// flush the buffer
		memset(buf, 0, sizeof(buf));
	}

	// check for read error
	if (rbytes == -1){
		memset(buf, 0, sizeof(buf));
		sprintf(buf, "Error: thread %d encountered a read error.\n", fargs->id);
		return NULL;
	}

	// free(pair);

	fprintf(stderr, "Thread %d finished\n", fargs->id);
	return NULL;
}

// perform parsing for RAO format
void parseRAO(struct interRegionPair *pair, char *buf, char *saveptr){
	char *tmp;
	int val;

	// skip first
	tmp = strtok_r(buf, " 	", &saveptr);
	val = atoi(strtok_r(NULL, " 	", &saveptr));
	if (val == 16)
		pair->end1.strand = 1;
	else
		pair->end1.strand = 0;

	pair->end1.chr = atoi(strtok_r(NULL, " 	", &saveptr));
	pair->end1.pos = atoi(strtok_r(NULL, " 	", &saveptr));

	// skip
	tmp = strtok_r(NULL, " 	", &saveptr);
	
	val = atoi(strtok_r(NULL, " 	", &saveptr));
	if (val == 16)
		pair->end2.strand = 1;
	else
		pair->end2.strand = 0;

	pair->end2.chr = atoi(strtok_r(NULL, " 	", &saveptr));
	pair->end2.pos = atoi(strtok_r(NULL, " 	", &saveptr));
	pair->end1.cuttingSite = 0;
	pair->end2.cuttingSite = 0;
}


// perform parsing for our own format
void parseOWN(struct interRegionPair *pair, char *buf, char *saveptr){
	pair->end1.chr = atoi(strtok_r(buf, "	 ", &saveptr));
	pair->end1.pos = atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end1.strand = (char)atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end1.cuttingSite = 0;

	pair->end2.chr = atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end2.pos = atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end2.strand = (char)atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end2.cuttingSite = 0;
}
