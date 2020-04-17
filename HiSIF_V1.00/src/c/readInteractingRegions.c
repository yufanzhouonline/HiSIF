/******************************************************************************
 * Purpose:
 * 	Given a file, read each line and extract the following:
 *
 * 	chr1#		strand(0/1)		pos1		chr2#		strand(0/1)		pos2
 *
 *		and store into a interRegionPair struct, then write to a fd.
 * Parameters:
 * 	const string& file						name of the file to read
 * 	int outfd									open file descriptor, could be a pipe
 *
 *
 * Returns:
 * 	0		success
 * 	-1		error
 *	Notes:
 *
 *****************************************************************************/
#if defined (__cplusplus)
extern "C"{
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include <regionStructs.h>

void parseRAO(struct interRegionPair *pair, char *buf, char *saveptr);
void parseOWN(struct interRegionPair *pair, char *buf, char *saveptr);


// option will be a t or an r, for rao or traditional
int readInteractingRegions(int infd, int pipes[25][2], int dist, char option)
{
	struct interRegionPair *pair = malloc(sizeof(struct interRegionPair));

	memset(pair, 0, sizeof(struct interRegionPair));
	int wbytes = -1, rbytes, index;

	char buf[1024], *saveptr;

	memset(buf, 0, sizeof(buf));

	while((rbytes = readline(infd, buf, sizeof(buf))) > 0){
		// fprintf(stderr, "read line!\n");
		// parse the string
		if (option == 't')
			parseOWN(pair, buf, saveptr);
		else
			parseRAO(pair, buf, saveptr);

		// fprintf(stderr, "chr# == %d\n", pair->end1.chr);
		if (pair->end1.chr != pair->end2.chr)
			index = 24;
		else
			index = pair->end1.chr-1;

		// fprintf(stderr, "%d\n", index);



		// distance filter
		if (abs(pair->end1.pos - pair->end2.pos) > dist){
			// while ((r_write(pipes[index][1], pair, sizeof(struct interRegionPair))) < sizeof(struct interRegionPair));
			dprintf(pipes[index][1], "%d\t%d\t%d\t%d\t%d\t%d\t\n", pair->end1.chr, pair->end1.pos, pair->end1.strand,
				pair->end2.chr, pair->end2.pos, pair->end2.strand);

			// check for write error
			if (rbytes == -1){
				perror("Error: could not write to the pipe.\n");
				return -1;
			}
			memset(buf, 0, sizeof(buf));
		}
	}

	if (rbytes == -1){
		funcErr("readInteractingRegions", "read error.", 0);
		return -1;
	}

	free(pair);
}

// perform parsing for RAO format
void parseRAO(struct interRegionPair *pair, char *buf, char *saveptr){
	char *tmp;
	int val;

	// skip first
	tmp = strtok_r(buf, " 	", &saveptr);
	
	// strand
	val = atoi(strtok_r(NULL, " 	", &saveptr));
	if (val == 16)
		pair->end1.strand = 1;
	else
		pair->end1.strand = 0;

	
	pair->end1.chr = getChrNum(strtok_r(NULL, " 	", &saveptr));
	pair->end1.pos = atoi(strtok_r(NULL, " 	", &saveptr));

	// skip
	tmp = strtok_r(NULL, " 	", &saveptr);
	
	val = atoi(strtok_r(NULL, " 	", &saveptr));
	if (val == 16)
		pair->end2.strand = 1;
	else
		pair->end2.strand = 0;

	pair->end2.chr = getChrNum(strtok_r(NULL, " 	", &saveptr));
	pair->end2.pos = atoi(strtok_r(NULL, " 	", &saveptr));
	pair->end1.cuttingSite = 0;
	pair->end2.cuttingSite = 0;
}


// perform parsing for our own format
void parseOWN(struct interRegionPair *pair, char *buf, char *saveptr){

	pair->end1.chr = getChrNum(strtok_r(buf, " 	", &saveptr));
	pair->end1.pos = atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end1.strand = (char)atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end1.cuttingSite = 0;


	pair->end2.chr = getChrNum(strtok_r(NULL, " 	", &saveptr));
	pair->end2.pos = atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end2.strand = (char)atoi(strtok_r(NULL, "	 ", &saveptr));
	pair->end2.cuttingSite = 0;
}


#if defined (__cplusplus)
}
#endif
