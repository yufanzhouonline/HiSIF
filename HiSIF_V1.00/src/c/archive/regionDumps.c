/******************************************************************************
 * Purpose:
 * 	Given the argument structure for the readInteractingRegions threaded
 * 	function, write each of the structure arrays to the files.
 *
 * Parameters:
 * 	struct regionIndex *index			pointer to array structure
 * 	int outfd								output file
 * 	int nbytes								bytes of entire array, 0 if chunk calc
 *
 * Notes:
 * 	if a size is passed, then simply write to the file.
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "lowstructs.h"

#if defined (__cplusplus)
extern "C"{
#endif

int regionDump(struct regionIndex *index, int outfd, int nbytes){
	// fprintf(stderr, "Entered regionDump\n");
	int i, bytes, bytes_to_write;

	// how much of this array do we need to write?
	if (nbytes == 0){
		bytes_to_write = countChunk(index) * sizeof(struct interRegionPair);
		// fprintf(stderr, "End dump, %d structs\n", bytes_to_write / sizeof(struct interRegionPair));
	}
	else{
		// fprintf(stderr, "Full dump\n");
		bytes_to_write = nbytes;
	}

	if (bytes_to_write == 0)
		return 0;

	// fprintf(stderr, "Dumping %d bytes\n", nbytes);

	// write to the file
	while ((bytes = r_write(outfd, index->start, bytes_to_write))
	< bytes_to_write);

	// fprintf(stderr, "Dumped %d bytes\n", bytes);

	if (bytes == -1){
		perror("Error regionDumps: error dumping one of the arrays to the files\n");
		return -1;
	}

	return 0;
}

// determine # of array elements to use
int countChunk(struct regionIndex *array){
	// fprintf(stderr, "CountChunk\n");
	int cnt = 0;
	
	struct interRegionPair *tmpPair;

	// fprintf(stderr, "start == %X\nend == %X\n", array->start, array->end);
	for (tmpPair = array->start; tmpPair != array->end; tmpPair++){
		if (tmpPair->end1.chr != 0)
			cnt++;
		else
			break;

	}
	return cnt;
}

#if defined (__cplusplus)
}
#endif
