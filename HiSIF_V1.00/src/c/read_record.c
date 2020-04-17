#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

// for r_read function
#include "usp.h"


/******************************************************************************
 * Purpose:
 * 	Given an open file fd, a record number, and a size of a record,
 * 	seek through fd to the record location and return the value.
 *
 * Parameters:
 * 	int fd							file descriptor of records
 * 	unsigned int recordNum		which record to access
 *		unsinged int recordSize 	how many bytes in 1 record?
 *		void* record					pointer to location to return
 *
 *
 *
 * Notes:
 * 	Assumes fd is open already for read
 * 	To use later, requires a cast to what you're reading
 *****************************************************************************/
int read_record(int fd, unsigned int recordNum, unsigned int recordSize, void* record){
	if (lseek(fd, recordNum * recordSize, SEEK_SET) < 0){
		fprintf(stderr, "read_record: failed to lseek the file.\n");
		return -1;
	}

	// now read the record into *record
	int nbytes = r_read(fd, record, recordSize);

	if (nbytes == -1){
		fprintf(stderr, "read_record: failed to read from fd.\n");
		return -1;
	}
}


/******************************************************************************
 * Purpose:
 * 	Given an open file fd, a record number, and a size of a record,
 * 	seek through fd to the record location and write the given value.
 *
 * Parameters:
 * 	int fd							file descriptor of records
 * 	unsigned int recordNum		which record to access
 *		unsinged int recordSize 	how many bytes in 1 record?
 *		void* record					pointer to location to return
 *
 *
 *
 * Notes:
 * 	Assumes fd is open already for write
 * 	To use later, requires a cast to what you're reading
 *****************************************************************************/
int write_record(int fd, unsigned int recordNum, unsigned int recordSize, void* record){
	if (lseek(fd, recordNum * recordSize, SEEK_SET) < 0){
		fprintf(stderr, "write_record: failed to lseek the file.\n");
		return -1;
	}

	// now read the record into *record
	int nbytes = r_write(fd, record, recordSize);

	if (nbytes == -1){
		fprintf(stderr, "write_record: failed to write from fd.\n");
		return -1;
	}
}
