/***************************************************
 * Purpose:
 * 	Code for a guard process, who reads and writes
 * 	to it's own directory. This ensures concurrent
 * 	file i/o with no overwriting
 *
 * Tasks:
 * 	1) Read pipe until empty
 * 	2) For each read cBundle:
 * 		a) open the corresponding file (based on length)
 * 		b) seek to the end (not to overwrite what was there)
 * 			if not empty
 * 		c) write it's cBundle
 * 		d) close the file
 *
 * Parameters:
 * 	int readpipe - reading data from
 *		char *dir - directory where files are stored
 *
 * Note:
 * 	Performance will take a hit with all the opening and
 * 	closing, could introduce a time-to-live here to keep
 * 	fd's open longer.
 *******************************************************/
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>

// read/write functions
#include <C_Utilities.hpp>

// cBundle
#include <CPP_Utilities.hpp>

int guard(int readpipe, char *dir){
	
	cBundle data;
	char file[2048];
	int fd;
	while ((r_read(readpipe, (char *)&data, sizeof(cBundle))) > 0){
		memset(file, 0, sizeof(file));
		sprintf(file, "%s/%d", dir, data.length);
		// printf("file %s\n", file);
		if ((fd = open(file, O_CREAT | O_WRONLY, 0777)) == -1){
			fprintf(stderr, "Error: guard could not open file.\n");
			return -1;
		}

		lseek(fd, 0, SEEK_END);
		r_write(fd, (char*)&data, sizeof(cBundle));
		close(fd);
	}

	return 0;
}
