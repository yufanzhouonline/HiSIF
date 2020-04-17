/************************************************
 * Purpose:
 * 	read and print a binary file
 ***********************************************/

#include "regionStructs.h"
#include "stdio.h"

int displayInteractingRegions(int fd){
	struct interRegionPair pair;
	int bytes, cnt = 1;
	while ((bytes = read(fd, &pair, sizeof(struct interRegionPair))) > 0){
		printf("%d\t%d\t%d\t%d\t%d\t%d\n",
			pair.end1.chr, pair.end1.pos, pair.end1.strand,
			pair.end2.chr, pair.end2.pos, pair.end2.strand);

	// printf("Bytes read == %d\nBytes wanted to read == %d\n\n", bytes, sizeof(struct interRegionPair));
		cnt++;
	}
	if (bytes == -1){
		perror("Error: read error.\n");
		return -1;
	}
	return 0;
}
