/********************************************************************
 * Purpose:
 * 	Given a list of filenames (chr1, ..., chr25),
 * 	Re-Arrange them in proper order.
 *
 * Parameters:
 * 	char **filelist					list of chrFiles
 *
 * Assumptions:
 * 	files are 'chr1.tmp'
 *******************************************************************/
#include <string.h>
#include <limits.h>

int extractChrNum(char *str);
void chrSort(char **filelist, int size){
	int i, val;

	char buf[NAME_MAX];
	for (i = 0; i < size; i++){
		val = extractChrNum(filelist[i])-1;
		// copy destination
		strcpy(buf, filelist[val]);
	
		// move to new destination
		strcpy(filelist[val], filelist[i]);

		// copy saved to old location
		strcpy(filelist[i], buf);
	}
}
