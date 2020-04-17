/********************************************************************
 * Purpose:
 * 	Check for '/' and add to the end of the string.
 * 	Used when concatenating directories and filenames.
 *
 * Parameters:
 * 	char *dirpath				directory string
 *
 *
 * Note:
 * 	Assumes that there is space to add the '/'
 ********************************************************************/
#if defined (__cplusplus)
extern "C"{
#endif
#include <string.h>

void addDirEnd(char *dirpath){
	int length = strlen(dirpath);
	if (dirpath[length-1] != '/'){
		dirpath[length] = '/';
		dirpath[length+1] = '\0';
	}
}

#if defined (__cplusplus)
}
#endif
