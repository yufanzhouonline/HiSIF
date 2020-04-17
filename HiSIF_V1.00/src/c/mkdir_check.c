#if defined (__cplusplus)
extern "C"{
#endif

// make a directory, if it does not already exist
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int mkdir_check(char *filepath){
	// does it exist?
	DIR *test = opendir(filepath);
	if (test){
		return 0;
	}

	// create directory
	return mkdir(filepath, 0777);
}


// remove this directory if it exists
int rmdir_check(char *filepath){
	// does it exist?
	DIR *test = opendir(filepath);
	if (test){
		char cmd[1024];
		memset(cmd, 0, sizeof(cmd));
		sprintf(cmd, "rm -r %s", filepath);
		system(cmd);
	}
	return 0;
}

#if defined (__cplusplus)
}
#endif
