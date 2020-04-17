#if defined (__cplusplus)
extern "C"{
#endif

#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

/**********************************************************************
Purpose: 
	returns a char pointer to an array of the filenames 

Error Checking:
	1) NULL parameter
	2) unable to open directory
	3) empty directory
	4) file not being a directory


Notes:
	note this does not return the <.> and <..> fd's in the directory
**********************************************************************/
static int isdirectory(char *path);


char **getfilelist(char *path, int *count) {
   struct dirent *direntp;
   DIR *dirp;

											// number of files in this directory */
	*count = 0;

	char **filelist;
																	// parameter check */
   if (path == NULL) {
		fprintf(stderr, "Error: null pathname provided\n");
      return NULL;
   }   

// is this a directory?

	if (!isdirectory(path)){
		fprintf(stderr, "Error: %s is not a directory\n", path);
		return NULL;
	}
		
// can we open this directory?
   if ((dirp = opendir(path)) == NULL) {
      perror ("Error: Failed to open directory\n");
      return NULL;
   }  

// count the number of files, to allocate proper space */

   while ((direntp = readdir(dirp)) != NULL){
      // printf("%s\n", direntp->d_name);
      (*count)++;
	}

//	(*count) -= 2;

						 				 // close the directory file to re-open */
   while ((closedir(dirp) == -1) && (errno == EINTR));


	if ((dirp = opendir(path)) == NULL) {
      perror ("Error: Failed to open directory\n");
      return NULL;
   } 

								 // allocate enough space in memory for array */
	filelist = (char**)malloc(sizeof(char*) * (*count));

	int i = 0;


														// read in the <.> and <..> */
//	direntp = readdir(dirp);	
//	direntp = readdir(dirp);	


															 // gather the filenames */
	char *tmp;
	while ((direntp = readdir(dirp)) != NULL){
		tmp = (char *)malloc(sizeof(char) * NAME_MAX);												

		strcpy(tmp, direntp->d_name);
		filelist[i] = tmp;
		// strcpy(filelist[i], tmp);
		i++;
	}

	while ((closedir(dirp) == -1) && (errno == EINTR));

   return filelist;
}

void freeFileList(char **filelist, int num){
	int i;
	for (i = 0; i < num; i++){
		free(filelist[i]);
	}

	free(filelist);
}

// is this filepath a directory?
static int isdirectory(char *path) {
   struct stat statbuf;

   if (stat(path, &statbuf) == -1)
      return 0;
   else
      return S_ISDIR(statbuf.st_mode);
}



#if defined (__cplusplus)
}
#endif
