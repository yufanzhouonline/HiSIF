/******************************************************************************
 * Purpose:
 *		Multi-processed version of the multi-threaded fragment processing.
 * 	1 fd for every output file, 1 process who reads and writes to it,
 *		and one process who reads from each input file.
 *	
 *		Using sort -u -k 2, remove duplicates and sort based on column 2
 *
 *  Parameters:
 *  	char *indirpath				path to input files
 *  	char *outdirpath				directory for output files
 *****************************************************************************/
#include <errno.h>
#include <limits.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include <utilities.h>
#include <usp.h>

void closeAll(int *array, int count);
void closeWrite(int array[25][2], int count);
void closeRead(int array[25][2], int count);


int frag_proc(char *indirpath, char *outdirpath, char *option){

		// pipes for sorting
	int pipes[25][2];
		// files to output sorted to
	int i, writefds[25], filecount, myId;
	
	char buf[1024];

	char **filelist;

	pid_t pid;

			// pid's of sorting processes
	pid_t sortingChildren[25];

	// pipe setup
	for (i = 0; i < 25; i++){
		if (pipe(pipes[i]) == -1){
			funcErr("frag_proc", "could not create pipe", 0);
			return -1;
		}
	}

	if (!isdirectory(indirpath) || !isdirectory(outdirpath)){
		funcErr("frag_proc", "one of these parameters is not a directory", 1);
		return -1;
	}

	filelist = getfilelist(indirpath, &filecount);

	int readfds[filecount];
	pid_t readingChildren[filecount];

	// check for '/'
	char newfilepath[NAME_MAX];
	
	// open the read files
	for (i = 0; i < filecount; i++){
		memset(newfilepath, 0, NAME_MAX);
		strcpy(newfilepath, indirpath);

		addDirEnd(newfilepath);
		/*
		if (newfilepath[strlen(indirpath)-1] != '/'){
			newfilepath[strlen(indirpath)] = '/';
			newfilepath[strlen(indirpath)+1] = '\0';
		}*/

		strcat(newfilepath, filelist[i]);

		if ((readfds[i] = open(newfilepath, O_RDONLY)) == -1)
			funcErr("frag_proc", "could not open read file.", 0);
	}

	freeFileList(filelist, filecount);
	// create the reading processes
	for (i = 0; i < filecount; i++){
		if ((pid = fork()) == -1){
			funcErr("frag_proc", "unable to fork", 0);
			return -1;
		}

		// children read the file and write to the pipes
		if (pid == 0){
			int retval, readfd;

			// save only what we need, close the rest
			// don't need the read ends
			closeRead(pipes, 25);
			retval = readInteractingRegions(readfds[i], pipes, 1000, option[1]);
			// cleanup
			closeAll(readfds, filecount);
			closeWrite(pipes, 25);
			return retval;
		}
		else		// save child's pid
			readingChildren[i] = pid;
	}

	// close the write ends for all the sort processes
	closeWrite(pipes, 25);
	// close all the read files, don't need
	closeAll(readfds, filecount);

	// open output files
	for (i = 0; i < 25; i++){
		memset(newfilepath, 0, NAME_MAX);
		memset(buf, 0, sizeof(buf));

		strcpy(newfilepath, outdirpath);

		if (newfilepath[strlen(outdirpath)-1] != '/'){
			newfilepath[strlen(outdirpath)] = '/';
			newfilepath[strlen(outdirpath)+1] = '\0';
		}

		sprintf(buf, "chr%d.tmp", i+1);
		strcat(newfilepath, buf);
		// fprintf(stderr, "%s\n", newfilepath);

		if ((writefds[i] = open(newfilepath, O_CREAT|O_WRONLY|O_TRUNC, 0666)) == -1){
			funcErr("frag_proc", "could not open write file.", 0);
			return -1;
		}
		// fprintf(stderr, "Parent! %d\n", writefds[i]);
	}

	// create the sorting processes
	for (i = 0; i < 25; i++){
		if ((pid = fork()) == -1){
			funcErr("frag_proc", "could not for for the sorting processes", 0);
			return -1;
		}
	
		// children
		if (pid == 0){
				
			// need to read from this pipe
			dup2(pipes[i][0], STDIN_FILENO);
			closeRead(pipes, 25);

			
			// fprintf(stderr, "Child! %d\n", writefds[i]);
			// need to write to this file
			dup2(writefds[i], STDOUT_FILENO);

			// closeWrite(pipes, 25);
			closeAll(writefds, 25);

			execlp("sort", "sort", "-u", "-k", "2", NULL);
			funcErr("frag_proc", "Unable to exec!!!", 0);
			return -1;
		}
		else
			sortingChildren[i] = pid;
	}

	// close all the things
	closeAll(writefds, 25);
	closeRead(pipes, 25);

	// wait for all processes
	while (waitpid(-1, NULL, 0)){
		if (errno == ECHILD)
			break;
	}
}


// close all file descriptors
void closeAll(int *array, int count){
	int i;
	for (i = 0; i < count; i++){
		if (close(array[i]) == -1)
			fprintf(stderr, "Error closing fd == %d\n", array[i]);
	}
}

// close the read ends
void closeRead(int array[25][2], int count){
	int i;
	for (i = 0; i < count; i++){
		// fprintf(stderr, "Closing array[%d][0] == %d\n", i, array[i][0]);
		if (close(array[i][0]) == -1)
			fprintf(stderr, "Error closing array[%d][0] == %d\n", i, array[i][0]);
	}
}

// close the write ends
void closeWrite(int array[25][2], int count){
	int i;
	for (i = 0; i < count; i++){
		// fprintf(stderr, "Closing array[%d][1] == %d\n", i, array[i][1]);
		if (close(array[i][1]) == -1)
			fprintf(stderr, "Error closing array[%d][1] == %d\n", i, array[i][1]);
	}
}
