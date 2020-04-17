/********************************************************************
 * Purpose:
 * 	Given an open .fa file, parse through it and determine the
 * 	length of the chromosome. This is accomplished by getting
 * 	the size of the file (skipping the header), and subtracting
 * 	the number of newline characters (minus 2 for the last line,
 * 	as well as the header).
 *
 *	Parameters:
 *		const char *filepath 	filepath to .fa file
 *
 *
 *	Algorithm:
 *		1) Fork child process to gain the line count of the file
 *			-used for removing the # of newline characters (-2)
 *		2) Main process will skip the first line, and lseek() to
 *			the end of the file, finding the file size
 *		3) When both processes are done, calculate the difference
 *
 * Notes:
 *		There will be a pipe that allows communication between the
 *		parent and child processes.
 *
 *		If the child cannot execl, it will write -1 to the pipe, and
 *		the parent can decide what to do then.
 *
 *******************************************************************/
#if defined (__cplusplus)
extern "C"{
#endif

#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <wait.h>

#include <usp.h>

int getchrlen(char *filepath){
	int length = 0, nl_count;
														  // .fa file descriptor */
	int fafd;
	int child_pid;


										  // used to determine the file size */
	off_t start, end;
												 // used to read in first line */
	char buf[1024];
	char *tmpbuf;
														// pipe of communication */
	int fd[2];

																	 // param check */																	
	if (filepath == NULL){
		fprintf(stderr, "Error: null filepath given\n");
		return -1;
	}
	
	if (pipe(fd) == -1){
		perror("Error: could not create pipe\n");
		return -1;
	}

	if ((fafd = open(filepath, O_RDONLY)) == -1){
			perror("Error: could not open file\n");
			return -1;
	}

	// fork here, child does wc -l filepath
	if ((child_pid = fork()) == -1){
		perror("Error: could not fork\n");
		return -1;
	}

															 // child code, wc -l */
	if (child_pid == 0){
											// close the read end, don't need */
		close(fd[0]);


																		 // read end */
		dup2(fafd, STDIN_FILENO);
		close(fafd);

																		// write end */
		dup2(fd[1], STDOUT_FILENO);
		close(fd[1]);

		execl("/usr/bin/wc", "wc", "-l", filepath, NULL);

		// COULD NOT EXEC!
		perror("Error: could not exec!\n");
		int tmp = -1;
		write(fd[1], &tmp, sizeof(int));
		close(fd[1]);
		close(fafd);
		return -1;

	} else{									 // parent code, get file size */
		
										  // close the write end, don't need */
		close(fd[1]);

		// read the first line!
		if (readline(fafd, buf, sizeof(buf)) < 0){
			perror("Error: could not read the first line\n");
			return -1;
		}
	
		start = lseek(fafd, 0, SEEK_CUR);
		end = lseek(fafd, 0, SEEK_END);

		length = end - start;

		
		// wait on child, and read from pipe
		while (wait(NULL) > 0 && errno != EINTR);

		// read from the pipe, if -1 there was an error, also return -1 */
		memset(buf, sizeof(buf), 0);
		while ((r_read(fd[0], buf, sizeof(buf)) < 0)){
			perror("Error: unable to read from the pipe\n");
			return -1;
		}

		//printf("Buffer got == %s\n", buf);
		// tokenize the read
		// tmpbuf = malloc(1024 * sizeof(char));

		// setup of strtok
		// fprintf(stderr, "before tokenize\n");
		tmpbuf = strtok(buf, "	 ");

		// printf("Value! %s\n", tmpbuf);
		
		nl_count = atoi(tmpbuf);

		if (nl_count == -1){
			fprintf(stderr, "Error: something failed in the child process\n");
			return -1;
		}

		length = length - nl_count + 1;
		
	}

																			// cleanup */
	close(fafd);
	close(fd[0]);
	return length;
}

#if defined (__cplusplus)
}
#endif
