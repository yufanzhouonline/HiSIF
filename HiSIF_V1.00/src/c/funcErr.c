/********************************************************************
 * Purpose:
 * 	Print an error message, system or otherwise.
 *
 * Parameters:
 * 	char *func						name of function where error occured
 * 	char *msg						error message
 * 	int cntl							was this a system error? 1 if yes
 *
 * Notes:
 * 	cntl == 2			this is an informational, surround with --- ---
 * 	cntl == 1			just print the error message
 * 	cntl == 0			call perror() as well as printing error message
 *******************************************************************/

#include <stdio.h>
#include <string.h>

void funcErr(char *func, char *msg, int cntl){
	char buf[1024];

	memset(&buf, 0, sizeof(buf));

	switch(cntl){
		case 2:
			sprintf(buf, "---%s:: %s---\n", func, msg);
			fprintf(stderr, buf);
			return;

		case 1:
			sprintf(buf, "Error %s:: %s\n", func, msg);
			fprintf(stderr, buf);
			return;
		case 0:
			sprintf(buf, "Error %s:: %s\n", func, msg);
			fprintf(stderr, buf);
			perror("\t>>>");
	}
	return;
}
