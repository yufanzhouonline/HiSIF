// close all filedescriptors in this array but the one specified
#if defined (__cplusplus)
extern "C"{
#endif

#include <unistd.h>

/*
 * Parameters:
 * 	int pipe[][2] - double array of pipes
 * 	int end - 0 for read, 1 for write, 2 for both
 * 	int size - total # of pipes
 * 	int notThis - the index to keep open
 */
int closeAllBut(int pipe[][2], int end, int size, int notThis){
	int i;
	for (i = 0; i < size; i++){
		if (i != notThis){
			switch (end){
				case 0:
					close(pipe[i][0]);
					break;
				case 1:
					close(pipe[i][1]);
					break;
				default:
					close(pipe[i][0]);
					close(pipe[i][1]);
					break;
			}
		}

	}

	return 0;
}

int closeAll(int pipe[][2], int end, int size){
	
	// close all but the first
	closeAllBut(pipe, end, size, 0);

	// close the first
	switch(end){
		case 0:
			close(pipe[0][0]);
			break;
		case 1:
			close(pipe[0][1]);
			break;
		default:
			close(pipe[0][0]);
			close(pipe[0][1]);
			break;
	}
	return 0;
}

#if defined (__cplusplus)
}
#endif
