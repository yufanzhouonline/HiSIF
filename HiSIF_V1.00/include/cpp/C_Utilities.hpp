#include <stdlib.h>
#include <stdio.h>
#include <utility>
#include <map>

using namespace std;

// C functions that are actually used in the C++ code
extern "C"{
	int getchrlen(char *filepath);
	char **getfilelist(char *path, int *count);
	void freeFileList(char **filelist, int num);
	void addDirEnd(char *dirpath);
	int readline(int fd, char *buf, int nbytes);
	ssize_t r_read(int fd, char *buf, size_t nbytes);
	ssize_t r_write(int fd, void *buf, size_t size);
	int mkdir_check(char *filepath);
	int rmdir_check(char *dirpath);

	int closeAllBut(int pipe[][2], int end, int size, int notThis);
	int closeAll(int pipe[][2], int end, int size);
}

// struct used to dump large dataset
struct dump{
	int one[3];
	int two[3];
	// used for the size of the vector in the map
	unsigned int size;
};
typedef struct dump prox_dump;

// struct to hold the values of the object: Frag_Info_With_Cutting_Site
struct fiwcs{
	int end1[4];
	int end2[4];
};
typedef struct fiwcs frag_info;


// used for guard processes
int guard(int, char *);
