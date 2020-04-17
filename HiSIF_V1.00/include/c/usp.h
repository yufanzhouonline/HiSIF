/******************************************************************************
 * Files included from the Unix Systems Programming Enviornment Textbook.
 *
 *	Source: http://usp.cs.utsa.edu/usp/
 *****************************************************************************/

#if defined (__cplusplus)
extern "C"{
#endif

int readline(int fd, char *buf, int nbytes);
ssize_t r_read(int fd, void *buf, size_t size);
ssize_t r_write(int fd, void *buf, size_t size);
int isdirectory(char *path);

#if defined (__cplusplus)
}
#endif
