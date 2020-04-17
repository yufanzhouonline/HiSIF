#if defined (__cplusplus)
extern "C"{
#endif

#include <errno.h>
#include <unistd.h>

// read amount, but account for interrupts
ssize_t r_read(int fd, void *buf, size_t size) {
   ssize_t retval;

   while (retval = read(fd, buf, size), retval == -1 && errno == EINTR) ;
   return retval;
}

#if defined (__cplusplus)
}
#endif
