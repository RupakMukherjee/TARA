#ifndef __CUFFT_UTILS_H__
#define __CUFFT_UTILS_H__

#include <stdio.h>
#include <stdlib.h>

#define CHECK_CUFFT(call)                                                      \
{                                                                              \
     cufftResult err;                                                          \
     if ( (err = (call)) != CUFFT_SUCCESS)                                      \
     {                                                                          \
         fprintf(stderr, "Got CUFFT error %d at %s:%d\n", err, __FILE__,        \
                 __LINE__);                                                     \
         exit(1);                                                               \
     }                                                                          \
}
#endif

