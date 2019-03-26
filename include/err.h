#pragma once
#include <stdio.h>

//error handle
#ifdef OMPI_MPI_H
#define FatalError(s)                                             \
  {                                                               \
    printf("Fatal error '%s' at %s:%d\n", s, __FILE__, __LINE__); \
    MPI_Abort(MPI_COMM_WORLD, 1);                                 \
  }
#else
#define FatalError(s)                                             \
  {                                                               \
    printf("Fatal error '%s' at %s:%d\n", s, __FILE__, __LINE__); \
    exit(1);                                                      \
  }
#endif