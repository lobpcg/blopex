/* This code was developed by Merico Argentati, Andrew Knyazev, Ilya Lashuk and Evgueni Ovtchinnikov */

#ifndef MULTIVECTOR_VOID_FUNCTION_PROTOTYPES
#define MULTIVECTOR_VOID_FUNCTION_PROTOTYPES

#include "interpreter.h"

#ifdef __cplusplus
extern "C" {
#endif

void*
CreateCopyMultiVector( void*, int copyValues );
void
DestroyMultiVector( void* );
int
MultiVectorWidth( void* v );
int
MultiVectorHeight( void* v );
void
MultiSetMask( void *vector, int *mask );
void
CopyMultiVector( void* src, void* dest );
void
ClearMultiVector( void* );
void
SetMultiVectorRandomValues( void* v, int seed );
void
MultiInnerProd( void*, void*,
                    int gh, int h, int w, void* v );
void
MultiInnerProdDiag( void* x, void* y,
                    int* mask, int n, void* diag );
void
MultiVectorByMatrix( void*,
                   int gh, int h, int w, void* v,
                   void* );
void MultiVectorByDiagonal( void* x,
                      int* mask, int n, void* diag,
                      void* y );
void
MultiVectorAxpy( double, void*, void* );
void MatMultiVec (void * data, void * x, void * y);
void
MultiVectorPrint(  void   *x, char* tag, int limit);
int
SerialSetupInterpreter( mv_InterfaceInterpreter *i );

/*
int
MultiVectorPrint( void* x, const char* fileName );
void*
MultiVectorRead( MPI_Comm comm, void*, const char* fileName );
*/

#ifdef __cplusplus
}
#endif

#endif /* MULTIVECTOR_VOID_FUNCTION_PROTOTYPES */
