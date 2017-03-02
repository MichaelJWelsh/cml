# CML
CML is a (fully C++ compatible) header-only C matrix library designed with a focus on portability, simplicity, and efficiency. Several computational-intense functions require/are enhanced by BLAS/LAPACK. The user is able to define optional flags prior to including this library in source code to do things such as removing BLAS/LAPACK as a dependency (while simultaneously removing several functions from API), or giving all functions static storage. CML operates on type 'MATRIX' for all processes. 'MATRIX' is column major so the user can easily interface with the most popular linear-algebra-related API's. 'errno' is used for error handling extensively. 'errno' is set only in the presence of an error. The user is free to disregard 'errno' if an error occurred, but the state of the resultant matrix is undefined. 


## Usage
CML can be used in either header-only mode or implementation mode. The header-only mode is used by default when included and allows including this header in other headers and does not contain the actual implementation. 

The implementation mode requires the user to define the preprocessor macro 'CML_IMPLEMENTATION' in one .c/.cpp file before #including this file, e.g.:
 ```C		
#define CML_IMPLEMENTATION
#include "cml.h"
```
IMPORTANT: Everytime you include "cml.h", you have to define the same optional flags. Not doing so will lead to compile-time errors.
