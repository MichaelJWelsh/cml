# CML
CML is a (fully C++ compatible) header-only C matrix library designed with a focus on portability, simplicity, and efficiency. Several computational-intense functions require/are enhanced by BLAS/LAPACK. The user is able to define optional flags prior to including this library in source code to do things such as removing BLAS/LAPACK as a dependency (while simultaneously removing several functions from API), or giving all functions static storage. CML operates on type 'MATRIX' for all processes. 'MATRIX' is column major so the user can easily interface with the most popular linear-algebra-related API's. 'errno' is used for error handling extensively. 'errno' is set only in the presence of an error. The user is free to disregard 'errno' if an error occurred, but the state of the resultant matrix is undefined. 


## Usage
CML can be used in either header-only mode or implementation mode. The header-only mode is used by default when included and allows including this header in other headers and does not contain the actual implementation. 

The implementation mode requires the user to define the preprocessor macro 'CML_IMPLEMENTATION' in one .c/.cpp file before ```#includ```ing this file, e.g.:
 ```C		
#define CML_IMPLEMENTATION
#include "cml.h"
```
IMPORTANT: Everytime you include "cml.h", you have to define the same optional flags. Not doing so will lead to compile-time errors.


## Example
```C
#include <stdlib.h>
#include <stdio.h>

#define CML_IMPLEMENTATION
#include "cml.h"

void print_matrix(MATRIX *m, int call_number) {
	printf("Call %d:\nrows: %zd\ncols: %zd\n", call_number, m->rows, m->cols);
	for (size_t i = 0; i < m->rows; ++i) {
		printf("\t");
		for (size_t j = 0; j < m->cols; ++j) {
			printf("%6.1f", cml_get(m, i, j));
		}
		printf("\n\n");
	}
	printf("\n\n");
}

int main(void) {
	MATRIX *A = cml_new(3,3);
	print_matrix(A, 1);
	cml_set_all(A, 5);
	print_matrix(A, 2);
	cml_set(A, 0, 1, 99);
	print_matrix(A, 3);
	MATRIX *B = cml_identity(A->rows);
	print_matrix(B, 4);
	cml_cpy_row(B, A, 2, 0);
	print_matrix(B, 5);
	MATRIX *result = cml_new(0, 0);
	cml_mul(A, B, result);
	print_matrix(A, 6);
	print_matrix(B, 7);
	print_matrix(result, 8);
	cml_mul(A, B, NULL);
	print_matrix(A, 9);
	print_matrix(B, 10);
	

	cml_free(result);
	cml_free(B);
	cml_free(A);
	return 0;
}
```
