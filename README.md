# QRDecomposition
This is a version of matrix QR decomposition and inversion parallelized with pthread written in C++.

To compile do the following in root directory of the project, where MATRIX_SIZE and NUMBER_OF_THREADS you should substitute with the specific values:
>cd tmp && cmake . && make && ./main MATRIX_SIZE NUMBER_OF_THREADS
