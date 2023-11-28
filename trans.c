/* 
 * Cache Lab 2 (Part B)
 * Author: Aiden Vandekerckhove
 * Date: 4/23/23
 * trans.c - Matrix transpose B = A^T
 * 
 * This program provides several functions for transposing a matrix of integers A into another matrix B.
 * Cache optimized versions of the transpose function are avaliable for (32 x 32), (64 x 64), and (61 x 67).
 * All the optimized transpose functions use variations on the block method to minimize misses from cache 
 * conflicts on a cache with parameters: (s = 5, E = 1, b = 5).
 *  
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */

#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);
void transpose32x32(int M, int N, int A[N][M], int B[M][N]);
void transpose64x64(int M, int N, int A[N][M], int B[M][N]);
void transpose61x67(int M, int N, int A[N][M], int B[M][N]);


/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]) {
    // Delegate to the correct helper function based on matrix size:
    if (M == 32 && N == 32) {
        transpose32x32(M, N, A, B);
    } else if (M == 64 && N == 64) {
        transpose64x64(M, N, A, B);
    } else if (M == 61 && N == 67) {
        transpose61x67(M, N, A, B);
    }
}

// Helper function to transpose a 32 x 32 matrix using the block method (size 8). 
void transpose32x32(int M, int N, int A[N][M], int B[M][N]) {
    // Loop over blocks:
    for (int blockRow = 0; blockRow < N; blockRow += 8) {
        for (int blockCol = 0; blockCol < M; blockCol += 8) {
            // Transpose 8 x 8 submatrix block:
            if (blockRow != blockCol) {
                // Normal block: inverse index transpose (B[j][i] = A[i][j]).
                for (int i = blockRow; i < (blockRow + 8); i++) {
                    for (int j = blockCol; j < (blockCol + 8); j++) {
                        B[j][i] = A[i][j];
                    }
                }
            } else {
                // Diagonal block: transpose using temp variables to optimize for cache.
                for (int i = blockRow; i < (blockRow + 8); i++) {
                    int temp0 =  A[i][blockCol];
                    int temp1 =  A[i][blockCol + 1];
                    int temp2 =  A[i][blockCol + 2];
                    int temp3 =  A[i][blockCol + 3];
                    int temp4 =  A[i][blockCol + 4];
                    int temp5 =  A[i][blockCol + 5];
                    int temp6 =  A[i][blockCol + 6];
                    int temp7 =  A[i][blockCol + 7];

                    B[blockCol][i] = temp0;
                    B[blockCol + 1][i] = temp1;
                    B[blockCol + 2][i] = temp2;
                    B[blockCol + 3][i] = temp3;
                    B[blockCol + 4][i] = temp4;
                    B[blockCol + 5][i] = temp5;
                    B[blockCol + 6][i] = temp6;
                    B[blockCol + 7][i] = temp7;
                }
            }
        }
    }
}

// Helper function to transpose a 64 x 64 matrix using the block method (size 4).
void transpose64x64(int M, int N, int A[N][M], int B[M][N]) {
    // This is the best I was able to get it :(

    for (int blockRow = 0; blockRow < N; blockRow += 4) {
        for (int blockCol = 0; blockCol < M; blockCol += 4) {
            // Transpose block (4 x 4): using temp variables to optimize for cache on diagonals.
            // Tested the normal way without temp variables on non-diagonals, but this was slightly better. 
            for (int iBlock = blockRow; iBlock < (blockRow + 4); iBlock++) {
                int temp0 =  A[iBlock][blockCol];
                int temp1 =  A[iBlock][blockCol + 1];
                int temp2 =  A[iBlock][blockCol + 2];
                int temp3 =  A[iBlock][blockCol + 3];

                B[blockCol][iBlock] = temp0;
                B[blockCol + 1][iBlock] = temp1;
                B[blockCol + 2][iBlock] = temp2;
                B[blockCol + 3][iBlock] = temp3;
            }
        }
    }
}

// Helper function to transpose a 61 x 67 matrix using the block method (size 16).
void transpose61x67(int M, int N, int A[N][M], int B[M][N]) {
    for (int blockRow = 0; blockRow < N; blockRow += 16) {
        for (int blockCol = 0; blockCol < M; blockCol += 16) {
            // Transpose 16 x 16 submatrix block: inverse index transpose (B[j][i] = A[i][j]).
            for (int i = blockRow; (i < (blockRow + 16)) && (i < N); i++) {
                for (int j = blockCol; (j < (blockCol + 16)) && (j < M); j++) {
                    B[j][i] = A[i][j];
                }
            }
        }
    }
}


/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N]) {
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions() {
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N]) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

