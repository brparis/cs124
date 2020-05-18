#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <vector> 
#include <cmath>

//fib1

int fib1(int n) {
    if (n == 0) {
        return 0;
    }
    else if (n == 1) {
        return 1;
    }
    else {
        return (fib1(n-1)% 65536+fib1(n-2)% 65536)% 65536;
    }
}


//fib2

int fib2(int n) {
    int fibarray[n+1];
    fibarray[0] = 0;
    fibarray[1] = 1;

    for (int i = 2; i <= n; i++) {
        fibarray[i] = (fibarray[i-1]% 65536 + fibarray[i-2]% 65536)% 65536;
    }
    return fibarray[n];
}

//fib3

struct matrix_holder
{
    int matrix[2][2];
};

// =============================================
// Implementation of naive matrix squaring
// resulting = mat1^2
// use pass by reference
// =============================================

void naivematrixsquaring(matrix_holder input, matrix_holder &result) {
//   for (int i = 0; i < 2; i++) {
//     for (int j = 0; j < 2; j++) {
//       for (int k = 0; k < 2; k++) {
//         // idk why this was so hard to think up
//         result.matrix[i][j] = (result.matrix[i][j] + (input.matrix[j][k] * input.matrix[k][j])%65536)%65536; 
//       }
//     }
//   }
    result.matrix[0][0] = (input.matrix[0][0] % 65536 *input.matrix[0][0] % 65536 + input.matrix[0][1] % 65536 * input.matrix[1][0] % 65536) % 65536;
    result.matrix[0][1] = (input.matrix[0][0] % 65536 *input.matrix[0][1] % 65536 + input.matrix[0][1] % 65536 * input.matrix[1][1] % 65536) % 65536;
    result.matrix[1][0] = (input.matrix[1][0] % 65536 *input.matrix[0][0] % 65536 + input.matrix[1][1] % 65536 * input.matrix[1][0] % 65536) % 65536;
    result.matrix[1][1] = (input.matrix[1][0] % 65536 *input.matrix[0][1] % 65536 + input.matrix[1][1] % 65536 * input.matrix[1][1] % 65536) % 65536;
}

void naivematrixmult(matrix_holder input1, matrix_holder input2, matrix_holder &result) {
//   for (int i = 0; i < 2; i++) {
//     for (int j = 0; j < 2; j++) {
//       for (int k = 0; k < 2; k++) {
//         // idk why this was so hard to think up
//         result.matrix[i][j] = (result.matrix[i][j] + (input1.matrix[j][k] * input2.matrix[k][j])%65536)%65536; 
//       }
//     }
//   }
    result.matrix[0][0] = (input1.matrix[0][0]% 65536 *input2.matrix[0][0]% 65536 + input1.matrix[0][1]% 65536 *input2.matrix[1][0] % 65536) % 65536;
    result.matrix[0][1] = (input1.matrix[0][0]% 65536 *input2.matrix[0][1]% 65536 + input1.matrix[0][1]% 65536 *input2.matrix[1][1] % 65536) % 65536;
    result.matrix[1][0] = (input1.matrix[1][0]% 65536 *input2.matrix[0][0]% 65536 + input1.matrix[1][1]% 65536 *input2.matrix[1][0] % 65536) % 65536;
    result.matrix[1][1] = (input1.matrix[1][0]% 65536 *input2.matrix[0][1]% 65536 + input1.matrix[1][1]% 65536 *input2.matrix[1][1] % 65536) % 65536;
}


int fib3(int n) {

    // don't want to deal with these base cases

    if (n == 0) {
        return 0;
    }

    if (n == 1) {
        return 1;
    }

    int arraysize = floor(log2(n)+2);
    //printf(" n is: %d\n and arraysize is: %d\n", n, arraysize);

    matrix_holder base_matrix;
    matrix_holder final_matrix;
    matrix_holder identity_matrix;
    matrix_holder intermediary_matrix;

    // initialize matrices
    for  (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            base_matrix.matrix[i][j] = 0;
            final_matrix.matrix[i][j] = 0;
            identity_matrix.matrix[i][j] = 0;
        }
    }

    // fix the base matrix such that it is [0 1; 1 1]
    base_matrix.matrix[0][0] = 0;
    base_matrix.matrix[0][1] = 1;
    base_matrix.matrix[1][0] = 1;
    base_matrix.matrix[1][1] = 1;
    //printf("BASE MATRIX selected \n %d   %d\n %d   %d\n", base_matrix.matrix[0][0], base_matrix.matrix[0][1], base_matrix.matrix[1][0], base_matrix.matrix[1][1]);

    final_matrix.matrix[0][0] = 1;
    final_matrix.matrix[1][1] = 1;
    identity_matrix.matrix[0][0] = 1;
    identity_matrix.matrix[1][1] = 1;
    //printf("FINAL MATRIX selected \n %d   %d\n %d   %d\n", final_matrix.matrix[0][0], final_matrix.matrix[0][1], final_matrix.matrix[1][0], final_matrix.matrix[1][1]);



    matrix_holder repeated_squares_arrray[arraysize];
    repeated_squares_arrray[0] = identity_matrix; // i.e., the identity (at this point)
    repeated_squares_arrray[1] = base_matrix;

    // populate all relevant repeating squares matrices
    for (int i = 2; i < arraysize; i++) {
        //printf("matrix to be squared vector %d matrix selected \n %d   %d\n %d   %d\n", i-1, repeated_squares_arrray[i-1].matrix[0][0], repeated_squares_arrray[i-1].matrix[0][1], repeated_squares_arrray[i-1].matrix[1][0], repeated_squares_arrray[i-1].matrix[1][1]);
        naivematrixsquaring(repeated_squares_arrray[i-1], repeated_squares_arrray[i]);
        // printf("repeated squares %d matrix selected \n %d   %d\n %d   %d\n", i, repeated_squares_arrray[i].matrix[0][0], repeated_squares_arrray[i].matrix[0][1], repeated_squares_arrray[i].matrix[1][0], repeated_squares_arrray[i].matrix[1][1]);
        // printf("we've arrived.\n");
    }

    for (int i = 1; i < arraysize+1; i++) {
        // and with one and then bit shift for comparison and determining whether we use the matrix
        //printf("here's n&1: %d\n", (n&1));
        if ((n & 1) == 1) {
            naivematrixmult(final_matrix, repeated_squares_arrray[i], intermediary_matrix);
            // printf("repeated squares matrix selected \n %d   %d\n %d   %d\n", repeated_squares_arrray[i].matrix[0][0], repeated_squares_arrray[i].matrix[0][1], repeated_squares_arrray[i].matrix[1][0], repeated_squares_arrray[i].matrix[1][1]);
            // printf("intermediary_matrix entry of interest %d\n", intermediary_matrix.matrix[0][1]);
            final_matrix = intermediary_matrix;
        }
        
        n = n >> 1;
        //printf("n is now %d\n", n);
    }

    return final_matrix.matrix[0][1]% 65536;// % 65536;
}


int main(int argc, char const *argv[])
{
    clock_t fib1_start_time, fib1_end_time, fib2_start_time, fib2_end_time, fib3_start_time, fib3_end_time;
    double fib1_cpu_time, fib2_cpu_time, fib3_cpu_time;
    
    int n = 50000;

    fib3_start_time = clock();
    int fib3ret = fib3(n);
    fib3_end_time = clock();
    fib3_cpu_time = ( (double) (fib3_end_time - fib3_start_time) ) / (double) CLOCKS_PER_SEC;
    printf("FIB3, %lf\n Fib3 returns: %i\n", fib3_cpu_time, fib3ret);

    fib2_start_time = clock();
    int fib2ret = fib2(n);
    fib2_end_time = clock();
    fib2_cpu_time = ( (double) (fib2_end_time - fib2_start_time) ) / (double) CLOCKS_PER_SEC;
    printf("FIB2, %lf\n Fib2 returns: %i\n", fib2_cpu_time, fib2ret);

    fib1_start_time = clock();
    int fib1ret = fib1(n);
    fib1_end_time = clock();
    fib1_cpu_time = ( (double) (fib1_end_time - fib1_start_time) ) / (double) CLOCKS_PER_SEC;
    printf("FIB1, %lf\n Fib1 returns: %i\n", fib1_cpu_time, fib1ret);

    return 0;
}