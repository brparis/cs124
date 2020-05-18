#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm> // only for personal calculations (for table)
 
using namespace std;


//////////////////////////////////////////////////////////////////////////////
// BEFORE SUBMITTING:                                                       //
// INPUT    ./strassen flag_of_choice dimension inputfile                   //
// OUTPUT   diagonal entries (one per line, including a trailing new line)  //
//////////////////////////////////////////////////////////////////////////////

/*

Let's consider the Strassen Setup.  First we divide our n x n dimension matrices into n/2 x n/2 submatrices:


|         |          |          |         |          |          |         |          |
|   C11   |    C12   |          |   A11   |    A12   |          |   B11   |    B12   |
|         |          |          |         |          |          |         |          |
----------------------    =     ----------------------    x     ----------------------
|         |          |          |         |          |          |         |          |
|   C21   |    C22   |          |   A21   |    A22   |          |   B21   |    B22   |
|         |          |          |         |          |          |         |          |


Using the method presented in CLRS, we create the following 10 matrices via addition/subtraction (10 n/2 x n/2 matrices ~O(n^2)):

S1 = B12 - B22
S2 = A11 + A12
S3 = A21 + A22
S4 = B21 - B11
S5 = A11 + A22
S6 = B11 + B22
S7 = A12 - A22
S8 = B21 + B22
S9 = A11 - A21
S10 = B11 + B12

Then, we recursively multiply n/2 x n/2 matrices 7 times to compute the following (sum or difference of products of A * B submatrices):

P1  =   A11 * S1    =   A11 * B12 - A11 * B22
P2  =   S2 * B22    =   A11 * B22 + A12 * B22
P3  =   S3 * B11    =   A21 * B11 + A22 * B11
P4  =   A22 * S4    =   A22 * B21 - A22 * B11
P5  =   S5 * S6     =   A11 * B11 + A11 * B22 + A22 * B11 + A22 * B22
P6  =   S7 * S8     =   A12 * B21 + A12 * B22 - A22 * B21 - A22 * B22
P7  =   S9 * S10    =   A11 * B11 + A11 * B12 - A21 * B11 - A21 * B12

Finally, we find the add and subract P matrices to construct the submatrices of C:

C11 = P5 + P4 - P2 + P6     =   A11 * B11 + A12 * B12 (inefficiently)
C12 = P1 + P2               =   A11 * B12 + A12 * B22 (inefficiently)
C21 = P3 + P4               =   A21 * B11 + A22 * B21 (inefficiently)
C22 = P5 + P1 - P3 - P7     =   A22 * B22 + A21 * B12 (inefficiently)


|         |          |          |                       |                       |
|   C11   |    C12   |          |   P5 + P4 - P2 + P6   |        P1 + P2        | 
|         |          |          |                       |                       |
----------------------    =     ----------------------- | -----------------------
|         |          |          |                       |                       |
|   C21   |    C22   |          |        P3 + P4        |   P5 + P1 - P3 - P7   | 
|         |          |          |                       |                       |

*/



//////////////////////
// helper functions //
//////////////////////

// ======================================================================================================================
// Implementation of NAIVE MATRIX MULTIPLICTION
// use this after the crossover point where recursively calling
// Strassen's to the base case would be less than ideal
//
// assumes square matrix (therefore, input dimension = output dimension)
// a * b = c
// don't necessarily have to pass by reference here
// (assumption crossover is low --> it might be faster to copy values to registers despite data > size pointer)
// i-k-j probably faster (caching) (machine dependent [pre-fetch streams, TLB], of course)
//
// ======================================================================================================================

void naivemm(vector< vector<int> > &a, vector< vector<int> > &b, vector< vector<int> > &c, int dimension) {
    // make sure matrix c is cleared out
    // optional revision: initialize a vector here with default value 0
    for (int i = 0; i < dimension;  i++) {
        for (int j = 0; j < dimension; j++) {
            c[i][j] = 0;
        }
    }

    // multiply in standard method (naive) ~O(n^3)
    for (int i = 0; i < dimension; i++) {
        for (int k = 0; k < dimension; k++) {
            for (int j = 0; j < dimension; j++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

// ======================================================================================================================
// Implementation of MATRIX ADDITION and SUBTRACTION
// for use in multiple stages of Strassen (calculating S matrices and final C matrices)
//
// assumes square matrix (therefore, input dimension = output dimension)
// a + b = c
// a - b = c
// pass by reference here
// (faster to send pointer than copy values to registers because data > size pointer)
//
// ======================================================================================================================
void addm(vector< vector<int> > &a, vector< vector<int> > &b, vector< vector<int> > &c, int dimension) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
}

void subtractm(vector< vector<int> > &a, vector< vector<int> > &b, vector< vector<int> > &c, int dimension) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
}



// ======================================================================================================================
// ZERO PADDING
// for use in initial stage of Strassen's if dimension is not a power of 2 or multiple of the crossOverSize
//
// assumes square matrix (therefore, row dimension = col dimension)
// enlarges to the nearest power of two or multiple of crossOverSize and pads unassigned elts with 0s
//
// pass by reference here
// (faster to send pointer than copy values to registers because data > size pointer)
//
// if this is a duplicate call to the function, an argument with expected new size is passed
//
// In the end, we need to remove the excess padding to return the final product
// An easy copying function limited by the original dimension can be employed here
// ======================================================================================================================

int zeroextendm(vector< vector<int> > &a, int dimension, int crossOverSize, int reSize) {
    int dimensionNew;
    if (reSize != 0) {
        dimensionNew = reSize;
    }
    else {
        dimensionNew = crossOverSize;
        while(dimensionNew < dimension) {
            dimensionNew *= 2;
        }

        int powerof2 = (int)(pow(2,(ceil(log2(dimension))))+0.5);  // added 0.5 to deal with possible rounding issues (int rounds down always)
        dimensionNew = min(dimensionNew, powerof2);
    }

    vector<int> newRow (dimensionNew, 0);

    // append zero elements to existing rows
    for (int i = 0; i < dimension; i++) {
        for (int j = dimension; j < dimensionNew; j++) {
            a[i].push_back(0);
        }
    }

    // add new rows of zeros
    for (int i = dimension; i < dimensionNew; i++) {
        a.push_back(newRow);
    }

    return dimensionNew;
}


void zeroremovem(vector< vector<int> > &c, int overSize, int destinationDimension) {
    // append zero elements to existing rows
    for (int i = 0; i < destinationDimension; i++) {
        for (int j = destinationDimension; j < overSize; j++ ) {
            c[i].pop_back();
        }
    }

    // add new rows of zeros
    for (int i = destinationDimension; i < overSize; i++) {
        c.pop_back();
    }
}





// ======================================================================================================================
// RANDOM MATRIX GENERATOR
// for use in benchmarking Strassen's performance
//
// creates square matrix of user-inputted dimension (therefore, row dimension = col dimension)
// enlarges to the nearest power of two and pads unassigned elts with 0s
//
// pass by reference here
// (faster to send pointer than copy empty values to registers because data > size pointer)
// ======================================================================================================================
void randm(vector< vector<int> > &a, int dimension) {
    random_device seed;                             // Code adapted from PA1
    static mt19937 generator(seed());               // Standard mersenne_twister_engine seeded with random_device seeder
    uniform_int_distribution<int> dis(-1, 1);   // distribution on which we randomly select

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            a[i][j] = dis(generator);
        }
    }
}





// ======================================================================================================================
// BERNOULLI ADJACENCY MATRIX GENERATOR
// for use in benchmarking triangles forumla (expected v. experimental)
//
// fills a square matrix with 1 or 0 based on Bernoulli probability
// diagonal = 0, A[i][j] = A[j][i] becausee symmetric matrix
//
// pass by reference here
// (faster to send pointer than copy empty values to registers because data > size pointer)
// ======================================================================================================================
void bernoullim(vector< vector<int> > &a, int dimension, float probability) {
    random_device seed;                             // Code adapted from PA1
    static mt19937 generator(seed());               // Standard mersenne_twister_engine seeded with random_device seeder
    bernoulli_distribution dis(probability);        // distribution on which we randomly select {0,1} with designated probability

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            // Diagonal = 0 (by definition for adj. mat)
            if (i == j) {
                a[i][j] = 0;
            }

            // cell and corresponding (complement) is unfilled, so fill
            if (i < j) {
                a[i][j] = int(dis(generator));  //bernoulli is bool, so cast to <int>
            }

            // corresponding cell has been filled, just copy
            if (i > j) {
                a[i][j] = a[j][i];
            }
            
        }
    }
}



// ============================================
// Comparing Two Matrices
// for use to validate output is correct
// ============================================
bool comparem(vector< vector<int> > &a, vector< vector<int> > &b, int dimension) {
  bool equal = true;
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      bool check = (a[i][j] == b[i][j]);
      equal = equal & check;
    }
  }
  return equal;
}


// ============================================
// Print Matrix
// for use to validate output is correct
// ============================================
void standardprintm(vector< vector<int> > &a, int dimension) {
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      printf("%d ", a[i][j]);
    }
    printf("\n");
  }
}


// ============================================
// Print Diagonal
// for auto-grader output (assignment spec)
// ============================================
void diagonalprintm(vector< vector<int> > &a, int dimension) {
  for (int i = 0; i < dimension; i++) {
    printf("%d ", a[i][i]);
    printf("\n");
  }
}






/////////////////////////////
// Strassen Implementation //
/////////////////////////////

// ======================================================================================================================
// Implementation of recursive Strassen's divide and conquer matrix multiplication algorithm
// 
//
// assumes square matrix (therefore, input dimension = output dimension)
// requires input of matrix dimension (currently must be a power of 2)
// matrices passed by reference (faster to send pointer than copy values to registers because data > size pointer)
// 
// recursively calls itself until the dimension of the input matrices are of the cross-over point size (crossOverSize)
// at which point naive matrix multiplication is employed to produce the output
//
// ======================================================================================================================

void strassenmm(vector< vector<int> > &a, vector< vector<int> > &b, vector< vector<int> > &c, int dimension, int crossOverSize) {

    // Check if dimension is of cross-over point size (crossOverSize)
    if (dimension <= crossOverSize) {
        naivemm(a, b, c, dimension);
        return;
    }

    ////////////////////////////////////////
    ////////////////////////////////////////
    // otherwise use Strassen's algorithm //
    ////////////////////////////////////////
    ////////////////////////////////////////

    dimension = dimension / 2;  // submatrix dimension (n/2) (DIVIDE & CONQUER!)



    /////////////////////////////////
    // construct submatrices (n/2) //
    /////////////////////////////////
    
    // first initialize column vector of ints with default value
    // vector<int> name(size, defaultValue)
    vector<int> column(dimension, 0);

    // now initialize the submatrices
    vector< vector<int> > A11(dimension, column), A12(dimension, column), A21(dimension, column), A22(dimension, column);
    vector< vector<int> > B11(dimension, column), B12(dimension, column), B21(dimension, column), B22(dimension, column);

    // fill the submatrices
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            A11[i][j] = a[i][j];                            // takes first n/2 rows and columns
            A12[i][j] = a[i][j + dimension];                // takes first n/2 rows and second n/2 columns
            A21[i][j] = a[i + dimension][j];                // takes second n/2 rows and first n/2 columns
            A22[i][j] = a[i + dimension][j + dimension];    // takes second n/2 rows and second n/2 columns
            B11[i][j] = b[i][j];                            // takes first n/2 rows and columns
            B12[i][j] = b[i][j + dimension];                // takes first n/2 rows and second n/2 columns
            B21[i][j] = b[i + dimension][j];                // takes second n/2 rows and first n/2 columns
            B22[i][j] = b[i + dimension][j + dimension];    // takes second n/2 rows and second n/2 columns

        }
    }


    ///////////////////////////////////////////////
    // calculate S matrices (+/- of submatrices) //
    ///////////////////////////////////////////////

    // S1 = B12 - B22
    // S2 = A11 + A12
    // S3 = A21 + A22
    // S4 = B21 - B11
    // S5 = A11 + A22
    // S6 = B11 + B22
    // S7 = A12 - A22
    // S8 = B21 + B22
    // S9 = A11 - A21
    // S10 = B11 + B12

    // first initialize column vector of ints with default value (use definition above)
    // vector<int> name(size, defaultValue)
    // vector<int> column(dimension, 0);

    // now initialize the S matrices
    vector< vector<int> >   S1(dimension, column), S2(dimension, column), S3(dimension, column), S4(dimension, column), S5(dimension, column), 
                            S6(dimension, column), S7(dimension, column), S8(dimension, column), S9(dimension, column), S10(dimension, column);

    // perform additions and subtractions
    subtractm(B12, B22, S1, dimension);
    addm(A11, A12, S2, dimension);
    addm(A21, A22, S3, dimension);
    subtractm(B21, B11, S4, dimension);
    addm(A11, A22, S5, dimension);
    addm(B11, B22, S6, dimension);
    subtractm(A12, A22, S7, dimension);
    addm(B21, B22, S8, dimension);
    subtractm(A11, A21, S9, dimension);
    addm(B11, B12, S10, dimension);



    //////////////////////////////////////////////////////////////////////////////
    // calculate P matrices (multiplication --> recursive call to strassenmm()) //
    //////////////////////////////////////////////////////////////////////////////

    // P1  =   A11 * S1
    // P2  =   S2 * B22
    // P3  =   S3 * B11
    // P4  =   A22 * S4
    // P5  =   S5 * S6
    // P6  =   S7 * S8
    // P7  =   S9 * S10

    // first initialize column vector of ints with default value (use definition above)
    // vector<int> name(size, defaultValue)
    // vector<int> column(dimension, 0);

    // now initialize the P matrices
    vector< vector<int> >   P1(dimension, column), P2(dimension, column), P3(dimension, column), P4(dimension, column), 
                            P5(dimension, column), P6(dimension, column), P7(dimension, column), P8(dimension, column);

    
    strassenmm(A11, S1, P1, dimension, crossOverSize);
    strassenmm(S2, B22, P2, dimension, crossOverSize);
    strassenmm(S3, B11, P3, dimension, crossOverSize);
    strassenmm(A22, S4, P4, dimension, crossOverSize);
    strassenmm(S5, S6, P5, dimension, crossOverSize);
    strassenmm(S7, S8, P6, dimension, crossOverSize);
    strassenmm(S9, S10, P7, dimension, crossOverSize);



    ////////////////////////////////////////////////////////
    // construct the final C submatrices (+/- P matrices) //
    ////////////////////////////////////////////////////////

    // C11 = P5 + P4 - P2 + P6
    // C12 = P1 + P2
    // C21 = P3 + P4
    // C22 = P5 + P1 - P3 - P7

    // first initialize column vector of ints with default value (use definition above)
    // vector<int> name(size, defaultValue)
    // vector<int> column(dimension, 0);

    // now initialize the P matrices
    vector< vector<int> > C11(dimension, column), C12(dimension, column), C21(dimension, column), C22(dimension, column);

    addm(P1, P2, C12, dimension);
    addm(P3, P4, C21, dimension);
    

    // ========================
    // C11 = ((P5 + P4) - P2) + P6
    // P4 := P5 + P4
    // P4 := P4 - P2
    // C11 = P4 + P6
    // ========================
    addm(P5, P4, P4, dimension);
    subtractm(P4, P2, P4, dimension);
    addm(P4, P6, C11, dimension);

    // ========================
    // C22 = ((P5 + P1) - P3) - P7
    // P5 := P5 + P1
    // P5 := P5 - P3
    // C22 = P5 - P7
    // ========================
    addm(P5, P1, P5, dimension);
    subtractm(P5, P3, P5, dimension);
    subtractm(P5, P7, C22, dimension);



    ///////////////////////////////////
    // Reconstruct PRODUCT Matrix, C //
    ///////////////////////////////////

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            c[i][j] = C11[i][j];                            // copies in first n/2 rows and columns
            c[i][j + dimension] = C12[i][j];                // copies in first n/2 rows and second n/2 columns
            c[i + dimension][j] = C21[i][j];                // copies in second n/2 rows and first n/2 columns
            c[i + dimension][j + dimension] = C22[i][j];    // copies in second n/2 rows and second n/2 columns
        }
    }

    // HOORAY!  You've finished this recursion step!
}






 
int main(int argc,char* argv[])
{

    //////////////////////////////////////////////////////////////////////////////
    // BEFORE SUBMITTING:                                                       //
    // INPUT    ./strassen flag_of_choice dimension inputfile                   //
    // OUTPUT   diagonal entries (one per line, including a trailing new line)  //
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////
    // command line parsing //
    //////////////////////////
    
    int crossOverSize = 32;
    int dimension = 0;
    bool is_this_part_three = false;
    float probability;

    // remember, Dhilan, argv[0] is "./randmst" //
    if(argc < 4) {
        printf("Too few command line arguments passed. You passed %d arguments!\n\n", argc);
        return 1;
    }
    
    else if (argc >= 4) {
        // convert arguments to integers //
        crossOverSize = atoi(argv[1]) ? atoi(argv[1]) : 32;
        dimension = atoi(argv[2]);
        
        
        // part three indicator check
        if (atoi(argv[1]) == -1) {
            crossOverSize = 140;
            dimension = 1024;
            is_this_part_three = true;
            probability = atof(argv[4]);
            cout << "Part three starting" << endl;
        }
    }

    
    vector<int> column1(dimension, 0);
    vector<int> column2(dimension, 0);
    vector<vector<int> >test1(dimension, column1), test2(dimension, column2);

    ////////////////////////////////
    // testing triangles (PART 3) //
    ////////////////////////////////

    if (is_this_part_three) {
        bernoullim(test1, dimension, probability);
        // cout << "Here's your adj mat: " << endl;
        // standardprintm(test1, dimension);
        // initialize resulting matrices
        vector<int> column(dimension, 0);
        vector< vector<int> > intermediaryResult(dimension, column), finalResult(dimension, column);
        strassenmm(test1, test1, intermediaryResult, dimension, crossOverSize);
        strassenmm(intermediaryResult, test1, finalResult, dimension, crossOverSize);

        int triangles = 0;
        for (int i = 0; i < dimension; i++) {
            triangles += finalResult[i][i];
        }
        triangles = triangles/6;

        cout << "Total triangle paths: " << triangles << " given probability p = " << probability << endl;

    }

    ///////////////////
    // RANDOM TRIALS //
    ///////////////////

    // fill with uniformly distributed random integers
    //randm(test1, dimension);
    //randm(test2, dimension);



    // FOR GRADESCOPE
    // open a file to read in values
    ifstream infile;
    infile.open(argv[3]);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            infile >> test1[i][j];
        }
    }
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            infile >> test2[i][j];
        }
    }
    infile.close();

    // cout << "input 1: " << endl;
    // standardprintm(test1, dimension);

    // cout << "\ninput 2: " << endl;
    // standardprintm(test2, dimension);

    

    // initialize a C matrix
    vector<int> column(dimension, 0);
    vector< vector<int> > C_strassen(dimension, column), C_naive(dimension, column);

    double compute_duration_strassen, compute_duration_naive;
    clock_t start_time, end_time;


    /////////////////////
    // Naive MM Timing //
    /////////////////////
    
    // cout << "Naive starting. " << endl;
    // start_time = clock();
    // naivemm(test1, test2, C_naive, dimension);
    // end_time = clock();

    // compute_duration_naive = ( (double) (end_time - start_time) ) / (double) CLOCKS_PER_SEC;
    // cout << "Naive done. " << "Duration: " << compute_duration << endl;




    /////////////////////
    // STRASSEN Timing //
    /////////////////////

    // cout << "Strassen starting. " << endl;

    // no padding necessary, naive will be called immediately
    if (dimension <= crossOverSize) {
        start_time = clock();
        strassenmm(test1, test2, C_strassen, dimension, crossOverSize);
        end_time = clock();
    }
    // power of 2, so run Strassen's
    else if (ceil(log2(dimension)) == floor(log2(dimension))) {
        start_time = clock();
        strassenmm(test1, test2, C_strassen, dimension, crossOverSize);
        end_time = clock();
    }
    // not power of 2, padding needed before applying Strassen's
    else {
        int newDimension;
        start_time = clock();
        newDimension = zeroextendm(test1, dimension, crossOverSize, 0);
        newDimension = zeroextendm(test2, dimension, crossOverSize, newDimension);
        newDimension = zeroextendm(C_strassen, dimension, crossOverSize, newDimension);
        strassenmm(test1, test2, C_strassen, newDimension, crossOverSize);
        zeroremovem(C_strassen, newDimension, dimension);
        end_time = clock();
    }

    
    compute_duration_strassen = ( (double) (end_time - start_time) ) / (double) CLOCKS_PER_SEC;

    // cout << "Strassen done. " << "Duration: " << compute_duration << endl;



   
    /////////////////////
    // MM VERIFICATION //
    /////////////////////

    // cout << "Strassen output: " << endl;
    // standardprintm(C_strassen, dimension);

    // cout << "Naive output: " << endl;
    // standardprintm(C_naive, dimension);

    // cout << "\n\nChecking correctness: " << endl;
    bool correct = comparem(C_naive, C_strassen, dimension);
    // cout << "Is it correct: " << correct << endl;

    
    
    
    ////////////////////////////////
    // FILE PRINTING (FOR TRIALS) //
    ////////////////////////////////

    // open a file in write/append mode
    // ofstream infile;
    // infile.open(argv[3]);

    // ofstream outfile;
    // outfile.open("trials.txt", ios::app);
    // outfile << compute_duration_naive << "," << compute_duration_strassen << "," << dimension << "," << crossOverSize << "," << correct << endl;
    // cout << compute_duration_naive << "," << compute_duration_strassen << "," << dimension << "," << crossOverSize << "," << correct << endl;
    
    //////////////////
    //   printing   //
    //////////////////

    // outfile << "///////////////////////////////////////////////////////////////////////////////////////" << endl;
    // outfile << "CPU Time (Avg), Max Dist, Total Distance (Avg), Numpoints, Numtrials, Dimension" << endl;
    // outfile << (double) cpu_time / (double) numtrials << endl;
    // outfile << max_on_trials <<  endl;
    // outfile << total_distance / (float) numtrials << endl;
    // outfile << numpoints << endl;
    // outfile << numtrials << endl;
    // outfile << dimension <<endl;
    // outfile << "///////////////////////////////////////////////////////////////////////////////////////" << endl;
    // outfile.close();

    // cout << total_distance / (float) numtrials << " " << numpoints << " " << numtrials << " " << dimension << endl;


    /////////////////////////////////////////////////////////////
    //   Output Diagonal Entries (1 per line) for GRADESCOPE   //
    /////////////////////////////////////////////////////////////
    diagonalprintm(C_strassen, dimension);

    return 0;
}
