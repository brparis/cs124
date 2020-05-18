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

void fillTracker(vector< vector<int> > &land, vector< vector<int> > &tracker, int dimension) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (land[i][j] == 1) {
                tracker[0][i]++; // increment row count
                tracker[1][j]++; // increment col count
            }
        }
    }
}

int findmin(vector< vector<int> > &tracker, int dimension) {
    int globalmin = tracker[0][0];

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < dimension; j++) {
            globalmin = (tracker[i][j] < globalmin) ? tracker[i][j] : globalmin;
        }
    }
    return globalmin;
}

void setmin(vector< vector<int> > &tracker, int dimension, int p){
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < dimension; j++) {
            tracker[i][j] = p;
        }
    }
}


void planters(vector< vector<int> > &land, vector< vector<int> > &tracker, int dimension, int p) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if ((tracker[0][i] > 0) && (tracker[0][j] > 0) && (land[i][j] == 1)){
                land[i][j] = 1;
                tracker[0][i]--; // decrement row count
                tracker[1][j]--; // decrement col count
            }
            else {
                land[i][j] = 0;
            }
        }
    }
}



// ======================================================================================================================
// RANDOM MATRIX GENERATOR
// ======================================================================================================================
void randm(vector< vector<int> > &land, int dimension) {
    random_device seed;           // Code adapted from PA1
    static mt19937 generator(seed());    // Standard mersenne_twister_engine seeded with random_device seeder
    uniform_int_distribution<int> dis(0, 1);    // distribution on which we randomly select

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            land[i][j] = dis(generator);
        }
    }
}



// ============================================
// Print Matrix
// for use to validate output is correct
// ============================================
void standardprintm(vector< vector<int> > &land, int dimension) {
  for (int i = 0; i < dimension; i++) {
    for (int j = 0; j < dimension; j++) {
      printf("%d ", land[i][j]);
    }
    printf("\n");
  }
}

void trackerprintm(vector< vector<int> > &tracker, int dimension) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < dimension; j++) {
      printf("%d ", tracker[i][j]);
    }
    printf("\n");
  }
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
    
    int dimension = 4;
    if (argc == 2){
        dimension = atoi(argv[1]);
    }

    vector<int> dimensionlong(dimension, 0);
    vector<vector<int> >land(dimension, dimensionlong), tracker(2, dimensionlong);

    // fill with uniformly distributed random integers
    randm(land, dimension);

    cout << "Original Land: " << endl;
    standardprintm(land, dimension);

    fillTracker(land, tracker, dimension);

    cout << "Tracker: " << endl;
    trackerprintm(tracker, dimension);

    int p = findmin(tracker, dimension);
    cout << "P: " << p << endl;

    setmin(tracker, dimension, p);
    cout << "Tracker, now: " << endl;
    trackerprintm(tracker, dimension);

    planters(land, tracker, dimension, p);
    cout << "Final Land: " << endl;
    standardprintm(land, dimension);
    return 0;
}
