#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>

using namespace std;


void binFighter(int selectedBin1, int selectedBin2, vector<int> &array){
    
    // find the smaller of selectedBin1 and selectedBin2
    // this is acceptable because the lower of the two will
    // correspond to a lower index in A by the search method
    // below

    int selectedBin = 0;
    if (selectedBin1 <= selectedBin2) {
        selectedBin = selectedBin1;
    }
    else {
        selectedBin = selectedBin2;
    }
    
    int iterSum = 0;
    int sizeOfA = array.size();
    int shiftIndex = 1;

    // for all non-empty bins, find "selectedBin"
    for (int i = 1; i < sizeOfA; i++) {
        iterSum += array[i];

        if (selectedBin <= iterSum) {
            shiftIndex = i;
            break;
        }
    }

    // now move a bin from this shiftIndex to the next
    if ((shiftIndex + 1) < sizeOfA) {
        array[shiftIndex] -= 1;
        array[shiftIndex + 1] += 1;
    }
    else {
        array[shiftIndex] -= 1;
        array.push_back(1);
    }

}

int balls_sim(int binsAndBalls) {
    
    // Let's begin by creating our A array having all bins with 0 balls in them
    // create 2 element array (binsAndBalls balls at index 0 and 0 at index 1 (no bins have 1 ball in it yet))
    // vector<int> name(size, defaultValue)
    vector<int> A (2, 0);
    A[0] = binsAndBalls;

    // Keep running sum over A from i = 1 to i = bins (don't include 0th index entries)
    // Could also do this by the following: (binsAndBalls - A[0])
    int nonEmptyBins = 0;

    static random_device seed;                                      // Code adapted from PA1
    static mt19937 generator(seed());                               // Standard mersenne_twister_engine seeded with random_device seeder
    static uniform_int_distribution<int> dis(1, binsAndBalls);      // distribution on which we randomly select


    for (int i = 0; i < binsAndBalls; i++) {
        int randomSelectedBin1 = dis(generator);
        int randomSelectedBin2 = dis(generator);

        // see if the Bin #s are strictly > nonEmptyBins
        // if so, place in a new bin; else, place in less loaded of corresponding existing bins

        if ((randomSelectedBin1 > nonEmptyBins) || (randomSelectedBin2 > nonEmptyBins)) {
            A[0] -= 1;
            A[1] += 1;
            nonEmptyBins += 1;
        }
        else {
            binFighter(randomSelectedBin1, randomSelectedBin2, A);
        }
    }

    int maxLoad = A.size() - 1;
    return maxLoad;
}



int main(int argc, char* argv[]) {

    if(argc != 2) {
        printf("Too few or too many command line arguments passed. You passed %d arguments!\n\n", argc);
        return 1;
    }

    int howManyBinsAndBalls = atoi(argv[1]);
    int foundMaxLoad = balls_sim(howManyBinsAndBalls);

    cout << "The maximum load given " << howManyBinsAndBalls << " bins and balls is " << foundMaxLoad << endl;
   
    ofstream outfile;
    outfile.open("trials.txt", ios::app);

    outfile << foundMaxLoad << endl;
    outfile.close();

}