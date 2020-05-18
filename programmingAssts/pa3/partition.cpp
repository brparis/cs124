/////////////////////////////////////////////////////////////
/*
    Various methods to solve the NUMBER PARTITION PROBLEM
    using a deterministic heuristic: Karmarkar-Karp

    Additionally, several rrandomization methods are
    employed to attempt to improve upon the deterministic
    Karmarkar-Karp residual output.

    Note:   implementation only outputs best found residue
            not the set assignments (though a simple
            extension of the program could achieve this).

    By: Dhilan Ramaprasad, Harvard EECS '21
    Date: April 2020

    ./partition FLAG ALGORITHM INPUTFILE
    Algorithm Keys:
    0:  Karmarkar-Karp
    1:  Repeated Random
    2:  Hill Climbing
    3:  Simulated Annealing
    11: Prepartitioned Repeated Random
    12: Prepartitioned Hill Climbing
    13: Prepartitioned Simulated Annealing
*/
/////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <cmath>
#include <random>
#include <cstring>
#include "maxheap.h"    // Homemade Binary Max Heap (Build, Insert, Extract, MaxHeapify)


using namespace std;
using namespace DhilHeap;
typedef long long llong;


// INPUTFILE LENGTH (pre-defined)
// ITERATIONS FOR RANDOM ALGS (pre-defined)
#define numLines 100
#define numIterations 25000


// RANDOM NUMBER GENERATOR
random_device seed;
static mt19937 generator(seed());   // Standard mersenne_twister_engine seeded with random_device seeder


// FUNCTION DECLARATIONS - DEFINITIONS BELOW main()
llong KarmarkarKarp(llong unsortedArray[]);
llong RepeatedRandom(llong unsortedArray[]);
llong HillClimbing(llong unsortedArray[]);
llong SimulatedAnnealing(llong unsortedArray[]);
llong PrepartitionedRepeatedRandom(llong unsortedArray[]);
llong PrepartitionedHillClimbing(llong unsortedArray[]);
llong PrepartitionedSimulatedAnnealing(llong unsortedArray[]);



int main(int argc,char* argv[]) { 
    // ./partition FLAG ALGORITHM INPUTFILE
    if (argc != 4) {
        cout << "USAGE ERR: ./partition FLAG ALGORITHM INPUTFILE" << endl;
        return 1;
    }

    int selectedAlgorithm = atoi(argv[2]);
    char* inFileName = argv[3];

    // Read in numLines from INPUTFILE into unsorted array
    ifstream infile;
    infile.open(inFileName);
    llong unsortedArray[numLines];
    for (int i = 0; i < numLines; i++) {
        infile >> unsortedArray[i];
    }
    infile.close();

    llong residue;

    switch(selectedAlgorithm) {
        case 0: 
            residue = KarmarkarKarp(unsortedArray);
            break;
        case 1:  
            residue = RepeatedRandom(unsortedArray);
            break;
        case 2:  
            residue = HillClimbing(unsortedArray);
            break;
        case 3:  
            residue = SimulatedAnnealing(unsortedArray);
            break;
        case 11: 
            residue = PrepartitionedRepeatedRandom(unsortedArray);
            break;
        case 12: 
            residue = PrepartitionedHillClimbing(unsortedArray);
            break;
        case 13: 
            residue = PrepartitionedSimulatedAnnealing(unsortedArray);
            break;
        default:
            cout << "USAGE ERR: Invalid Algorithm Selection {0,1,2,3,11,12,13}." << endl;
            return 1;
    }
    cout << residue << endl;
    return 0;
}




/////////////////////////////////////////////////////////////
// KARMARKAR-KARP:
//
// Repeatedly take largest two elements and find difference
// replace the larger element by the |difference| and the
// smaller by 0.  (Two coloring to reconstruct sets.)
// Continue until only 1 elt left in set := attainable residue
/////////////////////////////////////////////////////////////
llong KarmarkarKarp(llong unsortedArray[]) {
    MaxHeap karmarkarHeap(numLines);    // heap with capacity numLines
    karmarkarHeap.BuildHeap(unsortedArray, numLines);

    while (karmarkarHeap.size > 1) {
        llong a1 = karmarkarHeap.ExtractMax();
        llong a2 = karmarkarHeap.ExtractMax();

        llong differenceMagnitude = llabs(a1-a2);

        karmarkarHeap.Insert(differenceMagnitude);
    }

    if (karmarkarHeap.size != 1) {
        cout <<  "ERROR: Removed too many elements during KK." << endl;
    }

    llong residue = karmarkarHeap.ExtractMax();
    return residue;
}



/////////////////////////////////////////////////////////////
// REPEATED RANDOM:
//
// Randomly generate solutions to the NUMBER PARTITION PROB
// (i.e., assign list of input values into two sets)
//
// Repeat for numIterations and return best found residue
/////////////////////////////////////////////////////////////
llong RepeatedRandom(llong unsortedArray[]) {
    bernoulli_distribution bernoulliDis(0.5);    // distribution on which we randomly select {False,True} with fair chance of each
    llong residue = __LONG_LONG_MAX__;
    for (int i = 0; i < numIterations; i++) {
        llong set1 = 0;
        llong set2 = 0;
        for (int j = 0; j < numLines; j++) {
            if (bernoulliDis(generator)) {
                set1 += unsortedArray[j];
            }
            else {
                set2 += unsortedArray[j];
            }
        }
        llong iterationResidue = llabs(set1 - set2);
        residue = (iterationResidue < residue) ? iterationResidue : residue;
    }
    return residue;
}




/////////////////////////////////////////////////////////////
// HILL CLIMBING:
//
// Start with a randomly generated solution
// to the NUMBER PARTITION PROBLEM
// (i.e., assign list of input values into two sets)
//
// Repeat the following for numIterations 
// and return final residue: 
//
// With probability 1/2 change one element from a set
// to the other else swap two elements;
// thereby, finding a "neighbor" solution.
// 
// If better neighbor residue, take changes; else, revert.
/////////////////////////////////////////////////////////////

// ==========================================================
//                     HELPER FUNCTIONS
// ==========================================================
llong FindResidue(bool setAssignment[], llong unsortedArray[]) {
        llong set1 = 0;
        llong set2 = 0;
        for (int i = 0; i < numLines; i++) {
            if (setAssignment[i]) {
                set1 += unsortedArray[i];
            }
            else {
                set2 += unsortedArray[i];
            }
        }
        llong residue = llabs(set1 - set2);
        return residue;
}
// ==========================================================

llong HillClimbing(llong unsortedArray[]) {
    bernoulli_distribution bernoulliDis(0.5);                       // distribution on which we randomly select {False,True} with fair chance of each
    uniform_int_distribution <int> uniformDis(0,numLines - 1);      // distribution on which we randomly select an element to flip (create neighbor solution)

    // create initial random solution
    bool setAssignment[numLines];   
    for (int i = 0; i < numLines; i++) {
        setAssignment[i] = bernoulliDis(generator);
    }

    llong residue = FindResidue(setAssignment, unsortedArray);

    for (int i = 0; i < numIterations; i++) {
        /////////////////////////////////////
        // move to a neighboring solution  //
        // 0: change set assignment of elt //
        // 1: swap two elts' assignments   //
        /////////////////////////////////////
        int changingEltOne;
        int changingEltTwo;

        bool flipOrSwap = bernoulliDis(generator);  
        if (!flipOrSwap) {
            // flip one
            changingEltOne = uniformDis(generator);
            setAssignment[changingEltOne] = !setAssignment[changingEltOne];
        }
        else {
            // swap two
            changingEltOne = uniformDis(generator);
            changingEltTwo = uniformDis(generator);
            while (changingEltOne == changingEltTwo) {
                changingEltTwo = uniformDis(generator);
            }
            setAssignment[changingEltOne] = !setAssignment[changingEltOne];
            setAssignment[changingEltTwo] = !setAssignment[changingEltTwo];
        }


        llong iterationResidue = FindResidue(setAssignment, unsortedArray);


        if (iterationResidue < residue) {
            residue = iterationResidue;
        }
        else {
            // revert move to neighbor
            if (!flipOrSwap) {
                // undo flip one
                setAssignment[changingEltOne] = !setAssignment[changingEltOne];
            }
            else {
                // undo swap two
                setAssignment[changingEltOne] = !setAssignment[changingEltOne];
                setAssignment[changingEltTwo] = !setAssignment[changingEltTwo];
            }   
        }
    }
    return residue;
}




/////////////////////////////////////////////////////////////
// SIMULATED ANNEALING:
//
// Start with a randomly generated solution
// to the NUMBER PARTITION PROBLEM
// (i.e., assign list of input values into two sets)
//
// Change one element from a set to the other with p = 1/2;
// else, swap two elements' set assignments thereby, finding 
// a "neighbor" solution.
// Calculate residue of neighbor, move to neighbor
// if better residue.  Else, move to neighbor
// based on outcome of cooling schedule.  Else, revert.
//
// Continue moving for numIterations and return final residue
/////////////////////////////////////////////////////////////

// ==========================================================
//                     HELPER FUNCTIONS
// ==========================================================
double CoolingSchedule(int iter, llong residuePrime, llong residue) {
    double T = pow(10,10) * pow(0.8, floor(iter/300));
    double probability = exp(-(residuePrime - residue)/T);
    return probability;
}
// ==========================================================

llong SimulatedAnnealing(llong unsortedArray[]) {
    bernoulli_distribution bernoulliDis(0.5);                       // distribution on which we randomly select {False,True} with fair chance of each
    uniform_int_distribution <int> uniformDis(0,numLines - 1);      // distribution on which we randomly select an element to flip (create neighbor solution)
    
    bool initialSetAssignment[numLines];  
    llong initialResidue;
    llong finalResidue;

    // INITIALIZATION
    for (int i = 0; i < numLines; i++) {
        initialSetAssignment[i] = bernoulliDis(generator);
    }

    initialResidue = FindResidue(initialSetAssignment, unsortedArray);
    finalResidue = initialResidue;

    for (int i = 0; i < numIterations; i++) {
        /////////////////////////////////////////////
        //  temporarily move to neighbor solution  //
        //  0: change set assignment of elt        //
        //  1: swap two elts' assignments          //
        /////////////////////////////////////////////
        int changingEltOne;
        int changingEltTwo;

        bool flipOrSwap = bernoulliDis(generator);  

        if (!flipOrSwap) {
            // flip one elt
            changingEltOne = uniformDis(generator);
            initialSetAssignment[changingEltOne] = !initialSetAssignment[changingEltOne];
        }
        else {
            // swap two elts
            changingEltOne = uniformDis(generator);
            changingEltTwo = uniformDis(generator);
            while (changingEltOne == changingEltTwo) {
                changingEltTwo = uniformDis(generator);
            }
            initialSetAssignment[changingEltOne] = !initialSetAssignment[changingEltTwo];
            initialSetAssignment[changingEltTwo] = !initialSetAssignment[changingEltTwo];
        }


        llong intermediaryResidue = FindResidue(initialSetAssignment, unsortedArray);


        if (intermediaryResidue < initialResidue) {
            initialResidue = intermediaryResidue;
        }
        else {
            double probability = CoolingSchedule(i, intermediaryResidue, initialResidue);
            bernoulli_distribution coolingDis(probability);
            if (coolingDis(generator)) {
                initialResidue = intermediaryResidue;
            }
            else {
                // undo changes, revert move to neighbor
                if (!flipOrSwap) {
                    // undo flip one element
                    initialSetAssignment[changingEltOne] = !initialSetAssignment[changingEltOne];
                }
                else {
                    // undo swap two elements
                    initialSetAssignment[changingEltOne] = !initialSetAssignment[changingEltTwo];
                    initialSetAssignment[changingEltTwo] = !initialSetAssignment[changingEltTwo];
                }   
            }
        }

        finalResidue = (initialResidue < finalResidue) ? initialResidue : finalResidue;
    }
    return finalResidue;
}




/////////////////////////////////////////////////////////////
// PREPARTITIONED REPEATED RANDOM:
//
// Start with a randomly generated prepartition.
// Run Karmarkar-Karp to deterministically find residue
// given randomly generated prepartition.
//
// Repeat random prepartition generation and apply
// Karmarkar-Karp for numIterations and return best residue
/////////////////////////////////////////////////////////////

llong PrepartitionedRepeatedRandom(llong unsortedArray[]) {
    uniform_int_distribution <int> uniformDis(0,numLines - 1);      // distribution on which we randomly select an element to flip (create neighbor solution)
    llong residue = __LONG_LONG_MAX__;
    llong prepartitionedArray[numLines];


    for (int i = 0; i < numIterations; i++) {
        memset(prepartitionedArray, 0, sizeof(llong[numLines]));    // initialize to 0
        for (int j = 0; j < numLines; j++) {
            prepartitionedArray[uniformDis(generator)] += unsortedArray[j];
        }
        llong iterationResidue = KarmarkarKarp(prepartitionedArray);

        residue = (iterationResidue < residue) ? iterationResidue : residue;
    }
    return residue;
}





/////////////////////////////////////////////////////////////
// PREPARTITIONED HILL CLIMBING:
//
// Start with a randomly generated prepartition.
// Run Karmarkar-Karp to find its residue.
//
// Repeat the following for numIterations 
// and return final residue: 
//
// Change one assignment in the prepartition.
// This is a "neighbor" prepartition.
// Run Karmarkar-Karp to find its residue.
//
// If better neighbor residue, take the changes
// to the prepartition; else, revert.
/////////////////////////////////////////////////////////////

llong PrepartitionedHillClimbing(llong unsortedArray[]) {
    bernoulli_distribution bernoulliDis(0.5);                       // distribution on which we randomly select {False,True} with fair chance of each
    uniform_int_distribution <int> uniformDis(0,numLines - 1);      // distribution on which we randomly select an element to flip (create neighbor solution)

    ///////////////////////////////////////
    //    Initial Random Prepartition    //
    ///////////////////////////////////////
    int prepartitionAssignment[numLines];  // store prepartition assignments
    for (int i = 0; i < numLines; i++) {
        prepartitionAssignment[i] = uniformDis(generator);
    }

    llong prepartitionedArray[numLines];
    memset(prepartitionedArray, 0, sizeof(llong[numLines]));    // initialize to 0
    for (int j = 0; j < numLines; j++) {
        prepartitionedArray[prepartitionAssignment[j]] += unsortedArray[j];
    }

    llong residue = KarmarkarKarp(prepartitionedArray);

    ////////////////////////////////
    //   Neighbor  Prepartition   //
    ////////////////////////////////
    for (int i = 0; i < numIterations; i++) {
        int changingElt = uniformDis(generator);
        int priorAssignment = prepartitionAssignment[changingElt];
        int newPrepartitionLocation = uniformDis(generator);
        while (changingElt == newPrepartitionLocation) {
            newPrepartitionLocation = uniformDis(generator);
        }
        prepartitionAssignment[changingElt] = newPrepartitionLocation;

        memset(prepartitionedArray, 0, sizeof(llong[numLines]));    // reset to 0
        for (int j = 0; j < numLines; j++) {
            prepartitionedArray[prepartitionAssignment[j]] += unsortedArray[j];
        }

        llong neighborResidue = KarmarkarKarp(prepartitionedArray);
        if (neighborResidue < residue) {
            residue = neighborResidue;
        }
        else {
            // undo move to neighbor
            prepartitionAssignment[changingElt] = priorAssignment;
        }

    }
    return residue;
}





/////////////////////////////////////////////////////////////
// PREPARTITIONED SIMULATED ANNEALING:
//
// Start with a randomly generated prepartition.
// Run Karmarkar-Karp to find its residue.
//
// Repeat the following for numIterations 
// and return final residue: 
//
// Change one assignment in the prepartition.
// This is a "neighbor" prepartition.
// Run Karmarkar-Karp to find its residue.
//
// Calculate residue of neighbor, move to neighbor
// if better residue.  Else, move to neighbor
// based on outcome of cooling schedule.  Else, revert.
//
// Continue moving for numIterations and return final residue
/////////////////////////////////////////////////////////////

llong PrepartitionedSimulatedAnnealing(llong unsortedArray[]) {
    bernoulli_distribution bernoulliDis(0.5);                       // distribution on which we randomly select {False,True} with fair chance of each
    uniform_int_distribution <int> uniformDis(0,numLines - 1);      // distribution on which we randomly select an element to flip (create neighbor solution)

    ///////////////////////////////////////
    //    Initial Random Prepartition    //
    ///////////////////////////////////////
    int prepartitionAssignment[numLines];  // S
    // int finalPrepartitionAssignment[numLines];  // S"
    for (int i = 0; i < numLines; i++) {
        prepartitionAssignment[i] = uniformDis(generator);
    }

    llong prepartitionedArray[numLines];
    memset(prepartitionedArray, 0, sizeof(llong[numLines]));    // initialize to 0
    for (int j = 0; j < numLines; j++) {
        prepartitionedArray[prepartitionAssignment[j]] += unsortedArray[j];
    }

    llong initialResidue = KarmarkarKarp(prepartitionedArray);
    llong finalResidue = initialResidue;

    ////////////////////////////////
    //   Neighbor  Prepartition   //
    ////////////////////////////////
    for (int i = 0; i < numIterations; i++) {
        int changingElt = uniformDis(generator);
        int priorAssignment = prepartitionAssignment[changingElt];
        int newPrepartitionLocation = uniformDis(generator);
        while (changingElt == newPrepartitionLocation) {
            newPrepartitionLocation = uniformDis(generator);
        }
        prepartitionAssignment[changingElt] = newPrepartitionLocation;

        memset(prepartitionedArray, 0, sizeof(llong[numLines]));    // reset to 0
        for (int j = 0; j < numLines; j++) {
            prepartitionedArray[prepartitionAssignment[j]] += unsortedArray[j];
        }

        llong neighborResidue = KarmarkarKarp(prepartitionedArray);

        if (neighborResidue < initialResidue) {
            // move to neighbor
            initialResidue = neighborResidue;
        }
        else {
            double probability = CoolingSchedule(i, neighborResidue, initialResidue);
            bernoulli_distribution coolingDis(probability);
            if (coolingDis(generator)) {
                // move to neighbor
                initialResidue = neighborResidue;
            }
            else {
                // undo move to neighbor
                prepartitionAssignment[changingElt] = priorAssignment;
            }
        }
        
        finalResidue = (initialResidue < finalResidue) ? initialResidue : finalResidue;

    }
    return finalResidue;
}