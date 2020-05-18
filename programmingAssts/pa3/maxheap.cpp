/////////////////////////////////////////////////////////////
/*
    Here's a Max Heap array-based implementation good for: 
        - Build Heap                                       
        - Extract Max                                           
        - Insert           
    
    DECLARATIONS in maxheap.cpp

    By: Dhilan Ramaprasad, Harvard EECS '21
    Date: April 2020                                
*/
/////////////////////////////////////////////////////////////

#include "maxheap.h"
#include<iostream> 
#include<climits>
  
using namespace DhilHeap;
using namespace std;


// CONSTRUCTOR
MaxHeap::MaxHeap(int storageCapacity) {
    heapArray = new long long[storageCapacity];
    capacity = storageCapacity;
    size = 0;
}

// PARENT
// index of parent of node at given index i
int MaxHeap::Parent(int i) {
    return (i-1)/2;
}

// LEFT
// index of node left of node at given index i
int MaxHeap::Left(int i) {
    return (2*i + 1);
}

// RIGHT
// index of node right of node at given index i
int MaxHeap::Right(int i) {
    return (2*i + 2);
}

// Helper Function to swap arrray ELTS
void MaxHeap::swap(long long *x, long long *y) { 
    long long temp = *x; 
    *x = *y; 
    *y = temp; 
} 


///////////////////////////////////////////////////////////////
// MAX-HEAPIFY
// Given that the children of the node i in the Max-Heap H 
// are each the root of a Max-Heap, rearranges the tree rooted 
// at i to be a Max-Heap.  Call when heap property is broken.
//
// Ensure: i is the root of a Max-Heap
///////////////////////////////////////////////////////////////
void MaxHeap::MaxHeapify(int i) { 
    int l = this->Left(i); 
    int r = this->Right(i); 
    int largest; 
    if (l < this->size && this->heapArray[l] > this->heapArray[i]) {
        largest = l;
    }
    else {
        largest = i;
    }
    if (r < this->size && this->heapArray[r] > this->heapArray[largest]) {
        largest = r;
    }
    if (largest != i) { 
        this->swap(&this->heapArray[i], &this->heapArray[largest]); 
        this->MaxHeapify(largest);
    } 
} 


///////////////////////////////////////////////////////////////
// BUILD-HEAP
// Given an unordered array, makes it into a max heap.
// Begin at the leaf nodes and call max-heapify so that gradually, 
// all the elements, starting at the bottom, will follow the 
// heap property.
///////////////////////////////////////////////////////////////
void MaxHeap::BuildHeap(long long *unsortedArray, int arr_size) {
    if (arr_size > this->capacity) {
        cout << "OVERFLOW: Can't build heap.  Array too large." << endl;
        return;
    }
    for (int i = 0; i < arr_size; i++) {
        this->heapArray[i] = unsortedArray[i];
        this->size++;
    }
    for (int i = (this->size/2 - 1); i > -1; i--) {
        this->MaxHeapify(i);
    }
}


///////////////////////////////////////////////////////////////
// EXTRACT-MAX
// remove and return the maximal element (ROOT) of the
// Max-Heap H and then restore heap property via
// replacement of removed root with last leaf
// and a subsequent call to MaxHeapify()
//
// Require: H is non-empty Max-Heap
///////////////////////////////////////////////////////////////
long long MaxHeap::ExtractMax(){
    if (this->size <= 0) {
        return INT_MIN;
    }
    long long max = this->heapArray[0];
    this->heapArray[0] = this->heapArray[this->size - 1];
    this->size--;
    this->MaxHeapify(0);
    return max;
}


///////////////////////////////////////////////////////////////
// INSERT
// Add the node as a lowest leaf and then promote it upwards 
// until it is smaller than its parent.
///////////////////////////////////////////////////////////////
void MaxHeap::Insert(long long insertVal){

    if (this->size == this->capacity) {
        cout << "\nOVERFLOW: No more capacity in Max-Heap\n" << endl;
        return;
    }
    this->size++;
    this->heapArray[this->size - 1] = insertVal; // place in next available node
    int N = this->size - 1; // keep track of node currently containing new Val
    while (N != 0 && this->heapArray[this->Parent(N)] < this->heapArray[N]) {
        this->swap(&this->heapArray[this->Parent(N)], &this->heapArray[N]);
        N = this->Parent(N);
    }
}  

///////////////////////////////////////////////////////////////
// PRINT HEAP
// Prints heap in array format.  Nothing special.
// Explicitly limited by heap size.
///////////////////////////////////////////////////////////////
void MaxHeap::Print(){
    for (int i = 0; i < this->size; i++){
        cout << this->heapArray[i] << " ";
    }
    cout << endl;
}