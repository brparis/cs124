/////////////////////////////////////////////////////////////
/*
    Here's a Max Heap array-based implementation good for: 
        - Build Heap                                       
        - Extract Max                                           
        - Insert           
    
    METHODS in maxheap.cpp

    By: Dhilan Ramaprasad, Harvard EECS '21
    Date: April 2020                                
*/
/////////////////////////////////////////////////////////////

namespace DhilHeap {
    class MaxHeap {
    public:
        long long *heapArray; // pointer to int array of heap elts
        int capacity;   // originally allocated space (don't exceed)
        int size;       // current number of elts in MaxHeap
        
        // CONSTRUCTOR
        MaxHeap(int capacity);
        
        /////////////////////
        //     METHODS     //
        /////////////////////

        // PARENT
        // index of parent of node at given index i
        int Parent(int i);

        // LEFT
        // index of node left of node at given index i
        int Left(int i);

        // RIGHT
        // index of node right of node at given index i
        int Right(int i);

        // Helper Function to swap arrray ELTS
        void swap(long long *x, long long *y);

        ///////////////////////////////////////////////////////////////
        // MAX-HEAPIFY
        // Given that the children of the node i in the Max-Heap H 
        // are each the root of a Max-Heap, rearranges the tree rooted 
        // at i to be a Max-Heap.  Call when heap property is broken.
        ///////////////////////////////////////////////////////////////
        void MaxHeapify(int i);

        ///////////////////////////////////////////////////////////////
        // BUILD-HEAP
        // Given an unordered array, makes it into a max heap.
        // Begin at the leaf nodes and call max-heapify so that gradually, 
        // all the elements, starting at the bottom, will follow the 
        // heap property.
        ///////////////////////////////////////////////////////////////
        void BuildHeap(long long *unsortedArray, int size);

        ///////////////////////////////////////////////////////////////
        // EXTRACT-MAX
        // remove and return the maximal element (ROOT) of the
        // Max-Heap H and then restore heap property via
        // replacement of removed root with last leaf
        // and a subsequent call to MaxHeapify()
        ///////////////////////////////////////////////////////////////
        long long ExtractMax();

        ///////////////////////////////////////////////////////////////
        // INSERT
        // Add the node as a lowest leaf and then promote it upwards 
        // until it is smaller than its parent.
        ///////////////////////////////////////////////////////////////
        void Insert(long long insertVal);

        ///////////////////////////////////////////////////////////////
        // PRINT HEAP
        // Prints heap in array format.  Nothing special.
        // Explicitly limited by heap size.
        ///////////////////////////////////////////////////////////////
        void Print();
    };

}