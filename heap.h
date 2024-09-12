#ifndef HEAP_H
#define HEAP_H

// Structure for a node in the heap
typedef struct {
    void* o;
    int row;
    int col;
    int distance; //from PC
    int index;    //within heap array
} HeapNode;

// Structure for the min-heap
typedef struct {
    HeapNode** array;
    int capacity;
    int size;
} MinHeap;

// Function prototypes
MinHeap* createHeap(int capacity);
int isEmpty(MinHeap* heap);
void freeHeap(MinHeap* heap);
void enqueue(MinHeap* heap, HeapNode* node);
HeapNode* dequeue(MinHeap* heap);
void modifyDistance(MinHeap* heap, HeapNode* node, int distance);
void printHeap(MinHeap* heap);

#endif /* HEAP_H */