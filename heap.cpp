#include <stdio.h>
#include <stdlib.h>
#include <ncurses.h>
#include "heap.h"

// Function prototypes
MinHeap* createHeap(int capacity)
{
    MinHeap* heap = (MinHeap*)malloc(sizeof(MinHeap));
    HeapNode** arr = (HeapNode**)malloc(capacity*sizeof(HeapNode*));
    heap->size = 0;
    heap->capacity = capacity;
    heap->array = arr;

    return heap;
}

void freeHeap(MinHeap* heap)
{
    free(heap->array);
    free(heap);
}

void swap(MinHeap* heap, HeapNode* a, HeapNode* b)
{
    int aIndex = a->index;
    int bIndex = b->index;

    heap->array[aIndex] = b;
    heap->array[bIndex] = a;

    a->index = bIndex;
    b->index = aIndex;
}

void heapifyUp(MinHeap* heap, int index)
{
    if(index == 0) {return;}

    HeapNode* curNode = heap->array[index];
    HeapNode* parent = heap->array[(index-1)/2];
    while(curNode->distance < parent->distance)
    {
        swap(heap, curNode, parent);
        index = (index-1)/2;

        if(index == 0) {break;}

        curNode = heap->array[index];
        parent = heap->array[(index-1)/2];
    }
}

void heapifyDown(MinHeap* heap, int index)
{
    HeapNode* curNode = heap->array[index];
    int leftIndex = 2*index+1;
    int rightIndex = 2*index+2;
    if(rightIndex < heap->size)
    {
        HeapNode* lChild = heap->array[leftIndex];
        HeapNode* rChild = heap->array[rightIndex];
        HeapNode* smallChild = rChild;
        if(lChild->distance < rChild->distance)
        {
            smallChild = lChild;
        }

        if(smallChild->distance < curNode->distance)
        {
            swap(heap, smallChild, curNode);
            heapifyDown(heap, curNode->index);
        }
    }
    else if(leftIndex < heap->size)
    {
        HeapNode* lChild = heap->array[leftIndex];
        if(lChild->distance < curNode->distance)
        {
            swap(heap, curNode, lChild);
            heapifyDown(heap, curNode->index);
        }
    }
}

void enqueue(MinHeap* heap, HeapNode* node)
{
    node->index = heap->size;
    heap->array[heap->size] = node;
    heapifyUp(heap, heap->size);
    heap->size++;
}

HeapNode* dequeue(MinHeap* heap)
{
    HeapNode* dequeued = heap->array[0];
    heap->array[0] = heap->array[heap->size-1];
    heap->array[0]->index = 0;
    //heap->array[heap->size-1] == NULL;
    heap->size--;
    heapifyDown(heap, 0);

    return dequeued;
}

void modifyDistance(MinHeap* heap, HeapNode* node, int distance)
{
    node->distance = distance;
    heapifyUp(heap, node->index);
}

int isEmpty(MinHeap* heap)
{
    return heap->size;
}

void printHeap(MinHeap* heap)
{
    printw("Size: %d \nContents: ", heap->size);
    for(int i = 0; i < heap->size; i++)
    {
        //HeapNode* curNode = heap->array[i];
        printw("%d, ", heap->array[i]->distance);
    }
    printw("\n");
}