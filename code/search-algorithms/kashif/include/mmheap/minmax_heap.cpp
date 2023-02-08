
#include "minmax_heap.hpp"
#include "minmax_heap.h"
#include <vector>
#include <random>
#include <float.h>

#include <iostream> 
using namespace std;

// make min-max heap (heap has a fixed size)
void* mmheap_create(int heap_size, unsigned int query_pos)
{
    vector<struct query_result*>* heap = new vector<struct query_result*>;
    for(int i = 0; i < heap_size; i++)
    {
        struct query_result * qr = (struct query_result *) malloc(sizeof(struct query_result));
        if (qr == NULL)
        {
            fprintf(stderr, "Error in minmax_heap.cpp: Couldn't allocate memory for kNN heap element.");
            exit(1);
        }

        qr->node = NULL;
        qr->distance = FLT_MAX;
        qr->vector_id = (struct vid *) malloc(sizeof(struct vid));
        if (qr->vector_id == NULL)
        {
            fprintf(stderr, "Error in minmax_heap.cpp: Couldn't allocate memory for kNN heap element.");
            exit(1);
        }
        qr->vector_id->table_id = -1;
        qr->vector_id->set_id = -1;
        qr->vector_id->pos = -1;
        qr->query_vector_pos = 0;
        qr->time = 0;
        qr->num_checked_vectors = 0;
        qr->query_vector_pos = query_pos;
        qr->approx = 1;
        
        heap->push_back(qr);
    }
    return (void *) heap;
}

// insert a new element to the heap (delete max and insert new element)
void mmheap_insert(void * h, struct query_result *qr)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;

    pop_minmax_heap_max(heap->begin(), heap->end()); // put max at the end of the heap (array)

    struct query_result * max = heap->back(); // get max 
    free(max->vector_id);  // free allocated memory
    free(max);

    heap->back() = qr;
    push_minmax_heap(heap->begin(), heap->end());
}


void mmheap_print(void * h)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;
    int i = 0;
    for(struct query_result * num: *heap)
        printf("[k-%u]\t(%u, %u), d = %f, t = %f)\n", i++, num->vector_id->table_id, num->vector_id->set_id, num->distance, num->time);
}

// get min element from heap
struct query_result *mmheap_get_min(void * h)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;
    return heap->front();
}

// remove min element from heap
void mmheap_pop_min(void * h)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;
    pop_minmax_heap_min(heap->begin(), heap->end());
    heap->pop_back();
}

// get max element from heap
struct query_result *mmheap_get_max(void * h)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;
    if (heap->at(1)->distance > heap->at(2)->distance)
        return heap->at(1);
    else 
        return heap->at(2);
}

// remove max element from heap
void mm_heap_pop_max(void * h)
{
    vector<struct query_result*>* heap = (vector<struct query_result *>*) h;
    pop_minmax_heap_max(heap->begin(), heap->end());
    heap->pop_back();
}

void mmheap_destroy(void * h)
{
    vector<struct query_result *>* heap = (vector<struct query_result *>*) h;
    unsigned int heap_size = heap->size();
    struct query_result * qr = NULL;
    for (int i = 0; i < heap_size; i++)
    {
        qr = heap->at(i);
        free(qr->vector_id);
        free(qr);
    }
    delete(heap);
}