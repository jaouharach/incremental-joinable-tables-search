#include <iostream>
#include <vector>
#include <random>
#include<float.h>

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

#include "mmheap/minmax_heap.hpp"
#include "nn_struct.h"

using namespace std;
using namespace __gnu_pbds;


/* (START) ORDER STATISTICS TREE CODE  */

typedef
tree
    <
    pair<ts_type,unsigned long>, struct query_result*,
    less<pair<ts_type,unsigned long>>,
    rb_tree_tag,tree_order_statistics_node_update
    >
ostree;

// unsigned long insert_counter = 0; // to ensure a unique key for each nn in the ostree

struct query_result * init_query_result(unsigned int query_pos)
{
    struct query_result * qr = (struct query_result *) malloc(sizeof(struct query_result));
    if(qr == NULL)
    {
        fprintf(stderr, "Error in ostree.cpp: Couldn't allocate memory for new query result\n");
        exit(1);
    }
    qr->node = NULL;
    qr->vector_id = (struct vid *) malloc(sizeof(struct vid));
    if (qr->vector_id == NULL) {
        fprintf(stderr, "Error in ostree.cpp: Couldn't allocate "
                        "memory for kNN heap element.");
        exit(1);
    }
    qr->distance = FLT_MAX;
    qr->vector_id->table_id = -1;
    qr->vector_id->set_id = -1;
    qr->vector_id->pos = -1;
    qr->query_vector_pos = query_pos;
    qr->num_checked_vectors = 0;
    qr->approx = 1;

    return qr;
}

void *ostree_create(unsigned int tree_size, unsigned int query_pos, unsigned long *insert_counter)
{
    ostree *ost  = new ostree;
     for(int i = 0; i < tree_size; i++)
    {
        struct query_result *qr = init_query_result(query_pos);
        // qr->distance = 0.4 * (i + 0.43);
        ost->insert(make_pair(make_pair(qr->distance, *insert_counter), qr));
        (*insert_counter)++;
    }
    return (void *)ost;
}

struct query_result * ostree_get_max(void * tree)
{
    ostree * ost  = (ostree *) tree;
    unsigned long  size = ost->size();
    return ost->find_by_order(size - 1)->second;
}

struct query_result * ostree_get_min(void * tree)
{
    ostree * ost  = (ostree *) tree;
    unsigned long  size = ost->size();
    return ost->find_by_order(0)->second;
}

struct query_result * ostree_get(void * tree, int i)
{
    ostree * ost  = (ostree *) tree;
    unsigned long  size = ost->size();

    if(i >= size)
    {
        fprintf(stderr, "Error in ostree.cpp: cannot get get nn at position %d, position exceeds tree size!");
        exit(1);
    }
    return ost->find_by_order(i)->second;
}

void ostree_insert(void * tree, struct query_result * qr, unsigned long *insert_counter)
{
    ostree * ost  = (ostree *) tree;
    unsigned long  size = ost->size();

    // get the key for NN with the largest distance to the query
    pair<float,int> max_key =  ost->find_by_order(size - 1)->first;
    struct query_result * max_result = ost->find_by_order(size - 1)->second;

    
    qr->query_vector_pos = max_result->query_vector_pos; // copy query info
    free(max_result->vector_id);
    ost->erase(max_key);    // delete max from ostree  then insert the new element
    
    ost->insert(make_pair(make_pair(qr->distance, *insert_counter), qr));
    (*insert_counter)++;
}

void ostree_print(void * tree)
{
    ostree * ost  = (ostree *) tree;
    int tree_size = (int) ost->size();
    struct query_result *qr;
    printf("OSTREE (size = %d):\n", tree_size);
    for(int i = 0; i < tree_size; i++)
    {
        qr = ost->find_by_order(i)->second;
        printf("[e%d] q: %u,\tv:(%u, %u, %u),\td = %f,\tt = %f\n", i, qr->query_vector_pos, qr->vector_id->table_id, qr->vector_id->set_id, qr->vector_id->pos, qr->distance, qr->time);
    }
}

void ostree_destroy(void * tree)
{
    ostree * ost  = (ostree *) tree;
    int tree_size = (int) ost->size();
    struct query_result *qr;
    for(int i = 0; i < tree_size; i++)
    {
        qr = ost->find_by_order(i)->second;
        free(qr->vector_id);
        free(qr);
    }

    delete(ost);
}

/* (END) ORDER STATISTICS TREE CODE  */


/* (START) MIN MAX  HEAP CODE  */
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
    qr->query_vector_pos = max->query_vector_pos; // keep trach of query vector pos (in query column)
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
/* (END) MIN MAX  HEAP CODE  */