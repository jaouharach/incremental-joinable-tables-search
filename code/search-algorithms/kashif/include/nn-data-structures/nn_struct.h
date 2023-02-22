#include <stdio.h>
#include <stdint.h>

// #include "../../globals.h"
// #include "../dstree_node.h"
#include "../dstree_query_engine.h"

#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

// Oreder statistics tree
void * ostree_create(unsigned int size, unsigned int query_pos, unsigned long *insert_counter); 
struct query_result *  init_query_result(unsigned int query_pos);
void ostree_insert(void * ost, struct query_result *qr, unsigned long *insert_counter);
struct query_result * ostree_get_max(void * ost);
struct query_result * ostree_get_min(void * ost);
struct query_result * ostree_get(void * ost, int i);
void ostree_print(void * ost);
void ostree_destroy(void * ost);


// Min Max Heap
void* mmheap_create(unsigned long heap_size, unsigned int query_pos); 
struct query_result *mmheap_get_max(void * h);
void mmheap_pop_max(void * h);
struct query_result *mmheap_get_min(void * h);
void mmheap_pop_min(void * h);
void mmheap_insert(void * h, struct query_result *qr);
void mmheap_print(void * h);
void mmheap_destroy(void * h);


#ifdef __cplusplus
}
#endif
