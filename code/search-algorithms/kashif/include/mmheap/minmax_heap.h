#include <stdio.h>
#include <stdint.h>
#include "../../globals.h"
#include "../dstree_node.h"
#include "../dstree_query_engine.h"

// interface to link C++ to C code
#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

void* mmheap_create(int heap_size, unsigned int query_pos); 
struct query_result *mmheap_get_max(void * h);
struct query_result *mmheap_get_min(void * h);
void mmheap_insert(void * h, struct query_result *qr);
void mmheap_print(void * h);
void mmheap_destroy(void * h);

#ifdef __cplusplus
}
#endif
