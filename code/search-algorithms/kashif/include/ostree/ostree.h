#include <stdio.h>
#include <stdint.h>
#include "../../globals.h"
#include "../dstree_node.h"
#include "../dstree_query_engine.h"

#pragma once
#ifdef __cplusplus
extern "C"
{
#endif

void * ostree_create(unsigned int size, unsigned int query_pos, unsigned long *insert_counter); 
struct query_result *  init_query_result(unsigned int query_pos);

void ostree_insert(void * ost, struct query_result *qr, unsigned long *insert_counter);

struct query_result * ostree_get_max(void * ost);
struct query_result * ostree_get_min(void * ost);
struct query_result * ostree_get(void * ost, int i);

void ostree_print(void * ost);

void ostree_destroy(void * ost);

#ifdef __cplusplus
}
#endif
