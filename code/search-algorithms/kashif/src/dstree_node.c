//
//  dstree_node.c
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#include "../include/dstree_node.h"
#include "../config.h"
#include "../globals.h"
#include "../include/dstree_index.h"
#include "../include/dstree_query_engine.h"
#include "../include/pqueue.h"
#include "../include/nn-data-structures/nn_struct.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <pqueue.h>
#include <stdio.h>
#include <stdlib.h>

/**
 This function initializes a dstree root node.
 */

struct dstree_node *
dstree_root_node_init(struct dstree_index_settings *settings) {

  struct dstree_node *node = dstree_leaf_node_init(settings);
  int ts_length = settings->timeseries_size;
  int segment_size = settings->init_segments;

  node->node_segment_split_policies =
      malloc(sizeof(struct node_segment_split_policy) * 2);

  if (node->node_segment_split_policies == NULL) {
    fprintf(
        stderr,
        "Error in dstree_node.c: Could not allocate memory for root node segment \
                     split policies.\n");
    return FAILURE;
  }

  node->node_segment_split_policies[0].indicator_split_idx = 0; // Mean
  node->node_segment_split_policies[1].indicator_split_idx = 1; // Stdev
  node->num_node_segment_split_policies = 2; // Stdev and Mean

  // calc the split points by segmentSize
  short *split_points = NULL;
  split_points = malloc(sizeof(short) * segment_size);

  if (split_points == NULL) {
    fprintf(stderr, "Error in dstree_node.c: Could not \
                     allocate memory for root split points.\n");
    return FAILURE;
  }

  if (!calc_split_points(split_points, ts_length, segment_size)) {
    fprintf(stderr, "Error in dstree_node.c: Could not \
                     calculate the split points for the root.\n");
    return FAILURE;
  }

  if (!node_init_segments(node, split_points, segment_size)) {
    fprintf(stderr, "Error in dstree_node.c: Could not \
                     initialize the segments for the root.\n");
    return FAILURE;
  }

  if (!create_dstree_node_filename(settings, node, NULL)) {
    fprintf(stderr, "Error in dstree_node.c: Could not \
                     create a filename for the root node.\n");
    return FAILURE;
  }

  free(split_points);
  return node;
}

/**
 This function initalizes a dstree leaf node.
 */

struct dstree_node *
dstree_leaf_node_init(struct dstree_index_settings *settings) {
  COUNT_NEW_NODE

  struct dstree_node *node = malloc(sizeof(struct dstree_node));
  if (node == NULL) {
    fprintf(stderr, "error: could not allocate memory for new node.\n");
    return NULL;
  }
  node->right_child = NULL;
  node->left_child = NULL;
  node->parent = NULL;

  node->filename = NULL;

  node->node_segment_split_policies = NULL;
  node->num_node_segment_split_policies = 2;

  node->range = 0;

  node->level = 0;
  node->is_left = false;
  node->is_leaf = true;

  node->split_policy = NULL;
  node->node_points = NULL;
  node->hs_node_points = NULL;
  node->num_node_points = 0;
  node->num_hs_node_points = 0;

  node->node_segment_sketches = NULL;
  node->hs_node_segment_sketches = NULL;

  node->max_segment_length = 2;
  node->max_value_length = 10;

  node->file_buffer = NULL;

  node->node_size = 0;

  // node->gt = calloc(settings->max_leaf_size, sizeof(unsigned char));
  node->gt = calloc(settings->max_leaf_size, sizeof(label_type));
  node->fp = calloc(settings->max_leaf_size, sizeof(unsigned int));

  /* start kashif changes */
  node->vid = calloc(settings->max_leaf_size, sizeof(struct vid));

  if (pthread_mutex_init(&node->lock, NULL) != 0) 
  {
		fprintf(stderr, "error: could not initialize mutex lock.\n");
    return NULL;
	}
  
  /* end kashif changes */
  return node;
}

enum response node_init_segments(struct dstree_node *node, short *split_points,
                                 int segment_size) {

  node->num_node_points = segment_size;

  node->node_points = NULL;
  node->node_points = malloc(sizeof(short) * segment_size);
  if (node->node_points == NULL) {
    fprintf(stderr, "Error in node_init_segments(): Could not allocate memory "
                    "for node points.\n");
    return FAILURE;
  }

  for (int i = 0; i < segment_size; ++i) {
    node->node_points[i] = split_points[i];
  }

  node->hs_node_points = NULL;
  node->hs_node_points = malloc(sizeof(short) * segment_size * 2);
  if (node->hs_node_points == NULL) {
    fprintf(stderr, "Error in node_init_segments(): Could not allocate memory "
                    "for hs node points.\n");
    return FAILURE;
  }

  int min_length = 1; // mininum length of new segment = 1
  int num_hs_split_points = 0;

  calc_hs_split_points(node->hs_node_points, &node->num_hs_node_points,
                       node->node_points, segment_size, min_length);

  // allocate mem for array of sketches

  node->node_segment_sketches = NULL;
  node->node_segment_sketches =
      malloc(sizeof(struct segment_sketch) * segment_size);
  if (node->node_segment_sketches == NULL) {
    fprintf(stderr, "Error in node_init_segments(): Could not allocate memory "
                    "for node segment sketches.\n");
    return FAILURE;
  }

  node->hs_node_segment_sketches = NULL;
  node->hs_node_segment_sketches =
      malloc(sizeof(struct segment_sketch) * node->num_hs_node_points);
  if (node->hs_node_segment_sketches == NULL) {
    fprintf(stderr, "Error in node_init_segments(): Could not allocate memory "
                    "for hs node segment sketches.\n");
    return FAILURE;
  }

  // allocate memory for vertical indicators
  for (int i = 0; i < segment_size; ++i) {
    node->node_segment_sketches[i].indicators = NULL;
    node->node_segment_sketches[i].indicators =
        malloc(sizeof(ts_type) * 4); // node segment has 4 indicators

    if (node->node_segment_sketches[i].indicators == NULL) {
      fprintf(
          stderr,
          "Error in node_init_segments(): Could not allocate memory for node segment \
                     sketch indicator.\n");
      return FAILURE;
    }

    node->node_segment_sketches[i].indicators[0] = -FLT_MAX;
    node->node_segment_sketches[i].indicators[1] = FLT_MAX;
    node->node_segment_sketches[i].indicators[2] = -FLT_MAX;
    node->node_segment_sketches[i].indicators[3] = FLT_MAX;
    node->node_segment_sketches[i].num_indicators = 4;
  }

  // allocate memory for horizontal indicators
  for (int i = 0; i < node->num_hs_node_points; ++i) {
    node->hs_node_segment_sketches[i].indicators = NULL;
    node->hs_node_segment_sketches[i].indicators =
        malloc(sizeof(ts_type) * 4); // node segment has 4 indicators
    if (node->hs_node_segment_sketches[i].indicators == NULL) {
      fprintf(
          stderr,
          "Error in node_init_segments(): Could not allocate memory for hs node segment \
                     sketch indicator.\n");
      return FAILURE;
    }
    node->hs_node_segment_sketches[i].indicators[0] = -FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[1] = FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[2] = -FLT_MAX;
    node->hs_node_segment_sketches[i].indicators[3] = FLT_MAX;
    node->hs_node_segment_sketches[i].num_indicators = 4;
  }

  return SUCCESS;
}

enum response append_ts_to_node(struct dstree_index *index,
                                struct dstree_node *node, ts_type *timeseries) {

  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }
  /*
  node->file_buffer->buffered_list[idx] = NULL;
  node->file_buffer->buffered_list[idx] = malloc(sizeof(ts_type) * ts_length);
  */

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < ts_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }

  ++node->file_buffer->buffered_list_size;
  index->buffer_manager->current_count += ts_length;

  return SUCCESS;
}
/*
enum response append_ts_gt_to_node(struct dstree_index * index,
                                      struct dstree_node * node,
                                      ts_type * timeseries,
                                      unsigned char gt)
*/
enum response
append_ts_gt_to_node(struct dstree_index *index, struct dstree_node *node,
                     ts_type *timeseries, label_type gt, unsigned int fp)

{

  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }
  /*
  node->file_buffer->buffered_list[idx] = NULL;
  node->file_buffer->buffered_list[idx] = malloc(sizeof(ts_type) * ts_length);
  */

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < ts_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }

  // add the timeseries gt
  node->gt[(node->node_size - 1)] = gt;

  if (fp != -1)
    node->fp[node->node_size - 1] = fp;

  ++node->file_buffer->buffered_list_size;
  index->buffer_manager->current_count += ts_length;

  return SUCCESS;
}

enum response append_ts_to_child_node(struct dstree_index *index,
                                      struct dstree_node *node,
                                      ts_type *timeseries) {

  // fprintf(stderr, "IN APPEND TS TO CHILD NODE.\n");
  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < ts_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }

  ++node->file_buffer->buffered_list_size;

  return SUCCESS;
}
/*
enum response append_ts_gt_to_child_node(struct dstree_index * index,
                                         struct dstree_node * node,
                                         ts_type * timeseries,
                                         unsigned char gt)
*/
/*
enum response append_ts_gt_to_child_node(struct dstree_index * index,
                                         struct dstree_node * node,
                                         ts_type * timeseries,
                                         label_type gt)
*/
enum response
append_ts_gt_to_child_node(struct dstree_index *index, struct dstree_node *node,
                           ts_type *timeseries, struct dstree_node *parent,
                           int parent_idx, unsigned int fp)

{

  // fprintf(stderr, "IN APPEND TS TO CHILD NODE.\n");
  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;

  int ts_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * ts_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < ts_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = timeseries[i];
  }

  ++node->file_buffer->buffered_list_size;

  // node->gt[idx] = gt;
  node->gt[idx] = parent->gt[parent_idx];
  if (fp != -1)
    node->fp[idx] = parent->fp[parent_idx];

  return SUCCESS;
}

enum response
create_dstree_node_filename(struct dstree_index_settings *settings,
                            struct dstree_node *node,
                            struct dstree_node *parent_node) {
  int i;

  node->filename = malloc(sizeof(char) * (settings->max_filename_size));

  int l = 0;
  l += sprintf(node->filename + l, "%02d", node->num_node_points);

  // If this has level other than 0 then it is not a root node and as such it
  // does have some split data on its parent.

  if (node->level) {
    if (node->is_left) {
      l += sprintf(node->filename + l, "%s", "_L");
    } else {
      l += sprintf(node->filename + l, "%s", "_R");
    }

    // Add parent split policy, just the code 0 for mean and 1 for stdev
    struct node_split_policy *policy;
    policy = parent_node->split_policy;

    if (policy->indicator_split_idx) {
      l += sprintf(node->filename + l, "%s", "_1");
    } else {
      l += sprintf(node->filename + l, "%s", "_0");
    }

    l += sprintf(node->filename + l, "_(%d,%d,%g)", policy->split_from,
                 policy->split_to, policy->indicator_split_value);
  }

  // If its level is 0 then it is a root
  l += sprintf(node->filename + l, "_%d", node->level);

  return SUCCESS;
}

enum response update_node_statistics(struct dstree_node *node,
                                     ts_type *timeseries) {

  // update vertical node_segment_sketch
  for (int i = 0; i < node->num_node_points; i++) {
    int from = 0;
    int to = 0;

    from = get_segment_start(node->node_points, i);
    to = get_segment_end(node->node_points, i);

    if (!node_segment_sketch_update_sketch(&node->node_segment_sketches[i],
                                           timeseries, from, to)) {
      fprintf(stderr, "Error in dstree_index.c:  Could not update vertical "
                      "sketch for node segment.\n");
      return FAILURE;
    }
  }

  // update horizontal node_segment_sketch
  for (int i = 0; i < node->num_hs_node_points; i++) {
    if (!node_segment_sketch_update_sketch(
            &node->hs_node_segment_sketches[i], timeseries,
            get_segment_start(node->hs_node_points, i),
            get_segment_end(node->hs_node_points, i))) {
      fprintf(stderr, "Error in dstree_index.c:  Could not update horizontal "
                      "sketch for node segment.\n");
      return FAILURE;
    }
  }

  ++node->node_size;

  return SUCCESS;
}

/*

  This function calculates the euclidean distance of query to a
  given node.

  The dstree and isax both load all the time series in the node
  and compare each to the query.

  The function returns the smallest distance between the query
  and the time series in the node.

 */

/* start kashif changes */
int calculate_node_knn_distance_2(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, struct query_result *knn_results,
    // struct bsf_snapshot **bsf_snapshots, unsigned int *cur_bsf_snapshot,
    unsigned int *cur_size, float warping, struct vid * query_id, 
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int approx)
{
  // get the k-th distance from the results queue
  ts_type kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;
  int num_nn = 0;

  unsigned int ts_byte_size =
      sizeof(ts_type) * index->settings->timeseries_size;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_CHECKED_NODE
  COUNT_CHECKED_TS(node->node_size)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it
  if (node->file_buffer->buffered_list_size == 0) {
    COUNT_LOADED_NODE
    COUNT_LOADED_TS(node->node_size)
    COUNT_PARTIAL_LOAD_NODE_TIME_START

    node->file_buffer->buffered_list = get_all_time_series_in_node(index, node);
    node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr,
              "Error in dstree_index.c:  Could not retrieve all time series "
              "for node %s.\n",
              node->filename);
    }
    COUNT_PARTIAL_LOAD_NODE_TIME_END
  }
  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory
  double tS_bsf;
  double tE_bsf;
  struct timeval current_time_bsf;
  boolean update_snapshots = false;

  for (unsigned int idx = 0; idx < node->file_buffer->buffered_list_size;
       ++idx) {
      
      // skip query table if found in index
      if(index->vid_cache[(node->vid_pos) + idx].table_id == query_id->table_id)
        continue;

    struct query_result result;
    result = knn_results[k - 1];
    kth_bsf = result.distance;

    //(sid) warping distance here
    if (warping > 0) {
      distance = ts_warping_distance(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order, warping);
      // distance = get_dtw(query_ts_reordered, query_order,
      // node->file_buffer->buffered_list[idx],
      // warping*index->settings->timeseries_size,
      // index->settings->timeseries_size, kth_bsf);
    } else {
      distance = ts_euclidean_distance_reordered(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order);
    }

    if (distance <= kth_bsf) // (tmp change) <= instead of <
    {
      num_nn++;
      struct query_result object_result; // =  malloc(sizeof(struct query_result));
      object_result.node = node;
      object_result.distance = distance;
      if (index->settings->classify)
        object_result.label = index->gt_cache[(node->gt_pos) + idx];

      if (index->settings->track_file_pos)
      {
        object_result.file_pos = index->fp_cache[(node->fp_pos) + idx];
        // object_result.series = calloc (1, ts_byte_size);
        // mempcpy(object_result.series,node->file_buffer->buffered_list[idx],ts_byte_size);
      }
      if (index->settings->track_vector)
      {
        object_result.vector_id = (struct vid *) malloc(sizeof(struct vid));
        // object_result.vector_id->table_id = 101010;
        // object_result.vector_id->set_id = 10101;
        object_result.vector_id->table_id = index->vid_cache[(node->vid_pos) + idx].table_id;
        object_result.vector_id->set_id = index->vid_cache[(node->vid_pos) + idx].set_id;
        object_result.vector_id->pos = index->vid_cache[(node->vid_pos) + idx].pos;
        object_result.time = 0;
        object_result.num_checked_vectors = 0;
        strcpy(object_result.vector_id->raw_data_file, index->vid_cache[(node->vid_pos) + idx].raw_data_file);
      }
          
      int stored_at = queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);
      // printf("NN stored at %d, ", stored_at+1);
      update_snapshots = true;
      
      if (index->settings->track_vector)
        free(object_result.vector_id);
      

      if(approx == 0 && knn_results[stored_at].approx == 1)
      {
        knn_results[stored_at].approx = 0; // overwrite result that was found in approximate search
      }

      knn_results[stored_at].num_checked_vectors = checked_ts_count;
     

    }
  }

  if (node->file_buffer != NULL) {
    // clearing the data for this node
    for (int i = 0; i < index->settings->max_leaf_size; ++i) {
      free(node->file_buffer->buffered_list[i]);
    }
    free(node->file_buffer->buffered_list);
  }

  node->file_buffer->buffered_list = NULL;
  node->file_buffer->buffered_list_size = 0;
  //  printf("total added nn = %d\n", num_nn);
  return num_nn;
}

/* start kashif changes */
int thread_queue_bounded_sorted_insert(struct dstree_index * index, struct query_result *q, struct query_result d,
                                unsigned int *cur_size, unsigned int k, unsigned int thread_id) {
  struct query_result temp;
  temp.vector_id = malloc(sizeof(struct vid));
  if (temp.vector_id == NULL)
  {
      printf("Error in dstree_node.c: Couldn't allocate memory for temp query result.");
      exit(1);
  }
  size_t i;
  size_t newsize;
  
  /* the queue is full, ovewrite last element*/
  if (*cur_size == k)
  {
    q[k - 1].distance = d.distance;
    q[k - 1].approx = d.approx;
    q[k - 1].vector_id->table_id = d.vector_id->table_id;
    q[k - 1].vector_id->set_id = d.vector_id->set_id;
    q[k - 1].vector_id->pos = d.vector_id->pos;
    q[k - 1].query_vector_pos = d.query_vector_pos;
    strcpy(q[k - 1].vector_id->raw_data_file, d.vector_id->raw_data_file);
    

    COUNT_THREAD_PARTIAL_TIME_END(thread_id)
    update_thread_query_stats(index, thread_id);
    q[k - 1].time = index->stats->thread_query_total_cpu_time[thread_id];
    RESET_THREAD_QUERY_COUNTERS(thread_id)
    RESET_THREAD_PARTIAL_COUNTERS(thread_id)
    COUNT_THREAD_PARTIAL_TIME_START(thread_id)
  }
  else
  {
    q[*cur_size].distance = d.distance;
    q[*cur_size].approx = d.approx;
    q[*cur_size].vector_id->table_id = d.vector_id->table_id;
    q[*cur_size].vector_id->set_id = d.vector_id->set_id;
    q[*cur_size].vector_id->pos = d.vector_id->pos;
    q[*cur_size].query_vector_pos = d.query_vector_pos;
    strcpy(q[*cur_size].vector_id->raw_data_file, d.vector_id->raw_data_file);

    COUNT_THREAD_PARTIAL_TIME_END(thread_id)
    update_thread_query_stats(index, thread_id);
    q[*cur_size].time = index->stats->thread_query_total_cpu_time[thread_id];
    RESET_THREAD_QUERY_COUNTERS(thread_id)
    RESET_THREAD_PARTIAL_COUNTERS(thread_id)
    COUNT_THREAD_PARTIAL_TIME_START(thread_id)

    ++(*cur_size);
  }

  unsigned int idx, j;

  idx = 1;

  while (idx < *cur_size) {
    j = idx;
    while (j > 0 && ((q[j - 1]).distance > q[j].distance)) {
      /* start kashif changes */
      // temp = q[j];
      temp.distance = q[j].distance;
      // temp.time = q[j].time;
      temp.approx = q[j].approx;
      temp.vector_id->table_id = q[j].vector_id->table_id;
      temp.vector_id->set_id = q[j].vector_id->set_id;
      temp.vector_id->pos = q[j].vector_id->pos;
      temp.query_vector_pos = q[j].query_vector_pos;
      strcpy(temp.vector_id->raw_data_file, q[j].vector_id->raw_data_file);

      /* end kashif changes */


      q[j].distance = q[j - 1].distance;
      // q[j].time = q[j - 1].time;
      q[j].approx = q[j - 1].approx;
      q[j].vector_id->table_id = q[j - 1].vector_id->table_id;
      q[j].vector_id->set_id = q[j - 1].vector_id->set_id;
      q[j].vector_id->pos = q[j - 1].vector_id->pos;
      q[j].query_vector_pos = q[j - 1].query_vector_pos;
      strcpy(q[j].vector_id->raw_data_file, q[j - 1].vector_id->raw_data_file);
      
      // COUNT_THREAD_PARTIAL_TIME_END(thread_id)
      // update_thread_query_stats(index, thread_id);
      // q[j].time = index->stats->thread_query_total_cpu_time[thread_id];
      // RESET_THREAD_QUERY_COUNTERS(thread_id)
      // RESET_THREAD_PARTIAL_COUNTERS(thread_id)
      // COUNT_THREAD_PARTIAL_TIME_START(thread_id)

      q[j - 1].distance = temp.distance;
      // q[j - 1].time = temp.time;
      q[j - 1].approx = temp.approx;
      q[j - 1].vector_id->table_id = temp.vector_id->table_id;
      q[j - 1].vector_id->set_id = temp.vector_id->set_id;
      q[j - 1].vector_id->pos = temp.vector_id->pos;
      q[j - 1].query_vector_pos = temp.query_vector_pos;
      strcpy(q[j - 1].vector_id->raw_data_file, temp.vector_id->raw_data_file);

      // COUNT_THREAD_PARTIAL_TIME_END(thread_id)
      // update_thread_query_stats(index, thread_id);
      // q[j - 1].time = index->stats->thread_query_total_cpu_time[thread_id];
      // RESET_THREAD_QUERY_COUNTERS(thread_id)
      // RESET_THREAD_PARTIAL_COUNTERS(thread_id)
      // COUNT_THREAD_PARTIAL_TIME_START(thread_id)
      --j;
    }
    ++idx;
  }
  
  free(temp.vector_id);
  return 0;
}
/* end kashif changes */

int calculate_node_knn_distance_para_incr(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, struct query_result *knn_results,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx)
{
  // get the k-th distance from the results queue
  ts_type kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;
  unsigned int num_nn = 0;

  unsigned int ts_byte_size =
      sizeof(ts_type) * index->settings->timeseries_size;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_THREAD_CHECKED_NODE(thread_id)
  COUNT_THREAD_CHECKED_TS(node->node_size,thread_id)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it

  if (node->file_buffer->buffered_list_size == 0)
  {
    pthread_mutex_lock(&node->lock);
    if (node->file_buffer->buffered_list_size == 0) {
      COUNT_THREAD_LOADED_NODE(thread_id) 
      COUNT_THREAD_LOADED_TS(node->node_size,thread_id)
      // COUNT_PARTIAL_LOAD_NODE_TIME_START

      node->file_buffer->buffered_list = get_all_time_series_in_node_para_incr(index, node, thread_id);
      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

      if (node->file_buffer->buffered_list == NULL) {
        fprintf(stderr,
                "Error in dstree_index.c:  Could not retrieve all time series "
                "for node %s.\n",
                node->filename);
      }
      // COUNT_PARTIAL_LOAD_NODE_TIME_END
    }
    pthread_mutex_unlock(&node->lock);
  }

  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory
  double tS_bsf;
  double tE_bsf;
  struct timeval current_time_bsf;
  boolean update_snapshots = false;

  for (unsigned int idx = 0; idx < node->file_buffer->buffered_list_size;
       ++idx) {
        
      if(index->vid_cache[(node->vid_pos) + idx].table_id == query_id->table_id)
        continue;

    struct query_result result;
    result = knn_results[k - 1];
    kth_bsf = result.distance;

    //(sid) warping distance here
    if (warping > 0) {
      distance = ts_warping_distance(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order, warping);
      // distance = get_dtw(query_ts_reordered, query_order,
      // node->file_buffer->buffered_list[idx],
      // warping*index->settings->timeseries_size,
      // index->settings->timeseries_size, kth_bsf);
    } else {
      distance = ts_euclidean_distance_reordered(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order);
    }

    if (distance <= kth_bsf) // (tmp change) <= instead of <
    {
      num_nn++;
      struct query_result object_result; // =  malloc(sizeof(struct query_result));
      object_result.node = node;
      object_result.distance = distance;
      if (index->settings->classify)
        object_result.label = index->gt_cache[(node->gt_pos) + idx];

      if (index->settings->track_file_pos)
      {
        object_result.file_pos = index->fp_cache[(node->fp_pos) + idx];
        // object_result.series = calloc (1, ts_byte_size);
        // mempcpy(object_result.series,node->file_buffer->buffered_list[idx],ts_byte_size);
      }
      if (index->settings->track_vector)
      {
        object_result.vector_id = (struct vid *) malloc(sizeof(struct vid));
        object_result.vector_id->table_id = index->vid_cache[(node->vid_pos) + idx].table_id;
        object_result.vector_id->set_id = index->vid_cache[(node->vid_pos) + idx].set_id;
        object_result.vector_id->pos = index->vid_cache[(node->vid_pos) + idx].pos;
        object_result.time = 0;
        object_result.approx = approx;
        object_result.num_checked_vectors = 0;
        strcpy(object_result.vector_id->raw_data_file, index->vid_cache[(node->vid_pos) + idx].raw_data_file);
      }
          
      // int stored_at = thread_queue_bounded_sorted_insert(index, knn_results, object_result, cur_size, k, thread_id);
      int stored_at = queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);
      // printf("bsf stored at %d,\n", stored_at);
      update_snapshots = true;
      
      if (index->settings->track_vector)
        free(object_result.vector_id);

      // COUNT_THREAD_PARTIAL_TIME_END(thread_id)
      // update_thread_query_stats(index, thread_id);

      if(approx == 0 && knn_results[stored_at].approx == 1)
      {
        knn_results[stored_at].approx = 0; // overwrite result that was found in approximate search
      }
      
      // knn_results[stored_at].time = index->stats->thread_query_total_cpu_time[thread_id];
      // knn_results[stored_at].num_checked_vectors = thread_checked_ts_count[thread_id];
     
      // RESET_THREAD_QUERY_COUNTERS(thread_id)
      // RESET_THREAD_PARTIAL_COUNTERS(thread_id)
      // COUNT_THREAD_PARTIAL_TIME_START(thread_id)
      
    }
  }

  // if (node->file_buffer != NULL) {
  //   // clearing the data for this node
  //   for (int i = 0; i < index->settings->max_leaf_size; ++i) {
  //     free(node->file_buffer->buffered_list[i]);
  //   }
  //   free(node->file_buffer->buffered_list);
  // }

  // node->file_buffer->buffered_list = NULL;
  // node->file_buffer->buffered_list_size = 0;

  return num_nn;
}

int calculate_node_knn_distance_para_incr_mmheap(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, void *knn_heap,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx, unsigned long * insert_counter)
{
  // get the k-th distance from the results queue
  ts_type kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;
  unsigned int num_nn = 0;

  unsigned int ts_byte_size =
      sizeof(ts_type) * index->settings->timeseries_size;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_THREAD_CHECKED_NODE(thread_id)
  COUNT_THREAD_CHECKED_TS(node->node_size,thread_id)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it

  if (node->file_buffer->buffered_list_size == 0)
  {
    pthread_mutex_lock(&node->lock);
    if (node->file_buffer->buffered_list_size == 0) {
      COUNT_THREAD_LOADED_NODE(thread_id) 
      COUNT_THREAD_LOADED_TS(node->node_size,thread_id)
      // COUNT_PARTIAL_LOAD_NODE_TIME_START

      node->file_buffer->buffered_list = get_all_time_series_in_node_para_incr(index, node, thread_id);
      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

      if (node->file_buffer->buffered_list == NULL) {
        fprintf(stderr,
                "Error in dstree_index.c:  Could not retrieve all time series "
                "for node %s.\n",
                node->filename);
      }
      // COUNT_PARTIAL_LOAD_NODE_TIME_END
    }
    pthread_mutex_unlock(&node->lock);
  }

  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory
  double tS_bsf;
  double tE_bsf;
  struct timeval current_time_bsf;
  boolean update_snapshots = false;

  for (unsigned int idx = 0; idx < node->file_buffer->buffered_list_size;
       ++idx) {
        
      if(index->vid_cache[(node->vid_pos) + idx].table_id == query_id->table_id)
        continue;

    struct query_result * kth_result = mmheap_get_max(knn_heap);
    if(kth_result == NULL)// empty heap
      kth_bsf = FLT_MAX;
    else
      kth_bsf = kth_result->distance;

    //(sid) warping distance here
    if (warping > 0) {
      distance = ts_warping_distance(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order, warping);
      // distance = get_dtw(query_ts_reordered, query_order,
      // node->file_buffer->buffered_list[idx],
      // warping*index->settings->timeseries_size,
      // index->settings->timeseries_size, kth_bsf);
    } else {
      distance = ts_euclidean_distance_reordered(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order);
    }

    if (distance <= kth_bsf) // (tmp change) <= instead of <
    {
      num_nn++;
      struct query_result * object_result = malloc(sizeof(struct query_result)); // =  malloc(sizeof(struct query_result));
      object_result->node = node;
      object_result->distance = distance;
      if (index->settings->classify)
        object_result->label = index->gt_cache[(node->gt_pos) + idx];

      if (index->settings->track_file_pos)
      {
        object_result->file_pos = index->fp_cache[(node->fp_pos) + idx];
        // object_result.series = calloc (1, ts_byte_size);
        // mempcpy(object_result.series,node->file_buffer->buffered_list[idx],ts_byte_size);
      }
      if (index->settings->track_vector)
      {
        object_result->vector_id = (struct vid *) malloc(sizeof(struct vid));
        object_result->vector_id->table_id = index->vid_cache[(node->vid_pos) + idx].table_id;
        object_result->vector_id->set_id = index->vid_cache[(node->vid_pos) + idx].set_id;
        object_result->vector_id->pos = index->vid_cache[(node->vid_pos) + idx].pos;
        object_result->time = 0;
        object_result->approx = approx;
        object_result->num_checked_vectors = 0;
        strcpy(object_result->vector_id->raw_data_file, index->vid_cache[(node->vid_pos) + idx].raw_data_file);
      }
          

      mmheap_insert(knn_heap, object_result);
      update_snapshots = true;
    }
  }
  return num_nn;
}

int calculate_node_knn_distance_para_incr_ostree(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, void *knn_tree,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx, unsigned long * insert_counter)
{
  // get the k-th distance from the results queue
  ts_type kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;
  unsigned int num_nn = 0;

  unsigned int ts_byte_size =
      sizeof(ts_type) * index->settings->timeseries_size;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_THREAD_CHECKED_NODE(thread_id)
  COUNT_THREAD_CHECKED_TS(node->node_size,thread_id)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it

  if (node->file_buffer->buffered_list_size == 0)
  {
    pthread_mutex_lock(&node->lock);
    if (node->file_buffer->buffered_list_size == 0) {
      COUNT_THREAD_LOADED_NODE(thread_id) 
      COUNT_THREAD_LOADED_TS(node->node_size,thread_id)
      // COUNT_PARTIAL_LOAD_NODE_TIME_START

      node->file_buffer->buffered_list = get_all_time_series_in_node_para_incr(index, node, thread_id);
      node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

      if (node->file_buffer->buffered_list == NULL) {
        fprintf(stderr,
                "Error in dstree_index.c:  Could not retrieve all time series "
                "for node %s.\n",
                node->filename);
      }
      // COUNT_PARTIAL_LOAD_NODE_TIME_END
    }
    pthread_mutex_unlock(&node->lock);
  }

  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory
  double tS_bsf;
  double tE_bsf;
  struct timeval current_time_bsf;
  boolean update_snapshots = false;

  for (unsigned int idx = 0; idx < node->file_buffer->buffered_list_size;
       ++idx) {
        
      if(index->vid_cache[(node->vid_pos) + idx].table_id == query_id->table_id)
        continue;

    struct query_result * kth_result = ostree_get_max(knn_tree);
    kth_bsf = kth_result->distance;

    //(sid) warping distance here
    if (warping > 0) {
      distance = ts_warping_distance(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order, warping);
      // distance = get_dtw(query_ts_reordered, query_order,
      // node->file_buffer->buffered_list[idx],
      // warping*index->settings->timeseries_size,
      // index->settings->timeseries_size, kth_bsf);
    } else {
      distance = ts_euclidean_distance_reordered(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order);
    }

    if (distance <= kth_bsf) // (tmp change) <= instead of <
    {
      num_nn++;
      struct query_result * object_result = malloc(sizeof(struct query_result)); // =  malloc(sizeof(struct query_result));
      object_result->node = node;
      object_result->distance = distance;
      if (index->settings->classify)
        object_result->label = index->gt_cache[(node->gt_pos) + idx];

      if (index->settings->track_file_pos)
      {
        object_result->file_pos = index->fp_cache[(node->fp_pos) + idx];
        // object_result.series = calloc (1, ts_byte_size);
        // mempcpy(object_result.series,node->file_buffer->buffered_list[idx],ts_byte_size);
      }
      if (index->settings->track_vector)
      {
        object_result->vector_id = (struct vid *) malloc(sizeof(struct vid));
        object_result->vector_id->table_id = index->vid_cache[(node->vid_pos) + idx].table_id;
        object_result->vector_id->set_id = index->vid_cache[(node->vid_pos) + idx].set_id;
        object_result->vector_id->pos = index->vid_cache[(node->vid_pos) + idx].pos;
        object_result->time = 0;
        object_result->approx = approx;
        object_result->num_checked_vectors = 0;
        strcpy(object_result->vector_id->raw_data_file, index->vid_cache[(node->vid_pos) + idx].raw_data_file);
      }
          

      ostree_insert(knn_tree, object_result, insert_counter);
      // printf("bsf stored at %d,\n", stored_at);
      update_snapshots = true;
      
      // if (index->settings->track_vector)
      //   free(object_result.vector_id);

      // COUNT_THREAD_PARTIAL_TIME_END(thread_id)
      // update_thread_query_stats(index, thread_id);

      // if(approx == 0 && knn_results[stored_at].approx == 1)
      // {
      //   knn_results[stored_at].approx = 0; // overwrite result that was found in approximate search
      // }
      
      // knn_results[stored_at].time = index->stats->thread_query_total_cpu_time[thread_id];
      // knn_results[stored_at].num_checked_vectors = thread_checked_ts_count[thread_id];
     
      // RESET_THREAD_QUERY_COUNTERS(thread_id)
      // RESET_THREAD_PARTIAL_COUNTERS(thread_id)
      // COUNT_THREAD_PARTIAL_TIME_START(thread_id)
      
    }
  }

  // if (node->file_buffer != NULL) {
  //   // clearing the data for this node
  //   for (int i = 0; i < index->settings->max_leaf_size; ++i) {
  //     free(node->file_buffer->buffered_list[i]);
  //   }
  //   free(node->file_buffer->buffered_list);
  // }

  // node->file_buffer->buffered_list = NULL;
  // node->file_buffer->buffered_list_size = 0;

  return num_nn;
}


/* end kashif changes */

void calculate_node_knn_distance(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    ts_type bsf, unsigned int k, struct query_result *knn_results,
    struct bsf_snapshot **bsf_snapshots, unsigned int *cur_bsf_snapshot,
    unsigned int *cur_size, float warping)
{
  // get the k-th distance from the results queue
  ts_type kth_bsf = FLT_MAX;
  ts_type distance = FLT_MAX;

  unsigned int ts_byte_size =
      sizeof(ts_type) * index->settings->timeseries_size;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_CHECKED_NODE
  COUNT_CHECKED_TS(node->node_size)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it
  if (node->file_buffer->buffered_list_size == 0) {
    COUNT_LOADED_NODE
    COUNT_LOADED_TS(node->node_size)
    COUNT_PARTIAL_LOAD_NODE_TIME_START

    node->file_buffer->buffered_list = get_all_time_series_in_node(index, node);
    node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr,
              "Error in dstree_index.c:  Could not retrieve all time series "
              "for node %s.\n",
              node->filename);
    }
    COUNT_PARTIAL_LOAD_NODE_TIME_END
  }
  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory
  double tS_bsf;
  double tE_bsf;
  struct timeval current_time_bsf;
  boolean update_snapshots = false;

  for (unsigned int idx = 0; idx < node->file_buffer->buffered_list_size;
       ++idx) {
    struct query_result result;
    result = knn_results[k - 1];

    kth_bsf = result.distance;

    //(sid) warping distance here
    if (warping > 0) {
      distance = ts_warping_distance(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order, warping);
      // distance = get_dtw(query_ts_reordered, query_order,
      // node->file_buffer->buffered_list[idx],
      // warping*index->settings->timeseries_size,
      // index->settings->timeseries_size, kth_bsf);
    } else {
      distance = ts_euclidean_distance_reordered(
          query_ts_reordered, node->file_buffer->buffered_list[idx],
          offset, // offset is 0 for whole matching
          index->settings->timeseries_size, kth_bsf, query_order);
    }

    if (distance <= kth_bsf) // (tmp change) <= instead of <
    {
      struct query_result object_result; // =  malloc(sizeof(struct query_result));
      object_result.node = node;
      object_result.distance = distance;
      if (index->settings->classify)
        object_result.label = index->gt_cache[(node->gt_pos) + idx];

      if (index->settings->track_file_pos)
      {
        object_result.file_pos = index->fp_cache[(node->fp_pos) + idx];
        // object_result.series = calloc (1, ts_byte_size);
        // mempcpy(object_result.series,node->file_buffer->buffered_list[idx],ts_byte_size);
      }
      if (index->settings->track_vector)
      {
        object_result.vector_id = (struct vid *) malloc(sizeof(struct vid));
        // object_result.vector_id->table_id = 101010;
        // object_result.vector_id->set_id = 10101;
        object_result.vector_id->table_id = index->vid_cache[(node->vid_pos) + idx].table_id;
        object_result.vector_id->set_id = index->vid_cache[(node->vid_pos) + idx].set_id;
        object_result.vector_id->pos = index->vid_cache[(node->vid_pos) + idx].pos;
        strcpy(object_result.vector_id->raw_data_file, index->vid_cache[(node->vid_pos) + idx].raw_data_file);
      }

      queue_bounded_sorted_insert(knn_results, object_result, cur_size, k);
      update_snapshots = true;
      
      if (index->settings->track_vector)
        free(object_result.vector_id);
      
    }
  }
  // only print the snapshots after finished visiting leaf
  // if interested in the value of the actual neighbor, move this code inside
  // the for loop above
  if (cur_bsf_snapshot != NULL && update_snapshots) {
    gettimeofday(&current_time_bsf, NULL);
    tS_bsf = partial_time_start.tv_sec * 1000000 + (partial_time_start.tv_usec);
    tE_bsf = current_time_bsf.tv_sec * 1000000 + (current_time_bsf.tv_usec);

    for (int j = 0; j < k; ++j) {
      bsf_snapshots[j][*cur_bsf_snapshot].distance = knn_results[j].distance;
      bsf_snapshots[j][*cur_bsf_snapshot].time = tE_bsf - tS_bsf;
      bsf_snapshots[j][*cur_bsf_snapshot].checked_nodes = checked_nodes_count;

      if (index->settings->classify)
        bsf_snapshots[j][*cur_bsf_snapshot].label = knn_results[j].label;

      if (index->settings->track_file_pos) {
        bsf_snapshots[j][*cur_bsf_snapshot].file_pos = knn_results[j].file_pos;
        // bsf_snapshots[j][*cur_bsf_snapshot].series = calloc (1,
        // ts_byte_size);
        // mempcpy(bsf_snapshots[j][*cur_bsf_snapshot].series,knn_results[j].series,ts_byte_size);
      }

      if(index->settings->track_vector)
      {
        bsf_snapshots[j][*cur_bsf_snapshot].vector_id->table_id = knn_results[j].vector_id->table_id;
        bsf_snapshots[j][*cur_bsf_snapshot].vector_id->set_id = knn_results[j].vector_id->set_id;
        bsf_snapshots[j][*cur_bsf_snapshot].vector_id->pos = knn_results[j].vector_id->pos;
        bsf_snapshots[j][*cur_bsf_snapshot].query_vector_pos = knn_results[j].query_vector_pos;
        strcpy(bsf_snapshots[j][*cur_bsf_snapshot].vector_id->raw_data_file, 
              knn_results[j].vector_id->raw_data_file);
      }
    }
    ++(*cur_bsf_snapshot);
  }

  if (node->file_buffer != NULL) {
    // clearing the data for this node
    for (int i = 0; i < index->settings->max_leaf_size; ++i) {
      free(node->file_buffer->buffered_list[i]);
    }
    free(node->file_buffer->buffered_list);
  }

  node->file_buffer->buffered_list = NULL;
  node->file_buffer->buffered_list_size = 0;
}

ts_type calculate_node_distance(struct dstree_index *index,
                                struct dstree_node *node,
                                ts_type *query_ts_reordered, int *query_order,
                                unsigned int offset, ts_type bsf) {
  //    ts_type bsf = FLT_MAX;
  ts_type dist = FLT_MAX;

  // count the number of leaves and the number of time series that were checked
  // this is different from the count_loaded_node and counted_loaded_ts
  // which count the number of leaves and time series that were not found in
  // memory and had to be retrieved from disk
  // checked_nodes = loaded_nodes + nodes_in_memory
  COUNT_CHECKED_NODE
  COUNT_CHECKED_TS(node->node_size)

  // TEST THAT DATA IS FULLY IN MEM
  // If the leaf's data is in disk, load it
  if (node->file_buffer->buffered_list_size == 0) {
    COUNT_LOADED_NODE
    COUNT_LOADED_TS(node->node_size)
    COUNT_PARTIAL_LOAD_NODE_TIME_START

    node->file_buffer->buffered_list = get_all_time_series_in_node(index, node);
    node->file_buffer->buffered_list_size = node->file_buffer->disk_count;

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr,
              "Error in dstree_index.c:  Could not retrieve all time series "
              "for node %s.\n",
              node->filename);
    }
    COUNT_PARTIAL_LOAD_NODE_TIME_END
  }
  // If the leaf's data is in memory, proceed. A leaf's data is either fully in
  // disk or in memory

  dist = calculate_ts_in_node_distance(index, node, query_ts_reordered,
                                       query_order, offset, bsf);

  // clearing the data for this node
  for (int i = 0; i < index->settings->max_leaf_size; ++i) {
    free(node->file_buffer->buffered_list[i]);
  }

  free(node->file_buffer->buffered_list);

  node->file_buffer->buffered_list = NULL;
  node->file_buffer->buffered_list_size = 0;

  return dist;
}

ts_type calculate_ts_in_node_distance(struct dstree_index *index,
                                      struct dstree_node *node,
                                      ts_type *query_ts_reordered,
                                      int *query_order, unsigned int offset,
                                      ts_type bound) {
  ts_type bsf = bound;
  ts_type temp_dist;

  for (int idx = 0; idx < node->file_buffer->buffered_list_size; ++idx) {

    temp_dist = ts_euclidean_distance_reordered(
        query_ts_reordered, node->file_buffer->buffered_list[idx],
        offset, // offset is 0 for whole matching
        index->settings->timeseries_size, bsf, query_order);
    if (temp_dist < bsf) {
      bsf = temp_dist;
    }
  }

  return bsf;
}

/////// new code start

ts_type calculate_node_min_distance_warping(struct dstree_index *index,
                                            struct dstree_node *node,
                                            ts_type *query, float warping) {
  // Window fraction = 5%
  float window_frac = warping;

  ts_type sum = 0;
  short *points = node->node_points;
  int num_points = (int)node->num_node_points;
  int length = points[num_points - 1];

  // lower and upper envelope (raw)
  ts_type *le = malloc(sizeof(ts_type) * length);
  ts_type *ue = malloc(sizeof(ts_type) * length);

  // lower and upper envelope (segment)
  ts_type *l_hat = malloc(sizeof(ts_type) * num_points);
  ts_type *u_hat = malloc(sizeof(ts_type) * num_points);

  lower_upper_lemire(query, length, window_frac * length, le, ue);

  assign_max_segments(ue, u_hat, points, num_points);
  assign_min_segments(le, l_hat, points, num_points);

  for (int i = 0; i < num_points; ++i) {

    ts_type temp_dist = 0;
    // mu_min > u_hat
    if (node->node_segment_sketches[i].indicators[1] > u_hat[i]) {
      temp_dist +=
          pow(node->node_segment_sketches[i].indicators[1] - u_hat[i], 2);
    }
    // mu_max < l_hat
    if (node->node_segment_sketches[i].indicators[0] < l_hat[i]) {
      temp_dist +=
          pow(l_hat[i] - node->node_segment_sketches[i].indicators[0], 2);
    }

    sum += temp_dist * get_segment_length(points, i);
  }

  free(ue);
  free(le);
  free(u_hat);
  free(l_hat);

  return sum;
}

/////// new code end

ts_type calculate_node_min_distance(struct dstree_index *index,
                                    struct dstree_node *node, ts_type *query,
                                    float warping) {
  // (sid) LBKeogh lower bound here
  if (warping > 0)
    return calculate_node_min_distance_warping(index, node, query, warping);

  ts_type sum = 0;
  short *points = node->node_points;
  int num_points = (int)node->num_node_points;

  ts_type *mean_per_segment = malloc(sizeof(ts_type) * num_points);
  ts_type *stdev_per_segment = malloc(sizeof(ts_type) * num_points);

  calc_mean_per_segment(query, points, mean_per_segment, num_points);
  calc_stdev_per_segment(query, points, stdev_per_segment, num_points);

  for (int i = 0; i < num_points; ++i) {
    // use mean and standard deviation to estimate the distance
    ts_type temp_dist = 0;

    if ((stdev_per_segment[i] - node->node_segment_sketches[i].indicators[2]) *
            (stdev_per_segment[i] -
             node->node_segment_sketches[i].indicators[3]) >
        0) {
      temp_dist += pow(fmin(fabs(stdev_per_segment[i] -
                                 node->node_segment_sketches[i].indicators[2]),
                            fabs(stdev_per_segment[i] -
                                 node->node_segment_sketches[i].indicators[3])),
                       2);
    }

    if ((mean_per_segment[i] - node->node_segment_sketches[i].indicators[0]) *
            (mean_per_segment[i] -
             node->node_segment_sketches[i].indicators[1]) >
        0) {
      temp_dist += pow(fmin(fabs(mean_per_segment[i] -
                                 node->node_segment_sketches[i].indicators[0]),
                            fabs(mean_per_segment[i] -
                                 node->node_segment_sketches[i].indicators[1])),
                       2);
    }
    sum += temp_dist * get_segment_length(points, i);
  }

  // sum = sqrt(sum);

  free(mean_per_segment);
  free(stdev_per_segment);
  return sum;
}

ts_type calculate_node_max_distance(struct dstree_index *index,
                                    struct dstree_node *node, ts_type *query)
{
  ts_type sum = 0;
  short *points = node->node_points;
  int num_points = (int)node->num_node_points;

  ts_type *mean_per_segment = malloc(sizeof(ts_type) * num_points);
  ts_type *stdev_per_segment = malloc(sizeof(ts_type) * num_points);

  calc_mean_per_segment(query, points, mean_per_segment, num_points);
  calc_stdev_per_segment(query, points, stdev_per_segment, num_points);

  for (int i = 0; i < num_points; ++i) {
    // use mean and standard deviation to estimate the distance
    ts_type temp_dist = 0;

    temp_dist += pow(
        stdev_per_segment[i] + node->node_segment_sketches[i].indicators[2], 2);

    ts_type ub_threshold = (node->node_segment_sketches[i].indicators[0] +
                            node->node_segment_sketches[i].indicators[1]) /
                           2.0;

    if (mean_per_segment[i] <= ub_threshold) {
      temp_dist += pow(fabs(mean_per_segment[i] -
                            node->node_segment_sketches[i].indicators[0]),
                       2);
    } else {
      temp_dist += pow(fabs(mean_per_segment[i] -
                            node->node_segment_sketches[i].indicators[1]),
                       2);
    }
    sum += temp_dist * get_segment_length(points, i);
  }
  // sum = sqrt(sum);

  free(mean_per_segment);
  free(stdev_per_segment);

  if (sum == 0)
    return FLT_MAX;
  else
    return sum;
}

int queue_bounded_sorted_insert(struct query_result *q, struct query_result d,
                                unsigned int *cur_size, unsigned int k) {
  struct query_result temp;
  temp.vector_id = malloc(sizeof(struct vid));
  if (temp.vector_id == NULL)
  {
      printf("Error in dstree_node.c: Couldn't allocate memory for temp query result.");
      exit(1);
  }
  size_t i;
  size_t newsize;

  int stored_at = 0;
  /* the queue is full, ovewrite last element*/
  if (*cur_size == k)
  {
    stored_at = k - 1;
    q[k - 1].distance = d.distance;
    q[k - 1].node = d.node;
    q[k - 1].label = d.label;
    q[k - 1].file_pos = d.file_pos;
    q[k - 1].vector_id->table_id = d.vector_id->table_id;
    q[k - 1].vector_id->set_id = d.vector_id->set_id;
    q[k - 1].vector_id->pos = d.vector_id->pos;
    // q[k - 1].time = d.time;
    q[k - 1].approx = d.approx;
    q[k - 1].num_checked_vectors += d.num_checked_vectors;
    // q[k - 1].query_vector_pos = d.query_vector_pos; // because query pos is set at the beginning 
    strcpy(q[k - 1].vector_id->raw_data_file, d.vector_id->raw_data_file);
  }
  else
  {
    stored_at = *cur_size;
    q[*cur_size].distance = d.distance;
    q[*cur_size].node = d.node;
    q[*cur_size].label = d.label;
    q[*cur_size].file_pos = d.file_pos;
    q[*cur_size].vector_id->table_id = d.vector_id->table_id;
    q[*cur_size].vector_id->set_id = d.vector_id->set_id;
    q[*cur_size].vector_id->pos = d.vector_id->pos;
    // q[*cur_size].time = d.time;
    q[*cur_size].approx = d.approx;
    q[*cur_size].num_checked_vectors = d.num_checked_vectors;
    // q[*cur_size].query_vector_pos = d.query_vector_pos; // because query pos is set at the beginning 
    strcpy(q[*cur_size].vector_id->raw_data_file, d.vector_id->raw_data_file);

    ++(*cur_size);
  }

  unsigned int idx, j;

  idx = 1;

  while (idx < *cur_size) {
    j = idx;
    while (j > 0 && ((q[j - 1]).distance > q[j].distance)) {
      /* start kashif changes */
      // temp = q[j];

      if(stored_at == j)
        stored_at = j - 1;
      else if (stored_at = (j - 1))
        stored_at = j;

      temp.distance = q[j].distance;
      temp.node = q[j].node;
      temp.label = q[j].label;
      temp.file_pos = q[j].file_pos;
      temp.vector_id->table_id = q[j].vector_id->table_id;
      temp.vector_id->set_id = q[j].vector_id->set_id;
      temp.vector_id->pos = q[j].vector_id->pos;
      // temp.time = q[j].time;
      temp.approx = q[j].approx;
      temp.num_checked_vectors = q[j].num_checked_vectors;
      // temp.query_vector_pos = q[j].query_vector_pos; // because query pos is set at the beginning 
      strcpy(temp.vector_id->raw_data_file, q[j].vector_id->raw_data_file);
      
      /* end kashif changes */
      q[j].distance = q[j - 1].distance;
      q[j].node = q[j - 1].node;
      q[j].label = q[j - 1].label;
      q[j].file_pos = q[j - 1].file_pos;

      q[j].vector_id->table_id = q[j - 1].vector_id->table_id;
      q[j].vector_id->set_id = q[j - 1].vector_id->set_id;
      q[j].vector_id->pos = q[j - 1].vector_id->pos;
      // q[j].time = q[j - 1].time;
      q[j].approx = q[j - 1].approx;
      q[j].num_checked_vectors = q[j - 1].num_checked_vectors;
      // q[j].query_vector_pos = q[j - 1].query_vector_pos; // because query pos is set at the beginning 
      strcpy(q[j].vector_id->raw_data_file, q[j - 1].vector_id->raw_data_file);

      q[j - 1].distance = temp.distance;
      q[j - 1].node = temp.node;
      q[j - 1].label = temp.label;
      q[j - 1].file_pos = temp.file_pos;
      q[j - 1].vector_id->table_id = temp.vector_id->table_id;
      q[j - 1].vector_id->set_id = temp.vector_id->set_id;
      q[j - 1].vector_id->pos = temp.vector_id->pos;
      // q[j - 1].time = temp.time;
      q[j - 1].approx = temp.approx;
      q[j - 1].num_checked_vectors = temp.num_checked_vectors;
      // q[j - 1].query_vector_pos = temp.query_vector_pos; // because query pos is set at the beginning 
      strcpy(q[j - 1].vector_id->raw_data_file, temp.vector_id->raw_data_file);
      --j;
    }
    ++idx;
  }
  
  free(temp.vector_id);
  return stored_at;
}

/* start kashif changes */
// new function append vector and track table_id, set_id
enum response append_vector_to_node(struct dstree_index *index,
                                    struct dstree_node *node, ts_type *vector,
                                    unsigned int table_id,
                                    unsigned int set_id,
                                    unsigned int pos, 
                                    char * raw_data_file)
{
  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;
  int vector_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * vector_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the vector in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < vector_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = vector[i];
  }

  if (index->settings->track_vector) {
    node->vid[node->node_size - 1].table_id = table_id;
    node->vid[node->node_size - 1].set_id = set_id;
    node->vid[node->node_size - 1].pos = pos;
    strcpy(node->vid[node->node_size - 1].raw_data_file, raw_data_file);
  }


  ++node->file_buffer->buffered_list_size;
  index->buffer_manager->current_count += vector_length;

  return SUCCESS;
}

enum response append_vector_to_child_node(struct dstree_index *index,
              struct dstree_node *node, ts_type *vector, unsigned int table_id, 
              unsigned int set_id, unsigned int pos, char * raw_data_file)
{
  // fprintf(stderr, "IN APPEND TS TO CHILD NODE.\n");
  if (!get_file_buffer(index, node)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not get the \
                     file buffer for this node.\n");
    return FAILURE;
  }

  if (node->file_buffer == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Null file buffer for \
                     this node after creating it.\n");
    return FAILURE;
  }

  int idx = node->file_buffer->buffered_list_size;

  int vector_length = index->settings->timeseries_size;
  int max_leaf_size = index->settings->max_leaf_size;

  if (idx == 0) {
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list =
        malloc(sizeof(struct ts_type *) * max_leaf_size);

    if (node->file_buffer->buffered_list == NULL) {
      fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the buffered list. \n");
      return FAILURE;
    }
  }

  node->file_buffer->buffered_list[idx] =
      (ts_type *)index->buffer_manager->current_record;
  index->buffer_manager->current_record += sizeof(ts_type) * vector_length;
  index->buffer_manager->current_record_index++;

  if (node->file_buffer->buffered_list[idx] == NULL) {
    fprintf(stderr, "Error in dstree_index.c:  Could not \
                         allocate memory for the time series in\
                         the buffer.\n");
    return FAILURE;
  }

  for (int i = 0; i < vector_length; ++i) {
    node->file_buffer->buffered_list[idx][i] = vector[i];
  }

  if (index->settings->track_vector) {
    node->vid[node->node_size - 1].table_id = table_id;
    node->vid[node->node_size - 1].set_id = set_id;
    node->vid[node->node_size - 1].pos = pos;
    strcpy(node->vid[node->node_size - 1].raw_data_file, raw_data_file);
  }

  ++node->file_buffer->buffered_list_size;

  return SUCCESS;
}
/* end kashif changes */
