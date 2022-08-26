//  dstree_query_engine.c
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#include "../config.h"
#include "../globals.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/dstree_file_buffer.h"
#include "../include/dstree_file_buffer_manager.h"
#include "../include/dstree_index.h"
#include "../include/dstree_node.h"
#include "../include/dstree_query_engine.h"

#include "../include/pqueue.h"
#ifdef VALUES
#include <values.h>
#endif

struct query_result approximate_search(ts_type *query_ts,
                                       ts_type *query_ts_reordered,
                                       int *query_order, unsigned int offset,
                                       ts_type bsf,
                                       struct dstree_index *index) {

  struct query_result result;
  struct dstree_node *node = index->first_node;
  // ts_type bsf = FLT_MAX;  //no bsf known so far

  if (node != NULL) {
    // Traverse tree
    while (!node->is_leaf) {
      if (node_split_policy_route_to_left(node, query_ts)) {
        node = node->left_child;
      } else {
        node = node->right_child;
      }
    }

    result.distance = calculate_node_distance(index, node, query_ts_reordered,
                                              query_order, offset, bsf);
    result.node = node;
  } else {
    printf("Error: index is empty \n");
    result.node = NULL;
    result.distance = MAXFLOAT;
  }

  return result;
}

// get the k best neighbors from one leaf
void approximate_knn_search(ts_type *query_ts, ts_type *query_ts_reordered,
                            int *query_order, unsigned int offset, ts_type bsf,
                            struct dstree_index *index,
                            struct query_result *knn_results, unsigned int k,
                            struct bsf_snapshot **bsf_snapshots,
                            unsigned int *cur_bsf_snapshot,
                            unsigned int *curr_size, float warping) {

  struct query_result result;
  struct dstree_node *node = index->first_node;
  // ts_type bsf = FLT_MAX;  //no bsf known so far

  if (node != NULL) {
    // Traverse tree
    while (!node->is_leaf) {
      if (node_split_policy_route_to_left(node, query_ts)) {
        node = node->left_child;
      } else {
        node = node->right_child;
      }
    }

    calculate_node_knn_distance(index, node, query_ts_reordered, query_order,
                                offset, bsf, k, knn_results, bsf_snapshots,
                                cur_bsf_snapshot, curr_size, warping);
  } else {
    printf("Error in dstree_query_engine: null pointer to node.\n");
  }
}

struct query_result exact_search(ts_type *query_ts, ts_type *query_ts_reordered,
                                 int *query_order, unsigned int offset,
                                 struct dstree_index *index,
                                 ts_type minimum_distance, ts_type epsilon,
                                 ts_type delta) {

  ts_type bsf = FLT_MAX;

  struct query_result approximate_result = approximate_search(
      query_ts, query_ts_reordered, query_order, offset, bsf, index);

  struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    return approximate_result;
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  // struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  if (approximate_result.node != NULL) {
    // Insert approximate result in heap.
    pqueue_insert(pq, &approximate_result);
  }

  struct query_result *do_not_remove = &approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  while ((n = pqueue_pop(pq))) {
    if (n->distance > bsf_result.distance / (1 + epsilon)) {
      break;
    }
    if (n->node->is_leaf) // n is a leaf
    {

      ts_type distance =
          calculate_node_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance);

      if (distance < bsf_result.distance) {
        bsf_result.distance = distance;
        bsf_result.node = n->node;
      }
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      if (child_distance < bsf_result.distance / (1 + epsilon)) {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        // index->stats->query_lb_distance = sqrtf(child_distance);
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      if (child_distance < bsf_result.distance / (1 + epsilon)) {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        // index->stats->query_lb_distance = sqrtf(child_distance);
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    if (n != do_not_remove)
      free(n);
  }

  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    if (n != do_not_remove)
      free(n);
  }
  // Free the priority queue.

  pqueue_free(pq);

  COUNT_PARTIAL_TIME_END

  index->stats->query_refine_total_time = partial_time;

  index->stats->query_refine_input_time = partial_input_time;
  index->stats->query_refine_output_time = partial_output_time;
  index->stats->query_refine_load_node_time = partial_load_node_time;
  index->stats->query_refine_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_refine_seq_input_count = partial_seq_input_count;
  index->stats->query_refine_seq_output_count = partial_seq_output_count;
  index->stats->query_refine_rand_input_count = partial_rand_input_count;
  index->stats->query_refine_rand_output_count = partial_rand_output_count;

  index->stats->query_total_time =
      partial_time + index->stats->query_filter_total_time;
  index->stats->query_total_input_time =
      partial_input_time + index->stats->query_filter_input_time;
  index->stats->query_total_output_time =
      partial_output_time + index->stats->query_filter_output_time;
  index->stats->query_total_load_node_time =
      partial_load_node_time + index->stats->query_filter_load_node_time;
  index->stats->query_total_cpu_time = index->stats->query_total_time -
                                       index->stats->query_total_input_time -
                                       index->stats->query_total_output_time;

  index->stats->query_total_seq_input_count =
      partial_seq_input_count + index->stats->query_filter_seq_input_count;
  index->stats->query_total_seq_output_count =
      partial_seq_output_count + index->stats->query_filter_seq_output_count;
  index->stats->query_total_rand_input_count =
      partial_rand_input_count + index->stats->query_filter_rand_input_count;
  index->stats->query_total_rand_output_count =
      partial_rand_output_count + index->stats->query_filter_rand_output_count;

  index->stats->query_exact_distance = sqrtf(bsf_result.distance);
  index->stats->query_exact_node_filename = bsf_result.node->filename;
  index->stats->query_exact_node_size = bsf_result.node->node_size;
  ;
  index->stats->query_exact_node_level = bsf_result.node->level;

  index->stats->query_refine_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_refine_loaded_ts_count = loaded_ts_count;
  index->stats->query_refine_checked_nodes_count = checked_nodes_count;
  index->stats->query_refine_checked_ts_count = checked_ts_count;

  index->stats->query_total_loaded_nodes_count =
      loaded_nodes_count + index->stats->query_filter_loaded_nodes_count;
  index->stats->query_total_loaded_ts_count =
      loaded_ts_count + index->stats->query_filter_loaded_ts_count;
  index->stats->query_total_checked_nodes_count =
      checked_nodes_count + index->stats->query_filter_checked_nodes_count;
  index->stats->query_total_checked_ts_count =
      checked_ts_count + index->stats->query_filter_checked_ts_count;

  index->stats->queries_refine_total_time +=
      index->stats->query_refine_total_time;

  index->stats->queries_refine_input_time += partial_input_time;
  index->stats->queries_refine_output_time += partial_output_time;
  index->stats->queries_refine_load_node_time += partial_load_node_time;
  index->stats->queries_refine_cpu_time +=
      partial_time - partial_input_time - partial_output_time;
  index->stats->queries_refine_seq_input_count += partial_seq_input_count;
  index->stats->queries_refine_seq_output_count += partial_seq_output_count;
  index->stats->queries_refine_rand_input_count += partial_rand_input_count;
  index->stats->queries_refine_rand_output_count += partial_rand_output_count;

  index->stats->queries_total_input_time =
      index->stats->queries_refine_input_time +
      index->stats->queries_filter_input_time;
  index->stats->queries_total_output_time =
      index->stats->queries_refine_output_time +
      index->stats->queries_filter_output_time;
  index->stats->queries_total_load_node_time =
      index->stats->queries_refine_load_node_time +
      index->stats->queries_filter_load_node_time;
  index->stats->queries_total_cpu_time = index->stats->queries_refine_cpu_time +
                                         index->stats->queries_filter_cpu_time;

  index->stats->queries_total_time = index->stats->queries_refine_total_time +
                                     index->stats->queries_filter_total_time;

  index->stats->queries_total_seq_input_count =
      index->stats->queries_filter_seq_input_count +
      index->stats->queries_refine_seq_input_count;
  index->stats->queries_total_seq_output_count =
      index->stats->queries_filter_seq_output_count +
      index->stats->queries_refine_seq_output_count;
  index->stats->queries_total_rand_input_count =
      index->stats->queries_filter_rand_input_count +
      index->stats->queries_refine_rand_input_count;
  index->stats->queries_total_rand_output_count =
      index->stats->queries_filter_rand_output_count +
      index->stats->queries_refine_rand_output_count;

  // keep a running sum then divide by the total number of queries
  index->stats->queries_avg_checked_nodes_count +=
      index->stats->query_total_checked_nodes_count;
  index->stats->queries_avg_checked_ts_count +=
      index->stats->query_total_checked_ts_count;
  index->stats->queries_avg_loaded_nodes_count +=
      index->stats->query_total_loaded_nodes_count;
  index->stats->queries_avg_loaded_ts_count +=
      index->stats->query_total_loaded_ts_count;

  // COUNT_TOTAL_TIME_START
  return bsf_result;
}

struct query_result *
exact_de_knn_search(ts_type *query_ts, ts_type *query_ts_reordered,
                    int *query_order, unsigned int offset,
                    struct dstree_index *index, ts_type minimum_distance,
                    ts_type epsilon, ts_type r_delta, unsigned int k,
                    unsigned int q_id, char *qfilename) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary

  while ((n = pqueue_pop(pq))) {
    temp = knn_results[k - 1];
    kth_bsf = temp.distance;
    if (n->distance > kth_bsf / (1 + epsilon)) {
      break;
    }

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, kth_bsf, k, knn_results,
                                  NULL, NULL, &curr_size, 0);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  free(knn_results);
}

void exact_ng_knn_search(ts_type *query_ts, ts_type *query_ts_reordered,
                         int *query_order, unsigned int offset,
                         struct dstree_index *index, ts_type minimum_distance,
                         unsigned int k, unsigned int q_id, char *qfilename,
                         unsigned int nprobes) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  unsigned int cur_probes = 0;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  ++cur_probes;

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    update_query_stats(index, q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary
  // cur_probes is strictly less than nprobes because an approximate answer has
  // already been found
  while ((n = pqueue_pop(pq)) && (cur_probes < nprobes)) {

    if (n->distance > bsf_result.distance) {
      break;
    }

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, NULL, NULL, &curr_size, 0);

      // increase the number of visited leaves
      ++cur_probes;

    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  // free the results, eventually do something with them!!

  free(knn_results);
}

void exact_knn_search_max_policy(ts_type *query_ts, ts_type *query_ts_reordered,
                                 int *query_order, unsigned int offset,
                                 struct dstree_index *index,
                                 ts_type minimum_distance, ts_type epsilon,
                                 ts_type delta, unsigned int k,
                                 unsigned int q_id, char *qfilename) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // pqueue_t *knn_results= pqueue_init(k, cmp_pri,
  //			       get_pri, set_pri,
  //				       get_pos, set_pos);

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;
    update_query_stats(index, q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_max_pri,
                             set_max_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);
  root_pq_item->max_distance =
      calculate_node_max_distance(index, index->first_node, query_ts);
  ;

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary
  while ((n = pqueue_pop(pq))) {

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    kth_bsf = knn_results[k - 1].distance;

    if (n->distance > kth_bsf / (1 + epsilon)) // add epsilon+1
    {
      // we cannot break because a node with lower lb can still be in the queue
      continue;
    }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, NULL, NULL, &curr_size, 0);

    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      ts_type child_max_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);
      child_max_distance =
          calculate_node_max_distance(index, n->node->left_child, query_ts);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        mindist_result_left->max_distance = child_max_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);
      child_max_distance =
          calculate_node_max_distance(index, n->node->right_child, query_ts);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        mindist_result_right->max_distance = child_max_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }

  // report all the elements at once since algo cannot be incremental with
  // maxdist policy
  for (unsigned int pos = 1; pos <= k; ++pos) {
    bsf_result = knn_results[pos - 1];
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, pos, bsf_result);
    get_query_stats(index, pos);
    print_query_stats(index, q_id, pos, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  // free the results, eventually do something with them!!

  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);
  free(knn_results);
}

void exact_knn_search_track_pruning(ts_type *query_ts,
                                    ts_type *query_ts_reordered,
                                    int *query_order, unsigned int offset,
                                    struct dstree_index *index,
                                    ts_type minimum_distance, ts_type epsilon,
                                    ts_type delta, unsigned int k,
                                    unsigned int q_id, char *qfilename) {
  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;
    update_query_stats(index, q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary
  while ((n = pqueue_pop(pq))) {

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1
    for (unsigned int pos = found_knn; pos < k; ++pos) {
      bsf_result = knn_results[pos];

      if (n->distance > bsf_result.distance / (1 + epsilon)) // add epsilon+1
      {
        found_knn = pos + 1;
        COUNT_PARTIAL_TIME_END

        update_query_stats(index, q_id, found_knn, bsf_result);
        get_query_stats(index, found_knn);
        // print_query_stats(index, q_id, found_knn,qfilename);

        // reset the bsf for the next NN
        if (found_knn < k) {
          bsf_result = knn_results[found_knn];
        }

        RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
      }
    }

    if (found_knn == k) {
      // printf("found all kNN\n");
      break;
    }
    // report the pos-NN neighbors, then continue

    // if (n->distance > kth_bsf/(1 + epsilon))
    //{
    //   break;
    // }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, NULL, NULL, &curr_size, 0);

    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      } else if ((n->node->left_child != approximate_result.node)) {
        print_pruning_snapshots(n->node->left_child, kth_bsf, child_distance,
                                found_knn + 1, q_id, qfilename);
      }
      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      } else if ((n->node->right_child != approximate_result.node)) {
        print_pruning_snapshots(n->node->right_child, kth_bsf, child_distance,
                                found_knn + 1, q_id, qfilename);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  // free the results, eventually do something with them!!

  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.

  pqueue_free(pq);
  free(knn_results);
}

void dump_mindists(struct dstree_index *index, struct dstree_node *node,
                   ts_type *query_ts) {

  ts_type distance;
  ts_type QoS;

  distance = calculate_node_min_distance(index, node, query_ts, 0);
  QoS = get_node_QoS(index, node);

  printf("%*s%lf\t%d\t%d\t%lf\n", node->level, "", sqrtf(distance),
         node->num_node_points, node->level, sqrtf(QoS));

  if (!node->is_leaf) {
    dump_mindists(index, node->left_child, query_ts);
    dump_mindists(index, node->right_child, query_ts);
  }
}

void exact_knn_search_track_bsf(ts_type *query_ts, ts_type *query_ts_reordered,
                                int *query_order, unsigned int offset,
                                struct dstree_index *index,
                                ts_type minimum_distance, ts_type epsilon,
                                ts_type delta, unsigned int k,
                                unsigned int q_id, char *qfilename,
                                struct bsf_snapshot **bsf_snapshots,
                                unsigned int *cur_bsf_snapshot) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;
  unsigned int last_found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // tracking bsf steps for all kNNs

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, bsf_snapshots, cur_bsf_snapshot,
                         &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END

  for (int idx = 0; idx < k; ++idx) {
    bsf_snapshots[idx][*cur_bsf_snapshot].distance = knn_results[idx].distance;
    bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
    bsf_snapshots[idx][*cur_bsf_snapshot].series = NULL;
    bsf_snapshots[idx][*cur_bsf_snapshot].checked_nodes = checked_nodes_count;
  }
  ++(*cur_bsf_snapshot);

  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;
    update_query_stats(index, q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  while ((n = pqueue_pop(pq))) {

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1
    for (unsigned int pos = found_knn; pos < k; ++pos) {
      bsf_result = knn_results[pos];

      if (n->distance > bsf_result.distance / (1 + epsilon)) // add epsilon+1
      {
        last_found_knn = found_knn;
        found_knn = pos + 1;
        COUNT_PARTIAL_TIME_END

        update_query_stats(index, q_id, found_knn, bsf_result);
        get_query_stats(index, found_knn);
        print_query_stats(index, q_id, found_knn, qfilename);

        // get a snapshot of the bsfs at this point in time

        for (int idx = 0; idx < k; ++idx) {
          bsf_snapshots[idx][*cur_bsf_snapshot].distance =
              knn_results[idx].distance;
          bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
          bsf_snapshots[idx][*cur_bsf_snapshot].checked_nodes =
              checked_nodes_count;
        }
        ++(*cur_bsf_snapshot);

        // reset the bsf for the next NN
        if (found_knn < k) {
          bsf_result = knn_results[found_knn];
        }

        RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
        // report all results for found_knn - last_found_knn or print their
        // results
      }
    }

    if (found_knn == k) {
      // printf("found all kNN\n");
      break;
    }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, bsf_snapshots, cur_bsf_snapshot,
                                  &curr_size, 0);
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);

    print_query_stats(index, q_id, found_knn, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  for (int idx = 0; idx < k; ++idx) {
    bsf_snapshots[idx][*cur_bsf_snapshot].distance = knn_results[idx].distance;
    bsf_snapshots[idx][*cur_bsf_snapshot].time = partial_time;
    bsf_snapshots[idx][*cur_bsf_snapshot].checked_nodes = checked_nodes_count;
  }
  ++(*cur_bsf_snapshot);

  print_bsf_snapshots(index, q_id, k, qfilename, bsf_snapshots,
                      *cur_bsf_snapshot);

  // free the results, eventually do something with them!!

  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);
  free(knn_results);
}

void exact_de_progressive_knn_search(
    ts_type *query_ts, ts_type *query_ts_reordered, int *query_order,
    unsigned int offset, struct dstree_index *index, ts_type minimum_distance,
    ts_type epsilon, ts_type r_delta, unsigned int k, unsigned int q_id,
    char *qfilename, struct bsf_snapshot **bsf_snapshots,
    unsigned int *cur_bsf_snapshot) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, bsf_snapshots, cur_bsf_snapshot,
                         &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  // RESET_QUERY_COUNTERS()
  // RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  while ((n = pqueue_pop(pq))) {
    temp = knn_results[k - 1];
    kth_bsf = temp.distance;
    if (n->distance > kth_bsf / (1 + epsilon)) {
      break;
    }

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, bsf_snapshots, cur_bsf_snapshot,
                                  &curr_size, 0);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  COUNT_PARTIAL_TIME_END
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already

  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    // COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    print_perk_progressive_bsf_snapshots(index, q_id, found_knn, qfilename,
                                         bsf_snapshots, *cur_bsf_snapshot,
                                         bsf_result.distance, NULL, NULL, NULL);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    // COUNT_PARTIAL_TIME_START
  }
  /*
    for (int idx = 0; idx < k; ++idx)
    {
      bsf_snapshots[idx][*cur_bsf_snapshot].distance =
    knn_results[idx].distance; bsf_snapshots[idx][*cur_bsf_snapshot].time =
    partial_time;
    }
    ++(*cur_bsf_snapshot);
  */
  // print_progressive_bsf_snapshots(index, q_id,k,qfilename,bsf_snapshots,
  // *cur_bsf_snapshot);

  // free the results, eventually do something with them!!

  free(knn_results);
}
void exact_de_incr_progressive_knn_search(
    ts_type *query_ts, ts_type *query_ts_reordered, int *query_order,
    unsigned int offset, struct dstree_index *index, ts_type minimum_distance,
    ts_type epsilon, ts_type r_delta, unsigned int k, unsigned int q_id,
    char *qfilename, struct bsf_snapshot **bsf_snapshots,
    unsigned int *cur_bsf_snapshot, float warping, FILE *dataset_file,
    FILE *series_file) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  ts_type *series =
      calloc(1, sizeof(ts_type) * index->settings->timeseries_size);
  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, bsf_snapshots, cur_bsf_snapshot,
                         &curr_size, warping);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  // RESET_QUERY_COUNTERS()
  // RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, warping);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // FILE *series_file = fopen(filename, "a");
  // FILE *dataset_file = fopen(index->settings->dataset, "rb");

  while ((n = pqueue_pop(pq))) {

    // the first element of the queue is not used, thus pos-1
    for (unsigned int pos = found_knn; pos < k; ++pos) {
      bsf_result = knn_results[pos];
      // printf("n->distance = %g, bsf_result.distance = %g\n",
      // sqrt(n->distance), sqrt(bsf_result.distance));
      if (n->distance > bsf_result.distance) // add epsilon+1
      {
        found_knn = pos + 1;
        COUNT_PARTIAL_TIME_END

        update_query_stats(index, q_id, found_knn, bsf_result);
        get_query_stats(index, found_knn);
        print_query_stats(index, q_id, found_knn, qfilename);
        print_perk_progressive_bsf_snapshots(
            index, q_id, found_knn, qfilename, bsf_snapshots, *cur_bsf_snapshot,
            bsf_result.distance, dataset_file, series_file, series);
        // print_perk_progressive_bsf_snapshots(index,
        // q_id,found_knn,qfilename,bsf_snapshots, *cur_bsf_snapshot,
        // bsf_result.distance, NULL, NULL); printf("found NN = %u\n",
        // found_knn); fflush(stdout); reset the bsf for the next NN
        if (found_knn < k) {
          bsf_result = knn_results[found_knn];
        }

        // RESET_QUERY_COUNTERS()
        // RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
      }
    }

    if (found_knn >= k) {
      // printf("found all kNN\n");
      // fflush(stdout);
      break;
    }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, bsf_snapshots, cur_bsf_snapshot,
                                  &curr_size, warping);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance = calculate_node_min_distance(index, n->node->left_child,
                                                   query_ts, warping);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      // if ((child_distance < kth_bsf/(1+epsilon)) &&
      if ((child_distance < kth_bsf) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance = calculate_node_min_distance(index, n->node->right_child,
                                                   query_ts, warping);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      // if ((child_distance < kth_bsf/(1+epsilon))  &&
      if ((child_distance < kth_bsf) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    print_perk_progressive_bsf_snapshots(
        index, q_id, found_knn, qfilename, bsf_snapshots, *cur_bsf_snapshot,
        bsf_result.distance, dataset_file, series_file, series);
    // print_perk_progressive_bsf_snapshots(index,
    // q_id,found_knn,qfilename,bsf_snapshots, *cur_bsf_snapshot,
    // bsf_result.distance, NULL, NULL); report all results for found_knn -
    // last_found_knn or print their results RESET_QUERY_COUNTERS()
    // RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  /*
    for (int idx = 0; idx < k; ++idx)
    {
      bsf_snapshots[idx][*cur_bsf_snapshot].distance =
    knn_results[idx].distance; bsf_snapshots[idx][*cur_bsf_snapshot].time =
    partial_time;
    }
    ++(*cur_bsf_snapshot);
  */
  // print_progressive_bsf_snapshots(index, q_id,k,qfilename,bsf_snapshots,
  // *cur_bsf_snapshot);

  // free the results, eventually do something with them!!

  // fclose(series_file);
  // fclose(dataset_file);
  free(series);
  free(knn_results);
}

/* start kashif changes */
struct query_result *exact_de_incr_progressive_knn_search_2(
    ts_type *query_ts, ts_type *query_ts_reordered, int *query_order,
    unsigned int offset, struct dstree_index *index, ts_type minimum_distance,
    ts_type epsilon, ts_type r_delta, unsigned int k, unsigned int q_id,
    char *qfilename, double *total_query_set_time,
    unsigned int *total_checked_ts, struct bsf_snapshot **bsf_snapshots,
    unsigned int *cur_bsf_snapshot, float warping, FILE *dataset_file,
    FILE *series_file, unsigned int query_vector_pos)
{

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  ts_type *series =
      calloc(1, sizeof(ts_type) * index->settings->timeseries_size);
  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
    knn_results[idx].vector_id = malloc(sizeof(struct vid));
    knn_results[idx].vector_id->table_id = 1000000000000;
    knn_results[idx].vector_id->set_id = 1000000000000;
    knn_results[idx].vector_id->pos = 1000000000000;
    knn_results[idx].query_vector_pos = query_vector_pos;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, bsf_snapshots, cur_bsf_snapshot,
                         &curr_size, warping);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  // RESET_QUERY_COUNTERS()
  // RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, warping);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // FILE *series_file = fopen(filename, "a");
  // FILE *dataset_file = fopen(index->settings->dataset, "rb");

  while ((n = pqueue_pop(pq))) {

    // the first element of the queue is not used, thus pos-1
    for (unsigned int pos = found_knn; pos < k; ++pos) {
      bsf_result = knn_results[pos];
      // printf("n->distance = %g, bsf_result.distance = %g\n",
      // sqrt(n->distance), sqrt(bsf_result.distance));
      if (n->distance > bsf_result.distance) // add epsilon+1
      {
        found_knn = pos + 1;
        COUNT_PARTIAL_TIME_END

        update_query_stats(index, q_id, found_knn, bsf_result);
        // get_query_stats(index, found_knn);
        // print_query_stats(index, q_id, found_knn, qfilename);


        // printf("-- start knn -- -- --- -- -- -- -- -- -- -- -- --\n");

        // print_perk_progressive_bsf_snapshots(
        //     index, q_id, found_knn, qfilename, bsf_snapshots, *cur_bsf_snapshot,
        //     bsf_result.distance, dataset_file, series_file, series);

        // printf("-- end knn -- -- --- -- -- -- -- -- -- -- -- --\n");

        // print_perk_progressive_bsf_snapshots(index,
        // q_id,found_knn,qfilename,bsf_snapshots, *cur_bsf_snapshot,
        // bsf_result.distance, NULL, NULL); printf("found NN = %u\n",
        // found_knn); fflush(stdout); reset the bsf for the next NN
        if (found_knn < k) {
          bsf_result = knn_results[found_knn];
        }

        // RESET_QUERY_COUNTERS()
        // RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
      }
    }

    if (found_knn >= k) {
      // printf("found all kNN\n");
      // fflush(stdout);
      break;
    }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, bsf_snapshots, cur_bsf_snapshot,
                                  &curr_size, warping);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance = calculate_node_min_distance(index, n->node->left_child,
                                                   query_ts, warping);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      // if ((child_distance < kth_bsf/(1+epsilon)) &&
      if ((child_distance < kth_bsf) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance = calculate_node_min_distance(index, n->node->right_child,
                                                   query_ts, warping);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      // if ((child_distance < kth_bsf/(1+epsilon))  &&
      if ((child_distance < kth_bsf) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn, qfilename);

    // printf("-- start knn -- -- --- -- -- -- -- -- -- -- -- --\n");

    // print_perk_progressive_bsf_snapshots(
    //     index, q_id, found_knn, qfilename, bsf_snapshots, *cur_bsf_snapshot,
    //     bsf_result.distance, dataset_file, series_file, series);
    
    // printf("-- end knn -- -- --- -- -- -- -- -- -- -- -- --\n");
    // print_perk_progressive_bsf_snapshots(index,
    // q_id,found_knn,qfilename,bsf_snapshots, *cur_bsf_snapshot,
    // bsf_result.distance, NULL, NULL); report all results for found_knn -
    // last_found_knn or print their results RESET_QUERY_COUNTERS()
    // RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  /*
    for (int idx = 0; idx < k; ++idx)
    {
      bsf_snapshots[idx][*cur_bsf_snapshot].distance =
    knn_results[idx].distance; bsf_snapshots[idx][*cur_bsf_snapshot].time =
    partial_time;
    }
    ++(*cur_bsf_snapshot);
  */
  // print_progressive_bsf_snapshots(index, q_id,k,qfilename,bsf_snapshots,
  // *cur_bsf_snapshot);

  // free the results, eventually do something with them!!

  // fclose(series_file);
  // fclose(dataset_file);

  
  // add  query_vector time to total query_set time
  *total_query_set_time += index->stats->query_filter_total_time / 1000000;
  *total_checked_ts += index->stats->query_filter_checked_ts_count;
  free(series);

  // save results to later find top-k matches
  return knn_results;
}

struct query_result *exact_de_knn_search_2(
    ts_type *query_ts, ts_type *query_ts_reordered, int *query_order,
    unsigned int offset, struct dstree_index *index, ts_type minimum_distance,
    ts_type epsilon, ts_type r_delta, unsigned int k, unsigned int q_id,
    char *qfilename, double *total_query_set_time,
    unsigned int *total_checked_ts, unsigned int query_vector_pos)
{

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
    knn_results[idx].vector_id = malloc(sizeof(struct vid));
    knn_results[idx].vector_id->table_id = 1000000000000;
    knn_results[idx].vector_id->set_id = 1000000000000;
    knn_results[idx].vector_id->pos = 1000000000000;
    knn_results[idx].query_vector_pos = query_vector_pos;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  /* start kashif changes (to fix)*/
  // sig fault in this line
  // index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  // index->stats->query_approx_node_filename = approximate_result.node->filename;
  // index->stats->query_approx_node_size = approximate_result.node->node_size;
  // ;
  // index->stats->query_approx_node_level = approximate_result.node->level;
  /* end kashif changes */

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary

  while ((n = pqueue_pop(pq))) {
    temp = knn_results[k - 1];
    kth_bsf = temp.distance;
    if (n->distance > kth_bsf / (1 + epsilon)) {
      break;
    }

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, kth_bsf, k, knn_results,
                                  NULL, NULL, &curr_size, 0);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  // add  query_vector time to total query_set time
  *total_query_set_time += index->stats->query_filter_total_time / 1000000;
  *total_checked_ts += index->stats->query_filter_checked_ts_count;

  // save results to later find top-k matches
  return knn_results;
}

struct query_result *exact_de_progressive_knn_search_2(
    ts_type *query_ts, ts_type *query_ts_reordered, int *query_order,
    unsigned int offset, struct dstree_index *index, ts_type minimum_distance,
    ts_type epsilon, ts_type r_delta, unsigned int k, unsigned int q_id,
    char *qfilename, double *total_query_set_time, unsigned int *total_checked_ts,
    struct bsf_snapshot **bsf_snapshots,
    unsigned int *cur_bsf_snapshot, unsigned int query_vector_pos) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
    knn_results[idx].vector_id = malloc(sizeof(struct vid));
    knn_results[idx].vector_id->table_id = 1000000000000;
    knn_results[idx].vector_id->set_id = 1000000000000;
    knn_results[idx].vector_id->pos = 1000000000000;
    knn_results[idx].query_vector_pos = query_vector_pos;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, bsf_snapshots, cur_bsf_snapshot,
                         &curr_size, 0);

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    // update_query_stats(index,q_id, found_knn, approximate_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn,qfilename);
  }

  // RESET_QUERY_COUNTERS()
  // RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  while ((n = pqueue_pop(pq))) {
    temp = knn_results[k - 1];
    kth_bsf = temp.distance;
    if (n->distance > kth_bsf / (1 + epsilon)) {
      break;
    }

    // get the kth-distance from current bsfs
    // kth_bsf = pqueue_peek_last(knn_results)->distance;
    // temp_bsf = kth_bsf;

    // the first element of the queue is not used, thus pos-1

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, bsf_snapshots, cur_bsf_snapshot,
                                  &curr_size, 0);

      // if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      // epsilon)))
      //  break;

      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf / (1 + epsilon)) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }
  COUNT_PARTIAL_TIME_END
  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);

  // report the elements that were not reported already

  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    // COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    // get_query_stats(index, found_knn);
    // print_query_stats(index, q_id, found_knn, qfilename);
    print_perk_progressive_bsf_snapshots(index, q_id, found_knn, qfilename,
                                         bsf_snapshots, *cur_bsf_snapshot,
                                         bsf_result.distance, NULL, NULL, NULL);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    // COUNT_PARTIAL_TIME_START
  }
  /*
    for (int idx = 0; idx < k; ++idx)
    {
      bsf_snapshots[idx][*cur_bsf_snapshot].distance =
    knn_results[idx].distance; bsf_snapshots[idx][*cur_bsf_snapshot].time =
    partial_time;
    }
    ++(*cur_bsf_snapshot);
  */
  // print_progressive_bsf_snapshots(index, q_id,k,qfilename,bsf_snapshots,
  // *cur_bsf_snapshot);

  // add  query_vector time to total query_set time
  *total_query_set_time += index->stats->query_filter_total_time / 1000000;
  *total_checked_ts += index->stats->query_filter_checked_ts_count;

  // save results to later find top-k matches
  return knn_results;
}
/* end kashif changes */

void dstree_calc_tlb(ts_type *query_ts, struct dstree_index *index,
                     struct dstree_node *curr_node) {

  ts_type curr_lb_dist = 0;
  ts_type curr_exact_dist = 0;

  // curr_node = index->first_node
  // curr_dist = calculate_node_min_distance (index, index->first_node,
  // query_ts);

  // This is an internal node, traverse left and right children
  if (curr_node == NULL) {
    return;
  }

  if (!curr_node->is_leaf) {
    dstree_calc_tlb(query_ts, index, curr_node->left_child);
    dstree_calc_tlb(query_ts, index, curr_node->right_child);
  }
  // This is a leaf, calculate the tlb = lb_distance/exact_distance
  else {
    curr_lb_dist = calculate_node_min_distance(index, curr_node, query_ts, 0);
    if (curr_node->file_buffer->buffered_list_size == 0) {
      curr_node->file_buffer->buffered_list =
          get_all_time_series_in_node(index, curr_node);
      curr_node->file_buffer->buffered_list_size =
          curr_node->file_buffer->disk_count;

      if (curr_node->file_buffer->buffered_list == NULL) {
        fprintf(stderr,
                "Error in dstree_index.c:  Could not retrieve all time series "
                "for node %s.\n",
                curr_node->filename);
      }
    }
    total_ts_count =
        total_ts_count + curr_node->file_buffer->buffered_list_size;
    ++leaf_nodes_count;

    for (int idx = 0; idx < curr_node->file_buffer->buffered_list_size; ++idx) {
      curr_exact_dist = ts_euclidean_distance(
          query_ts, curr_node->file_buffer->buffered_list[idx],
          index->settings->timeseries_size);
      // printf("Leaf node: %ul exact_distance = %g lb_distance = %g\n",
      // leaf_nodes_count, curr_exact_dist, curr_lb_dist);
      total_tlb += sqrtf(curr_lb_dist / curr_exact_dist);
    }
  }
}

void get_query_stats(struct dstree_index *index, unsigned int found_knn) {

  if (total_ts_count != 0) {
    index->stats->query_pruning_ratio =
        1.0 - ((double)index->stats->query_total_checked_ts_count /
               index->stats->total_ts_count);
  }

  if (found_knn == 1) {
    if (index->stats->query_exact_distance != 0) {
      index->stats->query_eff_epsilon = (index->stats->query_approx_distance -
                                         index->stats->query_exact_distance) /
                                        index->stats->query_exact_distance;
    }
  }
}

void print_progressive_bsf_snapshots(struct dstree_index *index,
                                     unsigned int query_num, unsigned int k,
                                     char *queries,
                                     struct bsf_snapshot **bsf_snapshots,
                                     unsigned int cur_bsf_snapshot) {
  /*
  const char *filename = malloc(sizeof(char) *
  (strlen(index->settings->root_directory) + 18)); filename = strcpy(filename,
  index->settings->root_directory); filename = strcat(filename,
  "../raw_series.csv\0");

  FILE *file = fopen(filename, "a");
  */
  unsigned int ts_length = index->settings->timeseries_size;

  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < cur_bsf_snapshot; ++j) {
      /*
      fprintf(file,"%u,%u",query_num-1,j);
      for (unsigned int l = 0; l < ts_length; ++l)
        fprintf(file,",%g",bsf_snapshots[i][j].series[l]);
      fprintf(file,"\n");
      */
      if (j > 0) {
        if (bsf_snapshots[i][j - 1].distance == bsf_snapshots[i][j].distance) {
          break;
        }
      }
      printf("Query_bsf_snapshot_time_secs\t%lf\t%s\t%u\t%u\n",
             bsf_snapshots[i][j].time / 1000000, queries, query_num, i + 1);

      printf("Query_bsf_snapshot_distance\t%lf\t%s\t%u\t%u\n",
             sqrtf(bsf_snapshots[i][j].distance), queries, query_num, i + 1);
      printf("Query_bsf_snapshot_checked_nodes_count\t%lu\t%s\t%u\t%u\n",
             bsf_snapshots[i][j].checked_nodes, queries, query_num, i + 1);
    }
  }

  // fclose(file);
}
void print_perk_progressive_bsf_snapshots(
    struct dstree_index *index, unsigned int query_num, unsigned int found_knn,
    char *queries, struct bsf_snapshot **bsf_snapshots,
    unsigned int cur_bsf_snapshot, ts_type exact_distance, FILE *dataset_file,
    FILE *series_file, ts_type *series) {

  unsigned int ts_length = index->settings->timeseries_size;
  unsigned long file_pos_raw;

  for (unsigned int j = 0; j < cur_bsf_snapshot; ++j) {

    /*
    if (j > 0)
    {
      if ((bsf_snapshots [found_knn-1][j].distance == exact_distance) &&
    (bsf_snapshots[found_knn-1][j-1].distance == exact_distance))
        {
          break;
        }
    }
    */
    if (j > 0) {
      if ((bsf_snapshots[found_knn - 1][j].distance ==
           bsf_snapshots[found_knn - 1][j - 1].distance)) {
        continue;
      }
    }
    printf("(***)\n");
    printf("Query_bsf_snapshot_time_secs\t%lf\t%s\t%u\t%u\n",
           bsf_snapshots[found_knn - 1][j].time / 1000000, queries, query_num,
           found_knn);

    printf("Query_bsf_snapshot_distance\t%lf\t%s\t%u\t%u\n",
           sqrtf(bsf_snapshots[found_knn - 1][j].distance), queries, query_num,
           found_knn);
    printf("Query_bsf_snapshot_checked_nodes_count\t%lu\t%s\t%u\t%u\n",
           bsf_snapshots[found_knn - 1][j].checked_nodes, queries, query_num,
           found_knn);
    if (index->settings->classify) {
      printf("Query_bsf_snapshot_label\t%u\t%s\t%u\t%u\n",
             (unsigned int)bsf_snapshots[found_knn - 1][j].label, queries,
             query_num, found_knn);
    }
    if (index->settings->track_file_pos) {
      printf("Query_bsf_snapshot_file_pos\t%u\t%s\t%u\t%u\n",
             (unsigned int)bsf_snapshots[found_knn - 1][j].file_pos, queries,
             query_num, found_knn);
    
    if (index->settings->track_vector) {
      printf("Query_bsf_snapshot_vector_id\t(%u, %u)\t%s\t%u\t%u\n",
             (unsigned int)bsf_snapshots[found_knn - 1][j].vector_id->table_id, 
             (unsigned int)bsf_snapshots[found_knn - 1][j].vector_id->set_id, queries,
             query_num, found_knn);
    }
    printf("(***)\n");
      // printf("%u,%u,%u",query_num-1,found_knn-1,j);
      // printf("COUCOU");
      /*
        file_pos_raw = (((unsigned long) bsf_snapshots[found_knn-1][j].file_pos)
    - 1) * ts_length * sizeof(ts_type); fseek(dataset_file, 0L, SEEK_SET);
        fseek(dataset_file, file_pos_raw , SEEK_SET);
        fread(series, sizeof(ts_type), ts_length, dataset_file);
        fprintf(series_file,"%u,%u,%u,%u",query_num,found_knn,j,(unsigned int)
    bsf_snapshots[found_knn-1][j].file_pos);

    for (unsigned int l = 0; l < ts_length; ++l)
    {
        fprintf(series_file,",%g",series[l]);
       //printf(",%g",bsf_snapshots[found_knn-1][j].series[l]);
    }
    fprintf(series_file,"\n");
     */
      // printf("\n");
    }
  }
}

void print_bsf_snapshots(struct dstree_index *index, unsigned int query_num,
                         unsigned int k, char *queries,
                         struct bsf_snapshot **bsf_snapshots,
                         unsigned int cur_bsf_snapshot)
{

  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < cur_bsf_snapshot; ++j) {
      printf("Query_bsf_snapshot_time_secs\t%lf\t%s\t%u\t%u\n",
             bsf_snapshots[i][j].time / 1000000, queries, query_num, i + 1);
      printf("Query_bsf_snapshot_distance\t%lf\t%s\t%u\t%u\n",
             sqrtf(bsf_snapshots[i][j].distance), queries, query_num, i + 1);
    }
  }
}

void print_pruning_snapshots(struct dstree_node *node, ts_type node_bsf,
                             ts_type node_mindist, unsigned int k,
                             unsigned int query_num, char *queries) {

  printf("Query_pruning_snapshot_node_filename\t%s\t%s\t%u\t%u\n",
         node->filename, queries, query_num, k);

  printf("Query_pruning_snapshot_node_level\t%u\t%s\t%u\t%u\n", node->level,
         queries, query_num, k);

  printf("Query_pruning_snapshot_node_bsf\t%lf\t%s\t%u\t%u\n", node_bsf,
         queries, query_num, k);

  printf("Query_pruning_snapshot_node_mindist\t%lf\t%s\t%u\t%u\n", node_mindist,
         queries, query_num, k);
}

void print_query_stats(struct dstree_index *index, unsigned int query_num,
                       unsigned int found_knn, char *queries) {

  printf("Query_total_input_time_secs\t%lf\t%s\t%u\t%u\n",
         index->stats->query_total_input_time / 1000000, queries, query_num,
         found_knn);

  printf("Query_total_output_time_secs\t%lf\t%s\t%u\t%u\n",
         index->stats->query_total_output_time / 1000000, queries, query_num,
         found_knn);

  printf("Query_total_load_node_time_secs\t%lf\t%s\t%u\t%u\n",
         index->stats->query_total_load_node_time / 1000000, queries, query_num,
         found_knn);

  printf("Query_total_cpu_time_secs\t%lf\t%s\t%u\t%u\n",
         index->stats->query_total_cpu_time / 1000000, queries, query_num,
         found_knn);

  printf("Query_total_time_secs\t%lf\t%s\t%u\t%u\n",
         index->stats->query_total_time / 1000000, queries, query_num,
         found_knn);

  printf("Query_total_seq_input_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_seq_input_count, queries, query_num,
         found_knn);

  printf("Query_total_seq_output_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_seq_output_count, queries, query_num,
         found_knn);

  printf("Query_total_rand_input_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_rand_input_count, queries, query_num,
         found_knn);

  printf("Query_total_rand_output_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_rand_output_count, queries, query_num,
         found_knn);

  printf("Query_total_checked_nodes_count\t%u\t%s\t%u\t%u\n",
         index->stats->query_total_checked_nodes_count, queries, query_num,
         found_knn);

  printf("Query_total_checked_ts_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_checked_ts_count, queries, query_num,
         found_knn);

  printf("Query_total_loaded_nodes_count\t%u\t%s\t%u\t%u\n",
         index->stats->query_total_loaded_nodes_count, queries, query_num,
         found_knn);

  printf("Query_total_loaded_ts_count\t%llu\t%s\t%u\t%u\n",
         index->stats->query_total_loaded_ts_count, queries, query_num,
         found_knn);

  printf("Query_exact_distance\t%f\t%s\t%u\t%u\n",
         index->stats->query_exact_distance, queries, query_num, found_knn);
  if (index->settings->classify) {
    printf("Query_label\t%u\t%s\t%u\t%u\n",
           (unsigned int)index->stats->query_label, queries, query_num,
           found_knn);
  }
  if (index->settings->track_file_pos) {
    printf("Query_file_pos\t%u\t%s\t%u\t%u\n",
           (unsigned int)index->stats->query_file_pos, queries, query_num,
           found_knn);
  }
   if (index->settings->track_vector)
  {
    printf("Query_vector_id\t(%u, %u)\t%s\t%u\t%u\n",
          index->stats->query_vector_id->table_id, index->stats->query_vector_id->set_id,
          queries, query_num, found_knn);
  } 
  printf("Query_exact_node_filename\t%s\t%s\t%u\t%u\n",
         index->stats->query_exact_node_filename, queries, query_num,
         found_knn);

  printf("Query_exact_node_size\t%u\t%s\t%u\t%u\n",
         index->stats->query_exact_node_size, queries, query_num, found_knn);

  printf("Query_exact_node_level\t%u\t%s\t%u\t%u\n",
         index->stats->query_exact_node_level, queries, query_num, found_knn);

  printf("Query_pruning_ratio_level\t%f\t%s\t%u\t%u\n",
         index->stats->query_pruning_ratio, queries, query_num, found_knn);
}

void update_query_stats(struct dstree_index *index, unsigned int query_id,
                        unsigned int found_knn,
                        struct query_result bsf_result) {

  index->stats->query_total_time = partial_time;
  index->stats->query_total_input_time = partial_input_time;
  index->stats->query_total_output_time = partial_output_time;
  index->stats->query_total_load_node_time = partial_load_node_time;

  index->stats->query_total_cpu_time = index->stats->query_total_time -
                                       index->stats->query_total_input_time -
                                       index->stats->query_total_output_time;

  index->stats->query_total_seq_input_count = partial_seq_input_count;
  index->stats->query_total_seq_output_count = partial_seq_output_count;
  index->stats->query_total_rand_input_count = partial_rand_input_count;
  index->stats->query_total_rand_output_count = partial_rand_output_count;

  index->stats->query_total_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_total_loaded_ts_count = loaded_ts_count;
  index->stats->query_total_checked_nodes_count = checked_nodes_count;
  index->stats->query_total_checked_ts_count = checked_ts_count;

  index->stats->query_exact_distance = sqrtf(bsf_result.distance);
  if (index->settings->classify)
    index->stats->query_label = bsf_result.label;
  if (index->settings->track_file_pos)
    index->stats->query_file_pos = bsf_result.file_pos;
  if (index->settings->track_vector)
    index->stats->query_vector_id = bsf_result.vector_id;
  

  // index->stats->query_exact_node_filename = bsf_result.node->filename;
  // index->stats->query_exact_node_size = bsf_result.node->node_size;;
  // index->stats->query_exact_node_level = bsf_result.node->level;
}

ts_type get_node_QoS(struct dstree_index *index, struct dstree_node *node) {

  ts_type node_range_value = 0;
  for (int i = 0; i < node->num_node_points; ++i) {
    struct segment_sketch curr_node_segment_sketch =
        node->node_segment_sketches[i];

    // This is the QoS of this segment. QoS is the estimation quality evaluated
    // as = QoS = segment_length * (max_mean-min_mean) * ((max_mean-min_mean) +
    //     (max_stdev * max_stdev))
    // The smaller the QoS, the more effective the bounds are for similarity
    // estimation

    node_range_value += range_calc(curr_node_segment_sketch,
                                   get_segment_length(node->node_points, i));
  }

  return node_range_value;
}

void exact_incr_knn_search(ts_type *query_ts, ts_type *query_ts_reordered,
                           int *query_order, unsigned int offset,
                           struct dstree_index *index, ts_type minimum_distance,
                           ts_type epsilon, ts_type r_delta, unsigned int k,
                           unsigned int q_id, char *qfilename,
                           unsigned int nprobes) {

  unsigned int curr_size = 0;
  ts_type bsf = FLT_MAX;
  ts_type kth_bsf = FLT_MAX;
  ts_type temp_bsf = FLT_MAX;
  unsigned int cur_probes = 0;

  // the next NN found by incremental search
  unsigned int found_knn = 0;

  // queue containing kNN results
  struct query_result *knn_results = calloc(k, sizeof(struct query_result));
  for (int idx = 0; idx < k; ++idx) {
    knn_results[idx].node = NULL;
    knn_results[idx].distance = FLT_MAX;
  }

  // return k approximate results
  approximate_knn_search(query_ts, query_ts_reordered, query_order, offset, bsf,
                         index, knn_results, k, NULL, NULL, &curr_size, 0);

  ++cur_probes;

  // set the approximate result to be the first item in the queue
  struct query_result approximate_result = knn_results[0];
  // struct query_result bsf_result = approximate_result;

  COUNT_PARTIAL_TIME_END
  index->stats->query_filter_total_time = partial_time;

  index->stats->query_filter_input_time = partial_input_time;
  index->stats->query_filter_output_time = partial_output_time;
  index->stats->query_filter_load_node_time = partial_load_node_time;
  index->stats->query_filter_cpu_time =
      partial_time - partial_input_time - partial_output_time;
  index->stats->query_filter_seq_input_count = partial_seq_input_count;
  index->stats->query_filter_seq_output_count = partial_seq_output_count;
  index->stats->query_filter_rand_input_count = partial_rand_input_count;
  index->stats->query_filter_rand_output_count = partial_rand_output_count;

  index->stats->query_filter_loaded_nodes_count = loaded_nodes_count;
  index->stats->query_filter_loaded_ts_count = loaded_ts_count;
  index->stats->query_filter_checked_nodes_count = checked_nodes_count;
  index->stats->query_filter_checked_ts_count = checked_ts_count;

  index->stats->query_approx_distance = sqrtf(approximate_result.distance);
  index->stats->query_approx_node_filename = approximate_result.node->filename;
  index->stats->query_approx_node_size = approximate_result.node->node_size;
  ;
  index->stats->query_approx_node_level = approximate_result.node->level;

  index->stats->queries_filter_total_time +=
      index->stats->query_filter_total_time;

  index->stats->queries_filter_input_time +=
      index->stats->query_filter_input_time;
  index->stats->queries_filter_output_time +=
      index->stats->query_filter_output_time;
  index->stats->queries_filter_load_node_time +=
      index->stats->query_filter_load_node_time;
  index->stats->queries_filter_cpu_time += index->stats->query_filter_cpu_time;

  index->stats->queries_filter_seq_input_count +=
      index->stats->query_filter_seq_input_count;
  index->stats->queries_filter_seq_output_count +=
      index->stats->query_filter_seq_output_count;
  index->stats->queries_filter_rand_input_count +=
      index->stats->query_filter_rand_input_count;
  index->stats->queries_filter_rand_output_count +=
      index->stats->query_filter_rand_output_count;

  if (approximate_result.node != NULL) {
    index->stats->query_approx_distance = sqrtf(approximate_result.distance);
    index->stats->query_approx_node_filename =
        approximate_result.node->filename;
    index->stats->query_approx_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_approx_node_level = approximate_result.node->level;
  }

  // Early termination...
  if (approximate_result.distance == 0) {
    index->stats->query_exact_distance = sqrtf(approximate_result.distance);
    index->stats->query_exact_node_filename = approximate_result.node->filename;
    index->stats->query_exact_node_size = approximate_result.node->node_size;
    ;
    index->stats->query_exact_node_level = approximate_result.node->level;
    // return approximate_result;

    // IMPORTANT!!!!
    // fix this: increase found_knn and do not print until the end.
    update_query_stats(index, q_id, found_knn, approximate_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
  }

  RESET_QUERY_COUNTERS()
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  struct query_result bsf_result = approximate_result;

  pqueue_t *pq = pqueue_init(index->first_node->node_size, cmp_pri, get_pri,
                             set_pri, get_pos, set_pos);

  // if(approximate_result->node != NULL) {
  //     // Insert approximate result in heap.
  //    pqueue_insert(pq, approximate_result);
  // }

  // struct query_result *do_not_remove = approximate_result;

  // Add the root to the priority queue

  struct query_result *root_pq_item = malloc(sizeof(struct query_result));
  root_pq_item->node = index->first_node;
  root_pq_item->distance =
      calculate_node_min_distance(index, index->first_node, query_ts, 0);

  // initialize the lb distance to be the distance of the query to the root.

  pqueue_insert(pq, root_pq_item);

  struct query_result *n;
  struct query_result temp;

  // start off with 100 bsf steps, increase if necessary
  // cur_probes is strictly less than nprobes because an approximate answer has
  // already been found
  while ((n = pqueue_pop(pq))) {

    // the first element of the queue is not used, thus pos-1
    for (unsigned int pos = found_knn; pos < k; ++pos) {
      bsf_result = knn_results[pos];
      // printf("n->distance = %g, bsf_result.distance = %g\n",
      // sqrt(n->distance), sqrt(bsf_result.distance));
      if (n->distance > bsf_result.distance) // add epsilon+1
      {
        found_knn = pos + 1;
        COUNT_PARTIAL_TIME_END

        update_query_stats(index, q_id, found_knn, bsf_result);
        get_query_stats(index, found_knn);
        print_query_stats(index, q_id, found_knn, qfilename);

        // printf("found NN = %u\n", found_knn);
        // fflush(stdout);
        // reset the bsf for the next NN
        if (found_knn < k) {
          bsf_result = knn_results[found_knn];
        }

        RESET_QUERY_COUNTERS()
        RESET_PARTIAL_COUNTERS()
        COUNT_PARTIAL_TIME_START
      }
    }

    if (found_knn >= k) {
      // printf("found all kNN\n");
      // fflush(stdout);
      break;
    }
    // report the pos-NN neighbors, then continue

    // if (n->distance > kth_bsf/(1 + epsilon))
    //{
    //   break;
    // }

    if (n->node->is_leaf) // n is a leaf
    {
      // upon return, the queue will update the next best (k-foundkNN)th objects
      calculate_node_knn_distance(index, n->node, query_ts_reordered,
                                  query_order, offset, bsf_result.distance, k,
                                  knn_results, NULL, NULL, &curr_size, 0);
      /*
      if (r_delta != FLT_MAX && (knn_results[k-1].distance  <= r_delta * (1 +
      epsilon))) break;
      */
      // increase the number of visited leaves
    }
    // If it is an intermediate node calculate mindist for children
    // and push them in the queue
    else // n is an internal node
    {
      temp = knn_results[k - 1];
      kth_bsf = temp.distance;

      ts_type child_distance;
      child_distance =
          calculate_node_min_distance(index, n->node->left_child, query_ts, 0);

      // mindist_result_left->node->parent = n->node;
      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf) &&
          (n->node->left_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_left =
            malloc(sizeof(struct query_result));
        mindist_result_left->node = n->node->left_child;
        mindist_result_left->distance = child_distance;
        pqueue_insert(pq, mindist_result_left);
      }

      child_distance =
          calculate_node_min_distance(index, n->node->right_child, query_ts, 0);

      // if (child_distance < bsf_result.distance/(1 + epsilon) )
      if ((child_distance < kth_bsf) &&
          (n->node->right_child != approximate_result.node)) // add epsilon
      {
        struct query_result *mindist_result_right =
            malloc(sizeof(struct query_result));
        mindist_result_right->node = n->node->right_child;
        mindist_result_right->distance = child_distance;
        pqueue_insert(pq, mindist_result_right);
      }
    }
    // Free the node currently popped.
    // if(n != do_not_remove)
    free(n);
  }

  // report the elements that were not reported already
  for (unsigned int pos = found_knn; pos < k; ++pos) {
    bsf_result = knn_results[pos];
    found_knn = pos + 1;
    COUNT_PARTIAL_TIME_END
    update_query_stats(index, q_id, found_knn, bsf_result);
    get_query_stats(index, found_knn);
    print_query_stats(index, q_id, found_knn, qfilename);
    // report all results for found_knn - last_found_knn or print their results
    RESET_QUERY_COUNTERS()
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }

  // free the results, eventually do something with them!!

  // Free the nodes that were not popped.
  while ((n = pqueue_pop(pq))) {
    free(n);
  }
  // Free the priority queue.
  pqueue_free(pq);
  free(knn_results);
}

