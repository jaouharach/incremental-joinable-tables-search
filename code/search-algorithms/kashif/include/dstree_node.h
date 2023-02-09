//
//  dstree_node.h
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2012 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef dstreelib_dstree_node_h
#define dstreelib_dstree_node_h


#include "../config.h"
#include "../globals.h"    
#include "dstree_node_split.h"
#include "dstree_file_buffer.h"
#include "pqueue.h"
#include "dstree_query_engine.h"
#include <stdint.h>
#include <pthread.h>

/* start kashif changes */

// vector id
typedef struct vid { 
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos;
  char raw_data_file[300]; // name of the raw (bin) file where vector is store
} vid;


struct job {
  struct vid query_id;
  ts_type * query_vector;
  int * query_order;
  ts_type * query_vector_reordered;
  unsigned int worker_id; // worker who executed this job
} job;

struct pool {
	char cancelled;
	unsigned int remaining; // 0 is thread is not working 1 if thread has finished current job 2 if thread finishied last job (no more jobs to take)
	unsigned int num_threads;

	struct  job * job_array;
  unsigned int num_jobs;
	unsigned int job_counter;

	pthread_t *threads;
  pthread_cond_t *cond_thread_state;
  pthread_mutex_t *cond_mutex;
  
  unsigned int *executed_jobs_count;// nb of executed jobs for each thread
  void *(*function)(void *);

  struct worker_param *params;
  unsigned int num_working_threads;
  unsigned int *working;
} pool;

struct result_vid { 
  unsigned int table_id; // max = 4,294,967,295, must change type if dataset contains more that 4,294,967,295 tables
  uint16_t set_id; // max = 65535, must change type if dataset tables contain more than 65535 columns
  uint16_t pos; // max = 65535, must change type if dataset columns contain more than 65535 cells (vectors)
  uint8_t qpos; // max = 255, must change type if query column contain more than 255 cells (vectors)
  float time; 
  float distance;
  unsigned int num_checked_vectors;
} result_vid;

// result set id
struct result_sid { 
  unsigned int table_id;
  unsigned int set_id;
  char raw_data_file[300]; // name of the raw (json) file where set is store
  unsigned int overlap_size;
};

// result table id for keyword search
struct result_table { 
  unsigned int table_id;
  char raw_data_file[300]; // name of the raw (json) file where set is store
  ts_type min_distance;
  unsigned int num_min;
  unsigned int total_matches;
};

/* end kashif changes */
struct dstree_node {

  struct node_segment_split_policy * node_segment_split_policies;
  short * node_points;
  short * hs_node_points;
  struct segment_sketch * node_segment_sketches;
  struct segment_sketch * hs_node_segment_sketches;
  struct node_split_policy * split_policy; 

  struct dstree_node *left_child;
  struct dstree_node *right_child;
  struct dstree_node *parent;

  struct dstree_file_buffer * file_buffer;

  char * filename;

  mean_stdev_range range;  
 
  int num_node_segment_split_policies;  
  short num_node_points;  //number of vertical split points
  short num_hs_node_points;  //number of horizontal split points

  int max_segment_length; 
  int max_value_length; 
  
  unsigned int node_size;
  unsigned int level;

  unsigned char is_leaf;
  boolean is_left;

  unsigned int gt_pos;
  //unsigned char * gt;
  label_type *gt;
  unsigned int fp_pos;
  unsigned int  *fp;

  /* start kashif changes */
  unsigned int vid_pos;
  struct vid * vid;
  pthread_mutex_t lock;
  /* end kashif changes */
};
struct dstree_node * dstree_root_node_init(struct dstree_index_settings * settings) ;
struct dstree_node * dstree_leaf_node_init(struct dstree_index_settings * settings) ;
//struct dstree_node * dstree_root_node_init();
//struct dstree_node * dstree_leaf_node_init();
enum response node_init_segments(struct dstree_node * node, short * split_points, int segment_size);

enum response append_ts_to_node(struct dstree_index * index, struct dstree_node * node, ts_type * timeseries);  
enum response update_node_statistics(struct dstree_node * node, ts_type * timeseries);

enum response create_dstree_node_filename(struct dstree_index_settings *settings,
                                          struct dstree_node * node, struct dstree_node * parent_node);

//ts_type calculate_ts_in_node_distance (struct dstree_index *index, struct dstree_node *node, ts_type *query);
//ts_type calculate_node_distance (struct dstree_index *index, struct dstree_node *node, ts_type *query);
ts_type calculate_node_min_distance (struct dstree_index *index, struct dstree_node *node, ts_type *query, float warping);
ts_type calculate_ts_in_node_distance (struct dstree_index *index,
				       struct dstree_node *node,
				       ts_type *query_ts_reordered,
				       int * query_order,
				       unsigned int offset,
				       ts_type bound);
ts_type calculate_node_distance (struct dstree_index *index, struct dstree_node *node, ts_type *query_ts_reordered, int *query_order, unsigned int offset, ts_type bsf);
int queue_bounded_sorted_insert(struct  query_result *q, struct query_result d, unsigned int *cur_size, unsigned int k);
void calculate_node_knn_distance (struct dstree_index *index, struct dstree_node *node,
				  ts_type *query_ts_reordered, int *query_order,
				  unsigned int offset, ts_type bsf,
				  unsigned int k,
				  struct query_result  *knn_results,
				  struct bsf_snapshot ** bsf_snapshots,
				  unsigned int *cur_bsf_snapshot,
				  unsigned int * cur_size,
				  float warping);
/* start kashif changes */
int thread_queue_bounded_sorted_insert(struct dstree_index * index, struct query_result *q, struct query_result d,
                                unsigned int *cur_size, unsigned int k, unsigned int thread_id);
int calculate_node_knn_distance_2(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, struct query_result *knn_results,
    // struct bsf_snapshot **bsf_snapshots, unsigned int *cur_bsf_snapshot,
    unsigned int *cur_size, float warping, struct vid * query_id, 
    double * total_query_time, unsigned int * total_checked_vectors,
    unsigned int approx);

int calculate_node_knn_distance_para_incr(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, struct query_result *knn_results,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx);

int calculate_node_knn_distance_para_incr_ostree(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, void *knn_tree,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx, unsigned long *insert_counter);


int calculate_node_knn_distance_para_incr_mmheap(
    struct dstree_index *index, struct dstree_node *node,
    ts_type *query_ts_reordered, int *query_order, unsigned int offset,
    unsigned int k, void * knn_heap,
    unsigned int *cur_size, float warping, struct vid * query_id,
    double * total_query_set_time, unsigned int * total_checked_ts,
    unsigned int thread_id, unsigned int approx, unsigned long *insert_counter);
    
/* end kashif changes */
enum response append_ts_gt_to_child_node(struct dstree_index * index,
					 struct dstree_node * node,
					 ts_type * timeseries,
					 struct dstree_node * parent,
					 int parent_idx,
					 unsigned int fp);
enum response append_ts_gt_to_node(struct dstree_index * index,
				   struct dstree_node * node,
				   ts_type * timeseries,
				   label_type gt,
				   unsigned int fp);

/* start kashif changes */
enum response append_vector_to_node(struct dstree_index * index, struct dstree_node * node, 
            ts_type * vector, unsigned int table_id, unsigned int set_id, 
            unsigned int pos, char * raw_data_file);

enum response append_vector_to_child_node(struct dstree_index *index,
              struct dstree_node *node, ts_type *vector, unsigned int table_id, 
              unsigned int set_id, unsigned int pos, char * raw_data_file);
/* end kashif changes */
/*
  enum response append_ts_gt_to_node(struct dstree_index * index,
				   struct dstree_node * node,
				   ts_type * timeseries,
				   label_type gt);
enum response append_ts_gt_to_child_node(struct dstree_index * index,
					 struct dstree_node * node,
					 ts_type * timeseries,
					 label_type gt);
*/

#endif

