//
//  dstree_index.c
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
#include <stdio.h>
#include <stdlib.h>

#ifdef VALUES
#include <values.h>
#endif

#include "../include/calc_utils.h"
#include "../include/dstree_file_buffer.h"
#include "../include/dstree_file_buffer_manager.h"
#include "../include/dstree_index.h"
#include "../include/dstree_node.h"
#include "../include/dstree_node_split.h"
#include <errno.h>
#include <float.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>

/**
 This function initializes the settings of a dstree index
 */
struct dstree_index_settings *dstree_index_settings_init(
    const char *root_directory, unsigned int timeseries_size,
    unsigned int init_segments, unsigned int max_leaf_size,
    double buffered_memory_size, boolean is_index_new, const char *dataset) {

  if (is_index_new) {
    if (chdir(root_directory) == 0) {
      fprintf(stderr, "WARNING! Target index directory already exists. Please "
                      "delete or choose a new one.\n");
      exit(-1);
    }
  }
  mkdir(root_directory, 0777);

  struct dstree_index_settings *settings =
      malloc(sizeof(struct dstree_index_settings));
  if (settings == NULL) {
    fprintf(stderr, "Error in dstree_index.c: could not allocate memory for "
                    "index settings.\n");
    return NULL;
  }

  settings->root_directory = root_directory;
  settings->timeseries_size = timeseries_size;
  settings->init_segments = init_segments;
  settings->max_leaf_size = max_leaf_size;
  printf("buffered memory = %f MB\n", buffered_memory_size);
  settings->buffered_memory_size = buffered_memory_size;

  settings->dataset = dataset;

  /* Each leaf node has a file called:
     numVerticalSegments_LeftorRight_SplitIndicator_(idxFrom,idxTo,splitValue)_level

     numVerticalSegments cannot exceed 2 characters (max_segment_length)
     LeftorRight is one character : L or R
     SplitIndicator is 0 (mean) or 1(stdev)
     idxFrom: Start of segment that was split
     idxTo :  End of segment that was split
     splitValue: Value used to split segment, (cannot exceed max_value_length
     characters) level :  int number of punctuation marks (underscores:4,
     parentheses:2, commas:2): total 8 null character: 1
  */

  float segment_ends_size = ceil(log10(SHRT_MAX));
  float split_value_size = ceil(log10(INT_MAX) + 1);

  settings->max_filename_size =
      2 + 1 + 1 + 2 * (segment_ends_size) + 10 + split_value_size + 8 + 1;

  return settings;
}

/**
 This function initializes an isax index
 @param isax_index_settings *settings
 @return isax_index
 */
struct dstree_index *dstree_index_init(struct dstree_index_settings *settings) {
  struct dstree_index *index = malloc(sizeof(struct dstree_index));
  if (index == NULL) {
    fprintf(stderr, "Error in dstree_index.c: Could not allocate memory for "
                    "index structure.\n");
    return NULL;
  }

  index->settings = settings;
  index->total_records = 0;

  dstree_init_stats(index);

  if (!init_file_buffer_manager(index)) {
    fprintf(stderr, "Error in dstree_index.c:  Could not initialize the \
                       file buffer manager for this index.\n");
    return NULL;
  }

  return index;
}

/* start kashif changes */
enum response dstree_init_thread_stats(struct dstree_index *index, int num_threads) 
{
  index->stats->thread_query_total_input_time = calloc(num_threads, sizeof(double));
  index->stats->thread_query_total_output_time = calloc(num_threads, sizeof(double));
  index->stats->thread_query_total_load_node_time = calloc(num_threads, sizeof(double));
  index->stats->thread_query_total_cpu_time = calloc(num_threads, sizeof(double));
  index->stats->thread_query_total_time = calloc(num_threads, sizeof(double));

  index->stats->thread_query_total_loaded_nodes_count = calloc(num_threads, sizeof(unsigned int));
  index->stats->thread_query_total_checked_nodes_count = calloc(num_threads, sizeof(unsigned int));
  index->stats->thread_query_total_loaded_ts_count = calloc(num_threads, sizeof(unsigned long long));
  index->stats->thread_query_total_checked_ts_count = calloc(num_threads, sizeof(unsigned long long));

  return SUCCESS;
}
/* end kashif changes */

enum response dstree_init_stats(struct dstree_index *index) {
  index->stats = malloc(sizeof(struct stats_info));
  if (index->stats == NULL) {
    fprintf(stderr, "Error in dstree_index.c: Could not allocate memory for "
                    "stats structure.\n");
    return FAILURE;
  }

  /*INDEX STATISTICS*/
  index->stats->idx_traverse_tree_input_time = 0;
  index->stats->idx_traverse_tree_output_time = 0;
  index->stats->idx_traverse_tree_cpu_time = 0;
  index->stats->idx_traverse_tree_total_time = 0;

  index->stats->idx_traverse_tree_seq_input_count = 0;
  index->stats->idx_traverse_tree_seq_output_count = 0;
  index->stats->idx_traverse_tree_rand_input_count = 0;
  index->stats->idx_traverse_tree_rand_output_count = 0;

  index->stats->idx_append_ts_to_leaf_input_time = 0;
  index->stats->idx_append_ts_to_leaf_output_time = 0;
  index->stats->idx_append_ts_to_leaf_cpu_time = 0;
  index->stats->idx_append_ts_to_leaf_total_time = 0;

  index->stats->idx_append_ts_to_leaf_seq_input_count = 0;
  index->stats->idx_append_ts_to_leaf_seq_output_count = 0;
  index->stats->idx_append_ts_to_leaf_rand_input_count = 0;
  index->stats->idx_append_ts_to_leaf_rand_output_count = 0;

  index->stats->idx_evaluate_split_policies_input_time = 0;
  index->stats->idx_evaluate_split_policies_output_time = 0;
  index->stats->idx_evaluate_split_policies_cpu_time = 0;
  index->stats->idx_evaluate_split_policies_total_time = 0;

  index->stats->idx_evaluate_split_policies_seq_input_count = 0;
  index->stats->idx_evaluate_split_policies_seq_output_count = 0;
  index->stats->idx_evaluate_split_policies_rand_input_count = 0;
  index->stats->idx_evaluate_split_policies_rand_output_count = 0;

  index->stats->idx_split_node_input_time = 0;
  index->stats->idx_split_node_output_time = 0;
  index->stats->idx_split_node_cpu_time = 0;
  index->stats->idx_split_node_total_time = 0;

  index->stats->idx_split_node_seq_input_count = 0;
  index->stats->idx_split_node_seq_output_count = 0;
  index->stats->idx_split_node_rand_input_count = 0;
  index->stats->idx_split_node_rand_output_count = 0;

  index->stats->idx_building_input_time = 0;
  index->stats->idx_building_output_time = 0;
  index->stats->idx_building_cpu_time = 0;
  index->stats->idx_building_total_time = 0;

  index->stats->idx_building_seq_input_count = 0;
  index->stats->idx_building_seq_output_count = 0;
  index->stats->idx_building_rand_input_count = 0;
  index->stats->idx_building_rand_output_count = 0;

  index->stats->idx_writing_input_time = 0;
  index->stats->idx_writing_output_time = 0;
  index->stats->idx_writing_cpu_time = 0;
  index->stats->idx_writing_total_time = 0;

  index->stats->idx_writing_seq_input_count = 0;
  index->stats->idx_writing_seq_output_count = 0;
  index->stats->idx_writing_rand_input_count = 0;
  index->stats->idx_writing_rand_output_count = 0;

  index->stats->idx_reading_input_time = 0;
  index->stats->idx_reading_output_time = 0;
  index->stats->idx_reading_cpu_time = 0;
  index->stats->idx_reading_total_time = 0;

  index->stats->idx_reading_seq_input_count = 0;
  index->stats->idx_reading_seq_output_count = 0;
  index->stats->idx_reading_rand_input_count = 0;
  index->stats->idx_reading_rand_output_count = 0;

  index->stats->idx_total_input_time = 0;
  index->stats->idx_total_output_time = 0;
  index->stats->idx_total_cpu_time = 0;
  index->stats->idx_total_time = 0;

  index->stats->idx_total_seq_input_count = 0;
  index->stats->idx_total_seq_output_count = 0;
  index->stats->idx_total_rand_input_count = 0;
  index->stats->idx_total_rand_output_count = 0;

  index->stats->total_nodes_count = 0;
  index->stats->leaf_nodes_count = 0;
  index->stats->empty_leaf_nodes_count = 0;

  index->stats->idx_size_bytes = 0;
  index->stats->idx_size_blocks = 0;

  index->stats->min_fill_factor = FLT_MAX;
  index->stats->max_fill_factor = 0;
  index->stats->sum_fill_factor = 0;
  index->stats->sum_squares_fill_factor = 0;
  index->stats->avg_fill_factor = 0;
  index->stats->sd_fill_factor = 0;

  index->stats->min_height = FLT_MAX;
  index->stats->max_height = 0;
  index->stats->sum_height = 0;
  index->stats->sum_squares_height = 0;
  index->stats->avg_height = 0;
  index->stats->sd_height = 0;

  /*PER QUERY STATISTICS*/
  index->stats->query_filter_input_time = 0;
  index->stats->query_filter_output_time = 0;
  index->stats->query_filter_load_node_time = 0;
  index->stats->query_filter_cpu_time = 0;
  index->stats->query_filter_total_time = 0;

  index->stats->query_filter_seq_input_count = 0;
  index->stats->query_filter_seq_output_count = 0;
  index->stats->query_filter_rand_input_count = 0;
  index->stats->query_filter_rand_output_count = 0;

  index->stats->query_filter_loaded_nodes_count = 0;
  index->stats->query_filter_checked_nodes_count = 0;
  index->stats->query_filter_loaded_ts_count = 0;
  index->stats->query_filter_checked_ts_count = 0;

  index->stats->query_refine_input_time = 0;
  index->stats->query_refine_output_time = 0;
  index->stats->query_refine_load_node_time = 0;
  index->stats->query_refine_cpu_time = 0;
  index->stats->query_refine_total_time = 0;

  index->stats->query_refine_seq_input_count = 0;
  index->stats->query_refine_seq_output_count = 0;
  index->stats->query_refine_rand_input_count = 0;
  index->stats->query_refine_rand_output_count = 0;

  index->stats->query_refine_loaded_nodes_count = 0;
  index->stats->query_refine_checked_nodes_count = 0;
  index->stats->query_refine_loaded_ts_count = 0;
  index->stats->query_refine_checked_ts_count = 0;

  index->stats->query_total_input_time = 0;
  index->stats->query_total_output_time = 0;
  index->stats->query_total_load_node_time = 0;
  index->stats->query_total_cpu_time = 0;
  index->stats->query_total_time = 0;

  index->stats->query_total_seq_input_count = 0;
  index->stats->query_total_seq_output_count = 0;
  index->stats->query_total_rand_input_count = 0;
  index->stats->query_total_rand_output_count = 0;

  index->stats->query_total_loaded_nodes_count = 0;
  index->stats->query_total_checked_nodes_count = 0;
  index->stats->query_total_loaded_ts_count = 0;
  index->stats->query_total_checked_ts_count = 0;

  index->stats->query_exact_distance = 0;
  index->stats->query_exact_node_filename = NULL;
  index->stats->query_exact_node_size = 0;
  index->stats->query_exact_node_level = 0;
  index->stats->query_label = 0;

  index->stats->query_approx_distance = 0;
  index->stats->query_approx_node_filename = NULL;
  index->stats->query_approx_node_size = 0;
  index->stats->query_approx_node_level = 0;

  index->stats->query_lb_distance = 0;

  index->stats->query_tlb = 0;
  index->stats->query_eff_epsilon = 0;
  index->stats->query_pruning_ratio = 0;

  /*SUMMARY STATISTICS FOR ALL QUERIES*/
  index->stats->queries_refine_input_time = 0;
  index->stats->queries_refine_output_time = 0;
  index->stats->queries_refine_load_node_time = 0;
  index->stats->queries_refine_cpu_time = 0;
  index->stats->queries_refine_total_time = 0;

  index->stats->queries_refine_seq_input_count = 0;
  index->stats->queries_refine_seq_output_count = 0;
  index->stats->queries_refine_rand_input_count = 0;
  index->stats->queries_refine_rand_output_count = 0;

  index->stats->queries_filter_input_time = 0;
  index->stats->queries_filter_output_time = 0;
  index->stats->queries_filter_load_node_time = 0;
  index->stats->queries_filter_cpu_time = 0;
  index->stats->queries_filter_total_time = 0;

  index->stats->queries_filter_seq_input_count = 0;
  index->stats->queries_filter_seq_output_count = 0;
  index->stats->queries_filter_rand_input_count = 0;
  index->stats->queries_filter_rand_output_count = 0;

  index->stats->queries_total_input_time = 0;
  index->stats->queries_total_output_time = 0;
  index->stats->queries_total_load_node_time = 0;
  index->stats->queries_total_cpu_time = 0;
  index->stats->queries_total_time = 0;

  index->stats->queries_total_seq_input_count = 0;
  index->stats->queries_total_seq_output_count = 0;
  index->stats->queries_total_rand_input_count = 0;
  index->stats->queries_total_rand_output_count = 0;

  index->stats->queries_min_eff_epsilon = FLT_MAX;
  index->stats->queries_max_eff_epsilon = 0;
  index->stats->queries_sum_eff_epsilon = 0;
  index->stats->queries_sum_squares_eff_epsilon = 0;
  index->stats->queries_avg_eff_epsilon = 0;
  index->stats->queries_sd_eff_epsilon = 0;

  index->stats->queries_min_pruning_ratio = FLT_MAX;
  index->stats->queries_max_pruning_ratio = 0;
  index->stats->queries_sum_pruning_ratio = 0;
  index->stats->queries_sum_squares_pruning_ratio = 0;
  index->stats->queries_avg_pruning_ratio = 0;
  index->stats->queries_sd_pruning_ratio = 0;

  index->stats->queries_min_tlb = FLT_MAX;
  index->stats->queries_max_tlb = 0;
  index->stats->queries_sum_tlb = 0;
  index->stats->queries_sum_squares_tlb = 0;
  index->stats->queries_avg_tlb = 0;
  index->stats->queries_sd_tlb = 0;

  index->stats->tlb_ts_count = 0;
  index->stats->eff_epsilon_queries_count = 0;

  //    index->stats->total_queries_count = 0;

  /*COMBINED STATISTICS FOR INDEXING AND QUERY WORKLOAD*/
  index->stats->total_input_time = 0;
  index->stats->total_output_time = 0;
  index->stats->total_load_node_time = 0;
  index->stats->total_cpu_time = 0;
  index->stats->total_time = 0;
  index->stats->total_time_sanity_check = 0;

  index->stats->total_seq_input_count = 0;
  index->stats->total_seq_output_count = 0;
  index->stats->total_rand_input_count = 0;
  index->stats->total_rand_output_count = 0;

  index->stats->total_parse_time = 0;
  index->stats->total_ts_count = 0;

  return SUCCESS;
}

void dstree_get_index_stats(struct dstree_index *index) {
  index->stats->total_seq_input_count =
      index->stats->idx_building_seq_input_count +
      index->stats->idx_writing_seq_input_count +
      index->stats->idx_reading_seq_input_count;
  //    + index->stats->queries_total_seq_input_count;
  index->stats->total_seq_output_count =
      index->stats->idx_building_seq_output_count +
      index->stats->idx_writing_seq_output_count +
      index->stats->idx_reading_seq_output_count;
  //+ index->stats->queries_total_seq_output_count;
  index->stats->total_rand_input_count =
      index->stats->idx_building_rand_input_count +
      index->stats->idx_writing_rand_input_count +
      index->stats->idx_reading_rand_input_count;
  // + index->stats->queries_total_rand_input_count;
  index->stats->total_rand_output_count =
      index->stats->idx_building_rand_output_count +
      index->stats->idx_writing_rand_output_count +
      index->stats->idx_reading_rand_output_count;
  //+ index->stats->queries_total_rand_output_count;

  index->stats->total_input_time = index->stats->idx_building_input_time +
                                   index->stats->idx_writing_input_time +
                                   index->stats->idx_reading_input_time;
  //   + index->stats->queries_total_input_time;
  index->stats->total_output_time = index->stats->idx_building_output_time +
                                    index->stats->idx_writing_output_time +
                                    index->stats->idx_reading_output_time;
  //    + index->stats->queries_total_output_time;
  index->stats->total_cpu_time = index->stats->idx_building_cpu_time +
                                 index->stats->idx_writing_cpu_time +
                                 index->stats->idx_reading_cpu_time;
  //    + index->stats->queries_total_cpu_time;

  index->stats->total_time = index->stats->total_input_time +
                             index->stats->total_output_time +
                             index->stats->total_cpu_time;

  // index->stats->total_time_sanity_check = total_time;

  // index->stats->load_node_time = load_node_time;
  index->stats->total_parse_time = total_parse_time;

  // index->stats->loaded_nodes_count = loaded_nodes_count;
  index->stats->leaf_nodes_count = leaf_nodes_count;
  index->stats->empty_leaf_nodes_count = empty_leaf_nodes_count;

  // index->stats->checked_nodes_count = checked_nodes_count;
  index->stats->total_nodes_count = total_nodes_count;
  index->stats->total_ts_count = total_ts_count;

  dstree_get_index_footprint(index);
}

void dstree_get_index_footprint(struct dstree_index *index) {

  const char *filename =
      malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
  filename = strcpy(filename, index->settings->root_directory);
  filename = strcat(filename, "root.idx\0");

  struct stat st;
  unsigned int count_leaves;

  if (stat(filename, &st) == 0) {
    index->stats->idx_size_bytes = (long long)st.st_size;
    index->stats->idx_size_blocks = (long long)st.st_blksize;
  }

  count_leaves = index->stats->leaf_nodes_count;

  index->stats->avg_fill_factor =
      ((double)index->stats->sum_fill_factor) / count_leaves;
  index->stats->sum_squares_fill_factor -=
      (pow(index->stats->sum_fill_factor, 2) / count_leaves);
  index->stats->sd_fill_factor =
      sqrt(((double)index->stats->sum_squares_fill_factor) / count_leaves);

  index->stats->avg_height = ((double)index->stats->sum_height) / count_leaves;
  index->stats->sum_squares_height -=
      (pow(index->stats->sum_height, 2) / count_leaves);
  index->stats->sd_height =
      sqrt(((double)index->stats->sum_squares_height) / count_leaves);

  free(filename);
}

void dstree_print_index_stats(struct dstree_index *index, char *dataset) {
  /*
  printf("------------------------ \n");
  printf("INDEX SUMMARY STATISTICS \n");
  printf("------------------------ \n");
  */
  //  id = -1 for index and id = query_id for queries
  int id = -1;
  printf("Index_traverse_tree_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_traverse_tree_input_time / 1000000, dataset, id, id);
  printf("Index_traverse_tree_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_traverse_tree_output_time / 1000000, dataset, id,
         id);
  printf("Index_traverse_tree_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_traverse_tree_cpu_time / 1000000, dataset, id, id);
  printf("Index_traverse_tree_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_traverse_tree_total_time / 1000000, dataset, id, id);

  printf("Index_traverse_tree_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_traverse_tree_seq_input_count, dataset, id, id);

  printf("Index_traverse_tree_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_traverse_tree_seq_output_count, dataset, id, id);

  printf("Index_traverse_tree_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_traverse_tree_rand_input_count, dataset, id, id);

  printf("Index_traverse_tree_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_traverse_tree_rand_output_count, dataset, id, id);

  printf("Index_append_ts_to_leaf_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_input_time / 1000000, dataset, id,
         id);
  printf("Index_append_ts_to_leaf_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_output_time / 1000000, dataset, id,
         id);
  printf("Index_append_ts_to_leaf_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_cpu_time / 1000000, dataset, id,
         id);
  printf("Index_append_ts_to_leaf_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_total_time / 1000000, dataset, id,
         id);

  printf("Index_append_ts_to_leaf_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_seq_input_count, dataset, id, id);

  printf("Index_append_ts_to_leaf_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_seq_output_count, dataset, id, id);

  printf("Index_append_ts_to_leaf_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_rand_input_count, dataset, id, id);

  printf("Index_append_ts_to_leaf_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_append_ts_to_leaf_rand_output_count, dataset, id,
         id);

  printf("Index_evaluate_split_policies_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_input_time / 1000000,
         dataset, id, id);
  printf("Index_evaluate_split_policies_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_output_time / 1000000,
         dataset, id, id);
  printf("Index_evaluate_split_policies_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_cpu_time / 1000000, dataset,
         id, id);
  printf("Index_evaluate_split_policies_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_total_time / 1000000,
         dataset, id, id);

  printf("Index_evaluate_split_policies_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_seq_input_count, dataset, id,
         id);

  printf("Index_evaluate_split_policies_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_seq_output_count, dataset,
         id, id);

  printf("Index_evaluate_split_policies_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_rand_input_count, dataset,
         id, id);

  printf("Index_evaluate_split_policies_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_evaluate_split_policies_rand_output_count, dataset,
         id, id);

  printf("Index_split_node_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_split_node_input_time / 1000000, dataset, id, id);
  printf("Index_split_node_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_split_node_output_time / 1000000, dataset, id, id);
  printf("Index_split_node_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_split_node_cpu_time / 1000000, dataset, id, id);
  printf("Index_split_node_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_split_node_total_time / 1000000, dataset, id, id);

  printf("Index_split_node_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_split_node_seq_input_count, dataset, id, id);

  printf("Index_split_node_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_split_node_seq_output_count, dataset, id, id);

  printf("Index_split_node_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_split_node_rand_input_count, dataset, id, id);

  printf("Index_split_node_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_split_node_rand_output_count, dataset, id, id);

  printf("Index_building_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_building_input_time / 1000000, dataset, id, id);
  printf("Index_building_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_building_output_time / 1000000, dataset, id, id);
  printf("Index_building_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_building_cpu_time / 1000000, dataset, id, id);
  printf("Index_building_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_building_total_time / 1000000, dataset, id, id);

  printf("Index_building_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_building_seq_input_count, dataset, id, id);

  printf("Index_building_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_building_seq_output_count, dataset, id, id);

  printf("Index_building_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_building_rand_input_count, dataset, id, id);

  printf("Index_building_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_building_rand_output_count, dataset, id, id);

  printf("Index_writing_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_writing_input_time / 1000000, dataset, id, id);

  printf("Index_writing_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_writing_output_time / 1000000, dataset, id, id);

  printf("Index_writing_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_writing_cpu_time / 1000000, dataset, id, id);

  printf("Index_writing_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_writing_total_time / 1000000, dataset, id, id);

  printf("Index_writing_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_writing_seq_input_count, dataset, id, id);

  printf("Index_writing_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_writing_seq_output_count, dataset, id, id);

  printf("Index_writing_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_writing_rand_input_count, dataset, id, id);

  printf("Index_writing_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_writing_rand_output_count, dataset, id, id);

  printf("Index_reading_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_reading_input_time / 1000000, dataset, id, id);
  printf("Index_reading_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_reading_output_time / 1000000, dataset, id, id);
  printf("Index_reading_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_reading_cpu_time / 1000000, dataset, id, id);
  printf("Index_reading_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->idx_reading_total_time / 1000000, dataset, id, id);

  printf("Index_reading_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_reading_seq_input_count, dataset, id, id);

  printf("Index_reading_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_reading_seq_output_count, dataset, id, id);

  printf("Index_reading_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->idx_reading_rand_input_count, dataset, id, id);

  printf("Index_reading_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->idx_reading_rand_output_count, dataset, id, id);

  printf("Index_total_input_time_secs\t%lf\t%s\t%d\n",
         index->stats->total_input_time / 1000000, dataset, id, id);
  printf("Index_total_output_time_secs\t%lf\t%s\t%d\n",
         index->stats->total_output_time / 1000000, dataset, id, id);
  printf("Index_total_cpu_time_secs\t%lf\t%s\t%d\n",
         index->stats->total_cpu_time / 1000000, dataset, id, id);
  printf("Index_total_time_secs\t%lf\t%s\t%d\n",
         index->stats->total_time / 1000000, dataset, id, id);

  printf("Index_total_seq_input_count\t%llu\t%s\t%d\n",
         index->stats->total_seq_input_count, dataset, id, id);

  printf("Index_total_seq_output_count\t%llu\t%s\t%d\n",
         index->stats->total_seq_output_count, dataset, id, id);

  printf("Index_total_rand_input_count\t%llu\t%s\t%d\n",
         index->stats->total_rand_input_count, dataset, id, id);

  printf("Index_total_rand_output_count\t%llu\t%s\t%d\n",
         index->stats->total_rand_output_count, dataset, id, id);

  printf("Internal_nodes_count\t%lu\t%s\t%d\n",
         (index->stats->total_nodes_count - index->stats->leaf_nodes_count),
         dataset, id, id);

  printf("Leaf_nodes_count\t%lu\t%s\t%d\n", index->stats->leaf_nodes_count,
         dataset, id, id);

  printf("Empty_leaf_nodes_count\t%lu\t%s\t%d\n",
         index->stats->empty_leaf_nodes_count, dataset, id, id);

  printf("Total_nodes_count\t%lu\t%s\t%d\n", index->stats->total_nodes_count,
         dataset, id, id);

  double size_MB = (index->stats->idx_size_bytes) * 1.0 / (1024 * 1024);

  printf("Index_size_MB\t%lf\t%s\t%d\n", size_MB, dataset, id, id);

  printf("Minimum_fill_factor\t%f\t%s\t%d\n", index->stats->min_fill_factor,
         dataset, id, id);

  printf("Maximum_fill_factor\t%f\t%s\t%d\n", index->stats->max_fill_factor,
         dataset, id, id);

  printf("Average_fill_factor\t%f\t%s\t%d\n", index->stats->avg_fill_factor,
         dataset, id, id);

  printf("SD_height\t%f\t%s\t%d\n", index->stats->sd_height, dataset, id, id);

  printf("Minimum_height\t%u\t%s\t%d\n", index->stats->min_height, dataset, id,
         id);

  printf("Maximum_height\t%u\t%s\t%d\n", index->stats->max_height, dataset, id,
         id);

  printf("Average_height\t%f\t%s\t%d\n", index->stats->avg_height, dataset, id,
         id);

  printf("SD_height\t%f\t%s\t%d\n", index->stats->sd_height, dataset, id, id);

  printf("Total_ts_count\t%u\t%s\t%d\n", index->stats->total_ts_count, dataset,
         id, id);

  // for (int i = 0; i < index->stats->leaves_counter; ++i) {
  //   double fill_factor = ((double)index->stats->leaves_sizes[i]) /
  //                        index->settings->max_leaf_size;
  //   printf("Leaf_report_node_%d \t Height  %d  \t%s\t%d\n", (i + 1),
  //          index->stats->leaves_heights[i], dataset, id, id);
  //   printf("Leaf_report_node_%d \t Fill_Factor  %f \t%s\t%d\n", (i + 1),
  //          fill_factor, dataset, id, id);
  // }
}

void print_tlb_stats(struct dstree_index *index, unsigned int query_num,
                     char *queries) {

  printf("Query_avg_tlb\t%lf\t%s\t%u\n", total_tlb / total_ts_count, queries,
         query_num);
  printf("Leaf_nodes_count\t%u\t%s\t%d\n", leaf_nodes_count, queries,
         query_num);
  printf("Total_ts_count\t%u\t%s\t%d\n", total_ts_count, queries, query_num);
}

void dstree_index_destroy(struct dstree_index *index, struct dstree_node *node,
                          boolean is_index_new) {
  
  if (node->level == 0) // root
  {
    if (node->node_segment_split_policies != NULL)
      free(node->node_segment_split_policies);
    if (index->buffer_manager != NULL)
      destroy_buffer_manager(index);

    if (index->settings->classify)
    {
      if (index->gt_filename != NULL)
        free(index->gt_filename);
      if (index->gt_cache != NULL)
        free(index->gt_cache);
    }
    if(index->settings->track_file_pos)
    {
      if (index->fp_filename != NULL)
        free(index->fp_filename);
      if (index->fp_cache != NULL)
        free(index->fp_cache);
    }
    /* start kashif changes */
    if(index->settings->track_vector)
    {
      if (index->vid_filename != NULL)
        free(index->vid_filename);
      if (index->vid_cache != NULL)
        free(index->vid_cache);
    }
    if(index->settings->parallel)
    {
      if (index->stats->thread_query_total_input_time != NULL)
      {
        free(index->stats->thread_query_total_input_time);
        index->stats->thread_query_total_input_time = NULL;
      }
      if (index->stats->thread_query_total_output_time != NULL)
      {
        free(index->stats->thread_query_total_output_time);
        index->stats->thread_query_total_output_time = NULL;
      }
      if (index->stats->thread_query_total_load_node_time != NULL)
      {
        free(index->stats->thread_query_total_load_node_time);
        index->stats->thread_query_total_load_node_time = NULL;
      }
      if (index->stats->thread_query_total_cpu_time != NULL)
      {
        free(index->stats->thread_query_total_cpu_time);
        index->stats->thread_query_total_cpu_time = NULL;
      }
      if (index->stats->thread_query_total_time != NULL)
      {
        free(index->stats->thread_query_total_time);
        index->stats->thread_query_total_time = NULL;
      }

      if (index->stats->thread_query_total_loaded_nodes_count != NULL)
      {
        free(index->stats->thread_query_total_loaded_nodes_count);
        index->stats->thread_query_total_loaded_nodes_count = NULL;
      }
      if (index->stats->thread_query_total_checked_nodes_count != NULL)
      {
        free(index->stats->thread_query_total_checked_nodes_count);
        index->stats->thread_query_total_checked_nodes_count = NULL;
      }
      if (index->stats->thread_query_total_loaded_ts_count != NULL)
      {
        free(index->stats->thread_query_total_loaded_ts_count);
        index->stats->thread_query_total_loaded_ts_count = NULL;
      }
      if (index->stats->thread_query_total_checked_ts_count != NULL)
      {
        free(index->stats->thread_query_total_checked_ts_count);
        index->stats->thread_query_total_checked_ts_count = NULL;
      }
    }
    
    /* end kashif changes */
  }
  if (!node->is_leaf) {
    dstree_index_destroy(index, node->right_child, is_index_new);
    dstree_index_destroy(index, node->left_child, is_index_new);
  }

  if (node->split_policy != NULL) {
    free(node->split_policy);
  }

  if (node->filename != NULL) {
    free(node->filename);
  }

  if (node->gt != NULL) {
    free(node->gt);
  }

  if (node->fp != NULL) {
    free(node->fp);
  }

  /* start kashif changes */
  if (node->vid != NULL)
  {
    free(node->vid);
    node->vid = NULL;
  }
  if(node->is_leaf)
  {
    pthread_mutex_destroy(&node->lock);
  }
  /* end kashif changes */

  if (node->file_buffer != NULL) {
    if(index->settings->parallel)
    {
      // clearing the data for this node (alocated by some thread)
      for (int i = 0; i < node->file_buffer->buffered_list_size; ++i) {
        free(node->file_buffer->buffered_list[i]);
      }
    }
    free(node->file_buffer->buffered_list);
    node->file_buffer->buffered_list = NULL;
    node->file_buffer->buffered_list_size = 0;
    free(node->file_buffer);
  }

  if (node->node_segment_sketches != NULL) {
    // free indicators first
    for (int i = 0; i < node->num_node_points; ++i) {
      free(node->node_segment_sketches[i].indicators);
    }

    free(node->node_segment_sketches);
  }

  if (node->node_points != NULL) {
    free(node->node_points);
  }

  if (is_index_new) {
    if (node->hs_node_segment_sketches != NULL) {
      // free indicators first
      for (int i = 0; i < node->num_hs_node_points; ++i) {
        free(node->hs_node_segment_sketches[i].indicators);
      }

      free(node->hs_node_segment_sketches);
    }

    if (node->hs_node_points != NULL) {
      free(node->hs_node_points);
    }
  }
  free(node);
}

void destroy_buffer_manager(struct dstree_index *index) {

  if (index->buffer_manager != NULL) {
    struct dstree_file_map *currP;
    struct dstree_file_map *temp;

    temp = NULL;
    currP = index->buffer_manager->file_map;

    while (currP != NULL) {
      temp = currP;
      currP = currP->next;
      free(temp);
    }

    free(index->buffer_manager->mem_array);

    free(index->buffer_manager);
  }
}

enum response dstree_index_insert(struct dstree_index *index,
                                  ts_type *timeseries) 
{

  // traverse the index tree to find the appropriate node
  struct dstree_node *node = index->first_node;

  while (!node->is_leaf) {
    if (!update_node_statistics(node, timeseries)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (node_split_policy_route_to_left(node, timeseries))
      node = node->left_child;
    else
      node = node->right_child;
  }

  COUNT_PARTIAL_TIME_END
  index->stats->idx_traverse_tree_total_time += partial_time;
  index->stats->idx_traverse_tree_input_time += partial_input_time;
  index->stats->idx_traverse_tree_output_time += partial_output_time;
  index->stats->idx_traverse_tree_cpu_time +=
      partial_time - partial_input_time - partial_output_time;
  index->stats->idx_traverse_tree_seq_input_count += partial_seq_input_count;
  index->stats->idx_traverse_tree_seq_output_count += partial_seq_output_count;
  index->stats->idx_traverse_tree_rand_input_count += partial_rand_input_count;
  index->stats->idx_traverse_tree_rand_output_count +=
      partial_rand_output_count;

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  if (node->is_leaf) {

    if (!update_node_statistics(node, timeseries)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (!append_ts_to_node(index, node, timeseries)) {
      fprintf(stderr, "Error in dstree_index.c: could not append \
                        time series to node %s\n",
              node->filename);
      return FAILURE;
    }

    COUNT_PARTIAL_TIME_END
    index->stats->idx_append_ts_to_leaf_total_time += partial_time;
    index->stats->idx_append_ts_to_leaf_input_time += partial_input_time;
    index->stats->idx_append_ts_to_leaf_output_time += partial_output_time;
    index->stats->idx_append_ts_to_leaf_cpu_time +=
        partial_time - partial_input_time - partial_output_time;
    index->stats->idx_append_ts_to_leaf_seq_input_count +=
        partial_seq_input_count;
    index->stats->idx_append_ts_to_leaf_seq_output_count +=
        partial_seq_output_count;
    index->stats->idx_append_ts_to_leaf_rand_input_count +=
        partial_rand_input_count;
    index->stats->idx_append_ts_to_leaf_rand_output_count +=
        partial_rand_output_count;

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    // if split needed, split the node and refresh curr_node
    if (node->node_size >= index->settings->max_leaf_size) {

      struct node_split_policy curr_node_split_policy;
      ts_type max_diff_value = (FLT_MAX * (-1));
      ts_type avg_children_range_value = 0;
      short hs_split_point = -1;
      short *child_node_points;
      int num_child_node_points = 0;
      const int num_child_segments = 2; // by default split to two subsegments

      // we want to test every possible split policy for each segment

      // for each segment
      for (int i = 0; i < node->num_node_points; ++i) {
        struct segment_sketch curr_node_segment_sketch =
            node->node_segment_sketches[i];

        // This is the QoS of this segment. QoS is the estimation quality
        // evaluated as = QoS = segment_length * (max_mean_min_mean) *
        // ((max_mean_min_mean) +
        //     (max_stdev * max_stdev))
        // The smaller the QoS, the more effective the bounds are for similarity
        // estimation

        ts_type node_range_value = range_calc(
            curr_node_segment_sketch, get_segment_length(node->node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_node_segment_split_policy =
              node->node_segment_split_policies[j];
          // to hold the two child segments
          struct segment_sketch *child_node_segment_sketches = NULL;

          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);

          if (child_node_segment_sketches == NULL) {
            fprintf(stderr, "Error in dstree_index.c: could not allocate \
                            memory for the child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr, "Error in dstree_index.c: could not allocate\
                               memory for the child node segment sketches \
                               indicators for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else {
            fprintf(stderr, "Error in dstree_index.c: Split policy was not \
                            set properly for node %s\n",
                    node->filename);
            return FAILURE;
          }

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->node_points, i));
          }

          // diff_value represents the splitting benefit
          // B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
          // the higher the diff_value, the better is the splitting

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);
          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_node_segment_split_policy;
          }
          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }
          free(child_node_segment_sketches);
        }
      }

      // add trade-off for horizontal split
      max_diff_value = max_diff_value * 2;

      // we want to test every possible split policy for each horizontal segment
      for (int i = 0; i < node->num_hs_node_points; ++i) {
        struct segment_sketch curr_hs_node_segment_sketch =
            node->hs_node_segment_sketches[i];
        ts_type node_range_value =
            range_calc(curr_hs_node_segment_sketch,
                       get_segment_length(node->hs_node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_hs_node_segment_split_policy =
              node->node_segment_split_policies[j];

          struct segment_sketch *child_node_segment_sketches =
              NULL; // to hold the two child segments
          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);
          if (child_node_segment_sketches == NULL) {
            fprintf(stderr,
                    "Error in dstree_index.c: could not allocate memory \
                            for the horizontal child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_hs_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr,
                      "Error in dstree_index.c: could not allocate memory \
                             for the horizontal child node segment sketches indicators \
                             for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else
            printf("split policy not initialized properly\n");

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->hs_node_points, i));
          }

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);

          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_hs_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_hs_node_segment_split_policy;
            hs_split_point = get_hs_split_point(
                node->node_points, curr_node_split_policy.split_from,
                curr_node_split_policy.split_to, node->num_node_points);
          }

          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }

          free(child_node_segment_sketches);
        }
      }

      node->split_policy = NULL;
      node->split_policy = malloc(sizeof(struct node_split_policy));
      if (node->split_policy == NULL) {
        fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                        for the split policy of node  %s\n",
                node->filename);
        return FAILURE;
      }
      node->split_policy->split_from = curr_node_split_policy.split_from;
      node->split_policy->split_to = curr_node_split_policy.split_to;
      node->split_policy->indicator_split_idx =
          curr_node_split_policy.indicator_split_idx;
      node->split_policy->indicator_split_value =
          curr_node_split_policy.indicator_split_value;
      node->split_policy->curr_node_segment_split_policy =
          curr_node_split_policy.curr_node_segment_split_policy;

      COUNT_PARTIAL_TIME_END
      index->stats->idx_evaluate_split_policies_total_time += partial_time;
      index->stats->idx_evaluate_split_policies_input_time +=
          partial_input_time;
      index->stats->idx_evaluate_split_policies_output_time +=
          partial_output_time;
      index->stats->idx_evaluate_split_policies_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_evaluate_split_policies_seq_input_count +=
          partial_seq_input_count;
      index->stats->idx_evaluate_split_policies_seq_output_count +=
          partial_seq_output_count;
      index->stats->idx_evaluate_split_policies_rand_input_count +=
          partial_rand_input_count;
      index->stats->idx_evaluate_split_policies_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

      // when hs_split_point stays less than 0, it means that
      // considering splitting a vertical segment is not worth it
      // according to the QoS heuristic

      if (hs_split_point < 0) {
        num_child_node_points = node->num_node_points;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have the same number of segments as parent
        for (int i = 0; i < num_child_node_points; ++i) {
          child_node_points[i] = node->node_points[i];
        }
      } else {
        num_child_node_points = node->num_node_points + 1;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have one additional segment than the parent
        for (int i = 0; i < (num_child_node_points - 1); ++i) {
          child_node_points[i] = node->node_points[i];
        }
        child_node_points[num_child_node_points - 1] =
            hs_split_point; // initialize newly added point

        qsort(child_node_points, num_child_node_points, sizeof(short),
              compare_short);
      }

      // this will put the time series of this node in the
      // file_buffer->buffered_list aray it will include the time series in disk
      // and those in memory

      if (!split_node(index, node, child_node_points, num_child_node_points)) {
        fprintf(stderr, "Error in dstree_index.c: could not split node %s.\n",
                node->filename);
        return FAILURE;
      }

      free(child_node_points);

      node->file_buffer->do_not_flush = true;

      if (!get_file_buffer(index, node)) {
        fprintf(stderr, "Error in dstree_index.c: could not get the file \
                           buffer for node %s.\n",
                node->filename);
        return FAILURE;
      }

      ts_type **ts_list;
      ts_list = get_all_time_series_in_node(index, node);

      // copying the contents of the the node being split
      // in case it gets flushed from memory to disk

      // printf ("splitting node %s with size %d", node->filename,
      // node->file_buffer->buffered_list_size);
      for (int idx = 0; idx < index->settings->max_leaf_size; ++idx) {
        if (node_split_policy_route_to_left(node, ts_list[idx])) {
          if (!update_node_statistics(node->left_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at left child of\
                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          if (!append_ts_to_child_node(index, node->left_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to left child of \
                           node %s\n",
                    node->filename);
            return FAILURE;
          }

        } else {
          if (!update_node_statistics(node->right_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at right child of\
                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          if (!append_ts_to_child_node(index, node->right_child,
                                       ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to right child of \
                           node %s\n",
                    node->filename);
            return FAILURE;
          }
        }
      }

      for (int i = 0; i < index->settings->max_leaf_size; ++i) {
        free(ts_list[i]);
      }

      free(ts_list);

      if (!delete_file_buffer(index, node)) {
        fprintf(stderr,
                "Error in dstree_index.c: could not delete file buffer for \
                           node %s\n",
                node->filename);
        return FAILURE;
      }

      COUNT_PARTIAL_TIME_END
      index->stats->idx_split_node_total_time += partial_time;
      index->stats->idx_split_node_input_time += partial_input_time;
      index->stats->idx_split_node_output_time += partial_output_time;
      index->stats->idx_split_node_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_split_node_seq_input_count += partial_seq_input_count;
      index->stats->idx_split_node_seq_output_count += partial_seq_output_count;
      index->stats->idx_split_node_rand_input_count += partial_rand_input_count;
      index->stats->idx_split_node_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

    } // end if_split_node
  }   // end if_node_is_leaf

  return SUCCESS;
}

// enum response dstree_index_classify_insert(struct dstree_index *index,
// ts_type * timeseries, unsigned char gt)
enum response dstree_index_classify_insert(struct dstree_index *index,
                                           ts_type *timeseries, label_type gt,
                                           unsigned int fp) {

  // traverse the index tree to find the appropriate node
  struct dstree_node *node = index->first_node;

  while (!node->is_leaf) {
    if (!update_node_statistics(node, timeseries)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (node_split_policy_route_to_left(node, timeseries))
      node = node->left_child;
    else
      node = node->right_child;
  }

  COUNT_PARTIAL_TIME_END
  index->stats->idx_traverse_tree_total_time += partial_time;
  index->stats->idx_traverse_tree_input_time += partial_input_time;
  index->stats->idx_traverse_tree_output_time += partial_output_time;
  index->stats->idx_traverse_tree_cpu_time +=
      partial_time - partial_input_time - partial_output_time;
  index->stats->idx_traverse_tree_seq_input_count += partial_seq_input_count;
  index->stats->idx_traverse_tree_seq_output_count += partial_seq_output_count;
  index->stats->idx_traverse_tree_rand_input_count += partial_rand_input_count;
  index->stats->idx_traverse_tree_rand_output_count +=
      partial_rand_output_count;

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  if (node->is_leaf) {

    if (!update_node_statistics(node, timeseries)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (!append_ts_gt_to_node(index, node, timeseries, gt, fp)) {
      fprintf(stderr, "Error in dstree_index.c: could not append \
                        time series to node %s\n",
              node->filename);
      return FAILURE;
    }

    COUNT_PARTIAL_TIME_END
    index->stats->idx_append_ts_to_leaf_total_time += partial_time;
    index->stats->idx_append_ts_to_leaf_input_time += partial_input_time;
    index->stats->idx_append_ts_to_leaf_output_time += partial_output_time;
    index->stats->idx_append_ts_to_leaf_cpu_time +=
        partial_time - partial_input_time - partial_output_time;
    index->stats->idx_append_ts_to_leaf_seq_input_count +=
        partial_seq_input_count;
    index->stats->idx_append_ts_to_leaf_seq_output_count +=
        partial_seq_output_count;
    index->stats->idx_append_ts_to_leaf_rand_input_count +=
        partial_rand_input_count;
    index->stats->idx_append_ts_to_leaf_rand_output_count +=
        partial_rand_output_count;

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    // if split needed, split the node and refresh curr_node
    if (node->node_size >= index->settings->max_leaf_size) {

      struct node_split_policy curr_node_split_policy;
      ts_type max_diff_value = (FLT_MAX * (-1));
      ts_type avg_children_range_value = 0;
      short hs_split_point = -1;
      short *child_node_points;
      int num_child_node_points = 0;
      const int num_child_segments = 2; // by default split to two subsegments

      // we want to test every possible split policy for each segment

      // for each segment
      for (int i = 0; i < node->num_node_points; ++i) {
        struct segment_sketch curr_node_segment_sketch =
            node->node_segment_sketches[i];

        // This is the QoS of this segment. QoS is the estimation quality
        // evaluated as = QoS = segment_length * (max_mean_min_mean) *
        // ((max_mean_min_mean) +
        //     (max_stdev * max_stdev))
        // The smaller the QoS, the more effective the bounds are for similarity
        // estimation

        ts_type node_range_value = range_calc(
            curr_node_segment_sketch, get_segment_length(node->node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_node_segment_split_policy =
              node->node_segment_split_policies[j];
          // to hold the two child segments
          struct segment_sketch *child_node_segment_sketches = NULL;

          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);

          if (child_node_segment_sketches == NULL) {
            fprintf(stderr, "Error in dstree_index.c: could not allocate \
                            memory for the child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr, "Error in dstree_index.c: could not allocate\
                               memory for the child node segment sketches \
                               indicators for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else {
            fprintf(stderr, "Error in dstree_index.c: Split policy was not \
                            set properly for node %s\n",
                    node->filename);
            return FAILURE;
          }

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->node_points, i));
          }

          // diff_value represents the splitting benefit
          // B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
          // the higher the diff_value, the better is the splitting

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);
          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_node_segment_split_policy;
          }
          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }
          free(child_node_segment_sketches);
        }
      }

      // add trade-off for horizontal split
      max_diff_value = max_diff_value * 2;

      // we want to test every possible split policy for each horizontal segment
      for (int i = 0; i < node->num_hs_node_points; ++i) {
        struct segment_sketch curr_hs_node_segment_sketch =
            node->hs_node_segment_sketches[i];
        ts_type node_range_value =
            range_calc(curr_hs_node_segment_sketch,
                       get_segment_length(node->hs_node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_hs_node_segment_split_policy =
              node->node_segment_split_policies[j];

          struct segment_sketch *child_node_segment_sketches =
              NULL; // to hold the two child segments
          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);
          if (child_node_segment_sketches == NULL) {
            fprintf(stderr,
                    "Error in dstree_index.c: could not allocate memory \
                            for the horizontal child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_hs_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr,
                      "Error in dstree_index.c: could not allocate memory \
                             for the horizontal child node segment sketches indicators \
                             for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else
            printf("split policy not initialized properly\n");

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->hs_node_points, i));
          }

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);

          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_hs_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_hs_node_segment_split_policy;
            hs_split_point = get_hs_split_point(
                node->node_points, curr_node_split_policy.split_from,
                curr_node_split_policy.split_to, node->num_node_points);
          }

          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }

          free(child_node_segment_sketches);
        }
      }

      node->split_policy = NULL;
      node->split_policy = malloc(sizeof(struct node_split_policy));
      if (node->split_policy == NULL) {
        fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                        for the split policy of node  %s\n",
                node->filename);
        return FAILURE;
      }
      node->split_policy->split_from = curr_node_split_policy.split_from;
      node->split_policy->split_to = curr_node_split_policy.split_to;
      node->split_policy->indicator_split_idx =
          curr_node_split_policy.indicator_split_idx;
      node->split_policy->indicator_split_value =
          curr_node_split_policy.indicator_split_value;
      node->split_policy->curr_node_segment_split_policy =
          curr_node_split_policy.curr_node_segment_split_policy;

      COUNT_PARTIAL_TIME_END
      index->stats->idx_evaluate_split_policies_total_time += partial_time;
      index->stats->idx_evaluate_split_policies_input_time +=
          partial_input_time;
      index->stats->idx_evaluate_split_policies_output_time +=
          partial_output_time;
      index->stats->idx_evaluate_split_policies_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_evaluate_split_policies_seq_input_count +=
          partial_seq_input_count;
      index->stats->idx_evaluate_split_policies_seq_output_count +=
          partial_seq_output_count;
      index->stats->idx_evaluate_split_policies_rand_input_count +=
          partial_rand_input_count;
      index->stats->idx_evaluate_split_policies_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

      // when hs_split_point stays less than 0, it means that
      // considering splitting a vertical segment is not worth it
      // according to the QoS heuristic

      if (hs_split_point < 0) {
        num_child_node_points = node->num_node_points;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have the same number of segments as parent
        for (int i = 0; i < num_child_node_points; ++i) {
          child_node_points[i] = node->node_points[i];
        }
      } else {
        num_child_node_points = node->num_node_points + 1;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have one additional segment than the parent
        for (int i = 0; i < (num_child_node_points - 1); ++i) {
          child_node_points[i] = node->node_points[i];
        }
        child_node_points[num_child_node_points - 1] =
            hs_split_point; // initialize newly added point

        qsort(child_node_points, num_child_node_points, sizeof(short),
              compare_short);
      }

      // this will put the time series of this node in the
      // file_buffer->buffered_list aray it will include the time series in disk
      // and those in memory

      if (!split_node(index, node, child_node_points, num_child_node_points)) {
        fprintf(stderr, "Error in dstree_index.c: could not split node %s.\n",
                node->filename);
        return FAILURE;
      }

      free(child_node_points);

      node->file_buffer->do_not_flush = true;

      if (!get_file_buffer(index, node)) {
        fprintf(stderr, "Error in dstree_index.c: could not get the file \
                           buffer for node %s.\n",
                node->filename);
        return FAILURE;
      }

      ts_type **ts_list;
      ts_list = get_all_time_series_in_node(index, node);

      // copying the contents of the the node being split
      // in case it gets flushed from memory to disk

      // printf ("splitting node %s with size %d", node->filename,
      // node->file_buffer->buffered_list_size);
      for (int idx = 0; idx < index->settings->max_leaf_size; ++idx) {
        if (node_split_policy_route_to_left(node, ts_list[idx])) {
          if (!update_node_statistics(node->left_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at left child of\                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          // if(!append_ts_gt_to_child_node(index,
          // node->left_child,ts_list[idx],node->gt[idx], node->file_pos[idx]))
          if (!append_ts_gt_to_child_node(index, node->left_child, ts_list[idx],
                                          node, idx, fp)) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to left child of			\
                           node %s\n",
                    node->filename);
            return FAILURE;
          }
        } else {
          if (!update_node_statistics(node->right_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at right child of\
                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          // if(!append_ts_gt_to_child_node(index,
          // node->right_child,ts_list[idx],node->gt[idx],node->file_pos[idx]))
          if (!append_ts_gt_to_child_node(index, node->right_child,
                                          ts_list[idx], node, idx, fp)) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to right child of		\
                           node %s\n",
                    node->filename);
            return FAILURE;
          }
        }
      }

      for (int i = 0; i < index->settings->max_leaf_size; ++i) {
        free(ts_list[i]);
      }

      free(ts_list);

      free(node->gt);
      node->gt = NULL;

      if (fp != -1) {
        free(node->fp);
        node->fp = NULL;
      }

      if (!delete_file_buffer(index, node)) {
        fprintf(stderr,
                "Error in dstree_index.c: could not delete file buffer for \
                           node %s\n",
                node->filename);
        return FAILURE;
      }

      COUNT_PARTIAL_TIME_END
      index->stats->idx_split_node_total_time += partial_time;
      index->stats->idx_split_node_input_time += partial_input_time;
      index->stats->idx_split_node_output_time += partial_output_time;
      index->stats->idx_split_node_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_split_node_seq_input_count += partial_seq_input_count;
      index->stats->idx_split_node_seq_output_count += partial_seq_output_count;
      index->stats->idx_split_node_rand_input_count += partial_rand_input_count;
      index->stats->idx_split_node_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

    } // end if_split_node
  }   // end if_node_is_leaf

  return SUCCESS;
}

enum response dstree_index_write(struct dstree_index *index) {
  fprintf(stderr, ">>> Storing index: %s\n", index->settings->root_directory);
  const char *filename =
      malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
  filename = strcpy(filename, index->settings->root_directory);
  filename = strcat(filename, "root.idx\0");

  COUNT_PARTIAL_RAND_OUTPUT
  COUNT_PARTIAL_OUTPUT_TIME_START
  FILE *file = fopen(filename, "wb");
  COUNT_PARTIAL_OUTPUT_TIME_END

  free(filename);

  if (index->settings->classify) {
    const char *gt_filename =
        malloc(sizeof(char) * (strlen(index->settings->root_directory) + 8));
    strcpy(gt_filename, index->settings->root_directory);
    strcat(gt_filename, "gt.idx\0");
    index->gt_file = fopen(gt_filename, "wb");

    free(gt_filename);

    index->gt_pos_ctr = 0;
  }

  if (index->settings->track_file_pos)
  {
    const char *fp_filename =
        malloc(sizeof(char) * (strlen(index->settings->root_directory) + 8));
    strcpy(fp_filename, index->settings->root_directory);
    strcat(fp_filename, "fp.idx\0");
    index->fp_file = fopen(fp_filename, "wb");

    if (index->fp_file == NULL)
    {
      fprintf(stderr,
              "Error in dstree_index.c: Could not open"
              " the fp file. Reason = %s\n", strerror(errno));
      return FAILURE;
    }

    free(fp_filename);

    index->fp_pos_ctr = 0;
  }

  /* start kashif changes */
  if (index->settings->track_vector)
  {
    //open vid.edx
    const char *vid_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 9));
    strcpy(vid_filename, index->settings->root_directory);
    strcat(vid_filename, "vid.idx\0");
    index->vid_file = fopen(vid_filename, "wb");

    if (index->vid_file == NULL)
    {
      fprintf(stderr,
              "Error in dstree_index.c: Could not open"
              " the vid file. Reason = %s\n", strerror(errno));
      return FAILURE;
    }

    free(vid_filename);

    index->vid_pos_ctr = 0;
  }
  /* end kashif changes */

  if (file == NULL) {
    fprintf(stderr,
            "Error in dstree_index.c: Could not open"
            " the index file. Reason = %s\n",
            strerror(errno));
    return FAILURE;
  }

  unsigned int timeseries_size = index->settings->timeseries_size;
  unsigned int max_leaf_size = index->settings->max_leaf_size;
  unsigned int init_segments = index->settings->init_segments;
  double buffered_memory_size = index->settings->buffered_memory_size;
  unsigned char classify = index->settings->classify;
  unsigned char track_file_pos = index->settings->track_file_pos;
  unsigned char track_vector =  index->settings->track_vector;
  unsigned int dataset_size = index->first_node->node_size;

  // SETTINGS DATA
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_OUTPUT_TIME_START
  fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);
  fwrite(&buffered_memory_size, sizeof(double), 1, file);
  fwrite(&timeseries_size, sizeof(unsigned int), 1, file);
  fwrite(&init_segments, sizeof(unsigned int), 1, file);
  fwrite(&max_leaf_size, sizeof(unsigned int), 1, file);
  fwrite(&dataset_size, sizeof(unsigned int), 1, file);
  fwrite(&classify, sizeof(unsigned char), 1, file);
  fwrite(&track_file_pos, sizeof(unsigned char), 1, file);
  fwrite(&track_vector, sizeof(unsigned char), 1, file);
  COUNT_PARTIAL_OUTPUT_TIME_END

  // NODES AND FILE BUFFERS
  dstree_node_write(index, index->first_node, file);
  COUNT_PARTIAL_OUTPUT_TIME_START
  fseek(file, 0L, SEEK_SET);
  fwrite(&leaf_nodes_count, sizeof(unsigned long), 1, file);
  fclose(file);
  COUNT_PARTIAL_OUTPUT_TIME_END

  return SUCCESS;
}

struct dstree_index *dstree_index_read(const char *root_directory) {
  if (chdir(root_directory) != 0) {
    fprintf(stderr, "The index directory does not exist. "
                    "Please provide a valid directory.\n");
    exit(-1);
  }

  fprintf(stderr, ">>> Loading index: %s\n", root_directory);

  const char *filename = malloc(sizeof(char) * (strlen(root_directory) + 9));
  filename = strcpy(filename, root_directory);
  filename = strcat(filename, "root.idx\0");

  // printf ("Index file: %s\n",filename);
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  FILE *file = fopen(filename, "rb");
  COUNT_PARTIAL_INPUT_TIME_END

  free(filename);

  if (file == NULL) {
    fprintf(stderr,
            "Error in dstree_index.c: Could not open"
            " the index file. Reason = %s\n",
            strerror(errno));
    return FAILURE;
  }

  unsigned long count_leaves = 0;
  unsigned int timeseries_size = 0;
  unsigned int max_leaf_size = 0;
  unsigned int init_segments = 0;
  double buffered_memory_size = 0;
  boolean is_index_new = 0;
  unsigned char classify = 0;
  unsigned char track_file_pos = 0;
  unsigned char track_vector = 0;
  unsigned int dataset_size = 0;

  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fread(&count_leaves, sizeof(unsigned long), 1, file);
  fread(&buffered_memory_size, sizeof(double), 1, file);
  fread(&timeseries_size, sizeof(unsigned int), 1, file);
  fread(&init_segments, sizeof(unsigned int), 1, file);
  fread(&max_leaf_size, sizeof(unsigned int), 1, file);
  fread(&dataset_size, sizeof(unsigned int), 1, file);
  fread(&classify, sizeof(unsigned char), 1, file);
  fread(&track_file_pos, sizeof(unsigned char), 1, file);
  fread(&track_vector, sizeof(unsigned char), 1, file);
  COUNT_PARTIAL_INPUT_TIME_END

  struct dstree_index_settings *index_settings = dstree_index_settings_init(
      root_directory, timeseries_size, init_segments, max_leaf_size,
      buffered_memory_size, is_index_new, NULL);

  struct dstree_index *index = dstree_index_init(index_settings);

  index->stats->leaves_heights = calloc(count_leaves, sizeof(int));
  index->stats->leaves_sizes = calloc(count_leaves, sizeof(int));
  index->stats->leaves_counter = 0;

  index->settings->classify = classify;
  index->settings->track_file_pos = track_file_pos;
  index->settings->track_vector = track_vector;

  if (index->settings->classify) {
    const char *gt_filename =
        malloc(sizeof(char) * (strlen(index->settings->root_directory) + 7));
    strcpy(gt_filename, index->settings->root_directory);
    strcat(gt_filename, "gt.idx\0");
    index->gt_filename = gt_filename;
    index->gt_file = fopen(gt_filename, "rb");

    if (index->gt_file == NULL) {
      fprintf(stderr,
              "Error in dstree_index.c: Could not open"
              "the index file. Reason = %s\n",
              strerror(errno));
      return FAILURE;
    }
    index->gt_pos_ctr = 0;
    index->gt_cache = NULL;
    // index->gt_cache = malloc(sizeof(unsigned char) * dataset_size);
    // fread(index->gt_cache, sizeof(unsigned char), dataset_size,
    // index->gt_file);
    index->gt_cache = malloc(sizeof(label_type) * dataset_size);
    fread(index->gt_cache, sizeof(label_type), dataset_size, index->gt_file);
  }

  if (index->settings->track_file_pos) {
    const char *fp_filename =
        malloc(sizeof(char) * (strlen(index->settings->root_directory) + 7));
    strcpy(fp_filename, index->settings->root_directory);
    strcat(fp_filename, "fp.idx\0");
    index->fp_filename = fp_filename;
    index->fp_file = fopen(fp_filename, "rb");

    if (index->fp_file == NULL) {
      fprintf(stderr,
              "Error in dstree_index.c: Could not open"
              "the index file. Reason = %s\n",
              strerror(errno));
      return FAILURE;
    }
    index->fp_pos_ctr = 0;
    index->fp_cache = NULL;
    // index->gt_cache = malloc(sizeof(unsigned char) * dataset_size);
    // fread(index->gt_cache, sizeof(unsigned char), dataset_size,
    // index->gt_file);
    index->fp_cache = malloc(sizeof(unsigned int) * dataset_size);
    fread(index->fp_cache, sizeof(unsigned int), dataset_size, index->fp_file);
  }

  /* start kashif changes */
  if (index->settings->track_vector)
  {
    // open vid.idx
    const char *vid_filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 8));
    strcpy(vid_filename, index->settings->root_directory);
    strcat(vid_filename, "vid.idx\0");
    index->vid_filename = vid_filename;
    index->vid_file = fopen(vid_filename, "rb");

    if (index->vid_file == NULL)
    {
      fprintf(stderr, "Error in dstree_index.c: Could not open"
                      "the index file. Reason = %s\n",
              strerror(errno));
      return FAILURE;
    }
    index->vid_pos_ctr = 0;
    index->vid_cache = NULL;

    index->vid_cache = malloc(sizeof(struct vid) * dataset_size);
    fread(index->vid_cache, sizeof(struct vid), dataset_size, index->vid_file);

    // FILE * sc_file;
    // sc_file = fopen("sanity_check_2.txt", "a+");
    // for(int i = 0; i < dataset_size; i++){
    //   fprintf(sc_file, "(%u, %u, %u)", index->vid_cache[i].table_id, 
    //       index->vid_cache[i].set_id, 
    //       index->vid_cache[i].pos);
    // }
    // fclose(sc_file);
  }

  /* end kashif changes */

  index->first_node = dstree_node_read(index, file);
  COUNT_PARTIAL_INPUT_TIME_START
  fclose(file);
  COUNT_PARTIAL_INPUT_TIME_END

  return index;
}

enum response dstree_update_index_stats(struct dstree_index *index,
                                        struct dstree_node *node) {

  unsigned int height = node->level + 1;
  unsigned int threshold = index->settings->max_leaf_size;
  // at this point all time series are on disk
  unsigned int node_size = node->node_size;

  double node_fill_factor = (node_size * 100.0) / threshold;

  if (node_fill_factor < index->stats->min_fill_factor) {
    index->stats->min_fill_factor = node_fill_factor;
  }
  if (node_fill_factor > index->stats->max_fill_factor) {
    index->stats->max_fill_factor = node_fill_factor;
  }
  if (height < index->stats->min_height) {
    index->stats->min_height = height;
  }
  if (height > index->stats->max_height) {
    index->stats->max_height = height;
  }

  if (node_size == 0) {
    COUNT_EMPTY_LEAF_NODE
  }
  index->stats->sum_fill_factor += node_fill_factor;
  index->stats->sum_squares_fill_factor += pow(node_fill_factor, 2);

  index->stats->sum_height += height;
  index->stats->sum_squares_height += pow(height, 2);

  // printf("node size = %d\n", node_size);
  COUNT_TOTAL_TS(node_size)
}

enum response dstree_node_write(struct dstree_index *index,
                                struct dstree_node *node, FILE *file) {
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_SEQ_OUTPUT
  COUNT_PARTIAL_OUTPUT_TIME_START
  fwrite(&(node->is_leaf), sizeof(unsigned char), 1, file);
  fwrite(&(node->node_size), sizeof(unsigned int), 1, file);
  fwrite(&(node->level), sizeof(unsigned int), 1, file);
  COUNT_PARTIAL_OUTPUT_TIME_END

  if (node->is_leaf) {
    if (node->filename != NULL) {
      int filename_size = strlen(node->filename);
      COUNT_PARTIAL_SEQ_OUTPUT
      COUNT_PARTIAL_SEQ_OUTPUT
      COUNT_PARTIAL_OUTPUT_TIME_START
      fwrite(&filename_size, sizeof(int), 1, file);
      fwrite(node->filename, sizeof(char), filename_size, file);

      COUNT_PARTIAL_SEQ_OUTPUT
      COUNT_PARTIAL_SEQ_OUTPUT
      fwrite(&(node->num_node_points), sizeof(short), 1, file);
      fwrite(node->node_points, sizeof(short), node->num_node_points, file);

      for (int i = 0; i < node->num_node_points; ++i) {
        COUNT_PARTIAL_SEQ_OUTPUT
        COUNT_PARTIAL_SEQ_OUTPUT
        fwrite(&(node->node_segment_sketches[i].num_indicators), sizeof(int), 1,
               file);
        fwrite(node->node_segment_sketches[i].indicators, sizeof(ts_type),
               node->node_segment_sketches[i].num_indicators, file);
      }
      COUNT_PARTIAL_OUTPUT_TIME_END

      COUNT_PARTIAL_OUTPUT_TIME_END
      flush_buffer_to_disk(index, node);
      if (index->settings->classify) {
        node->gt_pos = index->gt_pos_ctr;

        COUNT_PARTIAL_OUTPUT_TIME_START
        fwrite(&(node->gt_pos), sizeof(unsigned int), 1, file);
        // fwrite(node->gt, sizeof(unsigned char), node->node_size,
        // index->gt_file);
        fwrite(node->gt, sizeof(label_type), node->node_size, index->gt_file);
        COUNT_PARTIAL_OUTPUT_TIME_END

        index->gt_pos_ctr += node->node_size;
      }
      if (index->settings->track_file_pos) {
        node->fp_pos = index->fp_pos_ctr;

        COUNT_PARTIAL_OUTPUT_TIME_START
        fwrite(&(node->fp_pos), sizeof(unsigned int), 1, file);
        // fwrite(node->gt, sizeof(unsigned char), node->node_size,
        // index->gt_file);
        fwrite(node->fp, sizeof(unsigned int), node->node_size, index->fp_file);
        COUNT_PARTIAL_OUTPUT_TIME_END

        index->fp_pos_ctr += node->node_size;
      }
      if (index->settings->track_vector)
      {
        node->vid_pos = index->vid_pos_ctr;

        // save vids to file
        COUNT_PARTIAL_OUTPUT_TIME_START
        fwrite(&(node->vid_pos), sizeof(unsigned int), 1, file);

        FILE * scf = fopen("sanity_check_1.txt", "a+");
        for(int w = 0; w < node->node_size; w++)
        {
          fprintf(scf, "(%u, %u, %u)", node->vid[w].table_id, node->vid[w].set_id, node->vid[w].pos);
        }
        fclose(scf);

        fwrite(node->vid, sizeof(struct vid), node->node_size, index->vid_file);
        COUNT_PARTIAL_OUTPUT_TIME_END

        index->vid_pos_ctr += node->node_size;
      }
      COUNT_LEAF_NODE
      // collect stats while traversing the index
      dstree_update_index_stats(index, node);
    } 
    else 
    {
      int filename_size = 0;
      COUNT_PARTIAL_SEQ_OUTPUT
      COUNT_PARTIAL_OUTPUT_TIME_START
      fwrite(&filename_size, sizeof(int), 1, file);
      COUNT_PARTIAL_OUTPUT_TIME_END
    }
  } else {
    COUNT_PARTIAL_SEQ_OUTPUT
    COUNT_PARTIAL_SEQ_OUTPUT
    COUNT_PARTIAL_SEQ_OUTPUT
    COUNT_PARTIAL_OUTPUT_TIME_START
    fwrite(node->split_policy, sizeof(struct node_split_policy), 1, file);
    fwrite(&(node->num_node_points), sizeof(short), 1, file);
    fwrite(node->node_points, sizeof(short), node->num_node_points, file);
    for (int i = 0; i < node->num_node_points; ++i) {
      COUNT_PARTIAL_SEQ_OUTPUT
      COUNT_PARTIAL_SEQ_OUTPUT
      fwrite(&(node->node_segment_sketches[i].num_indicators), sizeof(int), 1,
             file);
      fwrite(node->node_segment_sketches[i].indicators, sizeof(ts_type),
             node->node_segment_sketches[i].num_indicators, file);
    }
    COUNT_PARTIAL_OUTPUT_TIME_END
  }

  if (!node->is_leaf) {
    dstree_node_write(index, node->left_child, file);
    dstree_node_write(index, node->right_child, file);
  }

  return SUCCESS;
}

struct dstree_node *dstree_node_read(struct dstree_index *index, FILE *file) 
{

  struct dstree_node *node = NULL;
  // to initialize node values for leaf and internal nodes
  node = dstree_leaf_node_init(index->settings);

  unsigned char is_leaf = 0;
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_SEQ_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fread(&is_leaf, sizeof(unsigned char), 1, file);
  fread(&(node->node_size), sizeof(unsigned int), 1, file);
  fread(&(node->level), sizeof(unsigned int), 1, file);
  COUNT_PARTIAL_INPUT_TIME_END
  if (is_leaf) {
    node->is_leaf = 1;
    dstree_file_buffer_init(node);

    int filename_size = 0;

    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(&filename_size, sizeof(int), 1, file);
    COUNT_PARTIAL_INPUT_TIME_END

    node->file_buffer->disk_count = node->node_size;

    if (filename_size > 0) {
      node->filename = malloc(sizeof(char) * (filename_size + 1));

      COUNT_PARTIAL_SEQ_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      fread(node->filename, sizeof(char), filename_size, file);
      COUNT_PARTIAL_INPUT_TIME_END

      node->filename[filename_size] = '\0';

      COUNT_PARTIAL_SEQ_INPUT
      COUNT_PARTIAL_SEQ_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      fread(&(node->num_node_points), sizeof(short), 1, file);
      COUNT_PARTIAL_INPUT_TIME_END
      node->node_points = NULL;
      node->node_points = malloc(sizeof(short) * node->num_node_points);
      COUNT_PARTIAL_INPUT_TIME_START
      fread(node->node_points, sizeof(short), node->num_node_points, file);
      COUNT_PARTIAL_INPUT_TIME_END

      node->node_segment_sketches = NULL;
      node->node_segment_sketches =
          malloc(sizeof(struct segment_sketch) * node->num_node_points);

      for (int i = 0; i < node->num_node_points; ++i) {
        COUNT_PARTIAL_SEQ_INPUT
        COUNT_PARTIAL_SEQ_INPUT
        COUNT_PARTIAL_INPUT_TIME_START
        fread(&(node->node_segment_sketches[i].num_indicators), sizeof(int), 1,
              file);
        COUNT_PARTIAL_INPUT_TIME_END
        node->node_segment_sketches[i].indicators = NULL;
        node->node_segment_sketches[i].indicators = malloc(
            sizeof(ts_type) * node->node_segment_sketches[i].num_indicators);

        COUNT_PARTIAL_INPUT_TIME_START

        fread(node->node_segment_sketches[i].indicators, sizeof(ts_type),
              node->node_segment_sketches[i].num_indicators, file);
        COUNT_PARTIAL_INPUT_TIME_END
      }

      COUNT_LEAF_NODE
      
      if(index->settings->classify)
      {
        fread(&(node->gt_pos), sizeof(unsigned int), 1, file);
        fread(&(node->gt), sizeof(label_type), node->node_size, index->gt_file);
      }
      
      if(index->settings->track_file_pos)
      {
        fread(&(node->fp_pos), sizeof(unsigned int), 1, file);
      }

      if(index->settings->track_vector)
      {
        fread(&(node->vid_pos), sizeof(unsigned int), 1, file);
        // vids are read to vid_cache
        // fread(&(node->vid), sizeof(struct vid), node->node_size, index->vid_file);
      }

      index->stats->leaves_heights[index->stats->leaves_counter] =
          node->level + 1;
      index->stats->leaves_sizes[index->stats->leaves_counter] =
          node->node_size;
      ++(index->stats->leaves_counter);

      dstree_update_index_stats(index, node);
    } else {
      node->filename = NULL;
      node->node_size = 0;
    }
    // get all timeseries for this node into its this file buffer

    // get_all_time_series_in_node(index,node);
    // node->file_buffer->disk_count = 0;
  } else {
    node->split_policy = NULL;
    node->split_policy = malloc(sizeof(struct node_split_policy));
    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(node->split_policy, sizeof(struct node_split_policy), 1, file);
    fread(&(node->num_node_points), sizeof(short), 1, file);
    COUNT_PARTIAL_INPUT_TIME_END

    node->node_points = NULL;
    node->node_points = malloc(sizeof(short) * node->num_node_points);

    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(node->node_points, sizeof(short), node->num_node_points, file);
    COUNT_PARTIAL_INPUT_TIME_END

    node->node_segment_sketches = NULL;
    node->node_segment_sketches =
        malloc(sizeof(struct segment_sketch) * node->num_node_points);
    for (int i = 0; i < node->num_node_points; ++i) {
      COUNT_PARTIAL_SEQ_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      fread(&(node->node_segment_sketches[i].num_indicators), sizeof(int), 1,
            file);
      COUNT_PARTIAL_INPUT_TIME_END

      node->node_segment_sketches[i].indicators = NULL;
      node->node_segment_sketches[i].indicators = malloc(
          sizeof(ts_type) * node->node_segment_sketches[i].num_indicators);
      COUNT_PARTIAL_SEQ_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      fread(node->node_segment_sketches[i].indicators, sizeof(ts_type),
            node->node_segment_sketches[i].num_indicators, file);
      COUNT_PARTIAL_INPUT_TIME_END
    }

    node->file_buffer = NULL;
    node->is_leaf = 0;
    node->filename = NULL;
    // node->node_size = 0;
  }

  if (!is_leaf) {
    node->left_child = dstree_node_read(index, file);
    node->right_child = dstree_node_read(index, file);
  }

  return node;
}

/* start kashif changes */
enum response dstree_index_insert_vector(struct dstree_index *index,
                                  ts_type *vector, unsigned int table_id, 
                                  unsigned int set_id, unsigned int pos, 
                                  char * raw_data_file) 
{
  // traverse the index tree to find the appropriate node
  struct dstree_node *node = index->first_node;

  while (!node->is_leaf) {
    if (!update_node_statistics(node, vector)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (node_split_policy_route_to_left(node, vector))
      node = node->left_child;
    else
      node = node->right_child;
  }

  COUNT_PARTIAL_TIME_END
  index->stats->idx_traverse_tree_total_time += partial_time;
  index->stats->idx_traverse_tree_input_time += partial_input_time;
  index->stats->idx_traverse_tree_output_time += partial_output_time;
  index->stats->idx_traverse_tree_cpu_time +=
      partial_time - partial_input_time - partial_output_time;
  index->stats->idx_traverse_tree_seq_input_count += partial_seq_input_count;
  index->stats->idx_traverse_tree_seq_output_count += partial_seq_output_count;
  index->stats->idx_traverse_tree_rand_input_count += partial_rand_input_count;
  index->stats->idx_traverse_tree_rand_output_count +=
      partial_rand_output_count;

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  if (node->is_leaf) {
    // printf("curr node size = %d\n", node->node_size);
    // for(int z = 0; z < node->node_size; z++)
    // {
    //   printf("v(%u, %u, %u)", node->vid[z].table_id, node->vid[z].set_id, node->vid[z].pos);
    // }
    // printf("\n\n");

    if (!update_node_statistics(node, vector)) {
      fprintf(stderr, "Error in dstree_index.c: could not update \
                        statistics at node %s\n",
              node->filename);
      return FAILURE;
    }

    if (!append_vector_to_node(index, node, vector, table_id, set_id, pos, raw_data_file)) {
      fprintf(stderr, "Error in dstree_index.c: could not append \
                        time series to node %s\n",
              node->filename);
      return FAILURE;
    }

    COUNT_PARTIAL_TIME_END
    index->stats->idx_append_ts_to_leaf_total_time += partial_time;
    index->stats->idx_append_ts_to_leaf_input_time += partial_input_time;
    index->stats->idx_append_ts_to_leaf_output_time += partial_output_time;
    index->stats->idx_append_ts_to_leaf_cpu_time +=
        partial_time - partial_input_time - partial_output_time;
    index->stats->idx_append_ts_to_leaf_seq_input_count +=
        partial_seq_input_count;
    index->stats->idx_append_ts_to_leaf_seq_output_count +=
        partial_seq_output_count;
    index->stats->idx_append_ts_to_leaf_rand_input_count +=
        partial_rand_input_count;
    index->stats->idx_append_ts_to_leaf_rand_output_count +=
        partial_rand_output_count;

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    // if split needed, split the node and refresh curr_node
    if (node->node_size >= index->settings->max_leaf_size)
    {
      struct node_split_policy curr_node_split_policy;
      ts_type max_diff_value = (FLT_MAX * (-1));
      ts_type avg_children_range_value = 0;
      short hs_split_point = -1;
      short *child_node_points;
      int num_child_node_points = 0;
      const int num_child_segments = 2; // by default split to two subsegments

      // we want to test every possible split policy for each segment

      // for each segment
      for (int i = 0; i < node->num_node_points; ++i) 
      {
        struct segment_sketch curr_node_segment_sketch =
            node->node_segment_sketches[i];

        // This is the QoS of this segment. QoS is the estimation quality
        // evaluated as = QoS = segment_length * (max_mean_min_mean) *
        // ((max_mean_min_mean) +
        //     (max_stdev * max_stdev))
        // The smaller the QoS, the more effective the bounds are for similarity
        // estimation

        ts_type node_range_value = range_calc(
            curr_node_segment_sketch, get_segment_length(node->node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_node_segment_split_policy =
              node->node_segment_split_policies[j];
          // to hold the two child segments
          struct segment_sketch *child_node_segment_sketches = NULL;

          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);

          if (child_node_segment_sketches == NULL) {
            fprintf(stderr, "Error in dstree_index.c: could not allocate \
                            memory for the child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr, "Error in dstree_index.c: could not allocate\
                               memory for the child node segment sketches \
                               indicators for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_node_segment_split_policy, curr_node_segment_sketch,
                child_node_segment_sketches);
          else {
            fprintf(stderr, "Error in dstree_index.c: Split policy was not \
                            set properly for node %s\n",
                    node->filename);
            return FAILURE;
          }

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->node_points, i));
          }

          // diff_value represents the splitting benefit
          // B = QoS(N) - (QoS_leftNode + QoS_rightNode)/2
          // the higher the diff_value, the better is the splitting

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);
          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_node_segment_split_policy;
          }
          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }
          free(child_node_segment_sketches);
        }
      }

      // add trade-off for horizontal split
      max_diff_value = max_diff_value * 2;

      // we want to test every possible split policy for each horizontal segment
      for (int i = 0; i < node->num_hs_node_points; ++i) {
        struct segment_sketch curr_hs_node_segment_sketch =
            node->hs_node_segment_sketches[i];
        ts_type node_range_value =
            range_calc(curr_hs_node_segment_sketch,
                       get_segment_length(node->hs_node_points, i));

        // for every split policy
        for (int j = 0; j < node->num_node_segment_split_policies; ++j) {
          struct node_segment_split_policy curr_hs_node_segment_split_policy =
              node->node_segment_split_policies[j];

          struct segment_sketch *child_node_segment_sketches =
              NULL; // to hold the two child segments
          child_node_segment_sketches =
              malloc(sizeof(struct segment_sketch) * num_child_segments);
          if (child_node_segment_sketches == NULL) {
            fprintf(stderr,
                    "Error in dstree_index.c: could not allocate memory \
                            for the horizontal child node segment sketches for \
                            node  %s\n",
                    node->filename);
            return FAILURE;
          }

          for (int k = 0; k < num_child_segments; ++k) {
            child_node_segment_sketches[k].indicators = NULL;
            child_node_segment_sketches[k].indicators = malloc(
                sizeof(ts_type) * curr_hs_node_segment_sketch.num_indicators);
            if (child_node_segment_sketches[k].indicators == NULL) {
              fprintf(stderr,
                      "Error in dstree_index.c: could not allocate memory \
                             for the horizontal child node segment sketches indicators \
                             for node  %s\n",
                      node->filename);
              return FAILURE;
            }
          }

          if (is_split_policy_mean(curr_hs_node_segment_split_policy))
            mean_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else if (is_split_policy_stdev(curr_hs_node_segment_split_policy))
            stdev_node_segment_split_policy_split(
                &curr_hs_node_segment_split_policy, curr_hs_node_segment_sketch,
                child_node_segment_sketches);
          else
            printf("split policy not initialized properly\n");

          ts_type range_values[num_child_segments];
          for (int k = 0; k < num_child_segments; ++k) {
            struct segment_sketch child_node_segment_sketch =
                child_node_segment_sketches[k];
            range_values[k] =
                range_calc(child_node_segment_sketch,
                           get_segment_length(node->hs_node_points, i));
          }

          avg_children_range_value =
              calc_mean(range_values, 0, num_child_segments);

          ts_type diff_value = node_range_value - avg_children_range_value;

          if (diff_value > max_diff_value) {
            max_diff_value = diff_value;
            curr_node_split_policy.split_from =
                get_segment_start(node->hs_node_points, i);
            curr_node_split_policy.split_to =
                get_segment_end(node->hs_node_points, i);
            curr_node_split_policy.indicator_split_idx =
                curr_hs_node_segment_split_policy.indicator_split_idx;
            curr_node_split_policy.indicator_split_value =
                curr_hs_node_segment_split_policy.indicator_split_value;
            curr_node_split_policy.curr_node_segment_split_policy =
                curr_hs_node_segment_split_policy;
            hs_split_point = get_hs_split_point(
                node->node_points, curr_node_split_policy.split_from,
                curr_node_split_policy.split_to, node->num_node_points);
          }

          for (int k = 0; k < num_child_segments; ++k) {
            free(child_node_segment_sketches[k].indicators);
          }

          free(child_node_segment_sketches);
        }
      }

      node->split_policy = NULL;
      node->split_policy = malloc(sizeof(struct node_split_policy));
      if (node->split_policy == NULL) {
        fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                        for the split policy of node  %s\n",
                node->filename);
        return FAILURE;
      }
      node->split_policy->split_from = curr_node_split_policy.split_from;
      node->split_policy->split_to = curr_node_split_policy.split_to;
      node->split_policy->indicator_split_idx =
          curr_node_split_policy.indicator_split_idx;
      node->split_policy->indicator_split_value =
          curr_node_split_policy.indicator_split_value;
      node->split_policy->curr_node_segment_split_policy =
          curr_node_split_policy.curr_node_segment_split_policy;

      COUNT_PARTIAL_TIME_END
      index->stats->idx_evaluate_split_policies_total_time += partial_time;
      index->stats->idx_evaluate_split_policies_input_time +=
          partial_input_time;
      index->stats->idx_evaluate_split_policies_output_time +=
          partial_output_time;
      index->stats->idx_evaluate_split_policies_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_evaluate_split_policies_seq_input_count +=
          partial_seq_input_count;
      index->stats->idx_evaluate_split_policies_seq_output_count +=
          partial_seq_output_count;
      index->stats->idx_evaluate_split_policies_rand_input_count +=
          partial_rand_input_count;
      index->stats->idx_evaluate_split_policies_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

      // when hs_split_point stays less than 0, it means that
      // considering splitting a vertical segment is not worth it
      // according to the QoS heuristic

      if (hs_split_point < 0) 
      {
        num_child_node_points = node->num_node_points;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have the same number of segments as parent
        for (int i = 0; i < num_child_node_points; ++i) {
          child_node_points[i] = node->node_points[i];
        }
      }
      else 
      {
        num_child_node_points = node->num_node_points + 1;
        child_node_points = NULL;
        child_node_points = malloc(sizeof(short) * num_child_node_points);
        if (child_node_points == NULL) {
          fprintf(stderr, "Error in dstree_index.c: could not allocate memory \
                          for the child node segment points node  %s\n",
                  node->filename);
          return FAILURE;
        }
        // children will have one additional segment than the parent
        for (int i = 0; i < (num_child_node_points - 1); ++i) {
          child_node_points[i] = node->node_points[i];
        }
        child_node_points[num_child_node_points - 1] =
            hs_split_point; // initialize newly added point

        qsort(child_node_points, num_child_node_points, sizeof(short),
              compare_short);
      }

      // this will put the time series of this node in the
      // file_buffer->buffered_list aray it will include the time series in disk
      // and those in memory

      if (!split_node(index, node, child_node_points, num_child_node_points)) {
        fprintf(stderr, "Error in dstree_index.c: could not split node %s.\n",
                node->filename);
        return FAILURE;
      }

      free(child_node_points);

      node->file_buffer->do_not_flush = true;

      if (!get_file_buffer(index, node)) {
        fprintf(stderr, "Error in dstree_index.c: could not get the file \
                           buffer for node %s.\n",
                node->filename);
        return FAILURE;
      }

      ts_type **ts_list;
      ts_list = get_all_time_series_in_node(index, node);

      // copying the contents of the the node being split
      // in case it gets flushed from memory to disk

      // printf ("splitting node %s with size %d", node->filename,
      // node->file_buffer->buffered_list_size);
      for (int idx = 0; idx < index->settings->max_leaf_size; ++idx) {
        if (node_split_policy_route_to_left(node, ts_list[idx])) {
          if (!update_node_statistics(node->left_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at left child of\
                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          if (!append_vector_to_child_node(index, node->left_child, ts_list[idx],
                  node->vid[idx].table_id, node->vid[idx].set_id, node->vid[idx].pos, node->vid[idx].raw_data_file)) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to left child of \
                           node %s\n",
                    node->filename);
            return FAILURE;
          }

        } else {
          if (!update_node_statistics(node->right_child, ts_list[idx])) {
            fprintf(stderr, "Error in dstree_index.c: could not update \
                             statistics at right child of\
                             node %s\n",
                    node->filename);
            return FAILURE;
          }

          if (!append_vector_to_child_node(index, node->right_child,ts_list[idx],
                  node->vid[idx].table_id, node->vid[idx].set_id, node->vid[idx].pos, node->vid[idx].raw_data_file)) {
            fprintf(stderr, "Error in dstree_index.c: could not append \
                           time series to right child of \
                           node %s\n",
                    node->filename);
            return FAILURE;
          }
        }
      }

      for (int i = 0; i < index->settings->max_leaf_size; ++i) {
        free(ts_list[i]);
      }

      free(ts_list);

      if (!delete_file_buffer(index, node)) {
        fprintf(stderr,
                "Error in dstree_index.c: could not delete file buffer for \
                           node %s\n",
                node->filename);
        return FAILURE;
      }

      COUNT_PARTIAL_TIME_END
      index->stats->idx_split_node_total_time += partial_time;
      index->stats->idx_split_node_input_time += partial_input_time;
      index->stats->idx_split_node_output_time += partial_output_time;
      index->stats->idx_split_node_cpu_time +=
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_split_node_seq_input_count += partial_seq_input_count;
      index->stats->idx_split_node_seq_output_count += partial_seq_output_count;
      index->stats->idx_split_node_rand_input_count += partial_rand_input_count;
      index->stats->idx_split_node_rand_output_count +=
          partial_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

    } // end if_split_node
  }   // end if_node_is_leaf

  return SUCCESS;
}
/* end kashif changes */

/*
void cache_gt_file(struct dstree_index *index) {
    COUNT_INPUT_TIME_START
      //fread(index->gt_cache, sizeof(unsigned char),
index->first_node->node_size, index->gt_file); fread(index->gt_cache,
sizeof(label_type), index->first_node->node_size, index->gt_file);
    COUNT_INPUT_TIME_END
}
*/

