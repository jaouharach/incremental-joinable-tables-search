//
//  main.c
//  ds-tree C version
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#include "../config.h"
#include "../globals.h"
#include "../include/systemutils.h"

#include "../include/dstree_file_buffer.h"
#include "../include/dstree_file_buffer_manager.h"
#include "../include/dstree_file_loaders.h"
#include "../include/dstree_index.h"
#include "../include/dstree_node.h"
#include "../include/ts.h"
#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#ifdef VALUES
#include <values.h>
#endif

int main(int argc, char **argv) {
  INIT_STATS()
  COUNT_TOTAL_TIME_START

  static char *dataset = "/home/karima/myDisk/data/Cgenerator/data_current.txt";
  static char *queries =
      "/home/karima/myDisk/data/Cgenerator/query_current.txt";
  static char *dataset_hists =
      "/home/karima/myDisk/data/Cgenerator/data_current_hists.txt";

  static char *index_path = "out/";
  static unsigned int dataset_size = 0;
  static unsigned int queries_size = 5;
  static unsigned int time_series_size = 256;
  static unsigned int init_segments = 1;
  static unsigned int leaf_size = 100;
  static double buffered_memory_size = 6439.2;
  static int use_ascii_input = 0;
  static int mode = 0;
  static float minimum_distance = FLT_MAX;
  static float epsilon = 0; // by default perform exact search
  static float delta = 1; // by default perform exact search, delta = 0 means ng
  static unsigned int k = 1;
  static unsigned int nprobes = 0;
  struct dstree_index *index = NULL;
  boolean is_index_new = 1;
  boolean track_bsf = 0;
  boolean track_pruning = 0;
  boolean all_mindists = 0;
  boolean max_policy = 0;
  unsigned char incremental = 0;
  unsigned char classify = 0;
  unsigned char track_file_pos = 0;
  float warping = 0;

  /* start kashif changes */
  static char *result_dir = "";
  static char *raw_data_dir = "";
  static unsigned int total_data_files = 100; // number of datasets to be indexed
  static unsigned int total_columns = 0; // number of columns to be indexed
  unsigned char track_vector = 0;
  unsigned char keyword_search = 0;
  unsigned char parallel = 0; // parallel incremental query answering
  unsigned char store_results_in_disk = 0;
  unsigned int qset_num = 3;
  unsigned int min_qset_size = 5;
  unsigned int max_qset_size = 10;
  unsigned int num_top = 3;           // number top sets to be returned
  static unsigned int data_gb_size = 0; // datalake size in GB
  unsigned int num_k_values = 0;
  char * k_values_str = "1,3,5,10,30,50,100,300,500,1000,3000,5000";
  char * ground_truth_dir = "";
  unsigned int num_threads = 1;
  unsigned int stop_when_nn_dist_changes = 0; // = 0 return all knn, = 1 stop when nn distance changes, = 2 (similar to 1) but return all results in the last increment.
  
  char knn_data_structure [3][12] = {"sorted-arr\0", "minmax-heap\0", "ostree\0"};

  unsigned int nn_struct = 0; //  = 0 use sorted array, = 1 use min max heap, = 2 use ostree

  /* end kashif changes */

  // printf("new code\n");

  while (1) {
    static struct option long_options[] = {
        {"ascii-input", required_argument, 0, 'a'},
        {"buffer-size", required_argument, 0, 'b'},
        {"epsilon", required_argument, 0, 'c'},
        {"dataset", required_argument, 0, 'd'},
        {"dataset-hists", required_argument, 0, 'n'},
        {"delta", required_argument, 0, 'e'},
        {"queries-size", required_argument, 0, 'f'},
        {"track-bsf", no_argument, 0, 'g'},
        {"track-pruning", required_argument, 0, 'i'},
        {"all-mindists", required_argument, 0, 'j'},
        {"max-policy", required_argument, 0, 'm'},
        {"queries", required_argument, 0, 'q'},
        {"index-path", required_argument, 0, 'p'},
        {"dataset-size", required_argument, 0, 'z'},
        {"k", required_argument, 0, 'k'},
        {"mode", required_argument, 0, 'x'},
        {"minimum-distance", required_argument, 0, 's'},
        {"timeseries-size", required_argument, 0, 't'},
        {"leaf-size", required_argument, 0, 'l'},
        {"nprobes", required_argument, 0, 'o'},
        {"incremental", no_argument, 0, 'h'},
        {"classify", no_argument, 0, 'r'},
        {"track_file_pos", no_argument, 0, 'v'},
        {"warping", required_argument, 0, 'u'},
        {"help", no_argument, 0, '?'},

        /* start kashif changes */
        {"result-dir", required_argument, 0, '>'},
        {"raw-data-dir", required_argument, 0, '^'},
        {"dataset-GB-size", required_argument, 0, ']'},
        {"total-data-files", required_argument, 0,
         '|'}, // number of datasets to be indexed (tables)
        {"nq", required_argument, 0, '%'}, // number of query sets
        {"min-qset-size", required_argument, 0, 'y'}, // minimum query set set size (number of vectors)
        {"max-qset-size", required_argument, 0, 'w'}, // maxmum query set set size (number of vectors)
        {"top", required_argument, 0, '#'}, // maxmum query set set size (number of vectors)
        {"track-vector", no_argument, 0, '*'},
        {"keyword-search", no_argument, 0, ':'},
        {"k-values", required_argument, 0, ';'},
        {"ground-truth-dir", required_argument, 0, '&'},
        {"parallel", no_argument, 0, '-'},
        {"store-results-in-disk", required_argument, 0, '!'},
        {"num-threads", required_argument, 0, ')'},
        {"stop-when-nn-dist-changes", required_argument, 0, '$'}, // stop knn search when nn distance changes
        {"knn-data-structure", required_argument, 0, '('}, // in which data structure should we store knns 

        
        /* end kashif changes */
    };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long(argc, argv, "", long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 'q':
      queries = optarg;
      break;

    case 's':
      minimum_distance = atof(optarg);
      break;

    case 'c':
      epsilon = atof(optarg);
      break;

    case 'e':
      delta = atof(optarg);
      break;

    case 'b':
      buffered_memory_size = atof(optarg);
      break;

    case 'f':
      queries_size = atoi(optarg);
      if (queries_size < 1) {
        fprintf(stderr,
                "Please change the queries size to be greater than 0.\n");
        exit(-1);
      }
      break;

    case 'k':
      k = atoi(optarg);
      if (k < 1) {
        fprintf(stderr, "Please change the k to be greater than 0.\n");
        exit(-1);
      }
      break;

    case 'g':
      track_bsf = 1;
      break;

    case 'i':
      track_pruning = atoi(optarg);
      break;

    case 'j':
      all_mindists = atoi(optarg);
      break;

    case 'm':
      max_policy = atoi(optarg);
      break;

    case 'd':
      dataset = optarg;
      break;

    case 'n':
      dataset_hists = optarg;
      break;

    case 'p':
      index_path = optarg;
      break;

    case 'x':
      mode = atoi(optarg);
      break;

    case 'z':
      dataset_size = atoi(optarg);
      if (dataset_size < 1) {
        fprintf(stderr,
                "Please change the dataset size to be greater than 0.\n");
        exit(-1);
      }
      break;

    case 't':
      time_series_size = atoi(optarg);
      break;

    case 'o':
      nprobes = atoi(optarg);
      break;

    case 'r':
      classify = 1;
      break;

    case 'v':
      track_file_pos = 1;
      break;

    case 'u':
      warping = atof(optarg);
      if (warping < 0 || warping > 1) {
        fprintf(stderr, "Please enter valid warping fraction \n");
        exit(-1);
      }
      if (warping > 0)
        fprintf(stderr, "DTW\n");
      else
        fprintf(stderr, "Euclidean\n");
      break;

    case 'l':
      leaf_size = atoi(optarg);
      if (leaf_size <= 1) {
        fprintf(stderr, "Please change the leaf size to be greater than 1.\n");
        exit(-1);
      }
      break;

    case '?':
      printf("Usage:\n\
                        \t--Queries and Datasets should be single precision\n\
                        \t--floating points. They can be binary of ascii.\n\
                        \t--However performance is faster with binary files\n\
                        \t--dataset XX \t\t\tThe path to the dataset file\n\
                        \t--queries XX \t\t\tThe path to the queries file\n\
                        \t--dataset-size XX \t\tThe number of time series to load\n\
                        \t--queries-size XX \t\tThe number of queries to run\n\
                        \t--mode: 0=index, 1=query, 2=index & query  3=calc_tlb\t\t\n\
                        \t--index-path XX \t\tThe path of the output folder\n\
                        \t--buffer-size XX \t\tThe size of the buffer memory in MB\n\
                        \t--timeseries-size XX\t\tThe size of each time series\n\
                        \t--ascii-input X \t\t\0 for ascii files and 1 for binary files\n\
                        \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                        \t--help\n\n\
                        \t--**********************EXAMPLES**********************\n\n\
                        \t--*********************INDEX MODE*********************\n\n\
                        \t--bin/dstree --dataset XX --dataset-size XX             \n\n\
                        \t--          --index-path XX --timeseries-size XX --mode 0\n\n\
                        \t--*********************QUERY MODE*********************\n\n\
                        \t--bin/dstree --queries XX --queries-size XX             \n\n\
                        \t--           --index-path XX --mode 1                 \n\n\
                        \t--*****************INDEX AND QUERY MODE***************\n\n\
                        \t--bin/dstree --dataset XX --dataset-size XX             \n\n\
                        \t--          --timeseries-size XX --index-path XX      \n\n\
                        \t--           --queries XX --queries-size XX --mode 2  \n\n\
                        \t--****************************************************\n\n");

      return 0;
      break;
    case 'a':
      use_ascii_input = atoi(optarg);
      break;
    case 'h':
      incremental = 1;
      break;

    /* start kashif changes */
    case '%':
      qset_num = atoi(optarg);
      break;

    case 'y':
      min_qset_size = atoi(optarg);
      break;

    case 'w':
      max_qset_size = atoi(optarg);
      break;

    case '#':
      num_top = atoi(optarg);
      break;

    case '>':
      result_dir = optarg;
      break;

    case '^':
      raw_data_dir = optarg;
      break;

    case '|':
      total_data_files = atoi(optarg);
      break;
    case '*':
      track_vector = 1;
      break;

    case ':':
      keyword_search = 1;
      break;

    case ']':
      data_gb_size = atoi(optarg);
      break;

    case ';':
      k_values_str = optarg;
      break;

    case '&':
      ground_truth_dir = optarg;
      break;

    case '-':
      parallel = 1;
      break;

    case '!':
      store_results_in_disk = 1;
      break;

    case ')':
      num_threads = atoi(optarg);
      break;

    case '$':
      stop_when_nn_dist_changes = atoi(optarg);
      break;

    case '(':
      if (strcmp(optarg, knn_data_structure[0]) == 0)
      {
        // printf("sorted array\n");
        nn_struct = 0;
      }
      else if ((strcmp(optarg, knn_data_structure[1]) == 0))
      {
        // printf("minmax-heap\n");
        nn_struct = 1;
      }
      else if ((strcmp(optarg, knn_data_structure[2]) == 0))
      {
        // printf("os-tree\n");
        nn_struct = 2;
      }
      else
      {
        printf("please specify a a valid knn data structure! (sorted-arr, ostree, minmax-heap)\n");
        nn_struct = 0;
        exit(1);
      }

      break;

    /* end kashif changes */
    default:
      exit(-1);
      break;
    }
  }
  
  /* start kashif changes */
  // if (dataset_size < 1) {
  //   fprintf(stderr,"Current datasize is less that 1!. Please change dataset directory.\n");
  //   exit(-1);
  // }

  if (dataset_size == 0) // get nb of vectors in data lake (if not provided by the user)
    dataset_size = (unsigned int) get_total_data_vectors(dataset, total_data_files, &total_columns); // get total number of vectors in data repository
  
  if (data_gb_size == 0) // get data lake size in gb (if not provided by the user)
  {
    data_gb_size = get_data_gb_size(dataset, total_data_files);
  }
  printf("Start Experiment...\n\nTotal files: %d\nTotal columns %d\nTotal  vectors:\t%d\nData lake size in GB:\t%u\n\n\n", total_data_files, total_columns, dataset_size, data_gb_size);
  // printf("buffered memory = %f MB\n", buffered_memory_size);

  /* end kashif changes */

  

  minimum_distance = FLT_MAX;
  if (mode == 0) // only build and store the index
  {
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
    struct dstree_index_settings *index_settings = dstree_index_settings_init(
        index_path, time_series_size, init_segments, leaf_size,
        buffered_memory_size, is_index_new, dataset);

    if (index_settings == NULL) {
      fprintf(stderr,
              "Error main.c:  Could not initialize the index settings.\n");
      return -1;
    }
    
    index = dstree_index_init(index_settings);
    index->settings->track_vector = track_vector;

    index->first_node = dstree_root_node_init(index->settings);
    if (index == NULL) {
      fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
      return -1;
    }


    if (!use_ascii_input)
    {
      /* start kashif changes */
      if (classify)
      {        
        printf("Error in dstree.c: Function dstree_index_classify_multiple_binary_file() is not implemented!\n");
        exit(1);
      } 
      else 
      {
        dstree_index_multiple_binary_files(dataset, total_data_files, index);
        // dstree_index_binary_file(dataset, dataset_size, index);
      }
      /* end kashif changes */

      COUNT_PARTIAL_TIME_END
      index->stats->idx_building_total_time +=
          partial_time + index->stats->idx_traverse_tree_total_time +
          index->stats->idx_append_ts_to_leaf_total_time +
          index->stats->idx_evaluate_split_policies_total_time +
          index->stats->idx_split_node_total_time;
      index->stats->idx_building_input_time +=
          partial_input_time + index->stats->idx_traverse_tree_input_time +
          index->stats->idx_append_ts_to_leaf_input_time +
          index->stats->idx_evaluate_split_policies_input_time +
          index->stats->idx_split_node_input_time;
      index->stats->idx_building_output_time +=
          partial_output_time + index->stats->idx_traverse_tree_output_time +
          index->stats->idx_append_ts_to_leaf_output_time +
          index->stats->idx_evaluate_split_policies_output_time +
          index->stats->idx_split_node_output_time;
      index->stats->idx_building_cpu_time +=
          index->stats->idx_building_total_time -
          index->stats->idx_building_input_time -
          index->stats->idx_building_output_time;
      index->stats->idx_building_seq_input_count +=
          partial_seq_input_count +
          index->stats->idx_traverse_tree_seq_input_count +
          index->stats->idx_append_ts_to_leaf_seq_input_count +
          index->stats->idx_evaluate_split_policies_seq_input_count +
          index->stats->idx_split_node_seq_input_count;
      index->stats->idx_building_seq_output_count +=
          partial_seq_output_count +
          index->stats->idx_traverse_tree_seq_output_count +
          index->stats->idx_append_ts_to_leaf_seq_output_count +
          index->stats->idx_evaluate_split_policies_seq_output_count +
          index->stats->idx_split_node_seq_output_count;
      index->stats->idx_building_rand_input_count +=
          partial_rand_input_count +
          index->stats->idx_traverse_tree_rand_input_count +
          index->stats->idx_append_ts_to_leaf_rand_input_count +
          index->stats->idx_evaluate_split_policies_rand_input_count +
          index->stats->idx_split_node_rand_input_count;

      index->stats->idx_building_rand_output_count +=
          partial_rand_output_count +
          index->stats->idx_traverse_tree_rand_output_count +
          index->stats->idx_append_ts_to_leaf_rand_output_count +
          index->stats->idx_evaluate_split_policies_rand_output_count +
          index->stats->idx_split_node_rand_output_count;

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START
      if (!dstree_index_write(index)) {
        fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
        return -1;
      }
      COUNT_PARTIAL_TIME_END
      // COUNT_TOTAL_TIME_END
      index->stats->idx_writing_total_time = partial_time;
      index->stats->idx_writing_input_time = partial_input_time;
      index->stats->idx_writing_output_time = partial_output_time;
      index->stats->idx_writing_cpu_time =
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_writing_seq_input_count = partial_seq_input_count;
      index->stats->idx_writing_seq_output_count = partial_seq_output_count;
      index->stats->idx_writing_rand_input_count = partial_rand_input_count;
      index->stats->idx_writing_rand_output_count = partial_rand_output_count;

      dstree_get_index_stats(index);
      dstree_print_index_stats(index, dataset);
      // COUNT_TOTAL_TIME_START
    } 
    else 
    {
      /* start kashif changes */
      printf("Error in dstree.c: Function dstree_index_multiple_ascii_files() is not implemented!\n");
      exit(1);
      /* end kashif changes */
    }

  }
  else if (mode == 1) // read an existing index and execute queries
  {
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
    is_index_new = 0;
    index = dstree_index_read(index_path);
    fprintf(stderr, ">>> Index read successfully\n");
    index->settings->dataset = dataset;
    if (classify != 1) {
      index->settings->classify = 0;
      // if (index->gt_filename != NULL)
      //   free(index->gt_filename);
      // if (index->gt_cache != NULL)
      //   free(index->gt_cache);
    }
    if (!track_file_pos){
      index->settings->track_file_pos = 0;
      // if (index->fp_filename != NULL)
      //   free(index->fp_filename);
      // if (index->fp_cache != NULL)
      //   free(index->fp_cache);
    }
    /* start kashif changes */
    if (!track_vector)
    {
      index->settings->track_vector = 0;
      // if (index->vid_filename != NULL)
      //   free(index->vid_filename);
      // if (index->vid_cache != NULL)
      //   free(index->vid_cache);
    }
    /* end kashif changes */
    if (index == NULL) {
      fprintf(stderr, "Error main.c:  Could not read the index from disk.\n");
      return -1;
    }
    COUNT_PARTIAL_TIME_END
    // COUNT_TOTAL_TIME_END
    index->stats->idx_reading_total_time = partial_time;
    index->stats->idx_reading_input_time = partial_input_time;
    index->stats->idx_reading_output_time = partial_output_time;
    index->stats->idx_reading_cpu_time =
        partial_time - partial_input_time - partial_output_time;
    index->stats->idx_reading_seq_input_count = partial_seq_input_count;
    index->stats->idx_reading_seq_output_count = partial_seq_output_count;
    index->stats->idx_reading_rand_input_count = partial_rand_input_count;
    index->stats->idx_reading_rand_output_count = partial_rand_output_count;

    dstree_get_index_stats(index);
    // dstree_print_index_stats(index, dataset);
    // COUNT_TOTAL_TIME_START

    fprintf(stderr, ">>> Index loaded successfully from: %s\n", index_path);

    // calculate r_delta
    ts_type r_delta = FLT_MAX;

    if (delta > 0 && delta < 1) {
      FILE *fp;
      int bins;
      double *x, *y;
      int i, j;
      if ((fp = fopen(dataset_hists, "r")) != NULL) {
        fscanf(fp, "%d", &bins);

        x = (double *)calloc(bins, sizeof(double));
        y = (double *)calloc(bins, sizeof(double));

        for (int j = 0; j < bins; j++)
          fscanf(fp, "%lf%lf", &x[j], &y[j]);
        fclose(fp);
      } else
        printf("Error opening data distribution file.\n");

      for (j = 0; (j < bins) && (y[j] < delta); j++)
        ;
      j--;
      fprintf(stderr, "Using histogram bin %lf.\n", y[j]);

      r_delta =
          x[j] + (((delta - y[j]) * (x[j + 1] - x[j])) / (y[j + 1] - y[j]));
      fprintf(stderr, "Using r_delta = %lf.\n", r_delta);
    }

    index->settings->track_vector = track_vector;
    index->settings->parallel = parallel;
    /* start kashif changes */

    if(incremental && parallel)
    {
      printf("Parallel QA  (%d threads) ...\n", num_threads);
      dstree_multi_thread_variable_num_thread_parallel_incr_knn_query_multiple_binary_files(queries, qset_num, min_qset_size, max_qset_size, num_top, index,
                                    minimum_distance, epsilon, r_delta,k, track_bsf, track_pruning, all_mindists,
                                    max_policy, nprobes, incremental, result_dir, total_data_files, data_gb_size, 
                                    warping, keyword_search, k_values_str, ground_truth_dir, store_results_in_disk, num_threads, stop_when_nn_dist_changes, nn_struct);
    }
    else
    {
      printf("Sequential QA ...\n");
      dstree_knn_query_multiple_binary_files(queries, qset_num, min_qset_size, max_qset_size, num_top, index,
                                    minimum_distance, epsilon, r_delta,k, track_bsf, track_pruning, all_mindists,
                                    max_policy, nprobes, incremental, result_dir, total_data_files, data_gb_size, 
                                    warping, keyword_search, k_values_str, ground_truth_dir);
    }
    
    /* end kashif changes */
  } 
  else if (mode == 2) // build the index, execute queries and store the index
  {
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
    struct dstree_index_settings *index_settings = dstree_index_settings_init(
        index_path, time_series_size, init_segments, leaf_size,
        buffered_memory_size, is_index_new, dataset);

    if (index_settings == NULL) {
      fprintf(stderr,
              "Error main.c:  Could not initialize the index settings.\n");
      return -1;
    }

    index = dstree_index_init(index_settings);
    index->first_node = dstree_root_node_init(index->settings);
    if (index == NULL) {
      fprintf(stderr, "Error main.c:  Could not initialize the index.\n");
      return -1;
    }
    if (!use_ascii_input)
    {
      index->settings->track_vector = track_vector;

      /* start kashif changes */
      if (!dstree_index_multiple_binary_files(dataset, total_data_files, index))
      {
        fprintf(stderr, "Error main.c:  Could not build the index.\n");
        return -1;
      }
      /* end kashif changes */
      COUNT_PARTIAL_TIME_END
      // COUNT_TOTAL_TIME_END
      index->stats->idx_building_total_time = partial_time;
      index->stats->idx_building_input_time = partial_input_time;
      index->stats->idx_building_output_time = partial_output_time;
      index->stats->idx_building_cpu_time =
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_building_seq_input_count = partial_seq_input_count;
      index->stats->idx_building_seq_output_count = partial_seq_output_count;
      index->stats->idx_building_rand_input_count = partial_rand_input_count;
      index->stats->idx_building_rand_output_count = partial_rand_output_count;

      dstree_get_index_stats(index);
      dstree_print_index_stats(index, dataset);
      COUNT_TOTAL_TIME_START

      if (!dstree_knn_query_multiple_binary_files(queries, qset_num,
                                    min_qset_size, 
                                    max_qset_size, num_top, index,
                                    minimum_distance, epsilon, delta,
                                    k, track_bsf, track_pruning, all_mindists,
                                    max_policy, nprobes, incremental, result_dir, 
                                    total_data_files, data_gb_size, warping, keyword_search, k_values_str, ground_truth_dir)) {
        fprintf(stderr, "Error main.c:  Could not execute the query.\n");
        return -1;
      }

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START
      if (!dstree_index_write(index)) {
        fprintf(stderr, "Error main.c:  Could not save the index to disk.\n");
        return -1;
      }
      // COUNT_PARTIAL_TIME_END
      index->stats->idx_writing_total_time = partial_time;
      index->stats->idx_writing_input_time = partial_input_time;
      index->stats->idx_writing_output_time = partial_output_time;
      index->stats->idx_writing_cpu_time =
          partial_time - partial_input_time - partial_output_time;
      index->stats->idx_writing_seq_input_count = partial_seq_input_count;
      index->stats->idx_writing_seq_output_count = partial_seq_output_count;
      index->stats->idx_writing_rand_input_count = partial_rand_input_count;
      index->stats->idx_writing_rand_output_count = partial_rand_output_count;
      // ADD TIME COUNTERS
    }
    else
    {
      /* start kashif changes */
      fprintf(stderr, "Error in dstree.c: Kashif cannot take ascii input.\n");
      exit(1);
      /* end kashif changes */
    }
  }
  else if (mode == 3) // read an existing index and execute queries
  {
    /* start kashif changes */
      fprintf(stderr, "Error in dstree.c: Mode 3 not implemented for Kashif.\n");
      exit(1);
      /* end kashif changes */
  }
  else 
  {
    fprintf(
        stderr,
        "Please use a valid mode. run dstree --help for more information. \n");
    return -1;
  }

  COUNT_TOTAL_TIME_END
  fprintf(stderr,
          "Sanity check: combined indexing and querying times should be less "
          "than: %f secs \n",
          total_time / 1000000);

  dstree_index_destroy(index, index->first_node, is_index_new);

  if(index->stats->leaves_heights)
    free(index->stats->leaves_heights);
  if(index->stats->leaves_sizes)
    free(index->stats->leaves_sizes);
  free(index->stats);
  free(index->settings);
  free(index);

  // printf("\n");

  // malloc_stats_print(NULL, NULL, NULL);
  return 0;
}


