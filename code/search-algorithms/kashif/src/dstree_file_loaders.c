//
//  dstree_file_loaders.c
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//
#define _XOPEN_SOURCE 600 /* to avoid error with pthread_barried_t undefined */

#include "../include/dstree_file_loaders.h"
#include "../config.h"
#include "../globals.h"
#include "../include/dstree_query_engine.h"
#include "../include/kashif_utils.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

// CORRECT STATS FOR ASCII
enum response dstree_query_ascii_file(const char *ifilename, int q_num,
                                      const char delimiter,
                                      struct dstree_index *index,
                                      float minimum_distance, ts_type epsilon,
                                      ts_type r_delta) {
  FILE *ifile;
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "r");
  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n",
            ifilename);
    return FAILURE;
  }

  char *ts_str = NULL;
  size_t linecap = 0;
  ssize_t linelen;
  unsigned int q_loaded = 0;
  ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
  if (ts == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Querying..\
                         Could not allocate memory for time series!\n");
    return FAILURE;
  }

  COUNT_PARTIAL_INPUT_TIME_START
  while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 &&
         q_loaded < q_num) {
    COUNT_PARTIAL_INPUT_TIME_END
    COUNT_PARTIAL_SEQ_INPUT
    if (ts_str == NULL) {
      fprintf(stderr, "Error in dstree_file_loaders.c: Querying..\
                         Could not get the time series from file %s.\n",
              ifilename);
      return FAILURE;
    }
    if (!ts_parse_str(ts_str, ts, index->settings->timeseries_size,
                      &delimiter)) {
      fprintf(stderr, "Error in dstree_file_loaders.c:  Querying..Could not "
                      "parse the time series.\n");
      return FAILURE;
    }

    printf("\n\n");

    q_loaded++;
    COUNT_PARTIAL_INPUT_TIME_START
  }

  if (fclose(ifile)) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not close the query filename %s",
        ifilename);
    return FAILURE;
  }
  COUNT_PARTIAL_INPUT_TIME_END

  free(ts);
  free(ts_str);
  return SUCCESS;
}

enum response dstree_query_binary_file(const char *ifilename, int q_num,
                                       struct dstree_index *index,
                                       float minimum_distance, ts_type epsilon,
                                       ts_type r_delta) {
  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  FILE *ifile;
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "rb");
  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "File %s not found!\n", ifilename);
    exit(-1);
  }

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END
  unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz / ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    exit(-1);
  }

  unsigned int q_loaded = 0;
  ts_type *query_ts = calloc(1, sizeof(ts_type) * ts_length);
  ts_type *query_ts_reordered = calloc(1, sizeof(ts_type) * ts_length);
  int *query_order = calloc(1, sizeof(int) * ts_length);
  if (query_order == NULL)
    return FAILURE;

  while (q_loaded < q_num) {

    RESET_QUERY_COUNTERS()

    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(query_ts, sizeof(ts_type), ts_length, ifile);
    COUNT_PARTIAL_INPUT_TIME_END

    reorder_query(query_ts, query_ts_reordered, query_order, ts_length);

    struct query_result result =
        exact_search(query_ts, query_ts_reordered, query_order, offset, index,
                     minimum_distance, epsilon, r_delta);

    q_loaded++;

    get_query_stats(index, 1);
    print_query_stats(index, q_loaded, 1, ifilename);
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  COUNT_PARTIAL_TIME_END
  RESET_PARTIAL_COUNTERS()

  free(query_ts);
  free(query_ts_reordered);
  free(query_order);

  if (fclose(ifile)) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not close the query filename %s",
        ifilename);
    return FAILURE;
  }

  return SUCCESS;
}

enum response dstree_knn_query_binary_file(
    const char *ifilename, int q_num, struct dstree_index *index,
    float minimum_distance, ts_type epsilon, ts_type r_delta, unsigned int k,
    boolean track_bsf, boolean track_pruning, boolean all_mindists,
    boolean max_policy, unsigned int nprobes, unsigned char incremental,
    float warping) {

  struct bsf_snapshot **bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
  if (track_bsf) {
    max_bsf_snapshots = 10000;
    cur_bsf_snapshot = 0;

    bsf_snapshots = calloc(k, sizeof(struct bsf_snapshot *));
    for (unsigned int i = 0; i < k; ++i) {
      bsf_snapshots[i] = calloc(max_bsf_snapshots, sizeof(struct bsf_snapshot));
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
        bsf_snapshots[i][j].distance = FLT_MAX;
        bsf_snapshots[i][j].time = FLT_MAX;
        bsf_snapshots[i][j].series = NULL;
        bsf_snapshots[i][j].checked_nodes = -1;
        bsf_snapshots[i][j].label = 0;
      }
    }
  }

  RESET_PARTIAL_COUNTERS()

  COUNT_PARTIAL_TIME_START

  FILE *ifile;
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "rb");
  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "File %s not found!\n", ifilename);
    exit(-1);
  }

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END
  unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz / ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    exit(-1);
  }

  unsigned int q_loaded = 0;
  ts_type *query_ts = calloc(1, sizeof(ts_type) * ts_length);
  ts_type *query_ts_reordered = calloc(1, sizeof(ts_type) * ts_length);
  int *query_order = calloc(1, sizeof(int) * ts_length);
  if (query_order == NULL)
    return FAILURE;

  /*
      const char *filename = malloc(sizeof(char) *
    (strlen(index->settings->root_directory) + 18)); filename = strcpy(filename,
    index->settings->root_directory); filename = strcat(filename,
    "raw_series.csv\0");

    printf ("series_file = %s\n", filename);
    printf ("dataset_file = %s\n", index->settings->dataset);

     FILE *series_file = fopen(filename, "a");
     FILE *dataset_file = fopen(index->settings->dataset, "rb");
 */
  FILE *series_file = NULL;
  FILE *dataset_file = NULL;
  while (q_loaded < q_num) {

    RESET_QUERY_COUNTERS()

    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(query_ts, sizeof(ts_type), ts_length, ifile);
    COUNT_PARTIAL_INPUT_TIME_END

    reorder_query(query_ts, query_ts_reordered, query_order, ts_length);

    q_loaded++;

    if (track_bsf) {
      cur_bsf_snapshot = 0;
      if (incremental) {
        exact_de_incr_progressive_knn_search(
            query_ts, query_ts_reordered, query_order, offset, index,
            minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,
            bsf_snapshots, &cur_bsf_snapshot, warping, dataset_file,
            series_file);
      } else {
        exact_de_progressive_knn_search(
            query_ts, query_ts_reordered, query_order, offset, index,
            minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,
            bsf_snapshots, &cur_bsf_snapshot);
      }
      for (unsigned int i = 0; i < k; ++i) {
        for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
          bsf_snapshots[i][j].distance = FLT_MAX;
          bsf_snapshots[i][j].time = FLT_MAX;
          bsf_snapshots[i][j].series = NULL;
          bsf_snapshots[i][j].checked_nodes = -1;
          bsf_snapshots[i][j].label = 0;
          bsf_snapshots[i][j].vector_id->table_id = -1;
          bsf_snapshots[i][j].vector_id->set_id = -1;
        }
      }
    } else {
      exact_de_knn_search(query_ts, query_ts_reordered, query_order, offset,
                          index, minimum_distance, epsilon, r_delta, k,
                          q_loaded, ifilename);
    }

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  COUNT_PARTIAL_TIME_END
  RESET_PARTIAL_COUNTERS()

  free(query_ts);
  free(query_ts_reordered);
  free(query_order);

  if (fclose(ifile)) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not close the query filename %s",
        ifilename);
    return FAILURE;
  }

  if (track_bsf) {
    for (unsigned int i = 0; i < k; ++i) {
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
        if (bsf_snapshots[i][j].series != NULL)
          free(bsf_snapshots[i][j].series);
      }
      free(bsf_snapshots[i]);
    }
    free(bsf_snapshots);
  }
  // fclose(series_file);
  // fclose(dataset_file);

  return SUCCESS;
}

/* start kashif changes */
// new function query from multiple binary files
enum response dstree_knn_query_multiple_binary_files(
    const char *bin_files_directory, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top,
    struct dstree_index *index, float minimum_distance, ts_type epsilon,
    ts_type r_delta, unsigned int k, boolean track_bsf, boolean track_pruning,
    boolean all_mindists, boolean max_policy, unsigned int nprobes,
    unsigned char incremental, char *result_dir, unsigned int total_data_files,
    unsigned int dlsize, // total disk size of data files indexed in dstree
    float warping, unsigned char keyword_search, char * k_values_str, char * ground_truth_dir) {

  // struct bsf_snapshot **bsf_snapshots = NULL;
  // unsigned int max_bsf_snapshots;
  // unsigned int cur_bsf_snapshot;
  FILE *series_file = NULL;
  FILE *dataset_file = NULL;

  unsigned int * k_values = NULL;
  unsigned int num_k_values = 0;

  // extract k values from string "1,3,5,10" to [1, 3, 5, 10]
  k_values = get_k_values(k_values_str, &num_k_values);
  // printf("num k values = %u\n", num_k_values);

  if (k_values == NULL)
  {
    fprintf(stderr,
              "Error dstree.c:  Could not read set of k values.\n");
      return -1;
  }
  for(int u = 0; u < num_k_values; u++)
      printf(" k = %u\n", k_values[u]);

  int vector_length = index->settings->timeseries_size;
  int opened_files = 0, qvectors_loaded = 0;
  
  unsigned int offset = 0;

  // open source dir
  struct dirent *dfile;
  DIR *dir = opendir(bin_files_directory);

  if (!dir) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Unable to open directory stream! %s", bin_files_directory);
    exit(1);
  }

  // allocate memory for vector
  struct vector query_vector;
  query_vector.values = (ts_type *)malloc(sizeof(ts_type) * vector_length);

  ts_type *query_vector_reordered = calloc(1, sizeof(ts_type) * vector_length);
  int *query_order = calloc(1, sizeof(int) * vector_length);

  // query time
  double query_time = 0.0;
  unsigned int total_checked_ts = 0;
  unsigned int total_queries = qset_num;
  bool found_query = false; // throw error if no query set was found

  if (query_order == NULL)
    return FAILURE;

  /* Start experiment  - create experiment result directory */
  char *results_dir = make_result_directory(
      result_dir, total_data_files, qset_num, min_qset_size, max_qset_size);
  

  // initialize list of all knn results (from all query vectors in query set)
  struct result_vid **all_knn_results = NULL;
  struct query_result *curr_knn = NULL;
  struct vid *top_matches;
  struct vid query_id;

  // max k values for which we must save results
  unsigned int max_k = k_values[num_k_values - 1];

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  // for every file (table)
  while ((dfile = readdir(dir)) != NULL) {
    // skip directories
    // if (dfile->d_type != DT_REG)
    //     continue;
    if (qset_num == 0)
      break;

    // if file is binary file
    if (is_binaryfile(dfile->d_name)) {
      opened_files += 1;

      // get fill path of bin file
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, bin_files_directory);
      strcat(bin_file_path, "/");
      strcat(bin_file_path, dfile->d_name);

      // get binary table info
      int datasize, table_id, nsets, vector_length_in_filename;
      sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize,
             &table_id, &nsets, &vector_length_in_filename);

      // check if vector length in file name matches vector length passed as
      // argument
      if (vector_length_in_filename != vector_length) {
        fprintf(stderr,
                "Error in dstree_file_loaders.c:  Vector length passed in "
                "argumentes (--timeseries-size %d) does not match vector "
                "length in file (%d) %s.\n",
                vector_length_in_filename, vector_length, bin_file_path);
        return FAILURE;
      }

      /* read binary file */
      COUNT_PARTIAL_RAND_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      FILE *bin_file = fopen(bin_file_path, "rb");
      COUNT_PARTIAL_INPUT_TIME_END

      if (bin_file == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
                bin_file_path);
        return FAILURE;
      }

      /* Start processing file: read every vector in binary file*/
      int i = 0, j = 0, set_id = 0,
          total_bytes = (datasize * vector_length) + nsets;
      unsigned int nvec = 0u;
      ts_type val;

      while (total_bytes) // counts 4 bytes as one because every vector value is
                          // stored in 4 bytes
      {
        // beginning of a set of vectors
        if (i == 0) {
          if (qset_num == 0)
            break;

          // read first integer to check how many vactors in current set
          COUNT_PARTIAL_INPUT_TIME_START
          fread(&nvec, sizeof(nvec), 1, bin_file);
          COUNT_PARTIAL_INPUT_TIME_END

          total_bytes--;
          // query set does not fit requirments move to next set
          if ((unsigned int)nvec < min_qset_size ||
              (unsigned int)nvec > max_qset_size) {
            printf("skipping (%u, %u) size = %d\n", table_id, set_id, nvec);
            COUNT_PARTIAL_INPUT_TIME_START
            fseek(bin_file, nvec * 4 * vector_length, SEEK_CUR);
            COUNT_PARTIAL_INPUT_TIME_END
            i = 0;
            j = 0;
            total_bytes -= (nvec * vector_length);
            nvec = 0u;
            set_id += 1;
            
            continue;
          }

          found_query = true;

          query_time = 0.0;
          total_checked_ts = 0;
          qset_num--;

          all_knn_results = malloc(nvec * sizeof(struct result_vid *));
          for(int q = 0; q < nvec; q++)
          {
            all_knn_results[q] = calloc(k, sizeof(struct result_vid));
          }
          RESET_QUERY_COUNTERS()

          // set new set id
          query_vector.table_id = table_id;
          query_vector.set_id = set_id;
          query_vector.pos = 0;

          printf("\n+ Query (%d, %d) ...\n", table_id, set_id);

          set_id += 1;
          i++;
          j = 0;
        } else if (i <= (unsigned int)nvec * vector_length) {
          // end of vector but still in current set
          if (j > (vector_length - 1)) {
            j = 0;
            
            /*run query vector in dstree */
            reorder_query(query_vector.values, query_vector_reordered,
                          query_order, vector_length);
            qvectors_loaded += 1;

            // perform extact knn search
            // with incremental answering
            if (track_bsf) {
              query_id.table_id = query_vector.table_id;
              query_id.set_id = query_vector.set_id;
              query_id.pos = query_vector.pos;

              double query_vector_time = 0.0;
              unsigned int num_checked_vectors = 0;

              if (incremental) {
                curr_knn = exact_de_incr_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_vector_time,
                    &num_checked_vectors,
                    warping, dataset_file, series_file, &query_id, nvec, k_values, num_k_values);
              } 
              query_time += query_vector_time;
              total_checked_ts += num_checked_vectors;
            }
            // copy new knn(s) to knn_results array
            // printf("end of query for vector sent %u received %u\n", query_vector.pos, curr_knn[0].vector_id->pos);
            // printf("Retreiving %dnn results...\n", k);
            for (int t = 0; t < k; t++) 
            {
              all_knn_results[query_vector.pos][t].table_id = curr_knn[t].vector_id->table_id;
              all_knn_results[query_vector.pos][t].set_id = curr_knn[t].vector_id->set_id;
              all_knn_results[query_vector.pos][t].pos = curr_knn[t].vector_id->pos;
              all_knn_results[query_vector.pos][t].qpos = curr_knn[t].query_vector_pos;
              all_knn_results[query_vector.pos][t].time = curr_knn[t].time;
              all_knn_results[query_vector.pos][t].distance = curr_knn[t].distance;
            }
            query_vector.pos += 1;
            
            for (int t = 0; t < k; t++)
              free(curr_knn[t].vector_id);
            free(curr_knn);
          }

          COUNT_PARTIAL_SEQ_INPUT
          COUNT_PARTIAL_INPUT_TIME_START
          fread((void *)(&val), sizeof(val), 1, bin_file);
          COUNT_PARTIAL_INPUT_TIME_END

          total_bytes--;
          query_vector.values[j] = val;

          // end of last vector in current  set
          if (i == (unsigned int)nvec * vector_length) {
            /*run query vector in dstree */
            reorder_query(query_vector.values, query_vector_reordered,
                          query_order, vector_length);
            qvectors_loaded += 1;

            // perform extact knn search
            // with incremental answering
            if (track_bsf) {
              // cur_bsf_snapshot = 0;

              query_id.table_id = query_vector.table_id;
              query_id.set_id = query_vector.set_id;
              query_id.pos = query_vector.pos;
              
               double query_vector_time = 0.0;
              unsigned int num_checked_vectors = 0;

              if (incremental) {
                curr_knn = exact_de_incr_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_vector_time,
                    &num_checked_vectors,
                    warping, dataset_file, series_file, &query_id, nvec, k_values, num_k_values);
              } 
              query_time += query_vector_time;
              total_checked_ts += num_checked_vectors;
            }
            
            // append new knn(s) to knn_results array
            // printf("end of query for vector sent %u received %u\n", query_vector.pos, curr_knn[0].vector_id->pos);
            // printf("Retreiving %dnn results...\n", k);
            for (int t = 0; t < k; t++) {
              all_knn_results[query_vector.pos][t].table_id = curr_knn[t].vector_id->table_id;
              all_knn_results[query_vector.pos][t].set_id = curr_knn[t].vector_id->set_id;
              all_knn_results[query_vector.pos][t].pos = curr_knn[t].vector_id->pos;
              all_knn_results[query_vector.pos][t].qpos = curr_knn[t].query_vector_pos;
              all_knn_results[query_vector.pos][t].time = curr_knn[t].time;
              all_knn_results[query_vector.pos][t].distance = curr_knn[t].distance;
            }
            query_vector.pos = 0;

            for (int t = 0; t < k; t++)
              free(curr_knn[t].vector_id);
            free(curr_knn);
            
            // End of Query set
            /* Save query results to csv file */
            query_time /= 1000000;
            
            printf("Storing result to csv file...\n");

            for(int z = 0; z < num_k_values; z++)
            {
              unsigned int curr_k = k_values[z];
              
              char *query_result_file = make_file_path(results_dir, table_id, query_vector.set_id, nvec,
                                      total_data_files, dlsize, vector_length, curr_k);

              if(!save_to_query_result_file(query_result_file, table_id, query_vector.set_id,
                                      nvec, all_knn_results, curr_k))
              {
                fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't save query results to file %s.", query_result_file);
                exit(1);
              }
              free(query_result_file);
            }
            printf("end.\n");
            
            // printf("\nquery_time=%fsec\n", query_time);
            // printf("computing recall...\n");

            // compute recall
            // int k_idx;
            // unsigned int k = k_values[num_k_values - 1];
            // struct performance perf;
            // char ground_truth_file[255] = "";
            // int  num_gt_results = get_ground_truth_file(ground_truth_dir, query_vector.table_id, query_vector.set_id, ground_truth_file);
            // struct result_vid * ground_truth_results = get_ground_truth_results(ground_truth_file, num_gt_results);

            // printf("gt results: #resutls = %d\n", num_gt_results);
            // for(int gt = 0; gt < num_gt_results; gt++)
            // {
            //   printf("q = %u, s = (%u, %u, %u), d = %f\n", 
            //   ground_truth_results[gt].qpos, ground_truth_results[gt].table_id, 
            //   ground_truth_results[gt].set_id, ground_truth_results[gt].pos, ground_truth_results[gt].distance);
            // }
            // exit(1);

            // for(int q = 0; q < nvec; q++)
            // {
            //   printf("\n q = (%u, %u, %u) results:\n", query_vector.table_id, query_vector.set_id, q);
            //   for(int i = 0; i < k; i++)
            //   {
            //     perf = compute_recall_precision(ground_truth_results, num_gt_results, all_knn_results, i + 1, q, query_vector.table_id, query_vector.set_id);
            //     printf("k = %d,\tdistance = %f,\trecall = %f,\ttime = %f\n", i+1, all_knn_results[q][i].distance, perf.recall, all_knn_results[q][i].time/1000000);
            //   }

            //   exit(1);
              // for(int i = 0; i < num_k_values; i++)
              // {
              //   k_idx = k_values[i]-1;

              //   perf = compute_recall_precision(ground_truth_results, num_gt_results, all_knn_results, k_idx + 1, q, query_vector.table_id, query_vector.set_id);
              //   printf("k = %d,\tdistance = %f,\trecall = %f,\ttime = %f\n", k_idx+1, all_knn_results[q][k_idx].distance, perf.recall, all_knn_results[q][k_idx].time/1000000);
              // }

            //   printf("end results.\n\n");
            // }

            // free(ground_truth_results);
            // temp change
            // int k_values[] = {1, 3, 5, 7, 10, 30, 50, 70, 100, 300, 500, 700, 1000, 10000};
            // int num_k = 17;
            
            // int k_values[] = {1, 3, 5, 7, 10, 30, 50, 70, 100, 300, 500, 700, 1000};
            // int num_k = 13;
            // int curr_k;


            // printf("\nnvec : %d\n", nvec);
            // printf("k,\trecall,\tprecision,\ttime\n");

            // float time = 0.0;
            // struct performance perf;
            // for(int a = 0; a < num_k; a++)
            // {
            //   curr_k = k_values[a];
            //   time = 0.0;
            //   perf = compute_recall_precision(ground_truth_dir, all_knn_results, nvec, curr_k, query_vector.table_id, query_vector.set_id);
            //   for(int q = 0; q < nvec; q++)
            //     time += all_knn_results[q][curr_k-1].time;
            //   printf("%d,\t%f,\t%f,\t%f\n", curr_k, perf.recall, perf.precision, time/1000000);

            // }

            // if(keyword_search)
            // {
            //   // don't change these lines to allaow ui to fetch results
            //     struct result_table* top = get_top_tables_by_euclidean_distance(all_knn_results, knn_array_idx, num_top);
            //     for(int m = 0; m < num_top; m++)
            //     {
            //       printf("table-%u- in file @@%s$ min_distance=%.3f§ num_closest=%u# total_matches=%uµ\n", top[m].table_id, top[m].raw_data_file, top[m].min_distance, top[m].num_min, top[m].total_matches);
            //     }
            //     printf("\nquery_time=%fsec\n", query_time);
            //     free(top);
            //     // don't change these lines to allaow ui to fetch results
            // } 
            // else
            // {
            //     // don't change these lines to allaow ui to fetch results
            //     struct result_sid * top = get_top_sets(all_knn_results, knn_array_idx, num_top);
            //     for(int m = 0; m < num_top; m++)
            //     {
            //       printf("table-%u-column-%u- in file @@%s$ overlap=%u§\n", top[m].table_id, top[m].set_id, top[m].raw_data_file, top[m].overlap_size);
            //     }
            //     printf("\nquery_time=%fsec\n", query_time);
            //     free(top);
            //     // don't change these lines to allaow ui to fetch results
            // }

            // free memory
            for(int q = 0; q < nvec; q++)
            {
              free(all_knn_results[q]);
            }
            free(all_knn_results);
            
            RESET_PARTIAL_COUNTERS()
            COUNT_PARTIAL_TIME_START

            i = 0;
            j = 0;
            nvec = 0u;
            continue;
          }
          i++;
          j++;
        }
      }
      COUNT_PARTIAL_INPUT_TIME_START
      fclose(bin_file);
      COUNT_PARTIAL_INPUT_TIME_END
      
      /* end read binary file */
    }
  }

  if (found_query == false) {
    fprintf(stderr,
            "WARNING! Could not find any query in size range [%u - %u].\n",
            min_qset_size, max_qset_size);
    exit(-1);
  }
  if (opened_files == 0) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c:  Could not find any binary file "
            "in directory %s.\n",
            bin_files_directory);
    return FAILURE;
  }
  

  // free memory
  COUNT_INPUT_TIME_START
  closedir(dir);
  COUNT_INPUT_TIME_END
  
  COUNT_PARTIAL_TIME_END
  RESET_PARTIAL_COUNTERS()

  free(k_values);
  free(query_vector.values);
  free(query_vector_reordered);
  free(query_order);
  free(results_dir);
  
  return SUCCESS;
}

// parallel incremental query answering read the whole query column and submit all query vectors to teh search engine
enum response dstree_multi_thread_variable_num_thread_parallel_incr_knn_query_multiple_binary_files(
    const char *bin_files_directory, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top,
    struct dstree_index *index, float minimum_distance, ts_type epsilon,
    ts_type r_delta, unsigned int k, boolean track_bsf, boolean track_pruning,
    boolean all_mindists, boolean max_policy, unsigned int nprobes,
    unsigned char incremental, char *result_dir, unsigned int total_data_files,
    unsigned int dlsize, // total disk size of data files indexed in dstree
    float warping, unsigned char keyword_search, char * k_values_str, char * ground_truth_dir,
    unsigned char store_results_in_disk, unsigned int num_threads, unsigned int stop_when_nn_dist_changes, unsigned int nn_struct) 
{
  struct bsf_snapshot **bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;

  unsigned int * k_values = NULL;
  unsigned int num_k_values = 0;

  int vector_length = index->settings->timeseries_size;
  int opened_files = 0, qvectors_loaded = 0, curr_vector = 0;

  // query time (for cuurent query column)
  double total_query_time = 0.0;
  unsigned int total_checked_ts = 0;
  unsigned int total_queries = qset_num;
  bool found_query = false; // throw error if no query set was found
  unsigned int offset = 0;


   // initialize list of all knn results (from all query vectors in query set)
  struct result_vid **all_knn_results = NULL;
  int8_t ** recall_matrix = NULL; // recall matrix

  struct query_result *curr_knn = NULL;
  struct job * job_array = NULL;

  // extract k values from string "1,3,5,10" to [1, 3, 5, 10]
  k_values = get_k_values(k_values_str, &num_k_values);
  if (k_values == NULL)
  {
    fprintf(stderr,
              "Error dstree.c:  Could not read set of k values.\n");
      return -1;
  }
  for(int u = 0; u < num_k_values; u++)
      printf(" k = %u\n", k_values[u]);


  // open source dir
  struct dirent *dfile;
  DIR *dir = opendir(bin_files_directory);

  if (!dir) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Unable to open directory stream! %s", bin_files_directory);
    exit(1);
  }

  // create experiment result directory
  char *results_dir = make_result_directory(
      result_dir, total_data_files, qset_num, min_qset_size, max_qset_size);
  
  // max k values for which we must save results
  unsigned int max_k = k_values[num_k_values - 1];


  // for every file (table)
  while ((dfile = readdir(dir)) != NULL) {
    if (qset_num == 0)
      break;

    // if file is binary file
    if (is_binaryfile(dfile->d_name)) {
      opened_files += 1;

      // get fill path of bin file
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, bin_files_directory);
      strcat(bin_file_path, "/");
      strcat(bin_file_path, dfile->d_name);

      // get binary table info
      int datasize, table_id, nsets, vector_length_in_filename;
      sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize,
             &table_id, &nsets, &vector_length_in_filename);

      // check if vector length in file name matches vector length passed as
      // argument
      if (vector_length_in_filename != vector_length) {
        fprintf(stderr,
                "Error in dstree_file_loaders.c:  Vector length passed in "
                "argumentes (--timeseries-size %d) does not match vector "
                "length in file (%d) %s.\n",
                vector_length_in_filename, vector_length, bin_file_path);
        return FAILURE;
      }

      /* read binary file */
      FILE *bin_file = fopen(bin_file_path, "rb");
      if (bin_file == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
                bin_file_path);
        return FAILURE;
      }

      /* Start processing file: read every vector in binary file*/
      int i = 0, j = 0, set_id = 0,
          total_bytes = (datasize * vector_length) + nsets;
      unsigned int nvec = 0u;
      ts_type val;

      
      while (total_bytes) // counts 4 bytes as one because every vector value is
                          // stored in 4 bytes
      {
        // beginning of a set of vectors (column)
        if (i == 0)
        {
          if (qset_num == 0)
            break;

          // read first integer to check how many vactors in current set
          fread(&nvec, sizeof(nvec), 1, bin_file);

          total_bytes--;
          // if set does not fit query requirments move to next set
          if ((unsigned int)nvec < min_qset_size ||
              (unsigned int)nvec > max_qset_size) {
            fseek(bin_file, nvec * 4 * vector_length, SEEK_CUR);
            
            i = 0;
            j = 0;
            total_bytes -= (nvec * vector_length);
            nvec = 0u;
            set_id += 1;
            continue;
          }

          found_query = true;
          total_query_time = 0.0;
          total_checked_ts = 0;
          qset_num--;

          // found a new query set (query clumn): allocate memory for all query vectors and store them in a job array
          job_array = malloc(sizeof(struct job) * nvec);
          if(job_array == NULL)
          {
            fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for query vectors.");
            exit(1);
          }
          for(int q = 0; q < nvec; q++)
          {
            job_array[q].query_vector = calloc(1, sizeof(ts_type) * vector_length);
            job_array[q].query_vector_reordered = calloc(1, sizeof(ts_type) * vector_length);
            job_array[q].query_order = calloc(1, sizeof(int) * vector_length);
            
            if(job_array[q].query_order == NULL || job_array[q].query_vector_reordered == NULL || job_array[q].query_vector == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for query vectors.");
              exit(1);
            }

            job_array[q].query_id.table_id = table_id;
            job_array[q].query_id.set_id = set_id;
            job_array[q].query_id.pos = q;
          }

          // if the user requires the results to be stored in disk
          if(store_results_in_disk)
          {
            all_knn_results = malloc(k * nvec * sizeof(struct result_vid *)); // all knns of all query vectors
            if(all_knn_results == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for knn results.");
              exit(1);
            }
            for (int q = 0, i = 1; q < nvec; ++q) 
            { 
              all_knn_results[q] = calloc(k, sizeof(struct result_vid));
              if(all_knn_results[q] == NULL)
              {
                fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for parallel iqa results.");
                exit(1);
              }
              for(int idx = 0; idx < k; idx++)
              {
                all_knn_results[q][idx].table_id = 0;
                all_knn_results[q][idx].set_id = 0;
                all_knn_results[q][idx].pos = 0;
                all_knn_results[q][idx].qpos = job_array[q].query_id.pos;
                all_knn_results[q][idx].time = 0;
              }
            }
          }
          
          // recall_matrix = malloc(k * nvec * sizeof(int8_t *));
          // // interval where recall should be updated [start, end[ end not included
          // if(recall_matrix == NULL)
          // {
          //   fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for knn results.");
          //   exit(1);
          // }

          // for (int q = 0, i = 1; q < nvec; ++q) 
          // { 
          //   recall_matrix[q] = calloc(k, sizeof(int8_t));
          //   if(recall_matrix[q] == NULL)
          //   {
          //     fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for parallel iqa results.");
          //     exit(1);
          //   }
          // }

          curr_vector = 0;
          if(num_threads > nvec)
          {
            num_threads = nvec;
            fprintf(stderr, "Warning in dstree_file_loader.c: parameter num_thread greater than total query vectors.\n");
          }

          // initialise thread statistics (query time, etc.)
          INIT_THREAD_STATS(num_threads)
          printf("\nQuery (%d, %d), |Q| = %u ...\n", table_id, set_id, nvec);

          set_id += 1;
          i++;
          j = 0;
        } 
        else if (i <= (unsigned int)nvec * vector_length) 
        {
          // end of vector but still in current set
          if (j > (vector_length - 1)) {
            j = 0;
          
            reorder_query(job_array[curr_vector].query_vector, job_array[curr_vector].query_vector_reordered,
                          job_array[curr_vector].query_order, vector_length);
            
            curr_vector += 1;
            qvectors_loaded += 1;
            
          }        
          fread((void *)(&val), sizeof(val), 1, bin_file);

          total_bytes--;
          job_array[curr_vector].query_vector[j] = val;

          // end of last vector in current  set
          if (i == (unsigned int)nvec * vector_length) {
            // last vector, end of query set (column)
            reorder_query(job_array[curr_vector].query_vector, job_array[curr_vector].query_vector_reordered,
                          job_array[curr_vector].query_order, vector_length);
            
            qvectors_loaded += 1;
            curr_vector = 0;

            // load ground truth results (to measure recall)
            int  num_gt_results = 0;
            struct result_vid * ground_truth_results = NULL;
            // char ground_truth_file[255] = "";
            // int  num_gt_results = get_ground_truth_file(ground_truth_dir, job_array[0].query_id.table_id,job_array[0].query_id.set_id, ground_truth_file);
            // struct result_vid * ground_truth_results = get_ground_truth_results(ground_truth_file, num_gt_results);
            // printf("coordinator_thread:\t\tread ground truth results, num results = %d\n", num_gt_results);
            
            // create worker threads
            printf("coordinator_thread:\t\tinit thread pool with  %d threads, |Q| = %d.\n", num_threads, nvec);
            struct pool * thread_pool = (struct pool *) malloc(sizeof(struct pool));
            char data_struct [12] = "";

            if(nn_struct == 0)// use sorted array to store kNNs
            {
              strcpy(data_struct, "sorted-arr\0");
              init_thread_pool(thread_pool, index, epsilon, k,
                  exact_de_parallel_multi_thread_incr_knn_search, num_threads,  offset,
                  r_delta, total_checked_ts, &total_query_time, 
                  warping, all_knn_results, store_results_in_disk,
                  k_values, num_k_values, ground_truth_results,
                  num_gt_results, recall_matrix, nvec, job_array, vector_length, stop_when_nn_dist_changes);            

            }
            else if (nn_struct == 1) // use min max heap to store kNNs
            {
              strcpy(data_struct, "minmax-heap\0");
              init_thread_pool(thread_pool, index, epsilon, k,
                  exact_de_parallel_multi_thread_incr_knn_search_mmheap, num_threads,  offset,
                  r_delta, total_checked_ts, &total_query_time, 
                  warping, all_knn_results, store_results_in_disk,
                  k_values, num_k_values, ground_truth_results,
                  num_gt_results, recall_matrix, nvec, job_array, vector_length, stop_when_nn_dist_changes);            
            }
            else if (nn_struct == 2) // use ostree to store kNNs
            {
              strcpy(data_struct, "os-tree\0");
              init_thread_pool(thread_pool, index, epsilon, k,
                  exact_de_parallel_multi_thread_incr_knn_search_ostree, num_threads,  offset,
                  r_delta, total_checked_ts, &total_query_time, 
                  warping, all_knn_results, store_results_in_disk,
                  k_values, num_k_values, ground_truth_results,
                  num_gt_results, recall_matrix, nvec, job_array, vector_length, stop_when_nn_dist_changes);
            }
            else
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: unkown nn_struct value!\n");
              exit(1);
            }

            double max_cpu_time = 0, min_cpu_time = FLT_MAX; 
            int max_cpu_thread_id = -1, min_cpu_thread_id = -1;

            double ** thread_time = malloc(sizeof(double*) * num_threads);
            if(thread_time == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for thread time.\n");
              exit(1);
            }
            // sum up time for all the queryies executed by one thread
            for(int th = 0; th < num_threads; th++)
            {
              thread_time[th] = calloc(k, sizeof(double));
              if(thread_time[th] == NULL)
              {
                fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for thread time.\n");
                exit(1);
              }
            }
            
            // report all query time and query results
            char file_res [] = "query_results.csv\0";
            char file_time [] = "query_time.csv\0";
            FILE *fpt;
            fpt = fopen(file_time, "a");

            if (fpt == NULL) {
              fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s or file %s!\n", file_res, file_time);
              return FAILURE;
            }

            fprintf(fpt, "k,querytime,nb_threads,data_structure\n");
            
            for(int a = 0; a < nvec; a++)
            {
              int thread_id = job_array[a].worker_id;
              for(int b = 0; b < k; b++)
              {
                thread_time[thread_id][b] += all_knn_results[a][b].time;
              }
            }
            
            double * max_thread_time = malloc(k * sizeof(double));
            if(max_thread_time == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for max thread time array.\n");
              exit(1);
            }
            printf("k,\tt,\t#th,\tdata_structure\n");
            for(int b = 0; b < num_k_values; b++)
            {
              unsigned int curr_k = k_values[b];
              double max_cpu_time = 0, min_cpu_time = FLT_MAX; 
              int max_cpu_thread_id = -1, min_cpu_thread_id = -1;

              for(int th = 0; th < num_threads; th++)
              {
                if(thread_time[th][curr_k-1] > max_cpu_time)
                { 
                  max_cpu_time = thread_time[th][curr_k-1];
                  max_cpu_thread_id = th;
                }
              }
              max_thread_time[curr_k-1] = max_cpu_time;
              
              printf("%u,\t%f,\t%d,\t%s\n", curr_k, max_cpu_time/1000000, num_threads, data_struct);
              fprintf(fpt, "%u,%f,%d,%s\n", curr_k, max_cpu_time/1000000, num_threads, data_struct);
            }
            fclose(fpt);

            // measure recall at each nn, to visualize recall improvement/ degredattion
            // unsigned int num_gt_neighbors = 0;
            // float recall = 0.0;
            // float k_recall = 0.0;

            // for(int a = 0; a < 1; a++)
            // {
            //   // num_gt_neighbors = get_num_identical_results(ground_truth_results, num_gt_results, a);
            //   printf("recall report for query vector %d -- -- -- --\n");
            //   printf("k;\t\ts;\t\td;\t\trecall;\t\t\n");
              
            //   for(int b = 0; b < k; b++)
            //   {
            //     recall = compute_vector_recall(all_knn_results[a], ground_truth_results, a, b+1);
            //     printf("%d;\t\t(%u,%u,%u);\t\t%f;\t\t%f;\n", b+1, all_knn_results[a][b].table_id,
            //            all_knn_results[a][b].set_id, all_knn_results[a][b].pos, all_knn_results[a][b].distance, recall);
            //   }
            // }


            // printf("recall report for query column (%u, %u) -- -- -- --\n", job_array[0].query_id.table_id, job_array[0].query_id.set_id);
            //   printf("k;\t\trecall;\t\t\n");
            // for(int b = 0; b < k; b++)
            // {
            //   recall = compute_k_recall_from_matrix(all_knn_results, ground_truth_results, nvec, b+1, num_gt_results);
            //   printf("%d;\t\t%f;\n", b+1, recall);
            // }

            
            // for(int b = 0; b < k; b++)
            // {
            //   printf("recall report for k -- -- -- --\n");
            //   printf("k;\t\trecall;\t\t\n");
              
            //   for(int a = 1; a < nvec; a++)
            //   {
            //     recall = compute_vector_recall_identical_results(recall_matrix, a, b+1, b+1);
            //     printf("%d;\t\t(%u,%u,%u);\t\t%f;\t\t%f;\n", b+1, all_knn_results[a][b].table_id,
            //            all_knn_results[a][b].set_id, all_knn_results[a][b].pos, all_knn_results[a][b].distance, recall);
            //   }
            // }

            // store query results
            if(store_results_in_disk)
            {
              printf("start storing results to csv file...\n");

              for(int z = 0; z < num_k_values; z++)
              {
                unsigned int curr_k = k_values[z];
                
                char *query_result_file = make_file_path(results_dir, job_array[0].query_id.table_id, job_array[0].query_id.set_id, nvec,
                                        total_data_files, dlsize, vector_length, curr_k);

                if(!store_parallel_knn_results_in_disk(query_result_file, job_array[0].query_id.table_id, job_array[0].query_id.set_id,
                                        nvec, all_knn_results, curr_k, max_thread_time))
                {
                  fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't save query results to file %s.", query_result_file);
                  exit(1);
                }
                free(query_result_file);              
              }
              printf("end.\n");
            }

            // aggregate results and print joinable tables (disabled)
            if(0)  
            {
              if(keyword_search) // print results to be captured by ui
              {
                // don't change these lines to allow ui to fetch results
                struct result_table* top = get_top_tables_by_euclidean_distance(all_knn_results, k*nvec, num_top);
                for(int m = 0; m < num_top; m++)
                {
                  printf("table-%u- in file @@%s$ min_distance=%.3f§ num_closest=%u# total_matches=%uµ\n", top[m].table_id, top[m].raw_data_file, top[m].min_distance, top[m].num_min, top[m].total_matches);
                }
                free(top);
              } 
              else
              {
                unsigned int total_matching_columns = 0;
                struct result_sid * top = get_top_sets(all_knn_results, nvec, k, -1, &total_matching_columns);

                // don't change these lines to allaow ui to fetch results
                for(int m = 0; m < num_top; m++) // print results to be captured by ui
                {
                  printf("table-%u-column-%u- in file @@%s$ overlap=%u§\n", top[m].table_id, top[m].set_id, top[m].raw_data_file, top[m].overlap_size);
                }
                free(top);
              }
            }

            // free memory
            // kill worker thread and destroy barrier
            // pthread_cancel(worker_thread);
            // pthread_barrier_destroy(&knn_update_barrier);

            // free thread time
            free(max_thread_time);
            for(int th = 0; th < num_threads; th++)
            {
              free(thread_time[th]);
            }
            free(thread_time);
            
            // // free recall matrix
            // for (int q = 0, i = 1; q < nvec; ++q) 
            // { 
            //   free(recall_matrix[q]);
            // }
            // free(recall_matrix);
            
            // free ground truth results
            free(ground_truth_results);

            // free query vectors
            for(int q = 0; q < nvec; q++)
            {
              free(job_array[q].query_vector);
              free(job_array[q].query_vector_reordered);
              free(job_array[q].query_order);
            }
            free(job_array);
            job_array = NULL;

            // free knn results
            if(store_results_in_disk)
            {
              for (int q = 0, i = 1; q < nvec; ++q) 
              { 
                free(all_knn_results[q]);
              }
              free(all_knn_results);
            }
            
            FREE_THREAD_QUERY_COUNTERS
            FREE_THREAD_PARTIAL_COUNTERS

            i = 0;
            j = 0;
            nvec = 0u;
            continue;
          }
          i++;
          j++;
        }
      }
      fclose(bin_file);
      /* end read binary file */
    }
  }

  if (found_query == false) {
    fprintf(stderr,
            "WARNING! Could not find any query in size range [%u - %u].\n",
            min_qset_size, max_qset_size);
    exit(-1);
  }
  if (opened_files == 0) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c:  Could not find any binary file "
            "in directory %s.\n",
            bin_files_directory);
    return FAILURE;
  }
  
  // free memory
  closedir(dir);
  free(results_dir);
  free(k_values); 
  return SUCCESS;
}

enum response dstree_index_multiple_binary_files(const char *bin_files_directory,
                                   unsigned int total_data_files,
                                   struct dstree_index *index) {

  printf("\n+ must read %dfiles\n", total_data_files);

  int vector_length = index->settings->timeseries_size;
  int opened_files = 0;
  unsigned int total_vectors  = 0;
  unsigned int total_datasize  = 0;
  unsigned int total_read_files = 0;

  // allocate memory for vector
  struct vector v;
  v.values = (ts_type *)malloc(sizeof(ts_type) * vector_length);
  ts_type val;
  unsigned int nvec = 0u;

  // open source dir
  struct dirent *dfile;
  DIR *dir = opendir(bin_files_directory);
  if (!dir) {
    printf("Error in dstree_file_loader.c: Unable to open directory stream!");
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0) {
    if (is_binaryfile(dfile->d_name)) {
      total_data_files--;
      opened_files += 1;

      // get fill path of bin file
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, bin_files_directory);
      strcat(bin_file_path, "/");
      strcat(bin_file_path, dfile->d_name);

      // get binary table info
      int datasize, table_id, nsets, vector_length_in_filename;
      sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize,
             &table_id, &nsets, &vector_length_in_filename);
      total_datasize += datasize;
      // check if vector length in file name matches vector length passed as
      // argument
      if (vector_length_in_filename != vector_length) {
        fprintf(stderr,
                "Error in dstree_file_loaders.c:  Vector length passed in "
                "argumentes (--timeseries-size %d) does not match vector "
                "length in file (%d) %s.\n",
                vector_length_in_filename, vector_length, bin_file_path);
        return FAILURE;
      }

      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START

      COUNT_PARTIAL_RAND_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      FILE *bin_file = fopen(bin_file_path, "rb");
      total_read_files += 1;
      COUNT_PARTIAL_INPUT_TIME_END

      if (bin_file == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
                bin_file_path);
        return FAILURE;
      }

      /* Start processing file: read every vector in binary file*/
      unsigned int i = 0, j = 0, set_id = 0, vector_pos = 0,
                   total_bytes = (datasize * vector_length) + nsets;

      while (total_bytes) {
        if (i == 0) {
          i++;
          j = 0;

          // read first integer to check how many vactors in current set
          fread(&nvec, sizeof(nvec), 1, bin_file);
          total_bytes--;
          v.table_id = table_id;
          v.set_id = set_id;
          v.pos = 0;
          set_id += 1;

        } else if (i <= (unsigned int)nvec * vector_length) {
          // end of vector but still in current set
          if (j > (vector_length - 1)) {
            j = 0;
            /*Index v in dstree */
            if (!dstree_index_insert_vector(index, v.values, v.table_id,
                                            v.set_id, v.pos, dfile->d_name)) {
              fprintf(stderr, "Error in dstree_file_loaders.c:  Could not add "
                              "the time series to the index.\n");
              return FAILURE;
            }
            total_vectors += 1;
            // next vector position in set (/column)
            v.pos += 1;

            index->stats->idx_building_total_time += partial_time;
            index->stats->idx_building_input_time += partial_input_time;
            index->stats->idx_building_output_time += partial_output_time;

            index->stats->idx_building_seq_input_count +=
                partial_seq_input_count;
            index->stats->idx_building_seq_output_count +=
                partial_seq_output_count;
            index->stats->idx_building_rand_input_count +=
                partial_rand_input_count;
            index->stats->idx_building_rand_output_count +=
                partial_rand_output_count;
          }

          COUNT_PARTIAL_SEQ_INPUT
          COUNT_PARTIAL_INPUT_TIME_START
          fread((void *)(&val), sizeof(val), 1, bin_file);
          COUNT_PARTIAL_INPUT_TIME_END
          total_bytes--;
          v.values[j] = val;

          // end of last vector in current  set
          if (i == (unsigned int)nvec * vector_length) {
            /*Index v in dstree */
            if (!dstree_index_insert_vector(index, v.values, v.table_id,
                                            v.set_id, v.pos, dfile->d_name)) {
              fprintf(stderr, "Error in dstree_file_loaders.c:  Could not add "
                              "the time series to the index.\n");
              return FAILURE;
            }
            total_vectors += 1;
            // first vector position in next set (/column)
            v.pos = 0;

            index->stats->idx_building_total_time += partial_time;
            index->stats->idx_building_input_time += partial_input_time;
            index->stats->idx_building_output_time += partial_output_time;

            index->stats->idx_building_seq_input_count +=
                partial_seq_input_count;
            index->stats->idx_building_seq_output_count +=
                partial_seq_output_count;
            index->stats->idx_building_rand_input_count +=
                partial_rand_input_count;
            index->stats->idx_building_rand_output_count +=
                partial_rand_output_count;

            i = 0;
            j = 0;
            nvec = 0u;
            continue;
          }
          i++;
          j++;
        }
      }
      printf("(+) inserted %d vectors in table %u\n", datasize, table_id);

      COUNT_PARTIAL_INPUT_TIME_START
      if (fclose(bin_file)) {
        fprintf(
            stderr,
            "Error in dstree_file_loaders.c: Could not close the filename %s",
            bin_file_path);
        return FAILURE;
      }
      COUNT_PARTIAL_INPUT_TIME_END
      COUNT_PARTIAL_TIME_END
      /* end read processing file */
      index->stats->idx_building_total_time += partial_time;
      index->stats->idx_building_input_time += partial_input_time;
      index->stats->idx_building_output_time += partial_output_time;
      index->stats->idx_building_seq_input_count += partial_seq_input_count;
      index->stats->idx_building_seq_output_count += partial_seq_output_count;
      index->stats->idx_building_rand_input_count += partial_rand_input_count;
      index->stats->idx_building_rand_output_count += partial_rand_output_count;
      RESET_PARTIAL_COUNTERS()
      COUNT_PARTIAL_TIME_START
    }
  }
  COUNT_PARTIAL_INPUT_TIME_START
  closedir(dir);
  COUNT_PARTIAL_INPUT_TIME_END
  free(v.values);
  if (opened_files == 0) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c:  Could not find any binary file "
            "in directory %s.\n",
            bin_files_directory);
    return FAILURE;
  }

  printf("\ntotal_read_files = %d\n", total_read_files);
  printf("total_vectors = %d\n", total_vectors);
  printf("total_datasize = %d\n", total_datasize);
  return SUCCESS;
}
/* end kashif changes */
enum response dstree_knn_query_gt_binary_file(
    const char *ifilename, int q_num, struct dstree_index *index,
    float minimum_distance, ts_type epsilon, ts_type r_delta, unsigned int k,
    boolean track_bsf, boolean track_pruning, boolean all_mindists,
    boolean max_policy, unsigned int nprobes, unsigned char incremental) {

  struct bsf_snapshot **bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
  if (track_bsf) {
    max_bsf_snapshots = 10000;
    cur_bsf_snapshot = 0;

    bsf_snapshots = calloc(k, sizeof(struct bsf_snapshot *));
    for (unsigned int i = 0; i < k; ++i) {
      bsf_snapshots[i] = calloc(max_bsf_snapshots, sizeof(struct bsf_snapshot));
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
        bsf_snapshots[i][j].distance = FLT_MAX;
        bsf_snapshots[i][j].time = FLT_MAX;
        bsf_snapshots[i][j].series = NULL;
        bsf_snapshots[i][j].checked_nodes = -1;
      }
    }
  }
  /*
  const char *filename = malloc(sizeof(char) *
 (strlen(index->settings->root_directory) + 18)); filename = strcpy(filename,
 index->settings->root_directory); filename = strcat(filename,
 "raw_series.csv\0");

 printf ("series_file = %s\n", filename);
 printf ("dataset_file = %s\n", index->settings->dataset);

 FILE *series_file = fopen(filename, "w");
 FILE *dataset_file = fopen(index->settings->dataset, "rb");
 */
  FILE *series_file = NULL;
  FILE *dataset_file = NULL;

  RESET_PARTIAL_COUNTERS()

  COUNT_PARTIAL_TIME_START

  FILE *ifile;
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "rb");
  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "File %s not found!\n", ifilename);
    exit(-1);
  }

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END
  unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz / ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    exit(-1);
  }

  unsigned int q_loaded = 0;
  ts_type *query_ts = calloc(1, sizeof(ts_type) * ts_length);
  ts_type *query_ts_reordered = calloc(1, sizeof(ts_type) * ts_length);
  int *query_order = calloc(1, sizeof(int) * ts_length);
  if (query_order == NULL)
    return FAILURE;

  while (q_loaded < q_num) {

    RESET_QUERY_COUNTERS()

    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(query_ts, sizeof(ts_type), ts_length, ifile);
    COUNT_PARTIAL_INPUT_TIME_END

    reorder_query(query_ts, query_ts_reordered, query_order, ts_length);

    q_loaded++;

    if (track_bsf) {
      cur_bsf_snapshot = 0;
      if (incremental) {
        exact_de_incr_progressive_knn_search(
            query_ts, query_ts_reordered, query_order, offset, index,
            minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,
            bsf_snapshots, &cur_bsf_snapshot, 0, dataset_file, series_file);
      } else {
        exact_de_progressive_knn_search(
            query_ts, query_ts_reordered, query_order, offset, index,
            minimum_distance, epsilon, r_delta, k, q_loaded, ifilename,
            bsf_snapshots, &cur_bsf_snapshot);
      }
      for (unsigned int i = 0; i < k; ++i) {
        for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
          bsf_snapshots[i][j].distance = FLT_MAX;
          bsf_snapshots[i][j].time = FLT_MAX;
          bsf_snapshots[i][j].series = NULL;
          bsf_snapshots[i][j].checked_nodes = -1;
        }
      }
    } else {
      exact_de_knn_search(query_ts, query_ts_reordered, query_order, offset,
                          index, minimum_distance, epsilon, r_delta, k,
                          q_loaded, ifilename);
    }

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
  }
  COUNT_PARTIAL_TIME_END
  RESET_PARTIAL_COUNTERS()

  free(query_ts);
  free(query_ts_reordered);
  free(query_order);

  if (fclose(ifile)) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not close the query filename %s",
        ifilename);
    return FAILURE;
  }

  if (track_bsf) {
    for (unsigned int i = 0; i < k; ++i) {
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
        if (bsf_snapshots[i][j].series != NULL)
          free(bsf_snapshots[i][j].series);
      }
      free(bsf_snapshots[i]);
    }
    free(bsf_snapshots);
  }
  // fclose (dataset_file);
  // fclose (series_file);

  return SUCCESS;
}

enum response dstree_tlb_binary_file(const char *ifilename, int q_num,
                                     struct dstree_index *index,
                                     float minimum_distance) {

  FILE *ifile;

  ifile = fopen(ifilename, "rb");

  if (ifile == NULL) {
    fprintf(stderr, "File %s not found!\n", ifilename);
    exit(-1);
  }

  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  fseek(ifile, 0L, SEEK_SET);
  unsigned int ts_length = index->settings->timeseries_size;
  file_position_type total_records = sz / ts_length * sizeof(ts_type);
  unsigned int offset = 0;

  if (total_records < q_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    exit(-1);
  }

  unsigned int q_loaded = 0;
  ts_type *query_ts = calloc(1, sizeof(ts_type) * ts_length);

  while (q_loaded < q_num) {
    total_tlb = 0;
    total_ts_count = 0;
    leaf_nodes_count = 0;

    fread(query_ts, sizeof(ts_type), ts_length, ifile);

    dstree_calc_tlb(query_ts, index, index->first_node);

    q_loaded++;
    print_tlb_stats(index, q_loaded, ifilename);
  }

  free(query_ts);

  if (fclose(ifile)) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not close the query filename %s",
        ifilename);
    return FAILURE;
  }

  return SUCCESS;
}

enum response dstree_index_ascii_file(const char *ifilename,
                                      file_position_type ts_num,
                                      const char delimiter,
                                      struct dstree_index *index) {
  double parse_time = 0;
  int ts_length = index->settings->timeseries_size;
  ts_type *ts = NULL;

  ts = NULL;
  ts = malloc(sizeof(ts_type) * ts_length);

  if (ts == NULL) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
    return FAILURE;
  }

  FILE *ifile;
  COUNT_INPUT_TIME_START
  ifile = fopen(ifilename, "r");
  COUNT_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "File %s not found!\n", ifilename);
    exit(-1);
  }

  char *ts_str = NULL; //= malloc(sizeof(char) * 2000);
  size_t linecap = 0;
  ssize_t linelen;

  file_position_type ts_loaded = 0;

  int percentage = 100;
  if (percentage > 100) {
    percentage = (int)(ts_num / (file_position_type)100);
  }

  COUNT_INPUT_TIME_START
  while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 &&
         ts_loaded < ts_num) {
    COUNT_INPUT_TIME_END
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
    printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m", (ts_loaded + 1));
#endif
#endif
    COUNT_PARSE_TIME_START
    ts_parse_str(ts_str, ts, index->settings->timeseries_size, &delimiter);
    COUNT_PARSE_TIME_END

    if (!dstree_index_insert(index, ts)) {
      fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
      return FAILURE;
    }

    ts_loaded++;

    if (ts_loaded % percentage == 0) {
      float distance = 0;
    }

    COUNT_INPUT_TIME_START
  }

  free(ts_str);
  free(ts);
  COUNT_INPUT_TIME_START
  fclose(ifile);
  COUNT_INPUT_TIME_END

  return SUCCESS;
}

enum response dstree_index_binary_file(const char *ifilename,
                                       file_position_type ts_num,
                                       struct dstree_index *index) {
  double parse_time = 0;
  int ts_length = index->settings->timeseries_size;
  ts_type *ts = NULL;

  ts = NULL;
  ts = malloc(sizeof(ts_type) * ts_length);

  if (ts == NULL) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
    return FAILURE;
  }

  FILE *ifile;

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "rb");
  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
            ifilename);
    return FAILURE;
  }
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  COUNT_PARTIAL_INPUT_TIME_END
  file_position_type total_records =
      sz / index->settings->timeseries_size * sizeof(ts_type);

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END
  if (total_records < ts_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    return FAILURE;
  }

  file_position_type ts_loaded = 0;

  int percentage = 100;
  if (percentage > 100) {
    percentage = (int)(ts_num / (file_position_type)100);
  }

  while (ts_loaded < ts_num) {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
    printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m", (ts_loaded + 1));
#endif
#endif
    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
    COUNT_PARTIAL_INPUT_TIME_END

    if (!dstree_index_insert(index, ts)) {
      fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
      return FAILURE;
    }
    COUNT_PARTIAL_TIME_END

    index->stats->idx_building_total_time += partial_time;
    index->stats->idx_building_input_time += partial_input_time;
    index->stats->idx_building_output_time += partial_output_time;

    index->stats->idx_building_seq_input_count += partial_seq_input_count;
    index->stats->idx_building_seq_output_count += partial_seq_output_count;
    index->stats->idx_building_rand_input_count += partial_rand_input_count;
    index->stats->idx_building_rand_output_count += partial_rand_output_count;

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    ts_loaded++;

    if (ts_loaded % percentage == 0) {
      float distance = 0;
      // PRINT_STATS(distance)
    }
  }

  free(ts);
  COUNT_PARTIAL_INPUT_TIME_START
  if (fclose(ifile)) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c: Could not close the filename %s",
            ifilename);
    return FAILURE;
  }
  COUNT_PARTIAL_INPUT_TIME_END

  COUNT_PARTIAL_TIME_END
  index->stats->idx_building_total_time += partial_time;
  index->stats->idx_building_input_time += partial_input_time;
  index->stats->idx_building_output_time += partial_output_time;
  index->stats->idx_building_seq_input_count += partial_seq_input_count;
  index->stats->idx_building_seq_output_count += partial_seq_output_count;
  index->stats->idx_building_rand_input_count += partial_rand_input_count;
  index->stats->idx_building_rand_output_count += partial_rand_output_count;

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  return SUCCESS;
}

enum response dstree_index_classify_binary_file(const char *ifilename,
                                                file_position_type ts_num,
                                                struct dstree_index *index) {
  double parse_time = 0;
  int ts_length = index->settings->timeseries_size;
  ts_type *ts = NULL;

  int filename_size = strlen(ifilename) + 4;
  const char *gt_filename = malloc(sizeof(char) * filename_size);
  gt_filename = strcpy(gt_filename, "\0");
  gt_filename = strcat(gt_filename, ifilename);
  gt_filename = strcat(gt_filename, ".gt\0");

  ts = NULL;
  ts = malloc(sizeof(ts_type) * ts_length);

  // unsigned char ts_gt = 0;
  label_type ts_gt = 0;

  if (ts == NULL) {
    fprintf(
        stderr,
        "Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
    return FAILURE;
  }

  FILE *ifile;
  FILE *gt_file;

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  ifile = fopen(ifilename, "rb");
  gt_file = fopen(gt_filename, "rb");

  COUNT_PARTIAL_INPUT_TIME_END
  if (ifile == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
            ifilename);
    return FAILURE;
  }
  if (gt_file == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n",
            gt_filename);
    return FAILURE;
  }
  free(gt_filename);
  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_END);
  file_position_type sz = (file_position_type)ftell(ifile);
  COUNT_PARTIAL_INPUT_TIME_END
  file_position_type total_records =
      sz / index->settings->timeseries_size * sizeof(ts_type);

  COUNT_PARTIAL_RAND_INPUT
  COUNT_PARTIAL_INPUT_TIME_START
  fseek(ifile, 0L, SEEK_SET);
  COUNT_PARTIAL_INPUT_TIME_END
  if (total_records < ts_num) {
    fprintf(stderr, "File %s has only %llu records!\n", ifilename,
            total_records);
    return FAILURE;
  }

  // file_position_type ts_loaded = 0;
  unsigned int ts_loaded = 0;

  int percentage = 100;
  if (percentage > 100) {
    percentage = (int)(ts_num / (file_position_type)100);
  }

  while (ts_loaded < ts_num) {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
    printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m", (ts_loaded + 1));
#endif
#endif
    COUNT_PARTIAL_SEQ_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
    // fread(&ts_gt, sizeof(unsigned char), 1, gt_file);
    fread(&ts_gt, sizeof(label_type), 1, gt_file);
    COUNT_PARTIAL_INPUT_TIME_END

    if (!index->settings->track_file_pos) {
      if (!dstree_index_classify_insert(index, ts, ts_gt, -1)) {
        fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
        return FAILURE;
      }
    } else {
      if (!dstree_index_classify_insert(index, ts, ts_gt, (ts_loaded + 1))) {
        fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
        return FAILURE;
      }
    }

    COUNT_PARTIAL_TIME_END

    index->stats->idx_building_total_time += partial_time;
    index->stats->idx_building_input_time += partial_input_time;
    index->stats->idx_building_output_time += partial_output_time;

    index->stats->idx_building_seq_input_count += partial_seq_input_count;
    index->stats->idx_building_seq_output_count += partial_seq_output_count;
    index->stats->idx_building_rand_input_count += partial_rand_input_count;
    index->stats->idx_building_rand_output_count += partial_rand_output_count;

    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

    ts_loaded++;

    if (ts_loaded % percentage == 0) {
      float distance = 0;
      // PRINT_STATS(distance)
    }
  }

  free(ts);
  COUNT_PARTIAL_INPUT_TIME_START
  if (fclose(ifile)) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c: Could not close the filename %s",
            ifilename);
    return FAILURE;
  }

  if (fclose(gt_file)) {
    fprintf(stderr,
            "Error in dstree_file_loaders.c: Could not close the filename %s",
            gt_file);
    return FAILURE;
  }
  COUNT_PARTIAL_INPUT_TIME_END

  COUNT_PARTIAL_TIME_END
  index->stats->idx_building_total_time += partial_time;
  index->stats->idx_building_input_time += partial_input_time;
  index->stats->idx_building_output_time += partial_output_time;
  index->stats->idx_building_seq_input_count += partial_seq_input_count;
  index->stats->idx_building_seq_output_count += partial_seq_output_count;
  index->stats->idx_building_rand_input_count += partial_rand_input_count;
  index->stats->idx_building_rand_output_count += partial_rand_output_count;

  RESET_PARTIAL_COUNTERS()
  COUNT_PARTIAL_TIME_START

  return SUCCESS;
}

enum response reorder_query(ts_type *query_ts, ts_type *query_ts_reordered,
                            int *query_order, int ts_length) {

  q_index *q_tmp = malloc(sizeof(q_index) * ts_length);
  int i;

  if (q_tmp == NULL)
    return FAILURE;

  for (i = 0; i < ts_length; i++) {
    q_tmp[i].value = query_ts[i];
    q_tmp[i].index = i;
  }

  qsort(q_tmp, ts_length, sizeof(q_index), znorm_comp);

  for (i = 0; i < ts_length; i++) {

    query_ts_reordered[i] = q_tmp[i].value;
    query_order[i] = q_tmp[i].index;
  }
  free(q_tmp);

  return SUCCESS;
}

int znorm_comp(const void *a, const void *b) {
  q_index *x = (q_index *)a;
  q_index *y = (q_index *)b;

  if (fabsf(y->value) > fabsf(x->value))
    return 1;
  else if (fabsf(y->value) == fabsf(x->value))
    return 0;
  else
    return -1;
}

