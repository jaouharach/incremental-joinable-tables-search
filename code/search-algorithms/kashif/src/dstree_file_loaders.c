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

  struct bsf_snapshot **bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
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


  if (track_bsf) {
    max_bsf_snapshots = 0;
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
        bsf_snapshots[i][j].vector_id = malloc(sizeof(struct vid));
        bsf_snapshots[i][j].query_vector_pos = -1;
      }
    }
  }

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
  struct query_result *all_knn_results = NULL;
  struct query_result *curr_knn = NULL;
  struct vid *top_matches;
  struct vid query_id;
  int knn_array_idx = 0;

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

      while (total_bytes) // counts 4 bytes as one because every vector is
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
            COUNT_PARTIAL_INPUT_TIME_START
            fseek(bin_file, nvec * 4 * vector_length, SEEK_CUR);
            COUNT_PARTIAL_INPUT_TIME_END
            i = 0;
            j = 0;
            total_bytes -= (nvec * vector_length);
            nvec = 0u;
            continue;
          }
          found_query = true;

          query_time = 0.0;
          total_checked_ts = 0;
          qset_num--;

          all_knn_results = malloc(k * nvec * sizeof(struct query_result));
          for (int knn = 0; knn < (k * nvec); knn++) {
            all_knn_results[knn].vector_id = malloc(sizeof(struct vid));
          }

          knn_array_idx = 0;
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
            
            // printf("-- vector %u\n", query_vector.pos);

            // printf("before reorder\n");
            // for(int j = 0; j <index->settings->timeseries_size; j++)
            //   printf("%.3f - ", query_vector.values[j]);
            // printf("\n");
            /*run query vector in dstree */
            reorder_query(query_vector.values, query_vector_reordered,
                          query_order, vector_length);
            qvectors_loaded += 1;
            // printf("after reorder\n");
            // for(int j = 0; j <index->settings->timeseries_size; j++)
            //   printf("%.3f - ", query_vector.values[j]);
            // printf("\n");

            // perform extact knn search
            // with incremental answering
            if (track_bsf) {
              cur_bsf_snapshot = 0;

              query_id.table_id = query_vector.table_id;
              query_id.set_id = query_vector.set_id;
              query_id.pos = query_vector.pos;

              if (incremental) {
                curr_knn = exact_de_incr_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_time,
                    &total_checked_ts, bsf_snapshots, &cur_bsf_snapshot,
                    warping, dataset_file, series_file, &query_id);
                // printf("vector %d/%d, done.\n", query_vector.pos+1, nvec);
              } else {
                curr_knn = exact_de_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_time,
                    &total_checked_ts, bsf_snapshots, &cur_bsf_snapshot, query_vector.pos);
              }
              for (unsigned int b = 0; b < k; ++b) {
                for (unsigned int f = 0; f < max_bsf_snapshots; ++f) {
                  bsf_snapshots[b][f].distance = FLT_MAX;
                  bsf_snapshots[b][f].time = FLT_MAX;
                  bsf_snapshots[b][f].series = NULL;
                  bsf_snapshots[b][f].checked_nodes = -1;
                  bsf_snapshots[b][f].label = 0;
                  bsf_snapshots[b][f].vector_id->table_id = -1;
                  bsf_snapshots[b][f].vector_id->set_id = -1;
                  bsf_snapshots[b][f].vector_id->pos = -1;
                  bsf_snapshots[b][f].query_vector_pos = -1;
                  strcpy(bsf_snapshots[b][f].vector_id->raw_data_file, "");
                }
              }
            }
            // without incremental answering
            else {
              curr_knn = exact_de_knn_search_2(
                  query_vector.values, query_vector_reordered, query_order,
                  offset, index, minimum_distance, epsilon, r_delta, k,
                  qvectors_loaded, bin_file_path, &query_time,
                  &total_checked_ts, query_vector.pos);
            }

            // copy new knn(s) to knn_results array
            // printf("end of query for vector sent %u received %u\n", query_vector.pos, curr_knn[0].vector_id->pos);
            for (int t = 0; t < k; t++) 
            {
              all_knn_results[knn_array_idx].vector_id->table_id = curr_knn[t].vector_id->table_id;
              all_knn_results[knn_array_idx].vector_id->set_id = curr_knn[t].vector_id->set_id;
              all_knn_results[knn_array_idx].vector_id->pos = curr_knn[t].vector_id->pos;
              all_knn_results[knn_array_idx].query_vector_pos = curr_knn[t].query_vector_pos;
              strcpy(all_knn_results[knn_array_idx].vector_id->raw_data_file, curr_knn[t].vector_id->raw_data_file);
              // printf("match in file %s\n", (all_knn_results[knn_array_idx].vector_id->raw_data_file));
              all_knn_results[knn_array_idx].distance = curr_knn[t].distance;
              all_knn_results[knn_array_idx].time = curr_knn[t].time;
              knn_array_idx++;
            }
            query_vector.pos += 1;
            // printf("next query %u\n", query_vector.pos);
            if (knn_array_idx > (k * nvec)) {
              fprintf(stderr, "Error in dstree_file_loaders.c: Storing more results "
                     "that expected!");
              exit(1);
            }

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
              cur_bsf_snapshot = 0;

              query_id.table_id = query_vector.table_id;
              query_id.set_id = query_vector.set_id;
              query_id.pos = query_vector.pos;
              
              if (incremental) {
                curr_knn = exact_de_incr_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_time,
                    &total_checked_ts, bsf_snapshots, &cur_bsf_snapshot,
                    warping, dataset_file, series_file, &query_id);
                    // printf("vector %d/%d, done.\n", query_vector.pos+1, nvec);
              } else {
                curr_knn = exact_de_progressive_knn_search_2(
                    query_vector.values, query_vector_reordered, query_order,
                    offset, index, minimum_distance, epsilon, r_delta, k,
                    qvectors_loaded, bin_file_path, &query_time,
                    &total_checked_ts, bsf_snapshots, &cur_bsf_snapshot, query_vector.pos);
              }
              for (unsigned int b = 0; b < k; ++b) {
                for (unsigned int f = 0; f < max_bsf_snapshots; ++f) {
                  bsf_snapshots[b][f].distance = FLT_MAX;
                  bsf_snapshots[b][f].time = FLT_MAX;
                  bsf_snapshots[b][f].series = NULL;
                  bsf_snapshots[b][f].checked_nodes = -1;
                  bsf_snapshots[b][f].label = 0;
                  bsf_snapshots[b][f].vector_id->table_id = -1;
                  bsf_snapshots[b][f].vector_id->set_id = -1;
                  bsf_snapshots[b][f].vector_id->pos = -1;
                  bsf_snapshots[b][f].query_vector_pos = -1;
                  strcpy(bsf_snapshots[b][f].vector_id->raw_data_file, "");
                }
              }
            }
            // without incremental answering
            else {
              curr_knn = exact_de_knn_search_2(
                  query_vector.values, query_vector_reordered, query_order,
                  offset, index, minimum_distance, epsilon, r_delta, k,
                  qvectors_loaded, bin_file_path, &query_time,
                  &total_checked_ts, query_vector.pos);
            }
            // append new knn(s) to knn_results array
            // printf("end of query for vector sent %u received %u\n", query_vector.pos, curr_knn[0].vector_id->pos);
            for (int t = 0; t < k; t++) {
              all_knn_results[knn_array_idx].vector_id->table_id = curr_knn[t].vector_id->table_id;
              all_knn_results[knn_array_idx].vector_id->set_id = curr_knn[t].vector_id->set_id;
              all_knn_results[knn_array_idx].vector_id->pos = curr_knn[t].vector_id->pos;
              all_knn_results[knn_array_idx].query_vector_pos = curr_knn[t].query_vector_pos;
              strcpy(all_knn_results[knn_array_idx].vector_id->raw_data_file,  curr_knn[t].vector_id->raw_data_file);
              all_knn_results[knn_array_idx].distance = curr_knn[t].distance;
              all_knn_results[knn_array_idx].time = curr_knn[t].time;

              knn_array_idx++;
            }
            query_vector.pos = 0;

            for (int t = 0; t < k; t++)
              free(curr_knn[t].vector_id);
            free(curr_knn);
            
            if (knn_array_idx > (k * nvec)) {
              fprintf(stderr, "Error in dstree_file_loaders.c: Storing more results "
                     "that expected!");
              exit(1);
            }
            // End of Query set
            /* Save query results to csv file */
            query_time /= 1000000;
            printf("Storing result to csv file...");

            for(int z = 0; z < num_k_values; z++)
            {
              unsigned int curr_k = k_values[z];
              
              char *query_result_file = make_file_path(results_dir, table_id, query_vector.set_id, nvec,
                                      total_data_files, dlsize, vector_length, curr_k);

              if(!save_to_query_result_file(query_result_file, table_id, query_vector.set_id,
                                      k * nvec, all_knn_results, k, curr_k))
              {
                fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't save query results to file %s.", query_result_file);
                exit(1);
              }
              free(query_result_file);
            }
            
            printf("\nquery_time=%fsec\n", query_time);
            printf("computing recall...\n");

            // compute recall
            float recall = compute_recall(ground_truth_dir, all_knn_results, nvec, k, query_vector.table_id, query_vector.set_id);

            printf("\nrecall=%f\n", recall);

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
            for (int knn = 0; knn < (k * nvec); knn++)
              free(all_knn_results[knn].vector_id);
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

  free(query_vector.values);
  free(query_vector_reordered);
  free(query_order);
  free(results_dir);

  if (track_bsf) {
    for (unsigned int i = 0; i < k; ++i) {
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j) {
        if (bsf_snapshots[i][j].series != NULL)
          free(bsf_snapshots[i][j].series);

        if (bsf_snapshots[i][j].vector_id != NULL)
          free(bsf_snapshots[i][j].vector_id);
      }
      free(bsf_snapshots[i]);
    }
    free(bsf_snapshots);
  }

  return SUCCESS;
}

// parallel incremental query answering read the whole query column and submit all query vectors to teh search engine
enum response dstree_parallel_incr_knn_query_multiple_binary_files(
    const char *bin_files_directory, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top,
    struct dstree_index *index, float minimum_distance, ts_type epsilon,
    ts_type r_delta, unsigned int k, boolean track_bsf, boolean track_pruning,
    boolean all_mindists, boolean max_policy, unsigned int nprobes,
    unsigned char incremental, char *result_dir, unsigned int total_data_files,
    unsigned int dlsize, // total disk size of data files indexed in dstree
    float warping, unsigned char keyword_search, char * k_values_str, char * ground_truth_dir) {

  // worker thread variables, create worker thread
  #define THREAD_NUMS 2
  pthread_barrier_t knn_update_barrier;
  pthread_t worker_thread; 
  char finished = 0; 

  struct bsf_snapshot **bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
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
  int opened_files = 0, qvectors_loaded = 0, curr_vector = 0;
  
  unsigned int offset = 0;

  // open source dir
  struct dirent *dfile;
  DIR *dir = opendir(bin_files_directory);

  if (!dir) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Unable to open directory stream! %s", bin_files_directory);
    exit(1);
  }

  // query time (for cuurent query column)
  double total_query_time = 0.0;
  unsigned int total_checked_ts = 0;
  unsigned int total_queries = qset_num;
  bool found_query = false; // throw error if no query set was found


  // create experiment result directory
  char *results_dir = make_result_directory(
      result_dir, total_data_files, qset_num, min_qset_size, max_qset_size);
  

  // initialize list of all knn results (from all query vectors in query set)
  struct query_result **all_knn_results = NULL;
  struct query_result *curr_knn = NULL;
  ts_type ** query_vectors = NULL; // query column
  struct vid * query_id_arr = NULL;
  ts_type ** query_vectors_reordered = NULL;
  int **query_order_arr = NULL;

  int knn_array_idx = 0;
  // max k values for which we must save results
  unsigned int max_k = k_values[num_k_values - 1];

  // RESET_PARTIAL_COUNTERS()
  // COUNT_PARTIAL_TIME_START

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
      // COUNT_PARTIAL_RAND_INPUT
      // COUNT_PARTIAL_INPUT_TIME_START
      FILE *bin_file = fopen(bin_file_path, "rb");
      // COUNT_PARTIAL_INPUT_TIME_END

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

      while (total_bytes) // counts 4 bytes as one because every vector is
                          // stored in 4 bytes
      {
        // beginning of a set of vectors
        if (i == 0) {
          if (qset_num == 0)
            break;

          // read first integer to check how many vactors in current set
          // COUNT_PARTIAL_INPUT_TIME_START
          fread(&nvec, sizeof(nvec), 1, bin_file);
          // COUNT_PARTIAL_INPUT_TIME_END

          total_bytes--;
          // query set does not fit requirments move to next set
          if ((unsigned int)nvec < min_qset_size ||
              (unsigned int)nvec > max_qset_size) {
            // COUNT_PARTIAL_INPUT_TIME_START
            fseek(bin_file, nvec * 4 * vector_length, SEEK_CUR);
            // COUNT_PARTIAL_INPUT_TIME_END
            i = 0;
            j = 0;
            total_bytes -= (nvec * vector_length);
            nvec = 0u;
            continue;
          }

          found_query = true;
          total_query_time = 0.0;
          total_checked_ts = 0;
          qset_num--;

          // allocate memory for all query vectors
          query_vectors = malloc(sizeof(ts_type *) * nvec);
          query_id_arr = malloc(sizeof(struct vid) * nvec);
          query_vectors_reordered = malloc(sizeof(ts_type *) * nvec);
          query_order_arr = malloc(sizeof(int *) * nvec);
          
          if(query_vectors == NULL || query_id_arr == NULL|| query_vectors_reordered == NULL || query_order_arr == NULL)
          {
            fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for query vectors.");
            exit(1);
          }
          // set new set id
          for(int v = 0; v < nvec; v++)
          {
            query_vectors[v] = calloc(1, sizeof(ts_type) * vector_length);
            query_vectors_reordered[v] = calloc(1, sizeof(ts_type) * vector_length);
            query_order_arr[v] = calloc(1, sizeof(int) * vector_length);
            
            if(query_vectors_reordered[v] == NULL || query_order_arr[v] == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for query vectors.");
              exit(1);
            }

            query_id_arr[v].table_id = table_id;
            query_id_arr[v].set_id = set_id;
            query_id_arr[v].pos = v;
          }

          all_knn_results = malloc(k * nvec * sizeof(struct query_result *));
          if(all_knn_results == NULL)
          {
            fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for knn results.");
            exit(1);
          }

          for (int q = 0, i = 1; q < nvec; ++q) 
          { 
            all_knn_results[q] = calloc(k, sizeof(struct query_result));
            if(all_knn_results[q]  == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for parallel iqa results.");
              exit(1);
            }
            for(int idx = 0; idx < k; idx++)
            {
              all_knn_results[q][idx].node = NULL;
              all_knn_results[q][idx].distance = FLT_MAX;
              all_knn_results[q][idx].vector_id = malloc(sizeof(struct vid));
              if(all_knn_results[q][idx].vector_id  == NULL)
              {
                fprintf(stderr, "Error in dstree_file_loaders.c: Couldn't allocate memory for parallel iqa result ids.");
                exit(1);
              }
              all_knn_results[q][idx].vector_id->table_id = -1;
              all_knn_results[q][idx].vector_id->set_id = -1;
              all_knn_results[q][idx].vector_id->pos = -1;
              all_knn_results[q][idx].query_vector_pos = query_id_arr[q].pos;
              all_knn_results[q][idx].time = 0;
              all_knn_results[q][idx].num_checked_vectors = 0;
            }
          }

          curr_vector = 0;
          knn_array_idx = 0;
          // RESET_QUERY_COUNTERS()
          
          printf("\nQuery (%d, %d) ...\n", table_id, set_id);

          set_id += 1;
          i++;
          j = 0;
        } else if (i <= (unsigned int)nvec * vector_length) {
          // end of vector but still in current set
          if (j > (vector_length - 1)) {
            j = 0;
            
            // add vector
            reorder_query(query_vectors[curr_vector], query_vectors_reordered[curr_vector],
                          query_order_arr[curr_vector], vector_length);
            
            curr_vector += 1;
            qvectors_loaded += 1;
            
          }

          // COUNT_PARTIAL_SEQ_INPUT
          // COUNT_PARTIAL_INPUT_TIME_START
          fread((void *)(&val), sizeof(val), 1, bin_file);
          // COUNT_PARTIAL_INPUT_TIME_END

          total_bytes--;
          query_vectors[curr_vector][j]= val;

          // end of last vector in current  set
          if (i == (unsigned int)nvec * vector_length) {
            // last vector, end of query set (column)
            reorder_query(query_vectors[curr_vector], query_vectors_reordered[curr_vector],
                          query_order_arr[curr_vector], vector_length);
            
            qvectors_loaded += 1;
            curr_vector = 0;

            // (todo) perform extact parallel incremental knn search in 
            // setup param
            struct worker_param * param = malloc(sizeof(struct worker_param));
            if(param == NULL)
            {
              fprintf(stderr, "Error in dstree_file_loaders.c: Could't allocate memory for worker thread param.");
              exit(1);
            }
            param->bsf_snapshots_arr = bsf_snapshots;
            param->cur_bsf_snapshot_arr = cur_bsf_snapshot;
            param->dataset_file_arr = dataset_file; 
            param->epsilon = epsilon;
            param->index = index;
            param->k = k;
            param->minimum_distance = minimum_distance;
            param->num_query_vectors = nvec;
            param->offset = offset;
            param->qfilename_arr = NULL;
            param->query_id_arr = query_id_arr;
            param->query_order_arr = query_order_arr;
            param->query_ts_arr = query_vectors;
            param->query_ts_reordered_arr = query_vectors_reordered;
            param->r_delta = r_delta;
            param->total_checked_ts = &total_checked_ts;
            param->total_query_set_time = &total_query_time;
            param->warping = warping;
            param->global_knn_results = all_knn_results;
            param->finished = &finished;
            param->knn_update_barrier = &knn_update_barrier;

            // create worker thread
            pthread_create(&worker_thread, NULL, exact_de_parallel_incr_knn_search, (void *)param);
            // initialize barrier
            pthread_barrier_init(&knn_update_barrier, NULL, THREAD_NUMS);

            while(finished != 1)
            {
              printf("coordinator_thread:\t (Zzz)\twaiting...\n");
              pthread_barrier_wait(&knn_update_barrier);
              printf("coordinator_thread:\t (!)\treceived new knn results\n");
              // print approx results
              for(int q = 0; q < nvec; q++)
              {
                printf("knns results for vector: (%u, %u, %u)\n", query_id_arr[q].table_id, query_id_arr[q].set_id, query_id_arr[q].pos);
                for(int x = 0; x < k; x++)
                {
                  if(all_knn_results[q][x].distance != FLT_MAX)
                    printf("%dnn: v = (%u, %u, %u), d = %f\n", x, all_knn_results[q][x].vector_id->table_id, all_knn_results[q][x].vector_id->set_id, all_knn_results[q][x].vector_id->pos, all_knn_results[q][x].distance);
                  else
                    printf("%dnn: v = (-, -, -), d = FLT_MAX\n", x);

                }
              }
            }
            pthread_join(worker_thread, NULL); 

            // (todo) store query results
            // (todo) get query time
            // query_time /= 1000000;
            // printf("\nquery_time=%fsec\n", query_time);


            // (todo) compute recall
            // float recall = compute_recall(ground_truth_dir, all_knn_results, nvec, k, query_vector.table_id, query_vector.set_id);
            // printf("\nrecall=%f\n", recall);


            // free memory
            // free query vectors
            free(param);
            pthread_barrier_destroy(&knn_update_barrier);

            for(int v = 0; v < nvec; v++)
            {
              free(query_vectors[v]);
              free(query_vectors_reordered[v]);
              free(query_order_arr[v]);
            }
            free(query_id_arr);
            free(query_vectors);
            free(query_vectors_reordered);
            free(query_order_arr);

            query_id_arr = NULL;
            query_vectors = NULL;
            query_vectors_reordered = NULL;
            query_order_arr = NULL;

            // free knn results
            for (int q = 0, i = 1; q < nvec; ++q) 
            { 
              for(int idx = 0; idx < k; idx++)
              {
                free(all_knn_results[q][idx].vector_id);
              }
              free(all_knn_results[q]);
            }
            free(all_knn_results);
            

            // RESET_PARTIAL_COUNTERS()
            // COUNT_PARTIAL_TIME_START

            i = 0;
            j = 0;
            nvec = 0u;
            continue;
          }
          i++;
          j++;
        }
      }
      // COUNT_PARTIAL_INPUT_TIME_START
      fclose(bin_file);
      // COUNT_PARTIAL_INPUT_TIME_END
      
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
  // COUNT_INPUT_TIME_START
  closedir(dir);
  // COUNT_INPUT_TIME_END
  free(results_dir);
  free(k_values);
  // COUNT_PARTIAL_TIME_END
  // RESET_PARTIAL_COUNTERS()  
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

/* start kashif changes */
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
