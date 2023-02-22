/*SCRIPT TO CONVERT ONE BINARY FILE TO TABLE, SET AND VECTOR STRUCTURES*/
#include <dirent.h>
#include <libgen.h>
#include <linux/limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "globals.h"
#include <stdatomic.h>

#define NBITS 64
typedef float ts_type;

typedef struct vector {
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos;
  ts_type *values;
} vector;

struct pool_queue {
  struct vid *query_id;
  int *query_order;
  ts_type * query_vector;
  ts_type *query_vector_reordered;

	char free;
	struct pool_queue *next;
} pool_queue;

struct thread_param {
    struct pool *thread_pool;
    struct worker_param * work_param;
    uint8_t thread_id;
} thread_para;


struct performance {
  float recall;
  float precision;
} performance;
// struct pool {
// 	char cancelled;
// 	unsigned int remaining;
//   int8_t * thread_status; // 0 is thread is not working 1 if thread has finished current job 2 if thread finishied last job (no more jobs to take)
// 	unsigned int num_threads;

// 	struct pool_queue * query_queue;
// 	struct pool_queue *end_queue;

//   pthread_mutex_t * status_lock;
//   pthread_mutex_t queue_lock;
// 	pthread_cond_t queue_condition;

// 	pthread_t *threads;
//   unsigned int *task_count;
//   void *(*function)(void *);

//   struct worker_param *params;
// } pool;

// create result file name and path.
char *make_file_path(char *result_dir, unsigned int qtable_id,
                     unsigned int qset_id, unsigned int qsize, unsigned int l,
                     unsigned int dlsize, unsigned int vector_length, unsigned int k);

unsigned int * get_k_values(char * k_values_str, unsigned int * num_k_values);

// save query results to csv file
enum response save_to_query_result_file(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int max_k);

// store query results in disk, for multi thread knn search
enum response store_knn_results_in_disk(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int k);

// store query results in disk, for multi thread knn search (correct way)
enum response store_parallel_knn_results_in_disk(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int k, double *max_thread_time);

// create experiment results dir
char *make_result_directory(char *result_dir, unsigned int total_data_files,
                            unsigned int nq, unsigned int min_qset_size,
                            unsigned int max_qset_size);

// count digits in integer
int get_ndigits(unsigned int n);

/* convert IEEE double precision (64 bits) to float (/double) */
ts_type todecimal(char *s);
int toint(char *n);
char *get_set_bytes(FILE *f, int nvec, int start, int nbits, int vlen);
bool is_binaryfile(const char *filename);
int get_num_vec(FILE *f, int start, int nbits);
unsigned long get_total_data_vectors(char *bindir,
                unsigned int total_data_files, unsigned int * total_columns);
unsigned int get_data_gb_size(char *dl_dir, unsigned int total_data_files);
int get_ndigits(unsigned int n);
int delete_directory(char * dir_path);
struct result_sid *get_top_sets(struct result_vid **knn_results, int num_query_vectors, unsigned int k, 
                        int num_top, unsigned int * total_matching_columns);

struct result_table *get_top_tables_by_euclidean_distance(struct query_result *knn_results, int num_knn_results, 
                         unsigned int num_top);
int get_ground_truth_file(char * ground_truth_dir, int query_table_id, int query_set_id, char * ground_truth_file);
struct result_vid * get_ground_truth_results(char * ground_truth_file, int total_results);

struct performance compute_recall_precision(struct result_vid * ground_truth_results, int num_gt_results, struct result_vid ** knn_results, int k, int q, int query_table_id, int query_set_id);

float compute_recall(char * ground_truth_dir, struct result_vid ** knn_results, int num_query_vectors, int k, int query_table_id, int query_set_id);
void compute_one_query_vector_recall(struct result_vid * ground_truth_results, int  num_gt_results, struct query_result * knn_results, 
                              int first_nn, int last_nn, int8_t * recall_row);

float compute_recall_from_matrix(int8_t ** recall_matrix, int num_query_vectors, int k, int num_ground_truth_results);
float compute_vector_recall(struct result_vid * knn_results, struct result_vid * gt_results, int query_vector_pos, int k);
float compute_vector_recall_identical_results(int8_t ** recall_matrix, int query_vector_pos, int k, int num_identical_neighbors);

unsigned int get_num_identical_results(struct result_vid * ground_truth_results, int  num_gt_results, unsigned int query_vector_pos);

//  count elements in queue
unsigned int get_queue_size(struct pool_queue * queue, unsigned int thread_id);

// start thread in thread pool
static void * start_thread(void *arg);

// initialize thread pool
void init_thread_pool(struct pool* pool, struct dstree_index * index, ts_type epsilon, unsigned int k,
              void * (*thread_func)(void *), unsigned int num_threads, unsigned int offset,
              ts_type r_delta, unsigned int * total_checked_ts, double * total_query_time, 
              float warping, struct result_vid ** all_knn_results, unsigned char store_results_in_disk,
              unsigned int *k_values, unsigned int num_k_values, struct result_vid *ground_truth_results,
              unsigned int num_gt_results, int8_t * recall_matrix, unsigned int num_query_vectors,
              struct job * job_array, unsigned int vector_length, unsigned int stop_when_nn_dist_changes);

// fill job queue of the thread pool
struct pool_queue * pool_init_job_queue(struct vid *query_id_arr, int ** query_order_arr,
                    ts_type ** query_vectors, ts_type ** query_vectors_reordered,
                    unsigned int num_query_vectors, unsigned int vector_lengt);

// end pool of threads wait for all threads to finish.
void pool_end(struct pool *pool);

ts_type todecimal(char *s) {
  // printf("\ninput: %s\n", s);
  ts_type f;
  int sign, exp = 0;
  unsigned long long mant = 0;
  int i;

  sign = s[0] - '0';

  for (i = 1; i <= 11; i++)
    exp = exp * 2 + (s[i] - '0');

  if (exp > -1022) {
    mant = 1;
    exp -= 1023;
  } else {
    // Subnormal numbers
    mant = 0;
    exp = -1022;
  }
  for (i = 12; i < 64; i++) {
    mant = mant * 2 + (s[i] - '0');
  }
  /*

  printf("s = %d\n", sign);
  printf("exp = %d\n", exp);
  printf("m = %llu\n", mant);

  */

  f = pow(-1, sign) * pow(2, exp) * (mant / pow(2, 52));
  return f;
}

/* CONVERT BIN TO INT */
int toint(char *n) {
  int integer = 0;
  for (int i = 0; i < strlen(n); i++)
    integer = integer * 2 + (n[i] - '0');
  return integer;
}

bool is_binaryfile(const char *filename) {
  // check if filename has bin extesion.
  char *ext = ".bin";
  size_t nl = strlen(filename), el = strlen(ext);
  return nl >= el && !strcmp(filename + nl - el, ext);
}

// extract set bytes from bin file
char *get_set_bytes(FILE *f, int nvec, int start, int nbits, int vlen) {
  start += nbits;                  // skip bytes for num vectors
  int total = nvec * vlen * nbits; // total bits occupied by set values

  char *bytes = (char *)malloc(sizeof(char) * total);
  fseek(f, start, SEEK_SET);
  int i;
  for (i = 0; i < total; i++) {
    bytes[i] = fgetc(f);
  }
  bytes[i] = '\0';
  // //printf("\n%s\n", bytes);
  return bytes;
}

// get number of vectors in the current set (set that starts at 'start')
int get_num_vec(FILE *f, int start, int nbits) {
  char *bytes = (char *)malloc(sizeof(char) * nbits);
  int i;
  fseek(f, start, SEEK_SET);
  for (i = 0; i < nbits; i++) {
    bytes[i] = fgetc(f);
  }
  bytes[i] = '\0';
  // //printf("\n%s\n", bytes);
  return toint(bytes);
}

// read datasize from diffrent files and get total datasize of a datalake
unsigned long get_total_data_vectors(char *bindir,
                unsigned int total_data_files, unsigned int * total_columns) {
  struct dirent *dfile;
  DIR *dir = opendir(bindir);
  unsigned long total_vectors = 0;
  unsigned int datasize;
  unsigned int ncols;
  unsigned int table_id;
  if (!dir) {
    fprintf("Error in kashif_utils.c: Unable to open directory stream! %s", bindir);
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0) {
    // // skip directories
    // if (dfile->d_type != DT_REG)
    //   continue;
    if (is_binaryfile(dfile->d_name)) {
      total_data_files--;
      // ** get binary table info
      sscanf(dfile->d_name, "data_size%d_t%dc%d%*[^0123456789]", &datasize, &table_id, &ncols);
      total_vectors += datasize;
      *total_columns += ncols;
    }
  }
  closedir(dir);
  return total_vectors;
}

// get data lake size in GB
unsigned int get_data_gb_size(char *dl_dir, unsigned int total_data_files) {
  struct dirent *dfile;
  DIR *dir = opendir(dl_dir);
  float total_dlsize = 0.0;
  FILE *fp;

  if (!dir) {
    fprintf("Error in kashif_utils.c: Unable to open directory stream! %s", dl_dir);
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0) {
    // // skip directories
    // if (dfile->d_type != DT_REG)
    //   continue;

    if (is_binaryfile(dfile->d_name)) {
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, dl_dir);
      strcat(bin_file_path, "/");
      strcat(bin_file_path, dfile->d_name);
      total_data_files--;
      // ** get binary file size
      fp = fopen(bin_file_path, "r");
      if (fp == NULL) {
        printf("Could not get size of binary file %s. File not found!\n",
               dfile->d_name);
        exit(1);
      }
      fseek(fp, 0L, SEEK_END);
      total_dlsize += ftell(fp);
      fclose(fp);
    }
  }
  closedir(dir);
  return (unsigned int)round(total_dlsize /
                             1073741824); // return data size in GB
}

// extract k values from a string. ex: "1,3,5,10,30,50" to [1, 2, 3, 10, 30, 50]
unsigned int *  get_k_values(char * k_values_str, unsigned int * num_k_values)
{
  unsigned int k;
  unsigned int * k_values = NULL;

  char *pt;
  pt = strtok (k_values_str,",");
  while (pt != NULL) {
      k = atoi(pt);
      printf("k = %u - ", k);

      if(k > 0)
      {
        (*num_k_values) += 1;
        k_values = (unsigned int *) realloc(k_values, sizeof(*k_values) * (*num_k_values));
        if (k_values == NULL)
        {
          fprintf(stderr,"Error kashif_utils.c:  Could not allocate memory for set of k values.\n");
            return -1;
        }
        k_values[(*num_k_values) -1] = k;
      }
      pt = strtok (NULL, ",");
  }

  printf("\n");
  return k_values;
}

// count digits in integer
int get_ndigits(unsigned int n) {
  int total_digits = 0;
  while (n != 0) {
    n = n / 10;
    total_digits++;
  }
  return total_digits;
}

// create results dir
char *make_file_path(char *result_dir, unsigned int qtable_id,
                     unsigned int qset_id, unsigned int qsize,
                     unsigned int total_data_files, unsigned int dlsize,
                     unsigned int vector_length, unsigned int k) {
  COUNT_INPUT_TIME_START
  DIR *dir = opendir(result_dir);
  COUNT_INPUT_TIME_END
  if (!dir) {
    printf("WARNING! Experiment direstory '%s' does not exist!", result_dir);
    exit(1);
  }
  char *filepath = malloc(
      get_ndigits(qtable_id) + get_ndigits(qset_id) +
      get_ndigits(total_data_files) + get_ndigits(dlsize) +
      get_ndigits(vector_length) + get_ndigits(qsize) + get_ndigits(k) +
      strlen("TQ_Q_qsize_l_dlsize_len_k.csv") +
      strlen(result_dir));

  sprintf(filepath,
          "%s/"
          "TQ%u_Q%u_qsize%u_l%u_dlsize%u_len%u_k%u\0",
          result_dir, qtable_id, qset_id, qsize, total_data_files, dlsize,
          vector_length, k);


  COUNT_INPUT_TIME_START
  closedir(dir);
  COUNT_INPUT_TIME_END

  return filepath;
}

// save query results to csv file
enum response save_to_query_result_file(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int max_k) {
  FILE *fp;
  int i, j;
  COUNT_OUTPUT_TIME_START
  fp = fopen(csv_file, "w+");
  COUNT_OUTPUT_TIME_END
  if (fp == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n",
            csv_file);
    return FAILURE;
  }

  // write header
  COUNT_INPUT_TIME_START
  fprintf(fp, "TQ:Q, TS:S, q_pos, s_pos, q, s, d, time, k");
  double total_querytime = 0;
  unsigned int total_checked_vec = 0;
 
  // write results for a specific k value
  for (int q = 0; q < nqvec; q++)// vector counter
  {
    // printf("Collect values between %d and %d\n", i, (max_k + i));
    for(int s = 0; s < max_k; s++)// k counter
    {
      // printf("k = %d, time = %f -- %f\n", s, knn_results[s].time, knn_results[s].time/1000000);
      fprintf(fp, "\n");
      fprintf(fp, "%u:%u, %u:%u, %u, %u, [], [], %f, %.7f, %u", qtable_id, qset_id,
            knn_results[q][s].table_id, knn_results[q][s].set_id,
            knn_results[q][s].qpos, knn_results[q][s].pos, knn_results[q][s].distance, knn_results[q][s].time/1000000, max_k);
    
    
    // total_checked_vec += knn_results[s].num_checked_vectors;
    }
    total_querytime += knn_results[q][max_k-1].time;
  }
  fclose(fp);
  COUNT_OUTPUT_TIME_END

  total_checked_vec = 0;

  // add query time to file name and rename csv file
  char * new_csv_filename =  malloc(strlen(csv_file) + strlen("_runtime_ndistcalc_dataaccess.csv") + 20 + 1);
  sprintf(new_csv_filename, "%s_runtime%.4f_ndistcalc_dataaccess%u.csv\0", csv_file,  total_querytime/1000000, total_checked_vec);
  
  printf("%u, \t%f\n", max_k, total_querytime/1000000);
  int ret = rename(csv_file, new_csv_filename);
  
  free(new_csv_filename);
  return SUCCESS;
}

// store query results in disk, for multi thread knn search
enum response store_knn_results_in_disk(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int k)
{
  FILE *fp;
  int i, j;
  COUNT_OUTPUT_TIME_START
  fp = fopen(csv_file, "w+");
  COUNT_OUTPUT_TIME_END
  if (fp == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n",
            csv_file);
    return FAILURE;
  }

  // write header
  fprintf(fp, "TQ:Q, TS:S, q_pos, s_pos, q, s, d, time, k");
  double total_querytime = 0;
  unsigned int total_checked_vec = 0;
  double max_time = -1000;
  int max_time_query_vector = -1;

  // write results for a specific k value
  for (int q = 0; q < nqvec; q++)// vector counter
  {
    if(knn_results[q][k-1].time > max_time) // max total query time (total query time is time for the kth result)
    {
      max_time = knn_results[q][k-1].time;
      max_time_query_vector = q;
    }

    for(int s = 0; s < k; s++)// k counter
    {
      fprintf(fp, "\n");
      fprintf(fp, "%u:%u, %u:%u, %u, %u, [], [], %f, %.7f, %u", qtable_id, qset_id,
            knn_results[q][s].table_id, knn_results[q][s].set_id, knn_results[q][s].distance,
            knn_results[q][s].qpos, knn_results[q][s].pos, knn_results[q][s].time/1000000, k);
    total_checked_vec += knn_results[q][s].num_checked_vectors;
    }
    
  }
  
  total_querytime = max_time;
  fclose(fp);

  // add query time to file name and rename csv file
  
  char * new_csv_filename =  malloc(strlen(csv_file) + strlen("_runtime_ndistcalc_dataaccess.csv") + 50 + 1);
  sprintf(new_csv_filename, "%s_runtime%.4f_ndistcalc_dataaccess%u.csv\0", csv_file,  total_querytime/1000000, total_checked_vec);
  
  // printf("[k = %u] Combined total query time  = %f\n", k, total_querytime/1000000);
  int ret = rename(csv_file, new_csv_filename);
  
  free(new_csv_filename);
  return SUCCESS;
}

// store query results in disk, for multi thread knn search
enum response store_parallel_knn_results_in_disk(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int nqvec,
                                        struct result_vid **knn_results, unsigned int k, double *max_thread_time)
{
  FILE *fp;
  int i, j;
  COUNT_OUTPUT_TIME_START
  fp = fopen(csv_file, "w+");
  COUNT_OUTPUT_TIME_END
  if (fp == NULL) 
  {
    fprintf(stderr, "Error in kashif_utils.h: Could not open file %s, reason = %s!\n", csv_file, strerror(errno));
    return FAILURE;
  }

  // write header
  fprintf(fp, "TQ:Q, TS:S, q_pos, s_pos, q, s, d, time, k");
  unsigned int total_checked_vec = 0;
  
  // write results for a specific k value
  for (int q = 0; q < nqvec; q++)// vector counter
  {
    for(int s = 0; s < k; s++)// k counter
    {
      fprintf(fp, "\n");
      fprintf(fp, "%u:%u, %u:%u, %u, %u, [], [], %f, %.7f, %u", qtable_id, qset_id,
            knn_results[q][s].table_id, knn_results[q][s].set_id, knn_results[q][s].distance,
            knn_results[q][s].qpos, knn_results[q][s].pos, knn_results[q][s].time/1000000, k);
    total_checked_vec += knn_results[q][s].num_checked_vectors;
    }
    
  }
  
  fclose(fp);

  // add query time to file name and rename csv file
  char * new_csv_filename =  malloc(strlen(csv_file) + strlen("_runtime_ndistcalc_dataaccess.csv") + 20 + 1);
  sprintf(new_csv_filename, "%s_runtime%.4f_ndistcalc_dataaccess%u.csv\0", csv_file,  max_thread_time[k-1]/1000000, total_checked_vec);
  
  // printf("[k = %u] Combined total query time  = %f\n", k, max_thread_time[k-1]/1000000);
  int ret = rename(csv_file, new_csv_filename);
  
  free(new_csv_filename);
  return SUCCESS;
}

// delete non empty direcory
int delete_directory(char * dir_path)
{
    int ret = 0;
    ret = rmdir(dir_path);

    if (ret == 0)
        return 0;
    return -1;
}

// create results dir
char *make_result_directory(char *result_dir, unsigned int l, unsigned int nq,
                            unsigned int min_qset_size,
                            unsigned int max_qset_size) {
  char *result_dir_name =
      malloc(get_ndigits(l) + get_ndigits(nq) + get_ndigits(min_qset_size) +
             get_ndigits(max_qset_size) + strlen("/kashif_l_q_min_max") +
             strlen(result_dir) + 1);

  sprintf(result_dir_name, "%s/kashif_l%u_%uq_min%u_max%u", result_dir, l, nq,
          min_qset_size, max_qset_size);

  printf(">>> Result directory name: %s\n", result_dir_name);
  COUNT_OUTPUT_TIME_START
  DIR *dir = opendir(result_dir_name);
  if (dir) {
    delete_directory(result_dir_name);
  }
  mkdir(result_dir_name, 0777);
  closedir(dir);
  COUNT_OUTPUT_TIME_END
  
  return result_dir_name;
}

// get sets with the the largest number of matching vectors with the query (highest ovelap size)
struct result_sid *get_top_sets(struct result_vid **knn_results, int num_query_vectors, unsigned int k, 
                        int num_top, unsigned int *total_matching_columns)
{
  // frequency = number of matching vectors with the query set vectors
  int i, j;

  struct result_sid * distinct_sets = NULL;
  struct vid * distinct_match_vectors = NULL;

  *total_matching_columns = 0;
  unsigned int num_distinct_sets = 0;
  unsigned int num_distinct_vectors = 0;
  

  int found = 0;
  int new_column = -1, new_vector = -1;

  num_distinct_sets += 1;
  distinct_sets = realloc(distinct_sets, sizeof(struct result_sid));
  distinct_sets[0].table_id = knn_results[0][0].table_id;
  distinct_sets[0].set_id = knn_results[0][0].set_id;
  // strcpy(distinct_sets[0].raw_data_file, knn_results[0][0].raw_data_file);
  distinct_sets[0].overlap_size = 0;

  num_distinct_vectors += 1;
  distinct_match_vectors = realloc(distinct_match_vectors, sizeof(struct vid));
  distinct_match_vectors[0].table_id = knn_results[0][0].table_id;
  distinct_match_vectors[0].set_id = knn_results[0][0].set_id;
  distinct_match_vectors[0].pos = knn_results[0][0].pos;


  // get distinct sets & vectors
  for(int a = 0; a < num_query_vectors; a++)
  {
    // float first_dist = knn_results[a][0].distance;
    for(int b = 0; b < k; b++)
    {
      // // stop when distance changes
      // if(knn_results[a][b].distance > first_dist)
      //   break;

      new_column = 1;
      new_vector = 1;
      for(j = num_distinct_sets - 1; j >= 0; j--)
      {
        if(distinct_sets[j].table_id == knn_results[a][b].table_id)
          if(distinct_sets[j].set_id == knn_results[a][b].set_id)
          {
            new_column = 0;
            break;
          }
      }
      for(j = num_distinct_vectors - 1; j >= 0; j--)
      {
        if(distinct_match_vectors[j].table_id == knn_results[a][b].table_id)
          if(distinct_match_vectors[j].set_id == knn_results[a][b].set_id)
          if(distinct_match_vectors[j].pos == knn_results[a][b].pos)
          {
            new_vector = 0;
            break;
          }
      }
      // get distinct columns (sets)
      if(new_column == 1)
      {
        
        distinct_sets = realloc(distinct_sets, (sizeof(struct result_sid) * (num_distinct_sets+1)));
        distinct_sets[num_distinct_sets].table_id = knn_results[a][b].table_id;
        distinct_sets[num_distinct_sets].set_id = knn_results[a][b].set_id;
        distinct_sets[num_distinct_sets].overlap_size = 0;
        // strcpy(distinct_sets[num_distinct_sets].raw_data_file, knn_results[a][b].raw_data_file);
        num_distinct_sets += 1;
      }
      // get distinct vectors
      if(new_vector == 1)
      {
        
        distinct_match_vectors = realloc(distinct_match_vectors, (sizeof(struct vid) * (num_distinct_vectors+1)));
        distinct_match_vectors[num_distinct_vectors].table_id = knn_results[a][b].table_id;
        distinct_match_vectors[num_distinct_vectors].set_id = knn_results[a][b].set_id;
        distinct_match_vectors[num_distinct_vectors].pos = knn_results[a][b].pos;
        num_distinct_vectors += 1;
      }
    }
  }
  

  // print distinct column ids: 

  // printf("all knn results : \n ");
  // for(int a = 0; a < num_query_vectors; a++)
  // {
  //   for(int b = 0; b < k; b++)
  //     printf("%d, %d: (%u, %u, %u)\n", a, b+1, knn_results[a][b].table_id, knn_results[a][b].set_id, knn_results[a][b].pos);
  //   printf("\n");
  // }

  
  // measure column overlap with the query column
  for(i = num_distinct_sets - 1; i >= 0; i--)
  {
    for(j = num_distinct_vectors - 1; j >= 0; j--)
    {
      if(distinct_sets[i].table_id == distinct_match_vectors[j].table_id)
        if(distinct_sets[i].set_id == distinct_match_vectors[j].set_id)
        {
          distinct_sets[i].overlap_size += 1;
        }
    }
  }

  // fill the rest of 'num_top' with the last elemet
  if((num_top != -1) && (num_distinct_sets < num_top))
  {
    distinct_sets = realloc(distinct_sets, (sizeof(struct result_sid) * num_top));
    struct result_sid * last = &distinct_sets[num_distinct_sets - 1];

    for(i = num_distinct_sets; i < num_top; i++)
    {
      distinct_sets[i].table_id = last->table_id;
      distinct_sets[i].set_id = last->set_id;
      distinct_sets[i].overlap_size = last->overlap_size;
      strcpy(distinct_sets[i].raw_data_file, last->raw_data_file);
    }
  }

  free(distinct_match_vectors);

  *total_matching_columns = num_distinct_sets;
  return distinct_sets;
}

struct result_table *get_top_tables_by_euclidean_distance(struct query_result *knn_results, int num_knn_results, 
                         unsigned int num_top)
{
  int i, j, k;

  struct result_table * distinct_tables = NULL;
  unsigned int num_distinct_tables = 0;


  num_distinct_tables += 1;
  distinct_tables = realloc(distinct_tables, sizeof(struct result_table));
  distinct_tables[0].table_id = knn_results[0].vector_id->table_id;
  distinct_tables[0].min_distance = FLT_MAX;
  distinct_tables[0].num_min = 0;
  distinct_tables[0].total_matches = 0;
  strcpy(distinct_tables[0].raw_data_file, knn_results[0].vector_id->raw_data_file);

  int found;
  // get distinct sets
  for ( i = 1; i < num_knn_results; i++)
  {
    found = 0;
    for(j = num_distinct_tables - 1; j >= 0; j--)
    {
      if(distinct_tables[j].table_id == knn_results[i].vector_id->table_id)
      {
        found = 1;
        break;
      }
    }
    if(found == 0)
    {
      num_distinct_tables += 1;
      distinct_tables = realloc(distinct_tables, (sizeof(struct result_table) * num_distinct_tables));
      distinct_tables[num_distinct_tables - 1].table_id = knn_results[i].vector_id->table_id;
      distinct_tables[num_distinct_tables - 1].min_distance = FLT_MAX;
      distinct_tables[num_distinct_tables - 1].num_min = 0;
      distinct_tables[num_distinct_tables - 1].total_matches = 0;
      strcpy(distinct_tables[num_distinct_tables - 1].raw_data_file, knn_results[i].vector_id->raw_data_file);
    }
  }

  // get distinct vectors
  for(j = num_distinct_tables - 1; j >= 0; j--)
  {
    for (i = 0; i < num_knn_results; i++)
    {
      if(knn_results[i].vector_id->table_id == distinct_tables[j].table_id)
      {
        distinct_tables[j].total_matches += 1;
        if(knn_results[i].distance < distinct_tables[j].min_distance)
        {
          distinct_tables[j].min_distance = knn_results[i].distance;
          distinct_tables[j].num_min = 1;
        }
        else if (knn_results[i].distance == distinct_tables[j].min_distance)
        {
          distinct_tables[j].num_min += 1;
        }
      }
    }
  }

  // fill the rest withe the last element
  if(num_distinct_tables < num_top)
  {
    distinct_tables = realloc(distinct_tables, (sizeof(struct result_table) * num_top));
    struct result_table * last = &distinct_tables[num_distinct_tables - 1];

    for(i = num_distinct_tables; i < num_top; i++)
    {
      distinct_tables[i].table_id = last->table_id;
      distinct_tables[i].min_distance = last->min_distance;
      distinct_tables[i].num_min = last->num_min;
      distinct_tables[i].total_matches = last->total_matches;
      strcpy(distinct_tables[i].raw_data_file, last->raw_data_file);
    }
  }

  return distinct_tables;
}

int  get_ground_truth_file(char * ground_truth_dir, int query_table_id, int query_set_id, char * ground_truth_file)
{
  char queryid[255] = "";
  char *ret;
  char buffer[80];
  struct dirent *dfile;
  DIR *dir = opendir(ground_truth_dir);

  if (!dir)
  {
    fprintf(stderr, "Error in dstree_file_loaders.c: Unable to open directory stream! %s", ground_truth_dir);
    exit(1);
  }
  snprintf(queryid, 12, "TQ%d_Q%d", query_table_id, query_set_id);

  while ((dfile = readdir(dir)) != NULL)
  {
    ret = strstr(dfile->d_name, queryid);
    if(ret)
    {
      strcat(ground_truth_file, ground_truth_dir);
      strcat(ground_truth_file, "/");
      strcat(ground_truth_file, dfile->d_name);
      break;
    }
  }

  FILE * file = fopen(ground_truth_file, "r");

  if(file == NULL)
  {
    fprintf(stderr, "Error in kashif_utils: Cannot open gt file %s for column %s\n", ground_truth_file, queryid);
    exit(1);
  }
  int num_lines = 0;
  // read file and count number of lines
  while (fgets(buffer, 80, file)) {
      num_lines++;
  }
  closedir(dir);
  return num_lines - 1; // dont count reader
}

struct result_vid * get_ground_truth_results(char * ground_truth_file, int total_results)
{
  char buffer[80];
  FILE * file = fopen(ground_truth_file, "r");
  struct result_vid * ground_truth_results = malloc(sizeof(struct result_vid) * total_results);
  if(ground_truth_results == NULL)
  {
    fprintf(stderr, "Error in kashif_utile.h: Could'nt allocate memory for ground_truth_results");
    return 1;
  }
  

  int i = 0, header = 1;;
  while (fgets(buffer, 80, file)) {
      if(header == 1)
      {
        header = 0;
        continue;
      }

      // If you need all the values in a row
      char *token = strtok(buffer, ",");
      int col = 0;
      int table_id = 0, set_id = 0, vector_pos = 0, query_vector_pos = 0;
      float distance = 0;

      while (token) {
          // Just printing each integer here but handle as needed
          if(col == 1)
          {
              sscanf(token, " %d:%d", &table_id, &set_id);
          }
          else if (col == 2)
              query_vector_pos = atoi(token);
          else if (col == 3)
              vector_pos = atoi(token);
          else if (col == 6)
              distance = atof(token);

          token = strtok(NULL, ",");
          col++;
      }

      // printf("table: %d,  set: %d, vector: %d\n", table_id, set_id, vector_pos);
      ground_truth_results[i].table_id = table_id;
      ground_truth_results[i].set_id = set_id;
      ground_truth_results[i].pos = vector_pos;
      ground_truth_results[i].qpos = query_vector_pos;
      ground_truth_results[i].distance = distance;
      i++;
  }
  fclose(file);
  return ground_truth_results;
}

struct performance compute_recall_precision(struct result_vid * ground_truth_results, int num_gt_results, struct result_vid ** knn_results, int k, int q, int query_table_id, int query_set_id)
{
  struct performance perf;
  float num_matches = 0, kth_dist;
  int tp_fn = 0;


  kth_dist = knn_results[q][k - 1].distance;

  // count tp + fn
  for(int i = 0; i < num_gt_results; i++)
  {
    if(ground_truth_results[i].qpos == q && ground_truth_results[i].distance <= kth_dist)
      tp_fn++;
  }
  
  // count matches
  for(int x = 0; x < k; x++)
    for(int i = 0, j = 0; j < tp_fn && i < num_gt_results; i++)
    {
      if(ground_truth_results[i].qpos == knn_results[q][x].qpos)
      {
        j++;
        if(ground_truth_results[i].table_id == knn_results[q][x].table_id
        && ground_truth_results[i].set_id == knn_results[q][x].set_id
        && ground_truth_results[i].pos == knn_results[q][x].pos)
        {
          num_matches++;
          break;
        }
      }
    }

  // printf("k = %u, kthdistnace = %f, matches = %.1f, tp_fn = %d, #gt = %d\n", k, kth_dist, num_matches, tp_fn, num_gt_results);
  perf.recall = num_matches/(float)tp_fn;
  perf.precision = 0;

  return perf;
}

float compute_recall(char * ground_truth_dir, struct result_vid ** knn_results, int num_query_vectors, int k, int query_table_id, int query_set_id)
{
  char ground_truth_file[255] = "";
  int  num_gt_results = get_ground_truth_file(ground_truth_dir, query_table_id, query_set_id, ground_truth_file);
  struct result_vid * ground_truth_results = get_ground_truth_results(ground_truth_file, num_gt_results);

  int * exact_results_counter = calloc(num_query_vectors, sizeof(int));
  float num_matches = 0;
  int num_knn_results = num_query_vectors * k;

  // printf("total vectors: %d\n", num_query_vectors);
  for(int q = 0; q < num_query_vectors; q++)
  {
    for(int x = 0; x < k; x++)
      for(int i = 0; i < num_gt_results; i++)
      {
        if(ground_truth_results[i].table_id == knn_results[q][x].table_id
          && ground_truth_results[i].set_id == knn_results[q][x].set_id
          && ground_truth_results[i].qpos == knn_results[q][x].qpos
          && ground_truth_results[i].pos == knn_results[q][x].pos)
        {
          num_matches++;
          break;
        }
      }
  }

  for(int i = 0; i < num_gt_results; i++)
  {
    if(exact_results_counter[ground_truth_results[i].qpos] < k)
    {
      // printf("incr at %d\n", ground_truth_results[i].query_vector_pos);
      exact_results_counter[ground_truth_results[i].qpos]++;
    }
  }

  float exact_result_total_rows = 0.0;
  for(int i = 0; i < num_query_vectors; i++)
  {
    exact_result_total_rows += exact_results_counter[i];
  }
  
  // printf("%f / %f\n", num_matches, (float)num_gt_results);
  // printf("precision = %f\n", num_matches/exact_result_total_rows);

  free(exact_results_counter);
  free(ground_truth_results);
  return num_matches/(float)num_gt_results;
}

void compute_one_query_vector_recall(struct result_vid * ground_truth_results, int  num_gt_results, struct query_result * knn_results, 
                              int first_nn, int last_nn, int8_t * recall_row)
{
  int found_matches = 0;
  for(int x = first_nn; x < last_nn; x++) // iterate over the nn returned by incremental search
    for(int i = 0; i < num_gt_results; i++)
    {
      if(ground_truth_results[i].table_id == knn_results[x].vector_id->table_id
        && ground_truth_results[i].set_id == knn_results[x].vector_id->set_id
        && ground_truth_results[i].qpos == knn_results[x].query_vector_pos
        && ground_truth_results[i].pos == knn_results[x].vector_id->pos)
      {
        // update recall matrix
        recall_row[x] = 1;
        found_matches ++;
        break;
      }
    }
}


float compute_recall_from_matrix(int8_t ** recall_matrix, int num_query_vectors, int k, int num_ground_truth_results)
{
  float num_matches = 0.0;
  for(int q = 0; q < num_query_vectors; q++)
    for(int r = 0; r < k; r++)
    {
      num_matches += recall_matrix[q][r];
    }

  return num_matches / num_ground_truth_results;
}

float compute_k_recall_from_matrix(struct result_vid ** knn_results, struct result_vid * gt_results, int num_query_vectors, int k, int num_ground_truth_results)
{
  float num_matches = 0.0;
  for(int q = 0; q < num_query_vectors; q++)
    num_matches += compute_vector_recall(knn_results[q], gt_results, q, k);


  return num_matches / (k * num_query_vectors);
}

float compute_vector_recall(struct result_vid * knn_results, struct result_vid * gt_results, int query_vector_pos, int k)
{
  float num_matches = 0.0;
  int i = k;
  for(int x = 0, j = 0; x < k; x++)
  {
    while(i > 0)
    {
      if(gt_results[j].qpos == query_vector_pos)
      {
        if(gt_results[j].table_id == knn_results[x].table_id
          && gt_results[j].set_id == knn_results[x].set_id
          && gt_results[j].pos == knn_results[x].pos
          )
        {
          // update recall matrix
          num_matches ++;
          i--;
          break;
        }
        i--;
      }
      j++;
    }
  }

  return num_matches;
}

float compute_vector_recall_identical_results(int8_t ** recall_matrix, int query_vector_pos, int k, int num_identical_neighbors)
{
  float num_matches = 0.0;
  for(int r = 0; r < k; r++)
  {
    num_matches += recall_matrix[query_vector_pos][r];
  }

  return num_matches / num_identical_neighbors;
}

// count total identical nn for one query vector
unsigned int get_num_identical_results(struct result_vid * ground_truth_results, int  num_gt_results, unsigned int query_vector_pos)
{
  unsigned int num_identical_results = 0;
  for(int i = 0; i < num_gt_results; i++)
  {
    if(ground_truth_results[i].qpos == query_vector_pos)
    {
      num_identical_results ++;
    }
  }
  return num_identical_results;
}

//  count elements in queue
unsigned int get_queue_size(struct pool_queue * queue, unsigned int thread_id)
{
  if(queue == NULL)
  {
    fprintf(stderr, "Error in kashif_utils.c: cannot measure size of a null pointer to queue, call made by th%d.\n", thread_id);
    exit(1);
  }

  if(queue[0].next == NULL)
  {
    printf("Queue size = *1*\t");
    return 1;
  }
  int i = 0;
  do
  {
      i++;
  } while ((queue[i].next != NULL));
  

  printf("Queue size = %d\t", i+1);
  return i+1;
}

// start thread in thread pool
static void * start_thread(void *arg) 
{
  struct worker_param *param = (struct worker_param *) arg;
  struct pool * thread_pool = param->thread_pool;
  unsigned int thread_id = param->worker_id;

  int job_idx = -1, working, num_working_threads;
  unsigned int thread_done = 0, done = 1;

  while(1)
  {
    job_idx = __sync_fetch_and_add(&(thread_pool->job_counter), 1);
    if(job_idx >= thread_pool->num_jobs) // no more jobs
    {
      num_working_threads = __sync_sub_and_fetch(&(thread_pool->num_working_threads), 1);
      break;
    }

    struct job * curr_job = &(thread_pool->job_array[job_idx]);
    curr_job->worker_id = thread_id;

    // take query vector from job and add it to thread param
    param->query_id = &(curr_job->query_id);
    param->query_order = curr_job->query_order;
    param->query_ts = curr_job->query_vector;
    param->query_ts_reordered = curr_job->query_vector_reordered;
    thread_pool->function(param); // run knn search 

    thread_pool->executed_jobs_count[thread_id]++;
  }

	return NULL;
}

void init_thread_pool(struct pool* pool, struct dstree_index * index, ts_type epsilon, unsigned int k,
              void * (*thread_func)(void *), unsigned int num_threads, unsigned int offset,
              ts_type r_delta, unsigned int * total_checked_ts, double * total_query_time, 
              float warping, struct result_vid ** all_knn_results, unsigned char store_results_in_disk,
              unsigned int *k_values, unsigned int num_k_values, struct result_vid *ground_truth_results,
              unsigned int num_gt_results, int8_t * recall_matrix, unsigned int num_query_vectors,
              struct job * job_array, unsigned int vector_length, unsigned int stop_when_nn_dist_changes)
{
  // init query stats
  dstree_init_thread_stats(index, num_threads);
  
  pthread_barrier_t *knn_update_barrier = malloc((num_threads) * sizeof(pthread_barrier_t));
  char * finished = calloc(num_threads, sizeof(char));
  pool->threads = malloc((num_threads) * sizeof(pthread_t));
  pool->cond_thread_state = malloc((num_threads) * sizeof(pthread_cond_t));
  pool->cond_mutex = malloc((num_threads) * sizeof(pthread_mutex_t));
  pool->executed_jobs_count = calloc(num_threads, sizeof(unsigned int));
  pool->working = calloc(num_threads, sizeof(unsigned int));
  pool->params = malloc(sizeof(struct worker_param) * num_threads);
  
  pool->num_threads = num_threads;
  pool->function = thread_func;
  pool->job_array = job_array;
  pool->num_jobs = num_query_vectors;
  pool->job_counter = 0; // next job to be executed
  pool->num_working_threads = num_threads;
  pool->cancelled = 0;
  
  
  if(pool->threads == NULL || pool->working == NULL || finished == NULL || pool->params == NULL)
  {
    fprintf(stderr, "Error in kashif_utils.c: Couldn't allocate memory for thread pool.\n");
    exit(1);
  }
  
  // set thread parameters
	for (int i = 0; i < num_threads; i++) 
  {
    // init barriers
    pthread_barrier_init(&knn_update_barrier[i], NULL, 2);
    pthread_cond_init(&(pool->cond_thread_state[i]), NULL);
    pthread_mutex_init(&(pool->cond_mutex[i]), NULL);

    pool->params[i].worker_id = (uint8_t)i;
    pool->params[i].thread_pool = pool;

    pool->params[i].worker_id = i;
    pool->params[i].epsilon = epsilon;
    pool->params[i].index = index;
    pool->params[i].k = k;
    pool->params[i].offset = offset;
    
    pool->params[i].r_delta = r_delta;
    pool->params[i].total_checked_ts = total_checked_ts;
    pool->params[i].total_query_set_time = total_query_time;
    pool->params[i].warping = warping;
    pool->params[i].stop_when_nn_dist_changes = stop_when_nn_dist_changes;

    pool->params[i].global_knn_results = all_knn_results;
    pool->params[i].store_results_in_disk = store_results_in_disk;
    pool->params[i].k_values = k_values; // k values for which we want to record results
    pool->params[i].num_k_values = num_k_values;
    pool->params[i].ground_truth_results = ground_truth_results;
    pool->params[i].num_gt_results = num_gt_results;
    pool->params[i].global_recall_matrix = recall_matrix;
    pool->params[i].finished = &finished[i];
    pool->params[i].knn_update_barrier = &knn_update_barrier[i];

    RESET_THREAD_PARTIAL_COUNTERS(i)
    RESET_THREAD_QUERY_COUNTERS(i)

		pthread_create(&pool->threads[i], NULL, &start_thread, &pool->params[i]);
	}

  // get result incrementally
  unsigned int pool_done = 0;
  unsigned int thread_done = 0;
  unsigned int all_finished = 0;

  while(pool_done != 1)
  {
    pool_done = 0;
    all_finished = 0;
    for (int th = 0; th < num_threads; th++) 
    {
      __atomic_load(&(pool->working[th]), &thread_done, 1);
      thread_done = 0;

      if(!(atomic_compare_exchange_strong(&(pool->working[th]), &thread_done, thread_done)))
      {
        pthread_mutex_lock( &(pool->cond_mutex[th]));
        pthread_barrier_wait(&knn_update_barrier[th]); // wait for new results from thread
        
        pthread_cond_wait(&(pool->cond_thread_state[th]), &(pool->cond_mutex[th]));
        pthread_mutex_unlock( &(pool->cond_mutex[th]));

        // compute current recall
        // float recall_from_matrix = compute_recall_from_matrix(recall_matrix, num_query_vectors, k, num_gt_results);
        // printf("\ncoordinator_thread:\t (!)\tcurrent recall =%f\n", recall_from_matrix); 
      }
    }
    if((atomic_compare_exchange_strong(&(pool->num_working_threads), &all_finished, all_finished)))
    {
      pool_done = 1; // all theread are done no more jobs to do.
      break;
    }
  }  
  
  // free memory
  printf("\ncoordinator_thread:\t\tend, all workers have finished, freeing thread pool.\n");
  pool_end(pool);
  free(knn_update_barrier);
  free(finished);
}

// end thread pool, wait for all threads to finish
void pool_end(struct pool *pool) 
{
	struct pool *thread_pool = (struct pool *) pool;
	struct pool_queue *q;
	int i;

	thread_pool->cancelled = 1;
	for (i = 0; i < thread_pool->num_threads; i++) {
		pthread_join(thread_pool->threads[i], NULL);
	}
  free(thread_pool->threads);
  free(thread_pool->cond_thread_state);
  free(thread_pool->cond_mutex);
  free(thread_pool->executed_jobs_count);
  free(thread_pool->params);
  free(thread_pool->working);
  free(thread_pool);
}


struct pool_queue * pool_init_job_queue(struct vid *query_id_arr, int ** query_order_arr,
                    ts_type ** query_vectors, ts_type ** query_vectors_reordered,
                    unsigned int num_query_vectors, unsigned int vector_length)
{
  // fill the pools job queue with query vector, each query vector is a job waiting to be executed by a thread.
	struct pool_queue *query_queue = (struct pool_queue *) malloc(sizeof(struct pool_queue) * num_query_vectors);
	

  // last element in the queue
  int j = num_query_vectors - 1;
  query_queue[j].query_id = &query_id_arr[j];
  query_queue[j].query_order = query_order_arr[j];
  query_queue[j].query_vector = query_vectors[j];
  query_queue[j].query_vector_reordered = query_vectors_reordered[j];
  query_queue[j].next = NULL;
  query_queue[j].free = 0;

  // fill the rest of the queue starting from the end of the queue
  for(int i = num_query_vectors - 2; i >= 0; i--)
  {
      query_queue[i].query_id = &query_id_arr[i];
      query_queue[i].query_order = query_order_arr[i];
      query_queue[i].query_vector = query_vectors[i];
      query_queue[i].query_vector_reordered = query_vectors_reordered[i];
      query_queue[i].next =  &query_queue[i + 1];
      query_queue[i].free = 0;
  }

  return query_queue;
}
