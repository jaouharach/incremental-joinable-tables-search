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

#define NBITS 64
typedef float ts_type;

typedef struct vector {
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos;
  ts_type *values;
} vector;

// create result file name and path.
char *make_file_path(char *result_dir, unsigned int qtable_id,
                     unsigned int qset_id, unsigned int qsize, unsigned int l,
                     unsigned int dlsize, unsigned int vector_length,
                     float runtime, unsigned int total_checked_vec);

// save query results to csv file
enum response save_to_query_result_file(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int num_knns,
                                        struct query_result *knn_results);

// create experiment results dir
char *make_result_directory(char *result_dir, unsigned int total_data_files,
                            unsigned int nq, unsigned int min_qset_size,
                            unsigned int max_qset_size);

// count digits in integer
int get_ndigits(unsigned int n);

/* convert IEEE double precision (64 bits) to float (/double) */
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
                                    unsigned int total_data_files) {
  struct dirent *dfile;
  DIR *dir = opendir(bindir);
  unsigned long total_vectors = 0;
  unsigned int datasize;
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
      sscanf(dfile->d_name, "data_size%d%*[^0123456789]", &datasize);
      total_vectors += datasize;
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
                     unsigned int vector_length, float runtime,
                     unsigned int total_checked_vec) {
  DIR *dir = opendir(result_dir);
  if (!dir) {
    printf("WARNING! Experiment direstory '%s' does not exist!", result_dir);
    exit(1);
  }
  char *filepath = malloc(
      get_ndigits(qtable_id) + get_ndigits(qset_id) +
      get_ndigits(total_data_files) + get_ndigits(dlsize) +
      get_ndigits(vector_length) + get_ndigits((unsigned int)runtime) +
      get_ndigits(total_checked_vec) + get_ndigits(qsize) +
      strlen("TQ_Q_qsize_l_dlsize_len_runtime_ndistcalc_dataaccess.csv") +
      strlen(result_dir) +
      10 // float decimal precision for dlsize and runtime (.00)
      + 1);

  sprintf(filepath,
          "%s/"
          "TQ%u_Q%u_qsize%u_l%u_dlsize%u_len%u_runtime%.3f_ndistcalc_"
          "dataaccess%u.csv\0",
          result_dir, qtable_id, qset_id, qsize, total_data_files, dlsize,
          vector_length, runtime, total_checked_vec);

  closedir(dir);
  return filepath;
}

// save query results to csv file
enum response save_to_query_result_file(char *csv_file, unsigned int qtable_id,
                                        unsigned int qset_id, int num_knns,
                                        struct query_result *knn_results) {
  FILE *fp;
  int i, j;
  fp = fopen(csv_file, "w+");

  if (fp == NULL) {
    fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n",
            csv_file);
    return FAILURE;
  }

  // write header
  fprintf(fp, "TQ:Q, TS:S, qindex, sindex, q, s, d");

  // write results
  for (int i = 0; i < num_knns; i++) {
    fprintf(fp, "\n");
    fprintf(fp, "%u:%u, %u:%u, %u, %u, [], [], %.3f", qtable_id, qset_id,
            knn_results[i].vector_id->table_id, knn_results[i].vector_id->set_id,
            knn_results[i].query_vector_pos, knn_results[i].vector_id->pos,
            knn_results[i].distance);
  }
  fclose(fp);

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
  DIR *dir = opendir(result_dir_name);
  if (dir) {
    delete_directory(result_dir_name);
  }
  mkdir(result_dir_name, 0777);
  closedir(dir);
  
  return result_dir_name;
}

// get sets with the the largest number of matching vectors with the query.
struct result_sid *get_top_sets(struct query_result *knn_results, int num_knn_results, 
                         unsigned int num_top)
{
  // frequency = number of matching vectors with the query set vectors
  int i, j, k;

  struct result_sid * distinct_sets = NULL;
  struct vid * distinct_match_vectors = NULL;


  unsigned int num_distinct_sets = 0;
  int found = 0;

  num_distinct_sets += 1;
  distinct_sets = realloc(distinct_sets, sizeof(struct result_sid));
  distinct_sets[0].table_id = knn_results[0].vector_id->table_id;
  distinct_sets[0].set_id = knn_results[0].vector_id->set_id;
  distinct_sets[0].overlap_size = 0;
  strcpy(distinct_sets[0].raw_data_file, knn_results[0].vector_id->raw_data_file);

  // get distinct sets
  for ( i = 1; i < num_knn_results; i++)
  {
    found = 0;
    for(j = num_distinct_sets - 1; j >= 0; j--)
    {
      if(distinct_sets[j].table_id == knn_results[i].vector_id->table_id)
        if(distinct_sets[j].set_id == knn_results[i].vector_id->set_id)
        {
          found = 1;
          break;
        }
    }
    if(found == 0)
    {
      num_distinct_sets += 1;
      distinct_sets = realloc(distinct_sets, (sizeof(struct result_sid) * num_distinct_sets));
      distinct_sets[num_distinct_sets - 1].table_id = knn_results[i].vector_id->table_id;
      distinct_sets[num_distinct_sets - 1].set_id = knn_results[i].vector_id->set_id;
      distinct_sets[num_distinct_sets - 1].overlap_size = 0;
      strcpy(distinct_sets[num_distinct_sets - 1].raw_data_file, knn_results[i].vector_id->raw_data_file);
    }
  }

  unsigned int num_distinct_vectors = 0;
  
  num_distinct_vectors += 1;
  distinct_match_vectors = realloc(distinct_match_vectors, sizeof(struct vid));
  distinct_match_vectors[0].table_id = knn_results[0].vector_id->table_id;
  distinct_match_vectors[0].set_id = knn_results[0].vector_id->set_id;
  distinct_match_vectors[0].pos = knn_results[0].vector_id->pos;

  // get distinct vectors
  for ( i = 1; i < num_knn_results; i++)
  {
    found = 0;
    for(j = num_distinct_vectors - 1; j >= 0; j--)
    {
      if(distinct_match_vectors[j].table_id == knn_results[i].vector_id->table_id)
        if(distinct_match_vectors[j].set_id == knn_results[i].vector_id->set_id)
        if(distinct_match_vectors[j].pos == knn_results[i].vector_id->pos)
        {
          found = 1;
          break;
        }
    }
    if(found == 0)
    {
      num_distinct_vectors += 1;
      distinct_match_vectors = realloc(distinct_match_vectors, (sizeof(struct vid) * num_distinct_vectors));
      distinct_match_vectors[num_distinct_vectors - 1].table_id = knn_results[i].vector_id->table_id;
      distinct_match_vectors[num_distinct_vectors - 1].set_id = knn_results[i].vector_id->set_id;
      distinct_match_vectors[num_distinct_vectors - 1].pos = knn_results[i].vector_id->pos;
    }
  }

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

  // fill the rest withe the last elemet
  if(num_distinct_sets < num_top)
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

  return distinct_sets;
}