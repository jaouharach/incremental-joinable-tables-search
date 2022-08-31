#include <stdio.h>
#include <string.h>
#include<stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <dirent.h>
#include <getopt.h>
#include <math.h>
#include <sys/stat.h>
#include <float.h>
#include <time.h>
#include <linux/limits.h>
#include "../include/utils.h"

ts_type euclidian_distance(ts_type *q, ts_type *v, unsigned int len)
{
    ts_type sum = 0.0;
    for(int i = 0; i < len; ++i)
        sum = sum + ((q[i] - v[i]) * (q[i] - v[i]));
    return sum;
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


// new function get number of digits in integer
int get_ndigits(unsigned int n){
	int total_digits = 0;
	while(n!=0){
		//4
		n = n/10;
		++total_digits;
	}
	return total_digits;
}

// get data lake size in GB
unsigned int get_data_gb_size(char *dataset, unsigned int total_data_files) {
  struct dirent *dfile;
  DIR *dir = opendir(dataset);
  float total_dlsize = 0.0;
  FILE *fp;

  if (!dir) {
    printf("Here 21!\n");
    fprintf("Error in kashif_utils.c: Unable to open directory stream! %s", dataset);
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0) {
    if (is_binaryfile(dfile->d_name)) {
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, dataset);
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
  return (unsigned int)round(total_dlsize / 1073741824); // return data size in GB
}

bool is_binaryfile(const char *filename)
{
    // check if filename has bin extesion.
    char *ext = ".bin";
    size_t nl = strlen(filename), el = strlen(ext);
    return nl >= el && !strcmp(filename + nl - el, ext);
}

// new function save query results to csv file
void save_to_query_result_file(char * csv_file, unsigned int qtable_id, unsigned int qset_id, int num_knns, struct query_result * knn_results){
	FILE *fp;
	int i,j;
	fp = fopen(csv_file,"w+");

    if (fp == NULL) {
            printf("Error in bfed.c: Could not open file %s!\n", csv_file);
            exit(1);
    }

        // write header
        fprintf(fp, "TQ:Q, TS:S, qindex, sindex, q, s, d");
        

    // write results
    for(int i = 0; i < num_knns; i++){
        fprintf(fp, "\n");
        fprintf(fp, "%u:%u, %u:%u, %u, %u, [], [], %.3f", qtable_id, qset_id,
            knn_results[i].vector_id->table_id, knn_results[i].vector_id->set_id,
            knn_results[i].query_vector_pos, knn_results[i].vector_id->pos,
            knn_results[i].distance);
    }
	fclose(fp);
}


// new function make result file name and path.
char * make_file_path(char * result_dir, unsigned int qtable_id, unsigned int qset_id, unsigned int qsize, unsigned int l, unsigned int dlsize, unsigned int vector_length, float runtime, unsigned int total_checked_vec)
{
  runtime /= 1000000;
	DIR* dir = opendir(result_dir);
	if (!dir)
    {
		printf("WARNING! Experiment direstory '%s' does not exist!", result_dir);
		exit(1);
	}
    char * filepath = malloc(get_ndigits(qtable_id) + get_ndigits(qset_id) + get_ndigits(l)
							 + get_ndigits(dlsize) + get_ndigits(vector_length) + get_ndigits((unsigned int) runtime) + get_ndigits(total_checked_vec)
							 + get_ndigits(qsize) + strlen("TQ_Q_qsize_l_dlsize_len_runtime_ndistcalc_dataaccess.csv")
							 + strlen(result_dir)
							 + 6 // float decimal precision for dlsize and runtime (.00)
							 + 1);

	sprintf(filepath, "%s/TQ%u_Q%u_qsize%u_l%u_dlsize%u_len%u_runtime%.4f_ndistcalc_dataaccess%u.csv"
			, result_dir, qtable_id, qset_id, qsize, l, dlsize, vector_length, runtime, total_checked_vec);

	return filepath;
}


// new function make experiment results dir
char * make_result_directory(char* algorithm, char * result_dir, unsigned int l, unsigned int nq, unsigned int min_qset_size, unsigned int max_qset_size)
{
	char * result_dir_name = malloc(get_ndigits(l) + get_ndigits(nq)
									+ get_ndigits(min_qset_size) + get_ndigits(max_qset_size)
									+ strlen("/_l_q_min_max") + strlen(result_dir)+ strlen(algorithm) + 1);

	sprintf(result_dir_name, "%s/%s_l%u_%uq_min%u_max%u", result_dir, algorithm, l, nq, min_qset_size, max_qset_size);

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
struct result_sid *get_top_sets(struct query_result *knn_results, int num_knn_results, unsigned int num_top)
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
