#include<stdio.h>

typedef float ts_type;

typedef struct vector{
    unsigned int table_id;
    unsigned int set_id;
    unsigned int pos;
    ts_type * values;
} vector;

typedef struct query_result {
  ts_type distance;
  struct vid * vector_id; 
  unsigned int query_vector_pos; // position of the query vector in the query set 
};
// vector id
struct vid {
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos;
  char raw_data_file[300]; // binary file
};

// result set id
struct result_sid { 
  unsigned int table_id;
  unsigned int set_id;
  char raw_data_file[300]; // name of the raw (json) file where set is store
  unsigned int overlap_size;
};

int get_ndigits(unsigned int n);
unsigned int get_data_gb_size(char* dl_dir, unsigned int l);
bool is_binaryfile(const char *filename);
void bf_sequential_search( char * bin_files_directory, unsigned int vector_length, unsigned int qset_num, 
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top, char * result_dir,
    unsigned int total_data_files, unsigned int dataset_gb_size, unsigned int k);
void save_to_query_result_file(char * csv_file, unsigned int qtable_id, unsigned int qset_id, int num_knns, struct query_result * knn_results);
char * make_file_path(char * result_dir, unsigned int qtable_id, unsigned int qset_id, unsigned int qsize, unsigned int total_data_files, unsigned int dataset_gb_size, unsigned int vector_length, float runtime, unsigned int total_checked_ts);
char * make_result_directory(char* algorithm, char * result_dir, unsigned int total_data_files, unsigned int nq, unsigned int min_qset_size, unsigned int max_qset_size);
struct query_result * brute_force_knn_search(char * bin_files_dir, unsigned int total_data_files, unsigned int vector_length, struct vector q, unsigned int k, double * query_time, unsigned int *total_checked_ts, struct query_result *knn_results);
ts_type euclidian_distance(ts_type *q, ts_type *v, unsigned int len);
struct result_sid * get_top_sets(struct query_result * knn_results, int num_knn_results, unsigned int num_top);
struct query_result * brute_force_knn_search_optimized(char * bin_files_dir, unsigned int total_data_files, unsigned int vector_length, struct vector * qset, unsigned int qnvec, unsigned int k, double * query_time, unsigned int *total_checked_vec);
ts_type cosine_distance(ts_type *q, ts_type *v, unsigned int len);
ts_type euclidian_distance_with_early_abandoning(ts_type *q, ts_type *v, unsigned int len, ts_type bsf);
void bf_sequential_search_with_threshold(char * bin_files_directory, unsigned int vector_length, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top,
    char * result_dir, unsigned int total_data_files, //total num tables indexed in dstree
    unsigned int dataset_gb_size, unsigned int k, float theta);
struct query_result * brute_force_knn_search_optimized_with_sim_threshold(char * bin_files_dir, unsigned int total_data_files, unsigned int vector_length, struct vector * qset,unsigned int qnvec, unsigned int k, double * query_time, unsigned int *total_checked_vec, float theta);