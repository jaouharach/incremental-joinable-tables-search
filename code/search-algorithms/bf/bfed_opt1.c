/*SCRIPT TO CONVERT ONE BINARY FILE TO TABLE, SET AND VECTOR STRUCTURES*/
#include<stdio.h>
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
#include <valgrind/callgrind.h>
#include "include.h"

#define INT_MIN -999


int main(int argc, char const *argv[])
{   
    char *experiment_dir = "/home/jaouhara/Documents/Dissertation/dssdl/experiment-results/dstree/";
    char *queries_dir = "/home/jaouhara/Documents/Dissertation/dssdl/binary-data";
    char *data_lake_directory = "/home/jaouhara/Documents/Dissertation/dssdl/binary-data";


    unsigned int l = 0; //number of datasets to be indexed
    unsigned int vector_length = 0;
    unsigned int dlsize = 0; //datalake size in GB
    unsigned int k = 1;
    unsigned int qset_num = 0 ,min_qset_size = 0, max_qset_size = 0;
    unsigned int x = 0; //top x sets to be returned 

    while (1)
    {
        /*
            ./a.out --datalake '/home/jaouhara/Documents/Dissertation/dssdl/encode/binfiles/'  --vector-length 100  --queries '/home/jaouhara/Documents/Dissertation/dssdl/encode/binfiles'  --l 5 --nq 2 --min-qset-size 1 --max-qset-size 3 --k 3 --top 3  --experiment-dir '/home/jaouhara/Documents/Dissertation/dssdl/experiments-output/bf'
        */
         struct option long_options[] = {
            {"queries", required_argument, 0, 'q'}, // queries directory
            {"datalake", required_argument, 0, 'd'}, // data lake binary files directory
            {"nq", required_argument, 0, 'u'}, //number of query sets
            {"min-qset-size", required_argument, 0, 'y'}, //minimum query set set size (number of vectors)
            {"max-qset-size", required_argument, 0, 'w'}, //maxmum query set set size (number of vectors)
            {"top", required_argument, 0, '#'}, //maxmum query set set size (number of vectors)
            {"experiment-dir", required_argument, 0, '>'},
            {"l", required_argument, 0, '|'}, // number of datasets to be indexed (tables)
            {"k", required_argument, 0, 'k'},
            {"vector-length", required_argument, 0, 't'},
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long(argc, argv, "",
                            long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
        case 'q':
            queries_dir = optarg;
            break;

        case 'd':
            data_lake_directory = optarg;
            break;

        case 'k':
            k = atoi(optarg);
            if (k < 1)
            {
                fprintf(stderr, "Please change the k to be greater than 0.\n");
                exit(-1);
            }
            break;

        case 'u':
            qset_num = atoi(optarg);
            break;
        
        case 'y':
            min_qset_size = atoi(optarg);
            break;
            
        case 'w':
            max_qset_size = atoi(optarg);
            break;

        case '#':
            x = atoi(optarg);
            break;

        case '>':
            experiment_dir = optarg;
            break;

        case '|':
            l = atoi(optarg);
            break;;

        case 't':
            vector_length = atoi(optarg);
            break;

        default:
            exit(-1);
            break;
        }
    }
    dlsize = get_dlsize(data_lake_directory, l);
    printf("Pretty Data lake size in gb = %u\n", dlsize);    

    CALLGRIND_START_INSTRUMENTATION;
    bf_sequential_search(
    data_lake_directory, 
    vector_length,
    qset_num,
    min_qset_size, 
    max_qset_size,
    x,
    experiment_dir,
    l,
    dlsize,
    k
    );
    CALLGRIND_STOP_INSTRUMENTATION;
    CALLGRIND_DUMP_STATS;
    printf("\nCongrat! End of experiment.\n"); 
    printf("Datalake size:\t\t\t%uGB\n\nTotal table:\t\t\t%u\n\nNumber of queries:\t\t%u\n\nMin query size:\t\t\t%u\n\nMax query size:\t\t\t%u\n", dlsize, l, qset_num, min_qset_size, max_qset_size);
    return 0;
}


void bf_sequential_search(
    char * bin_files_directory, 
    unsigned int vector_length,
    unsigned int qset_num,
    unsigned int min_qset_size, 
    unsigned int max_qset_size,
    unsigned int x,
    char * experiment_dir,
    unsigned int l, //total num tables indexed in dstree
    unsigned int dlsize, // //total disk size of tables indexed in dstree
    unsigned int k
    )
{
    int opened_files = 0; 
    char * algorithm = "bfed_opt1";
    // open source dir
    struct dirent *dfile;
    DIR *dir = opendir(bin_files_directory);

    if (!dir)
    {
        printf("Unable to open directory stream %s!", bin_files_directory);
        exit(1);
    }

    // allocate memory for vector
    struct Vector * query_set = malloc(0);
    struct Vector query_vector; 
    query_vector.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
 
	double query_time = 0.0; 
    unsigned int total_checked_vec = 0;
    unsigned int total_queries = qset_num;
    bool found_query = false;

    char * results_dir =  make_result_directory(algorithm, experiment_dir, l, qset_num, min_qset_size, max_qset_size);
    //initialize list of all knn results (from all query vectors in query set)
    struct query_result * all_knn_results = NULL;
    struct vid * top_matches;
    int next_vec;

    // for every file (table)
    while ((dfile = readdir(dir)) != NULL && qset_num != 0)
    {
        //if file is binary file
        if(is_binaryfile(dfile->d_name))
        {
            opened_files += 1;

            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = ""; 
            strcat(bin_file_path, bin_files_directory);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bfed_opt1.c:  Vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            FILE * bin_file = fopen (bin_file_path, "rb");
            if (bin_file == NULL)
            {
                printf("Error in bfed_opt1.c: File %s not found!\n", bin_file_path);
                exit(1);
            }


            /* Start processing file: read every vector in binary file*/
            int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            uint32_t nvec = 0u;
            ts_type val;

            while(total_bytes)
            {
                // beginning of a set of vectors
                if(i == 0)
                { 
                    if(qset_num == 0)
                        break;

                    //read first integer to check how many vactors in current set
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    total_bytes--;

                    // query set does not fit requirments move to next set
                    if((unsigned int)nvec < min_qset_size || (unsigned int)nvec > max_qset_size)
                    {
                        fseek(bin_file, nvec * 4 *vector_length, SEEK_CUR);
                        i = 0;
                        j = 0;
                        total_bytes -= nvec * vector_length;
                        continue;
                    }
                    found_query = true;
                    query_set = (struct Vector *) realloc(query_set, sizeof(struct Vector) * nvec);
                    printf("\n----------------------------------------------------------------\n");
                    printf("Query %u/%u.\nRun Query: Q = (%d, %d) |Q| = %d", (total_queries-qset_num)+1, total_queries , table_id, set_id, nvec);
                    query_time = 0.0;
                    total_checked_vec = 0;
                    next_vec = 0;
                    qset_num--;

                    query_vector.table_id = table_id;
                    query_vector.set_id = set_id;

                    set_id = set_id + 1;
                    i++;
                    j = 0;
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    // end of vector but still in current set
                    if(j > 99){
                        j = 0; 
                        
                        // Append vector to query set
                        query_set[next_vec] = query_vector;
                        next_vec += 1;
                    }

                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    total_bytes--;
                    query_vector.values[j] = val;

                    // end of last vector in current  set
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        query_vector.values[j] = val;
                        query_set[next_vec] = query_vector;
                        next_vec = 0;
                        
                        /*run query set */
                        // Perform brute force knn search for all query vectors
                        all_knn_results =  brute_force_knn_search_optimized(bin_files_directory, l, vector_length, query_set, nvec, k, &query_time, &total_checked_vec);
                        

                        /* End of Query set */
                        printf("\nTotal Query time = %.20f\nTotal checked vectors = %u\n", query_time, total_checked_vec);
                        /* Save query results to csv file */
                        char * query_result_file = make_file_path(results_dir, query_vector.table_id, query_vector.set_id, nvec, l, dlsize, vector_length, query_time, total_checked_vec);

                        printf("\nQuery result file path: %s\n", query_result_file);
                        save_to_query_result_file(query_result_file, table_id, query_vector.set_id, k * nvec, all_knn_results);

                        struct vid * top_matches =  get_top_x(k * nvec,  all_knn_results, x);
                        for(int k =0; k < x; k++)
                        printf("TOP %d:\t (%u, %u)\n", k+1,top_matches[k].table_id, top_matches[k].set_id);
                        printf("----------------------------------------------------------------");

                        i = 0; j = 0;
                        nvec = 0u;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
        }
    }
    free(all_knn_results);
}



struct query_result * brute_force_knn_search_optimized(char * bin_files_dir, unsigned int l, unsigned int vector_length, struct Vector * qset, unsigned int qnvec, unsigned int k, double * query_time, unsigned int *total_checked_vec)
{
    ts_type * bsf = (ts_type *) malloc(sizeof(ts_type) * qnvec);
    unsigned int * next_knn = (unsigned int *) malloc(sizeof(unsigned int) * qnvec);
    //queue containing kNN results for every vector
    struct query_result * curr_knn = (struct query_result *) malloc(sizeof(struct query_result) * k*qnvec);
    

    struct Vector v; 
    v.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
    clock_t query_vec_time = 0;


    for (int i = 0; i < qnvec; ++i)
    {
        bsf[i] = FLT_MAX;
        next_knn[i] = 0;
    }

    for (int j = 0; j < qnvec*k; ++j)
    {
        curr_knn[j].vector_id = malloc(sizeof(struct vid));
        if(curr_knn[j].vector_id == NULL)
            exit(1);
        curr_knn[j].vector_id->set_id = 0;
        curr_knn[j].vector_id->table_id = 0;
        curr_knn[j].distance = FLT_MAX;
    }

    // Open binary files  dir
    struct dirent *dfile;
    DIR *dir = opendir(bin_files_dir);

    if (!dir)
    {
      printf("Unable to open directory stream %s!", bin_files_dir);
      exit(1);
    }
    while ((dfile = readdir(dir)) != NULL && l > 0)
    {

        if(is_binaryfile(dfile->d_name))
        {
            --l;
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = ""; 
            strcat(bin_file_path, bin_files_dir);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bfed_opt1.c:  Vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            FILE * bin_file = fopen (bin_file_path, "rb");
            if (bin_file == NULL)
            {
                printf("Error in bfed_opt1.c: File %s not found!\n", bin_file_path);
                exit(1);
            }

            ts_type val; 
            uint32_t nvec = 0u;
            ts_type d = 0.0;
            unsigned int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            
            //read every candidate set and candidate vector in binary file
            clock_t begin_query_in_curr_file = clock();
            while(total_bytes)
            {
                if(i == 0)
                { 
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    total_bytes--;

                    v.table_id = table_id;
                    v.set_id = set_id;
                    set_id = set_id +  1;

                    if(v.table_id == qset[0].table_id && v.set_id == qset[0].set_id)
                    {
                        // do not match query set to itself
                        fseek(bin_file, nvec * 4 *vector_length, SEEK_CUR);
                        i = 0;
                        j = 0;
                        total_bytes -= nvec * vector_length;
                        
                        continue;
                    }
                    
                    
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    if(j > 99)
                    {    
                        //printf("CAN VECTOR (%d, %d)-------------------------------\n", v.table_id, v.set_id);                    
                        j = 0;
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidian_distance_with_early_abandoning(qset[h].values, v.values, vector_length, bsf[h]);
                            // if current vector is closer to query vector
                            if (d < bsf[h])
                            {
                                bsf[h] = d;
                                // move previous bsfs to put curr bsf as first
                                for(int m = 1; m < k; ++m)
                                {
                                    curr_knn[h*k+m].distance = curr_knn[(h*k+m)-1].distance;
                                    curr_knn[h*k+m].vector_id->set_id = curr_knn[(h*k+m)-1].vector_id->set_id;
                                    curr_knn[h*k+m].vector_id->table_id = curr_knn[(h*k+m)-1].vector_id->table_id;
                                }

                                curr_knn[h*k].distance = bsf[h];
                                curr_knn[h*k].vector_id->set_id = v.set_id;
                                curr_knn[h*k].vector_id->table_id = v.table_id;
                                next_knn[h] = 1;

                            }
                            // if its another vector with same bsf
                            else if(d == bsf[h] && next_knn[h] < k)
                            {
                                curr_knn[(h*k)+next_knn[h]].distance = bsf[h];
                                curr_knn[(h*k)+next_knn[h]].vector_id->set_id = v.set_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->table_id = v.table_id;
                                next_knn[h] = next_knn[h] + 1;
                            }
                        }    
                    }
                    // read new float value
                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    total_bytes--;
                    v.values[j] = val;
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        if(j != 99){
                            printf("ERROR !j= %d", j);
                            exit(1);
                        }
                            
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidian_distance_with_early_abandoning(qset[h].values, v.values, vector_length, bsf[h]);
                            
                            if (d < bsf[h])
                            {
                                bsf[h] = d;
                                // move previous bsfs to put curr bsf as first
                                for(int m = 1; m < k; ++m)
                                {
                                    curr_knn[h*k+m].distance = curr_knn[(h*k+m)-1].distance;
                                    curr_knn[h*k+m].vector_id->set_id = curr_knn[(h*k+m)-1].vector_id->set_id;
                                    curr_knn[h*k+m].vector_id->table_id = curr_knn[(h*k+m)-1].vector_id->table_id;
                                }

                                curr_knn[h*k].distance = bsf[h];
                                curr_knn[h*k].vector_id->set_id = v.set_id;
                                curr_knn[h*k].vector_id->table_id = v.table_id;
                                next_knn[h] = 1;

                                
                            }
                            // if its another vector with same bsf
                            else if(d == bsf[h] && next_knn[h] < k)
                            {
                                curr_knn[(h*k)+next_knn[h]].distance = bsf[h];
                                curr_knn[(h*k)+next_knn[h]].vector_id->set_id = v.set_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->table_id = v.table_id;
                                next_knn[h] = next_knn[h] + 1;
                            }
                        }

                        i = 0; j = 0;
                        nvec = 0u;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            query_vec_time += clock() - begin_query_in_curr_file;
        fclose(bin_file);
        }
    }
    *query_time += (double)query_vec_time/CLOCKS_PER_SEC;
    free(bsf);
    free(next_knn);
    return curr_knn;
}


ts_type euclidian_distance_with_early_abandoning(ts_type *q, ts_type *v, unsigned int len, ts_type bsf)
{
    int s = 0; 
    ts_type d = 0.0;
    while(s < len)
    {
        d += (q[s] - v[s]) * (q[s] - v[s]);
        
        if(d > bsf)
        {
            break;
        }
        s++;  
    }
    return d;
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
unsigned int get_dlsize(char* dl_dir, unsigned int l){
	struct dirent *dfile;
    DIR *dir = opendir(dl_dir);
    float total_dlsize = 0.0;
	FILE* fp;

    if (!dir)
    {
        printf("Unable to open directory stream!");
        exit(1);
    }

    while ((dfile = readdir(dir)) != NULL && l > 0)
    {
		if(is_binaryfile(dfile->d_name)){
			char bin_file_path[PATH_MAX + 1] = ""; 
        	strcat(bin_file_path, dl_dir);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);
			l--;
			// ** get binary file size 
	    	fp = fopen(bin_file_path, "r");
			if (fp == NULL) {
				printf("Could not get size of binary file %s. File not found!\n", dfile->d_name);
				exit(1);
			}
			fseek(fp, 0L, SEEK_END);
			total_dlsize += ftell(fp);
			fclose(fp);
		}        
    }
    printf("Data lake size in gb = %.6f\n", total_dlsize / 1073741824);
   
	return (unsigned int) round(total_dlsize/1073741824);
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
            printf("Error in bfed_opt1.c: Could not open file %s!\n", csv_file);
            exit(1);
    }

        // write header
        fprintf(fp, "TQ:Q, TS:S, qindex, sindex, q, s, d");
        

    // write results
    for(int i = 0; i < num_knns; i++){
        fprintf(fp, "\n");
        fprintf(fp,"%u:%u, %u:%u, 0, 0, [], [], %.3f", qtable_id, qset_id, knn_results[i].vector_id->table_id, knn_results[i].vector_id->set_id, knn_results[i].distance);
    }
	fclose(fp);
}


// new function make result file name and path.
char * make_file_path(char * experiment_dir, unsigned int qtable_id, unsigned int qset_id, unsigned int qsize, unsigned int l, unsigned int dlsize, unsigned int vector_length, float runtime, unsigned int total_checked_vec)
{
	DIR* dir = opendir(experiment_dir);
	if (!dir)
    {
		printf("WARNING! Experiment direstory '%s' does not exist!", experiment_dir);
		exit(1);
	}
    char * filepath = malloc(get_ndigits(qtable_id) + get_ndigits(qset_id) + get_ndigits(l)
							 + get_ndigits(dlsize) + get_ndigits(vector_length) + get_ndigits((unsigned int) runtime) + get_ndigits(total_checked_vec)
							 + get_ndigits(qsize) + strlen("TQ_Q_qsize_l_dlsize_len_runtime_ndistcalc_dataaccess.csv")
							 + strlen(experiment_dir)
							 + 6 // float decimal precision for dlsize and runtime (.00)
							 + 1);

	sprintf(filepath, "%s/TQ%u_Q%u_qsize%u_l%u_dlsize%u_len%u_runtime%.4f_ndistcalc_dataaccess%u.csv"
			, experiment_dir, qtable_id, qset_id, qsize, l, dlsize, vector_length, runtime, total_checked_vec);

	return filepath;
}


// new function make experiment results dir
char * make_result_directory(char* algorithm, char * experiment_dir, unsigned int l, unsigned int nq, unsigned int min_qset_size, unsigned int max_qset_size)
{
	char * result_dir_name = malloc(get_ndigits(l) + get_ndigits(nq)
									+ get_ndigits(min_qset_size) + get_ndigits(max_qset_size)
									+ strlen("/_l_q_min_max") + strlen(experiment_dir)+ strlen(algorithm) + 1);

	sprintf(result_dir_name, "%s/%s_l%u_%uq_min%u_max%u", experiment_dir, algorithm, l, nq, min_qset_size, max_qset_size);

	printf("result directory name: %s\n", result_dir_name);
	DIR* dir = opendir(result_dir_name);
	if (dir)
  {
      printf("WARNING! Results directory already exists. Please delete directory : %s.\n", result_dir_name);
      exit(-1);
  }
  mkdir(result_dir_name, 0777);
  
  return result_dir_name;
}

// new function get top x matching sets from array of knn results
struct vid * get_top_x(int num_knn_results, struct query_result * knn_results,unsigned int x)
{
	// if number of results is already equal to x  
	//if(x == num_knn_results)
  //	exit(1);

	// frequency = number of matching vectors with the query set vectors
	int i, j, k;
	int * frequency_array = (int *) malloc(num_knn_results * sizeof(int));


	// array of top x sets and 
	int * max_freqs = (int *) malloc(x * sizeof(int));
	struct vid * top = (struct vid *) calloc(x, sizeof(struct vid));


	for(i = 0; i<num_knn_results; i++)
		frequency_array[i] = -1;

	for(i = 0; i<x; i++)
		max_freqs[i] = INT_MIN;

	//count occurence of each (table_id, set_id)
	for (i = 0; i < num_knn_results; i++)
	{
		int Count = 1;
		for(j = i + 1; j < num_knn_results; j++)
		{
    		if(knn_results[i].vector_id->table_id == knn_results[j].vector_id->table_id)
    		{
    			if((knn_results[i].vector_id->set_id == knn_results[j].vector_id->set_id))
    			{
    				Count++;
    				frequency_array[j] = 0;
    			}
    		}
    	}
    	if(frequency_array[i] != 0)
    	{
    		frequency_array[i] = Count;
		}
	}


	//use frequency array to find  top x most frequent ids 
	for (i = 0; i < num_knn_results; i++)
  	{
  		for(j = 0; j < x; j++)
  		{
  			if(frequency_array[i] > max_freqs[j])
  			{
  				//move curr top k to be top k+1
	  			for(k = x-1; k > j; k--)
	  			{
	  				//printf("replacing %d\n with ", max_freqs[k]);
					max_freqs[k] = max_freqs[k-1];
					top[k].table_id = top[k-1].table_id;
	  				top[k].set_id = top[k-1].set_id;

				}
				//replace curr top k
				max_freqs[j] = frequency_array[i];
	  			top[j].table_id = knn_results[i].vector_id->table_id;
  				top[j].set_id = knn_results[i].vector_id->set_id;

  				break;
  			}
  		}
  	}
	return top;
}