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
#include "../include/utils.h"
#include "../include/stats.h"
#include <errno.h>

int main(int argc, char const *argv[])
{   
    INIT_STATS()
    COUNT_TOTAL_TIME_START

    char *result_dir = "";
    char *queries_dir = "";
    char *dataset = "";


    unsigned int total_data_files = 0; //number of data files be indexed
    unsigned int vector_length = 0;
    unsigned int data_gb_size = 0; //datalake size in GB
    unsigned int k = 1;
    unsigned int qset_num = 0 ,min_qset_size = 0, max_qset_size = 0;
    unsigned int top = 0; //top x sets to be returned 
    unsigned char identical_knn_search = 0;
    while (1)
    {
        /*
            ./a.out --datalake '/home/jaouhara/Documents/Dissertation/dssdl/encode/binfiles/'  --vector-length 100  --queries '/home/jaouhara/Documents/Dissertation/dssdl/encode/binfiles'  --l 5 --nq 2 --min-qset-size 1 --max-qset-size 3 --k 3 --top 3  --experiment-dir '/home/jaouhara/Documents/Dissertation/dssdl/experiments-output/bf'
        */
         struct option long_options[] = {
            {"queries", required_argument, 0, 'q'}, // queries directory
            {"dataset", required_argument, 0, 'd'}, // data lake binary files directory
            {"nq", required_argument, 0, 'u'}, //number of query sets
            {"min-qset-size", required_argument, 0, 'y'}, //minimum query set set size (number of vectors)
            {"max-qset-size", required_argument, 0, 'w'}, //maxmum query set set size (number of vectors)
            {"top", required_argument, 0, '#'}, //maxmum query set set size (number of vectors)
            {"result-dir", required_argument, 0, '>'}, // where query results will be stored
            {"total-data-files", required_argument, 0, '|'}, // number of datasets to be indexed (tables)
            {"k", required_argument, 0, 'k'},
            {"vector-length", required_argument, 0, 't'},
            {"identical-knn-search", no_argument, 0, 'i'},
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
            dataset = optarg;
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
            top = atoi(optarg);
            break;

        case '>':
            result_dir = optarg;
            break;

        case '|':
            total_data_files = atoi(optarg);
            break;

        case 't':
            vector_length = atoi(optarg);
            break;

        case 'i':
            identical_knn_search = 1;
            break;

        default:
            exit(-1);
            break;
        }
    }
    data_gb_size = get_data_gb_size(dataset, total_data_files);
    printf("Total data size in gb = %u\n", data_gb_size);    

    bf_sequential_search(queries_dir, dataset, vector_length, qset_num, min_qset_size, max_qset_size,
                            top, result_dir, total_data_files, data_gb_size, k, identical_knn_search);

    printf("\n>>>> Congrat! End of experiment."); 
    printf("\n>>>> Results are stored in %s", result_dir); 
    COUNT_TOTAL_TIME_END
    printf("\nSanity check: combined indexing and querying times should be less "
            "than: %f secs \n",
            total_time / 1000000);
    // printf("\n>>>> Summary:\nDatalake size:\t\t\t%uGB\n\nTotal table:\t\t\t%u\n\nNumber of queries:\t\t%u\n\nMin query size:\t\t\t%u\n\nMax query size:\t\t\t%u\n", data_gb_size, total_data_files, qset_num, min_qset_size, max_qset_size);
    return 0;
}


void bf_sequential_search(char * queries, char * dataset, unsigned int vector_length, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top, char * result_dir,
    unsigned int total_data_file, unsigned int data_gb_size, unsigned int k, unsigned char identical_knn_search)
{
    printf("bf_sequential_search();\n");
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START
    int opened_files = 0; 
    char * algorithm = "bf";
    // open source dir
    struct dirent *dfile;
    COUNT_PARTIAL_INPUT_TIME_START
    DIR *dir = opendir(queries);
    COUNT_PARTIAL_INPUT_TIME_END

    unsigned int total_knns = 0;

    if (!dir)
    {
        printf("Unable to open directory stream %s!", queries);
        exit(1);
    }

    // allocate memory for vector
    struct vector * query_set = malloc(0);
    struct vector query_vector; 
    query_vector.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
 
    unsigned int total_checked_vec = 0;
    unsigned int total_queries = qset_num;
    bool found_query = false;

    char * results_dir =  make_result_directory(algorithm, result_dir, total_data_file, qset_num, min_qset_size, max_qset_size);
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
            strcat(bin_file_path, queries);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bf.c:  vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }
            
            /* read binary file */
            COUNT_PARTIAL_INPUT_TIME_START
            FILE * bin_file = fopen (bin_file_path, "rb");
            COUNT_PARTIAL_INPUT_TIME_END
            if (bin_file == NULL)
            {
                printf("Error in bf.c: File %s not found! because of error: %s\n", bin_file_path, strerror(errno));
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
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;

                    // query set does not fit requirments move to next set
                    if((unsigned int)nvec < min_qset_size || (unsigned int)nvec > max_qset_size)
                    {
                        COUNT_PARTIAL_INPUT_TIME_START
                        fseek(bin_file, nvec * 4 *vector_length, SEEK_CUR);
                        COUNT_PARTIAL_INPUT_TIME_END
                        i = 0;
                        j = 0;
                        total_bytes -= nvec * vector_length;
                        //printf("size = %u, skip();\n", nvec);
                        continue;
                    }
                    found_query = true;
                    query_set = (struct vector *) realloc(query_set, sizeof(struct vector) * nvec);
                    // (todo: free values memory)
                    for(int t = 0; t < nvec; t++)
                    {
                        query_set[t].values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
                        if(query_set[t].values == NULL)
                        {
                            printf("Error in bf.c: couldn't allocate memory for query values.");
                            exit(1);
                        }
                    }
                    printf("\nQuery %u/%u. id:(%u, %u), size = %u\n\n", (total_queries-qset_num)+1, total_queries, table_id, set_id, nvec);
                    
                    total_checked_vec = 0;
                    next_vec = 0;
                    qset_num--;

                    query_vector.table_id = table_id;
                    query_vector.set_id = set_id;
                    query_vector.pos = 0;
                    set_id = set_id + 1;
                    total_knns = 0; 

                    i++;
                    j = 0;
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    // end of vector but still in current set
                    if(j > (vector_length - 1)){
                        j = 0; 
                        
                        // Append vector to query set
                        query_set[next_vec].table_id = query_vector.table_id;
                        query_set[next_vec].set_id = query_vector.set_id;
                        query_set[next_vec].pos = query_vector.pos;
                        for(int x = 0; x < vector_length; x++)
                        {
                            query_set[next_vec].values[x] = query_vector.values[x];
                        }
                        
                        query_vector.pos += 1;
                        next_vec += 1;
                    }
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;
                    query_vector.values[j] = val;

                    // end of last vector in current  set
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        query_vector.values[j] = val;
                        query_set[next_vec].table_id = query_vector.table_id;
                        query_set[next_vec].set_id = query_vector.set_id;
                        query_set[next_vec].pos = query_vector.pos;
                        for(int x = 0; x < vector_length; x++)
                        {
                            query_set[next_vec].values[x] = query_vector.values[x];
                        }
                        next_vec = 0;
                        
                        // only return vectors with distance 0
                        if(identical_knn_search == 1)
                        {
                            printf("\n>>> performing identical knn search ...\n");
                           /*run query set */
                            // Perform brute force knn search for all query vectors
                            all_knn_results =  brute_force_identical_knn_search_optimized(dataset, total_data_file, vector_length, query_set, nvec, &total_checked_vec, &total_knns);

                            // printf("\nDone. total knns = %d\n.", total_knns);
                            COUNT_PARTIAL_TIME_END
                            query_time = partial_time - (partial_input_time + partial_output_time);
                            
                            /* End of Query set */
                            /* Save query results to csv file */
                            char * query_result_file = make_file_path(results_dir, query_vector.table_id, query_vector.set_id, nvec, total_data_file, data_gb_size, vector_length, query_time, total_checked_vec);
                            save_to_query_result_file(query_result_file, table_id, query_vector.set_id, total_knns, all_knn_results);


                            if(total_knns != 0)
                            {
                                struct result_sid * top = get_top_sets(all_knn_results, total_knns, num_top);
                                for(int m = 0; m < num_top; m++)
                                {
                                printf("table-%u-column-%u- in @@%s$ overlap=%u§\n", top[m].table_id, top[m].set_id, top[m].raw_data_file, top[m].overlap_size);
                                }
                                printf("\nquery_time=%fsec\n", query_time/1000000);
                                free(top);
                            }

                            for (int knn = 0; knn < (total_knns); knn++)
                                free(all_knn_results[knn].vector_id);
                            free(all_knn_results);
                            free(query_result_file);
                        }
                        else
                        {
                            printf("\n>>> performing knn search ...\n");
                            // Perform brute force knn search for all query vectors
                            all_knn_results =  brute_force_knn_search_optimized(dataset, total_data_file, vector_length, query_set, nvec, k, &total_checked_vec);

                            COUNT_PARTIAL_TIME_END
                            query_time = partial_time - (partial_input_time + partial_output_time);
                            // printf("\nquerytime  = %f - (%f + %f) = %f\n", partial_time, partial_input_time, partial_output_time, query_time);
                            
                            /* End of Query set */
                            /* Save query results to csv file */
                            char * query_result_file = make_file_path(results_dir, query_vector.table_id, query_vector.set_id, nvec, total_data_file, data_gb_size, vector_length, query_time, total_checked_vec);
                            save_to_query_result_file(query_result_file, table_id, query_vector.set_id, k * nvec, all_knn_results);

                            struct result_sid * top = get_top_sets(all_knn_results, k * nvec, num_top);
                            for(int m = 0; m < num_top; m++)
                            {
                                printf("table-%ucolumn-%u- in @@%s$ overlap=%u§\n", top[m].table_id, top[m].set_id, top[m].raw_data_file, top[m].overlap_size);
                            }
                            printf("\nquery_time=%fsec\n", query_time/1000000);
                            

                            for (int knn = 0; knn < (k*nvec); knn++)
                                free(all_knn_results[knn].vector_id);
                            free(all_knn_results);
                            free(top);
                            free(query_result_file);
                        }

                        RESET_PARTIAL_COUNTERS()
                        COUNT_PARTIAL_TIME_START

                        for(int t = 0; t < nvec; t++)
                        {
                            free(query_set[t].values);
                        }

                        i = 0; j = 0;
                        nvec = 0u;
                        query_vector.pos = 0;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            fclose(bin_file);
        }
    }
    COUNT_PARTIAL_INPUT_TIME_START
    closedir(dir);
    COUNT_PARTIAL_INPUT_TIME_END
    free(query_set);
    free(query_vector.values);
    free(results_dir);

    COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
}



struct query_result * brute_force_knn_search_optimized(char * dataset, unsigned int total_data_files, unsigned int vector_length, 
    struct vector * qset, unsigned int qnvec, unsigned int k, unsigned int *total_checked_vec)
{
    
    float * bsf = (float *) malloc(sizeof(float) * qnvec);
    unsigned int * curr_size = (unsigned int *) malloc(sizeof(unsigned int) * qnvec);
    //queue containing kNN results for every vector
    struct query_result * all_knn_results = (struct query_result *) malloc(sizeof(struct query_result) * k*qnvec);
    

    struct vector v; 
    v.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
    if (v.values == NULL)
    {
        printf("Error in bf.c: Couldn't allocate memory for temp vector values.");
        exit(1);
    }
    struct query_result * temp = malloc(sizeof(struct query_result));
    temp->vector_id = malloc(sizeof(struct vid));
    struct query_result * repl = malloc(sizeof(struct query_result));
    repl->vector_id = malloc(sizeof(struct vid));
    temp->distance = FLT_MAX;
    repl->distance = FLT_MAX;

    if (temp->vector_id == NULL || repl->vector_id == NULL)
    {
        printf("Error in bf.c: Couldn't allocate memory for temp query result.");
        exit(1);
    }

    for (int i = 0; i < qnvec; ++i)
    {
        bsf[i] = FLT_MAX;
        curr_size[i] = 0;
    }

    for (int j = 0; j < qnvec*k; ++j)
    {
        all_knn_results[j].vector_id = malloc(sizeof(struct vid));
        if(all_knn_results[j].vector_id == NULL)
        {
            printf("Error in bf.c: Couldn't allocate memory for vector id in knn results.");
            exit(1);
        }
        all_knn_results[j].vector_id->set_id = -1;
        all_knn_results[j].vector_id->table_id = -1;
        all_knn_results[j].vector_id->pos = -1;
        all_knn_results[j].query_vector_pos = -1;
        strcpy(all_knn_results[j].vector_id->raw_data_file, "");
        all_knn_results[j].distance = FLT_MAX;
    }

    // Open binary files  dir
    struct dirent *dfile;
    COUNT_PARTIAL_INPUT_TIME_START
    DIR *dir = opendir(dataset);
    COUNT_PARTIAL_INPUT_TIME_END
    if (!dir)
    {
      printf("Unable to open directory stream %s!", dataset);
      exit(1);
    }

    while ((dfile = readdir(dir)) != NULL && total_data_files > 0)
    {

        if(is_binaryfile(dfile->d_name))
        {
            char * raw_file_name = dfile->d_name;
            total_data_files--;
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = ""; 
            strcat(bin_file_path, dataset);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bf.c:  vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            COUNT_PARTIAL_INPUT_TIME_START
            FILE * bin_file = fopen (bin_file_path, "rb");
            COUNT_PARTIAL_INPUT_TIME_END
            if (bin_file == NULL)
            {
                printf("Error in bf.c: File %s not found!\n", bin_file_path);
                exit(1);
            }

            ts_type val; 
            uint32_t nvec = 0u;
            ts_type d = 0.0;
            unsigned int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            
            //read every candidate set and candidate vector in binary file
            while(total_bytes)
            {
                if(i == 0)
                { 
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;

                    v.table_id = table_id;
                    v.set_id = set_id;
                    v.pos = 0;
                    set_id = set_id +  1;

                    if(v.table_id == qset[0].table_id && v.set_id == qset[0].set_id)
                    {
                        // do not match query set to itself
                        COUNT_PARTIAL_INPUT_TIME_START
                        fseek(bin_file, nvec * 4 *vector_length, SEEK_CUR);
                        COUNT_PARTIAL_INPUT_TIME_END
                        i = 0;
                        j = 0;
                        total_bytes -= nvec * vector_length;
                        continue;
                    }
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    if(j > (vector_length - 1))
                    {    
                        j = 0;
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidean_distance(qset[h].values, v.values, vector_length);
                            // d = fabs(d);
                            ts_type kth_bsf_dist = all_knn_results[h*k+(k-1)].distance;
                            if (d <= kth_bsf_dist)
                            {
                                query_result_cpy_vector(repl, &v, qset[h].pos, d, raw_file_name);
                                queue_bounded_sorted_insert(&all_knn_results[h*k], repl, &curr_size[h], k);
                            }
                        }
                        v.pos += 1;
                    }
                    // read new float value
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;
                    v.values[j] = val;
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidean_distance(qset[h].values, v.values, vector_length);                               
                            // d = fabs(d);
                            ts_type kth_bsf_dist = all_knn_results[h*k+(k-1)].distance;
                            if (d <= kth_bsf_dist)
                            {
                                query_result_cpy_vector(repl, &v, qset[h].pos, d, raw_file_name);
                                queue_bounded_sorted_insert(&all_knn_results[h*k], repl, &curr_size[h], k);
                            }
                        }
                        
                        i = 0; j = 0;
                        nvec = 0u;
                        v.pos = 0;
                        continue;
                    }
                    i++;
                    j++;
                    
                }
            }
        COUNT_PARTIAL_INPUT_TIME_START
        fclose(bin_file);
        COUNT_PARTIAL_INPUT_TIME_END
        }
    }
    COUNT_PARTIAL_INPUT_TIME_START
    closedir(dir);
    COUNT_PARTIAL_INPUT_TIME_END
    free(temp->vector_id);
    free(repl->vector_id);
    free(temp);
    free(repl);
    free(bsf);
    free(v.values);
    free(curr_size);
    
    return all_knn_results;
}

// no k is specified, get all vectors with distance 0 to the the query vectors in the the query column.
struct query_result * brute_force_identical_knn_search_optimized(char * dataset, unsigned int total_data_files, unsigned int vector_length, 
    struct vector * qset, unsigned int qnvec, unsigned int *total_checked_vec, unsigned int * total_knns)
{
    
    float * bsf = (float *) malloc(sizeof(float) * qnvec);
    unsigned int * curr_size = (unsigned int *) malloc(sizeof(unsigned int) * qnvec);

    //queue containing kNN results for every vector
    struct query_result * all_knn_results = NULL;

    struct vector v; 
    v.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
    if (v.values == NULL)
    {
        printf("Error in bf.c: Couldn't allocate memory for temp vector values.");
        exit(1);
    }
    struct query_result * temp = malloc(sizeof(struct query_result));
    temp->vector_id = malloc(sizeof(struct vid));
    struct query_result * repl = malloc(sizeof(struct query_result));
    repl->vector_id = malloc(sizeof(struct vid));
    temp->distance = FLT_MAX;
    repl->distance = FLT_MAX;

    if (temp->vector_id == NULL || repl->vector_id == NULL)
    {
        printf("Error in bf.c: Couldn't allocate memory for temp query result.");
        exit(1);
    }

    for (int i = 0; i < qnvec; ++i)
    {
        bsf[i] = FLT_MAX;
        curr_size[i] = 0;
    }
    // Open binary files  dir
    struct dirent *dfile;
    COUNT_PARTIAL_INPUT_TIME_START
    DIR *dir = opendir(dataset);
    COUNT_PARTIAL_INPUT_TIME_END
    if (!dir)
    {
      printf("Unable to open directory stream %s!", dataset);
      exit(1);
    }

    while ((dfile = readdir(dir)) != NULL && total_data_files > 0)
    {

        if(is_binaryfile(dfile->d_name))
        {
            char * raw_file_name = dfile->d_name;
            total_data_files--;
            // get fill path of bin file
            char bin_file_path[PATH_MAX + 1] = ""; 
            strcat(bin_file_path, dataset);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // skip query table if found in data lake
            if(table_id == qset[0].table_id)
                continue;
            
            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bf.c:  vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            COUNT_PARTIAL_INPUT_TIME_START
            FILE * bin_file = fopen (bin_file_path, "rb");
            COUNT_PARTIAL_INPUT_TIME_END
            if (bin_file == NULL)
            {
                printf("Error in bf.c: File %s not found!\n", bin_file_path);
                exit(1);
            }

            ts_type val; 
            uint32_t nvec = 0u;
            ts_type d = 0.0;
            unsigned int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            
            //read every candidate set and candidate vector in binary file
            while(total_bytes)
            {
                if(i == 0)
                { 
                    i++;
                    j = 0;
                    //read first integer to check how many vactors in current set
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;

                    v.table_id = table_id;
                    v.set_id = set_id;
                    v.pos = 0;
                    set_id = set_id +  1;

                    // if(v.table_id == qset[0].table_id && v.set_id == qset[0].set_id)
                    // {
                    //     // do not match query set to itself
                    //     COUNT_PARTIAL_INPUT_TIME_START
                    //     fseek(bin_file, nvec * 4 *vector_length, SEEK_CUR);
                    //     COUNT_PARTIAL_INPUT_TIME_END
                    //     i = 0;
                    //     j = 0;
                    //     total_bytes -= nvec * vector_length;
                    //     continue;
                    // }
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    if(j > (vector_length - 1))
                    {    
                        j = 0;
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidean_distance(qset[h].values, v.values, vector_length);
                            if (d == 0.0)
                            {
                                printf("+ Exact match: (%u, %u, %u)\tmatches\t(%u, %u, %u)\n", 
                                        v.table_id, v.set_id, v.pos, qset[h].table_id, qset[h].set_id, qset[h].pos);
                                        
                                *total_knns += 1;
                                all_knn_results = (struct query_result *) realloc(all_knn_results, sizeof(struct query_result) * (*total_knns));
                            
                                all_knn_results[*total_knns - 1].vector_id = malloc(sizeof(struct vid));
                                if(all_knn_results[*total_knns - 1].vector_id == NULL)
                                {
                                    printf("Error in bf.c: Couldn't allocate memory for vector id in knn results.");
                                    exit(1);
                                }
                                all_knn_results[*total_knns - 1].vector_id->table_id = v.table_id;
                                all_knn_results[*total_knns - 1].vector_id->set_id = v.set_id;
                                all_knn_results[*total_knns - 1].vector_id->pos = v.pos;
                                all_knn_results[*total_knns - 1].query_vector_pos = qset[h].pos;
                                strcpy(all_knn_results[*total_knns - 1].vector_id->raw_data_file, raw_file_name);
                                all_knn_results[*total_knns - 1].distance = d;
                            }
                        }
                        v.pos += 1;
                    }
                    // read new float value
                    COUNT_PARTIAL_INPUT_TIME_START
                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    COUNT_PARTIAL_INPUT_TIME_END
                    total_bytes--;
                    v.values[j] = val;
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidean_distance(qset[h].values, v.values, vector_length);                               
                            if (d == 0.0)
                            {
                                // printf("++ new exact match.");
                                *total_knns += 1;
                                all_knn_results = (struct query_result *) realloc(all_knn_results, sizeof(struct query_result) * (*total_knns));
                            
                                all_knn_results[*total_knns - 1].vector_id = malloc(sizeof(struct vid));
                                if(all_knn_results[*total_knns - 1].vector_id == NULL)
                                {
                                    printf("Error in bf.c: Couldn't allocate memory for vector id in knn results.");
                                    exit(1);
                                }
                                all_knn_results[*total_knns - 1].vector_id->table_id = v.table_id;
                                all_knn_results[*total_knns - 1].vector_id->set_id = v.set_id;
                                all_knn_results[*total_knns - 1].vector_id->pos = v.pos;
                                all_knn_results[*total_knns - 1].query_vector_pos = qset[h].pos;
                                strcpy(all_knn_results[*total_knns - 1].vector_id->raw_data_file, raw_file_name);
                                all_knn_results[*total_knns - 1].distance = d;
                            }
                        }
                        
                        i = 0; j = 0;
                        nvec = 0u;
                        v.pos = 0;
                        continue;
                    }
                    i++;
                    j++;
                    
                }
            }
        COUNT_PARTIAL_INPUT_TIME_START
        fclose(bin_file);
        COUNT_PARTIAL_INPUT_TIME_END
        }
    }
    COUNT_PARTIAL_INPUT_TIME_START
    closedir(dir);
    COUNT_PARTIAL_INPUT_TIME_END
    free(temp->vector_id);
    free(repl->vector_id);
    free(temp);
    free(repl);
    free(bsf);
    free(v.values);
    free(curr_size);
    
    return all_knn_results;
}