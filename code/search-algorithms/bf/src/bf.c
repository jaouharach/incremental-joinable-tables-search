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
#include "../include/stats.h"

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

        default:
            exit(-1);
            break;
        }
    }
    data_gb_size = get_data_gb_size(dataset, total_data_files);
    printf("Total data size in gb = %u\n", data_gb_size);    

    bf_sequential_search(dataset, vector_length, qset_num, min_qset_size, max_qset_size,
                            top, result_dir, total_data_files, data_gb_size, k);

    printf("\n>>>> Congrat! End of experiment."); 
    printf("\n>>>> Results are stored in %s", result_dir); 
    COUNT_TOTAL_TIME_END
    printf("\nSanity check: combined indexing and querying times should be less "
            "than: %f secs \n",
            total_time / 1000000);
    // printf("\n>>>> Summary:\nDatalake size:\t\t\t%uGB\n\nTotal table:\t\t\t%u\n\nNumber of queries:\t\t%u\n\nMin query size:\t\t\t%u\n\nMax query size:\t\t\t%u\n", data_gb_size, total_data_files, qset_num, min_qset_size, max_qset_size);
    return 0;
}


void bf_sequential_search(char * dataset, unsigned int vector_length, unsigned int qset_num,
    unsigned int min_qset_size, unsigned int max_qset_size, unsigned int num_top, char * result_dir,
    unsigned int total_data_file, unsigned int data_gb_size, unsigned int k)
{
    int opened_files = 0; 
    char * algorithm = "bfed";
    // open source dir
    struct dirent *dfile;
    DIR *dir = opendir(dataset);

    if (!dir)
    {
        printf("Unable to open directory stream %s!", dataset);
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
            strcat(bin_file_path, dataset);strcat(bin_file_path, "/");strcat(bin_file_path, dfile->d_name);

            // get binary table info 
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name,"data_size%d_t%dc%d_len%d_noznorm.bin",&datasize,&table_id,&nsets,&vector_length_in_filename);

            // check if vector length in file name matches vector length passed as argument
            if(vector_length_in_filename != vector_length)
            {
                printf("Error in bfed.c:  vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            FILE * bin_file = fopen (bin_file_path, "rb");
            if (bin_file == NULL)
            {
                printf("Error in bfed.c: File %s not found!\n", bin_file_path);
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
                    query_set = (struct vector *) realloc(query_set, sizeof(struct vector) * nvec);
                    printf("Query %u/%u.\n\n", (total_queries-qset_num)+1, total_queries);
                    
                    RESET_QUERY_TIME()
                    
                    total_checked_vec = 0;
                    next_vec = 0;
                    qset_num--;

                    query_vector.table_id = table_id;
                    query_vector.set_id = set_id;
                    query_vector.pos = 0;
                    set_id = set_id + 1;
                    i++;
                    j = 0;
                }
                else if(i <= (unsigned int)nvec*vector_length)
                {
                    // end of vector but still in current set
                    if(j > vector_length - 1){
                        j = 0; 
                        
                        // Append vector to query set
                        query_set[next_vec] = query_vector;
                        query_vector.pos += 1;
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
                        all_knn_results =  brute_force_knn_search_optimized(dataset, total_data_file, vector_length, query_set, nvec, k, &total_checked_vec);

                        /* End of Query set */
                        /* Save query results to csv file */
                        char * query_result_file = make_file_path(results_dir, query_vector.table_id, query_vector.set_id, nvec, total_data_file, data_gb_size, vector_length, query_time, total_checked_vec);
                        save_to_query_result_file(query_result_file, table_id, query_vector.set_id, k * nvec, all_knn_results);

                        struct result_sid * top = get_top_sets(all_knn_results, k * nvec, num_top);
                        for(int m = 0; m < num_top; m++)
                        {
                        printf("column-%u- in @@%s$ overlap=%u§\n", top[m].set_id, top[m].raw_data_file, top[m].overlap_size);
                        }
                        printf("\nquery_time=%fsec\n", query_time/1000000);
                        

                        for (int knn = 0; knn < (k*nvec); knn++)
                            free(all_knn_results[knn].vector_id);
                        free(all_knn_results);
                        free(top);
                        free(query_result_file);

                        i = 0; j = 0;
                        nvec = 0u;
                        query_vector.pos = 0;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
        }
    }
    closedir(dir);
    free(query_set);
    free(query_vector.values);
    free(results_dir);
}



struct query_result * brute_force_knn_search_optimized(char * dataset, unsigned int total_data_files, unsigned int vector_length, 
    struct vector * qset, unsigned int qnvec, unsigned int k, unsigned int *total_checked_vec)
{
    ts_type * bsf = (ts_type *) malloc(sizeof(ts_type) * qnvec);
    unsigned int * next_knn = (unsigned int *) malloc(sizeof(unsigned int) * qnvec);
    //queue containing kNN results for every vector
    struct query_result * curr_knn = (struct query_result *) malloc(sizeof(struct query_result) * k*qnvec);
    

    struct vector v; 
    v.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);


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
        curr_knn[j].vector_id->set_id = -1;
        curr_knn[j].vector_id->table_id = -1;
        curr_knn[j].vector_id->pos = -1;
        curr_knn[j].query_vector_pos = -1;
        strcpy(curr_knn[j].vector_id->raw_data_file, "");
        curr_knn[j].distance = FLT_MAX;
    }

    // Open binary files  dir
    struct dirent *dfile;
    DIR *dir = opendir(dataset);

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
                printf("Error in bfed.c:  vector length passed in argumentes (--len %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
                exit(1);
            }

            /* read binary file */
            FILE * bin_file = fopen (bin_file_path, "rb");
            if (bin_file == NULL)
            {
                printf("Error in bfed.c: File %s not found!\n", bin_file_path);
                exit(1);
            }

            ts_type val; 
            uint32_t nvec = 0u;
            ts_type d = 0.0;
            unsigned int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            
            //read every candidate set and candidate vector in binary file
            RESET_QUERY_TIME()
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
                    v.pos = 0;
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
                    COUNT_QUERY_TIME_START
                    if(j > vector_length - 1)
                    {    
                        j = 0;
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidian_distance(qset[h].values, v.values, vector_length);

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
                                    curr_knn[h*k+m].vector_id->pos = curr_knn[(h*k+m)-1].vector_id->pos;
                                    curr_knn[h*k+m].query_vector_pos = curr_knn[(h*k+m)-1].query_vector_pos;
                                    strcpy(curr_knn[h*k+m].vector_id->raw_data_file, curr_knn[(h*k+m)-1].vector_id->raw_data_file);
                                }

                                curr_knn[h*k].distance = bsf[h];
                                curr_knn[h*k].vector_id->set_id = v.set_id;
                                curr_knn[h*k].vector_id->table_id = v.table_id;
                                curr_knn[h*k].vector_id->pos = v.pos;
                                curr_knn[h*k].query_vector_pos = qset[h].pos;
                                strcpy(curr_knn[h*k].vector_id->raw_data_file, raw_file_name);
                                
                                next_knn[h] = 1;

                                
                            }
                            // if its another vector with same bsf
                            else if(d == bsf[h] && next_knn[h] < k)
                            {
                                curr_knn[(h*k)+next_knn[h]].distance = bsf[h];
                                curr_knn[(h*k)+next_knn[h]].vector_id->set_id = v.set_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->table_id = v.table_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->pos = v.pos;
                                curr_knn[(h*k)+next_knn[h]].query_vector_pos = qset[h].pos;
                                strcpy(curr_knn[(h*k)+next_knn[h]].vector_id->raw_data_file, raw_file_name);
                                next_knn[h] = next_knn[h] + 1;
                            }
                        }
                        v.pos += 1;
                    }
                    COUNT_QUERY_TIME_END
                    // read new float value
                    fread((void*)(&val), sizeof(val), 1, bin_file);
                    COUNT_QUERY_TIME_START
                    total_bytes--;
                    v.values[j] = val;
                    if(i == (unsigned int)nvec*vector_length)
                    {   
                        *total_checked_vec += qnvec;
                        for(int h = 0; h < qnvec; h++)
                        {
                            d = euclidian_distance(qset[h].values, v.values, vector_length);

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
                                    curr_knn[h*k+m].vector_id->pos = curr_knn[(h*k+m)-1].vector_id->pos;
                                    curr_knn[h*k+m].query_vector_pos = curr_knn[(h*k+m)-1].query_vector_pos;
                                    strcpy(curr_knn[h*k+m].vector_id->raw_data_file, curr_knn[(h*k+m)-1].vector_id->raw_data_file);
                                }

                                curr_knn[h*k].distance = bsf[h];
                                curr_knn[h*k].vector_id->set_id = v.set_id;
                                curr_knn[h*k].vector_id->table_id = v.table_id;
                                curr_knn[h*k].vector_id->pos = v.pos;
                                curr_knn[h*k].query_vector_pos = qset[h].pos;
                                strcpy(curr_knn[h*k].vector_id->raw_data_file, raw_file_name);
                                next_knn[h] = 1;

                                
                            }
                            // if its another vector with same bsf
                            else if(d == bsf[h] && next_knn[h] < k)
                            {
                                curr_knn[(h*k)+next_knn[h]].distance = bsf[h];
                                curr_knn[(h*k)+next_knn[h]].vector_id->set_id = v.set_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->table_id = v.table_id;
                                curr_knn[(h*k)+next_knn[h]].vector_id->pos = v.pos;
                                curr_knn[(h*k)+next_knn[h]].query_vector_pos = qset[h].pos;
                                strcpy(curr_knn[(h*k)+next_knn[h]].vector_id->raw_data_file, raw_file_name);
                                next_knn[h] = next_knn[h] + 1;
                            }
                        }
                        

                        i = 0; j = 0;
                        nvec = 0u;
                        v.pos = 0;
                        continue;
                    }
                    i++;
                    j++;
                    COUNT_QUERY_TIME_END
                }
            }
        fclose(bin_file);
        }
    }
    closedir(dir);
    free(bsf);
    free(v.values);
    free(next_knn);
    return curr_knn;
}
