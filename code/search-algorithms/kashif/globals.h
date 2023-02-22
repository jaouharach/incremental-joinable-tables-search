//
//  globals.h
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2012 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef dstree_globals_h
#define dstree_globals_h

#include "config.h"
#include "stdlib.h"
#include "stdbool.h"
#include "errno.h"
#include "string.h"
#include "float.h"
#include <math.h>
//#include <jemalloc/jemalloc.h>


///// TYPES /////
typedef unsigned short segment_type;
typedef float ts_type;
typedef unsigned char label_type;
//typedef unsigned short label_type;

typedef unsigned long long file_position_type;
typedef unsigned long long root_mask_type;
typedef long long mean_stdev_range;

enum response {FAILURE=0, SUCCESS=1};
enum insertion_mode {PARTIAL = 1, 
                     TMP = 2, 
                     FULL = 4,
                     NO_TMP = 8};

enum buffer_cleaning_mode {FULL_CLEAN, TMP_ONLY_CLEAN, TMP_AND_TS_CLEAN};
enum node_cleaning_mode {DO_NOT_INCLUDE_CHILDREN = 0,
                         INCLUDE_CHILDREN = 1};

typedef bool boolean;

///// DEFINITIONS /////
#define MINVAL -2000000
#define MAXVAL 2000000
#define DELIMITER ' '
#define TRUE 1
#define FALSE 0
#define BUFFER_REALLOCATION_RATE  2 

///// GLOBAL VARIABLES /////
int FLUSHES;

///// BNECHMARKING /////
//#ifdef BENCHMARK
		#include <time.h>
		#include <sys/time.h>
 
        double tS;
        double tE;

        struct timeval total_time_start;
        struct timeval parse_time_start;
        struct timeval input_time_start;
        struct timeval output_time_start;
        struct timeval load_node_start;
        struct timeval current_time;
        
        struct timeval fetch_start;
        struct timeval fetch_check_start;
        double total_input_time;
        double total_output_time;
        double total_parse_time;
        double load_node_time;
        double total_time;


        struct timeval partial_time_start;
        struct timeval partial_input_time_start;
        struct timeval partial_output_time_start;
        struct timeval partial_load_node_time_start;         

        /* start kashif changes */
        struct timeval * now;
        struct timeval *thread_partial_time_start;
        struct timeval *thread_partial_input_time_start;
        struct timeval *thread_partial_output_time_start;

        double *thread_partial_time;
        double *thread_partial_input_time;
        double *thread_partial_output_time;

        /* end kashif changes */

        double partial_time;
        double partial_input_time;
        double partial_output_time;
        double partial_load_node_time;

        unsigned long long partial_seq_input_count;
        unsigned long long partial_seq_output_count;
        unsigned long long partial_rand_input_count;
        unsigned long long partial_rand_output_count;

        unsigned long total_nodes_count;
        unsigned long leaf_nodes_count;
        unsigned long empty_leaf_nodes_count;
        unsigned long loaded_nodes_count;
        unsigned long checked_nodes_count;
        unsigned long loaded_ts_count;
        unsigned long checked_ts_count;
        unsigned long total_ts_count;
        unsigned long total_queries_count;

        /* start kashif changes */
        unsigned long *thread_loaded_nodes_count;
        unsigned long *thread_checked_nodes_count;
        unsigned long *thread_loaded_ts_count;
        unsigned long *thread_checked_ts_count;
        /* start kashif changes */

        ts_type total_tlb;

        #define INIT_STATS() total_input_time = 0;\
                             total_output_time = 0;\
                             total_time = 0;\
                             total_parse_time = 0;\
                             load_node_time= 0;\
			     partial_time = 0;\			     
			     partial_input_time = 0;\
			     partial_output_time = 0;\
			     partial_load_node_time = 0;\
                             partial_seq_input_count = 0;\
                             partial_seq_output_count = 0;\
                             partial_rand_input_count = 0;\
                             partial_rand_output_count = 0;\
			     total_nodes_count = 0;\
			     leaf_nodes_count = 0;\
			     empty_leaf_nodes_count = 0;\
			     loaded_nodes_count = 0;\
			     loaded_ts_count = 0;\			     
			     checked_ts_count = 0;\
                             checked_nodes_count = 0;\			    
			     total_ts_count = 0;

        #define INIT_THREAD_STATS(num_threads) thread_partial_time = calloc(num_threads, sizeof(double));\			     
			     thread_partial_input_time = calloc(num_threads, sizeof(double));\
			     thread_partial_output_time = calloc(num_threads, sizeof(double));\
			     thread_loaded_nodes_count = calloc(num_threads, sizeof(unsigned long));\
			     thread_loaded_ts_count = calloc(num_threads, sizeof(unsigned long));\			     
			     thread_checked_ts_count = calloc(num_threads, sizeof(unsigned long));\
                             thread_checked_nodes_count = calloc(num_threads, sizeof(unsigned long));\
                             thread_partial_time_start = calloc(num_threads, sizeof(struct timeval));\
                             thread_partial_input_time_start = calloc(num_threads, sizeof(struct timeval));\
                             thread_partial_output_time_start = calloc(num_threads, sizeof(struct timeval));\
                             now = calloc(num_threads, sizeof(struct timeval));


//                             filtered_ts_disk_count = 0;	\
//                            filtered_ts_mem_count = 0;	\
//                             total_filtered_ts_count = 0;	\
//                            refined_ts_disk_count = 0;	\
//                             refined_ts_mem_count = 0;	\
//                             total_refined_ts_count = 0;


/*
        printf("input\t output\t parse\t nodes\t checked_nodes\t loaded_nodes\t loaded_records\t distance\t load_node_time\t total\n");
        #define PRINT_STATS(result_distance) printf("%lf\t %lf\t %lf\t %d\t %d\t %d\t %lld\t %lf\t %lf\t %lf\n", \
        total_input_time, total_output_time, \
        total_parse_time, total_tree_nodes, \
        checked_nodes, loaded_nodes, \
        loaded_records, result_distance, load_node_time, total_time);
*/
/*
#define PRINT_STATS(idx_size, \
		    avg_fill, \
		    min_fill, \
		    max_fill, \
		    avg_ht, \
		    min_ht, \
		    max_ht,\
		    true_distance, \
		    lb_distance,  \
                    loaded_ts_count, \
                    total_ts_count); \
*/
  
/*
        #define COUNT_NEW_NODE() total_nodes_count++; 
        #define COUNT_LEAF_NODE() leaf_nodes_count++;
        #define COUNT_LOADED_NODE() loaded_nodes_count++;
        #define COUNT_CHECKED_NODE() checked_nodes_count++;
        #define COUNT_LOADED_RECORD() loaded_records_count++;
*/
        #define COUNT_NEW_NODE ++total_nodes_count; 
        #define COUNT_LEAF_NODE ++leaf_nodes_count;
        #define COUNT_EMPTY_LEAF_NODE ++empty_leaf_nodes_count;
        #define COUNT_TOTAL_TS(num_ts) total_ts_count+=num_ts; //actual ts inserted in index

        #define COUNT_CHECKED_NODE ++checked_nodes_count;
        #define COUNT_LOADED_NODE ++loaded_nodes_count;
        #define COUNT_LOADED_TS(num_ts) loaded_ts_count +=num_ts; //ts loaded to answer query
        #define COUNT_CHECKED_TS(num_ts) checked_ts_count +=num_ts; //ts loaded to answer query


        /* start kashif changes */
        #define COUNT_THREAD_CHECKED_NODE(thread_id) ++thread_checked_nodes_count[thread_id];
        #define COUNT_THREAD_LOADED_NODE(thread_id) ++thread_loaded_nodes_count[thread_id];
        #define COUNT_THREAD_LOADED_TS(num_ts,thread_id) thread_loaded_ts_count[thread_id] +=num_ts; //ts loaded to answer query
        #define COUNT_THREAD_CHECKED_TS(num_ts,thread_id) thread_checked_ts_count[thread_id] +=num_ts; //ts loaded to answer query
        /* end kashif changes */

/*
        #define COUNT_FILTERED_TS_DISK(num_ts) filtered_ts_disk_count +=num_ts;
        #define COUNT_FILTERED_TS_MEM(num_ts) filtered_ts_mem_count +=num_ts;
        #define COUNT_REFINED_TS_DISK(num_ts) refined_ts_disk_count +=num_ts;
        #define COUNT_REFINED_TS_MEM(num_ts) refined_ts_mem_count +=num_ts;
*/


/*
        #define RESET_QUERY_COUNTERS loaded_nodes_count = 0;\
                                     checked_nodes_count = 0;\ 
                                     filtered_ts_disk_count = 0;\
                                     filtered_ts_mem_count = 0;\
                                     refined_ts_disk_count = 0;	\
                                     refined_ts_mem_count = 0;\
                                     exact_distance = 0;\
                                     approx_distance = 0;				     
*/

      #define RESET_QUERY_COUNTERS() loaded_nodes_count = 0;\
                                     loaded_ts_count = 0;\
                                     checked_nodes_count = 0;\
                                     checked_ts_count = 0;
				     
      #define RESET_PARTIAL_COUNTERS() partial_seq_output_count = 0;\
                                       partial_seq_input_count = 0;\
                                       partial_rand_output_count = 0;\
                                       partial_rand_input_count = 0;\
				       partial_input_time = 0;\
				       partial_output_time = 0;\
				       partial_load_node_time = 0;\
				       partial_time = 0;

        /* start kashif changes */
        #define RESET_THREAD_QUERY_COUNTERS(thread_id) thread_loaded_nodes_count[thread_id] = 0;\
                                        thread_loaded_ts_count[thread_id] = 0;\
                                        thread_checked_nodes_count[thread_id] = 0;\
                                        thread_checked_ts_count[thread_id] = 0;

        #define FREE_THREAD_QUERY_COUNTERS free(thread_loaded_nodes_count);\
                                        free(thread_loaded_ts_count);\
                                        free(thread_checked_nodes_count);\
                                        free(thread_checked_ts_count);\
                                        thread_loaded_nodes_count = NULL;\
                                        thread_loaded_ts_count = NULL;\
                                        thread_checked_nodes_count = NULL;\
                                        thread_checked_ts_count = NULL;


        #define RESET_THREAD_PARTIAL_COUNTERS(thread_id) thread_partial_input_time[thread_id] = 0;\
                                        thread_partial_output_time[thread_id] = 0;\
                                        thread_partial_time[thread_id] = 0;\

        #define FREE_THREAD_PARTIAL_COUNTERS free(thread_partial_input_time_start);\
                                        free(thread_partial_output_time_start);\
                                        free(thread_partial_time_start);\
                                        free(thread_partial_input_time);\
                                        free(thread_partial_output_time);\
                                        free(thread_partial_time);\
                                        free(now);\
                                        thread_partial_input_time_start = NULL;\
                                        thread_partial_output_time_start = NULL;\
                                        thread_partial_time_start = NULL;\
                                        thread_partial_input_time = NULL;\
                                        thread_partial_output_time = NULL;\
                                        thread_partial_time = NULL;\
                                        now = NULL;     
        /* end kashif changes */
       #define PRINT_QUERY_COUNTERS() printf("loaded_nodes and checked_ts = %lu and %lu\n", loaded_nodes_count, checked_ts_count);
       #define PRINT_PARTIAL_COUNTERS() printf("seq_output and partial_time = %llu and %f\n",partial_seq_output_count,  partial_time);

        #define COUNT_PARTIAL_SEQ_INPUT ++partial_seq_input_count;
        #define COUNT_PARTIAL_SEQ_OUTPUT ++partial_seq_output_count;
        #define COUNT_PARTIAL_RAND_INPUT ++partial_rand_input_count;
        #define COUNT_PARTIAL_RAND_OUTPUT ++partial_rand_output_count;


        #define COUNT_INPUT_TIME_START gettimeofday(&input_time_start, NULL);   
        #define COUNT_OUTPUT_TIME_START gettimeofday(&output_time_start, NULL); 
        #define COUNT_TOTAL_TIME_START gettimeofday(&total_time_start, NULL);

        #define COUNT_PARTIAL_TIME_START gettimeofday(&partial_time_start, NULL);
        #define COUNT_PARTIAL_INPUT_TIME_START gettimeofday(&partial_input_time_start, NULL);   
        #define COUNT_PARTIAL_OUTPUT_TIME_START gettimeofday(&partial_output_time_start, NULL);
        #define COUNT_PARTIAL_LOAD_NODE_TIME_START gettimeofday(&partial_load_node_time_start, NULL);

        /* start kashif changes */

        #define COUNT_THREAD_PARTIAL_TIME_START(thread_id) gettimeofday(&thread_partial_time_start[thread_id], NULL);
        #define COUNT_THREAD_PARTIAL_INPUT_TIME_START(thread_id) gettimeofday(&thread_partial_input_time_start[thread_id], NULL);   
        #define COUNT_THREAD_PARTIAL_OUTPUT_TIME_START(thread_id) gettimeofday(&thread_partial_output_time_start[thread_id], NULL);

        /* end kashif changes */
        #define COUNT_PARSE_TIME_START gettimeofday(&parse_time_start, NULL);   
        #define COUNT_LOAD_NODE_START gettimeofday(&load_node_start, NULL);
        #define COUNT_INPUT_TIME_END  gettimeofday(&current_time, NULL);\
                                      tS = input_time_start.tv_sec*1000000 + (input_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      total_input_time += (tE - tS); 
        #define COUNT_OUTPUT_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = output_time_start.tv_sec*1000000 + (output_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_output_time += (tE - tS); 
        #define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = total_time_start.tv_sec*1000000 + (total_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_time += (tE - tS);
        #define COUNT_PARTIAL_INPUT_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = partial_input_time_start.tv_sec*1000000 + (partial_input_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      partial_input_time += (tE - tS);
        #define COUNT_PARTIAL_LOAD_NODE_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = partial_load_node_time_start.tv_sec*1000000 + (partial_load_node_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                      partial_load_node_time += (tE - tS); 
        #define COUNT_PARTIAL_OUTPUT_TIME_END gettimeofday(&current_time, NULL); \
                                      tS = partial_output_time_start.tv_sec*1000000 + (partial_output_time_start.tv_usec); \
				      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      partial_output_time += (tE - tS); 
        #define COUNT_PARTIAL_TIME_END  gettimeofday(&current_time, NULL); \
                                      tS = partial_time_start.tv_sec*1000000 + (partial_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      partial_time += (tE - tS); 

        /* start kashif changes */
        #define COUNT_THREAD_PARTIAL_INPUT_TIME_END(thread_id)  gettimeofday(&now[thread_id], NULL); \
                                      tS = thread_partial_input_time_start[thread_id].tv_sec*1000000 + (thread_partial_input_time_start[thread_id].tv_usec); \
                                      tE = now[thread_id].tv_sec*1000000 + (now[thread_id].tv_usec); \
                                      thread_partial_input_time[thread_id] += (tE - tS);
        #define COUNT_THREAD_PARTIAL_LOAD_NODE_TIME_END(thread_id)  gettimeofday(&now[thread_id], NULL); \
                                      tS = thread_partial_load_node_time_start[thread_id].tv_sec*1000000 + (thread_partial_load_node_time_start[thread_id].tv_usec); \
                                      tE = now[thread_id].tv_sec*1000000 + (now[thread_id].tv_usec); \
                                      thread_partial_load_node_time[thread_id] += (tE - tS); 
        #define COUNT_THREAD_PARTIAL_OUTPUT_TIME_END(thread_id) gettimeofday(&now[thread_id], NULL); \
                                      tS = thread_partial_output_time_start[thread_id].tv_sec*1000000 + (thread_partial_output_time_start[thread_id].tv_usec); \
				      tE = now[thread_id].tv_sec*1000000  + (now[thread_id].tv_usec); \
                                      thread_partial_output_time[thread_id] += (tE - tS); 
        #define COUNT_THREAD_PARTIAL_TIME_END(thread_id) gettimeofday(&now[thread_id], NULL); \
                                      tS = thread_partial_time_start[thread_id].tv_sec*1000000 + (thread_partial_time_start[thread_id].tv_usec); \
                                      tE = now[thread_id].tv_sec*1000000  + (now[thread_id].tv_usec); \
                                      thread_partial_time[thread_id] += (tE - tS); \
        /* end kashif changes */

        #define COUNT_PARSE_TIME_END  gettimeofday(&current_time, NULL);  \
                                      tS = parse_time_start.tv_sec*1000000 + (parse_time_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      total_parse_time += (tE - tS); 
        #define COUNT_LOAD_NODE_END   gettimeofday(&current_time, NULL);  \
                                      tS = load_node_start.tv_sec*1000000 + (load_node_start.tv_usec); \
                                      tE = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                      load_node_time += (tE - tS); 
/*
   #else
        #define INIT_STATS() ;
        #define PRINT_STATS() ;
        #define COUNT_INPUT_TIME_START ;
        #define COUNT_INPUT_TIME_END ;
        #define COUNT_OUTPUT_TIME_START ;
        #define COUNT_OUTPUT_TIME_END ;
        #define COUNT_TOTAL_TIME_START ;
        #define COUNT_TOTAL_TIME_END ;

        #define COUNT_CHECKED_NODE;
        #define COUNT_LOADED_NODE;
        #define COUNT_LOADED_TS(num_ts);
        #define COUNT_CHECKED_TS(num_ts);

        #define COUNT_NEW_NODE; 
        #define COUNT_LEAF_NODE;
        #define COUNT_EMPTY_LEAF_NODE;

        #define COUNT_PARTIAL_TIME_START ;
        #define COUNT_PARTIAL_TIME_END ;
        #define COUNT_PARTIAL_INPUT_TIME_START ;
        #define COUNT_PARTIAL_INPUT_TIME_END ;
        #define COUNT_PARTIAL_OUTPUT_TIME_START ;
        #define COUNT_PARTIAL_OUTPUT_TIME_END ;
        #define COUNT_PARTIAL_LOAD_NODE_TIME_START ;
        #define COUNT_PARTIAL_LOAD_NODE_TIME_END ;
        #define RESET_QUERY_COUNTERS;
        #define RESET_PARTIAL_COUNTERS;
        #define COUNT_PARSE_TIME_START ;
        #define COUNT_PARSE_TIME_END ;
        #define COUNT_PARTIAL_SEQ_INPUT ;
        #define COUNT_PARTIAL_SEQ_OUTPUT;
        #define COUNT_PARTIAL_RAND_INPUT;
        #define COUNT_PARTIAL_RAND_OUTPUT;

    #endif
*/
#endif
