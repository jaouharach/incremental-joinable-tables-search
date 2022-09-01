// benchmark
#include <time.h>
#include <sys/time.h>

// timers 
double start;
double end;

struct timeval total_time_start; // query time and index creation
struct timeval total_query_time_start; // total query time
struct timeval query_time_start; // query time

struct timeval partial_time_start;
struct timeval partial_input_time_start;
struct timeval partial_output_time_start;

struct timeval current_time;


double total_time;
double total_query_time; // for all queries
double query_time; // for one query column (query set)

double partial_time;
double partial_input_time;
double partial_output_time;

#define INIT_STATS() total_time = 0;\
                    total_query_time = 0;\
                    query_time = 0;\
                    partial_input_time = 0;\
                    partial_output_time = 0;\
                    partial_time = 0;
            

#define COUNT_NEW_QUERY_TIME(qt) total_query_time += qt;

#define RESET_QUERY_TIME() query_time = 0;

#define RESET_PARTIAL_COUNTERS() partial_input_time = 0;\
                                partial_output_time = 0;\
                                partial_time = 0;

#define COUNT_PARTIAL_INPUT_TIME_START gettimeofday(&partial_input_time_start, NULL);   
#define COUNT_PARTIAL_OUTPUT_TIME_START gettimeofday(&partial_output_time_start, NULL); 
#define COUNT_TOTAL_TIME_START gettimeofday(&total_time_start, NULL);
#define COUNT_PARTIAL_TIME_START gettimeofday(&partial_time_start, NULL);
#define COUNT_QUERY_TIME_START gettimeofday(&query_time_start, NULL);


#define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = total_time_start.tv_sec * 1000000 + (total_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                total_time += (end - start);

#define COUNT_PARTIAL_INPUT_TIME_END  gettimeofday(&current_time, NULL);\
                                start = partial_input_time_start.tv_sec*1000000 + (partial_input_time_start.tv_usec); \
                                end = current_time.tv_sec*1000000 + (current_time.tv_usec); \
                                partial_input_time += (end - start); 
#define COUNT_PARTIAL_OUTPUT_TIME_END gettimeofday(&current_time, NULL); \
                                start = partial_output_time_start.tv_sec*1000000 + (partial_output_time_start.tv_usec); \
                                end = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                partial_output_time += (end - start); 
#define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = total_time_start.tv_sec*1000000 + (total_time_start.tv_usec); \
                                end = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                total_time += (end - start);
#define COUNT_PARTIAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = partial_time_start.tv_sec*1000000 + (partial_time_start.tv_usec); \
                                end = current_time.tv_sec*1000000  + (current_time.tv_usec); \
                                partial_time += (end - start); 

#define COUNT_QUERY_TIME_END  gettimeofday(&current_time, NULL); \
                                start = query_time_start.tv_sec * 1000000 + (query_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                query_time += (end - start);
