// benchmark
#include <time.h>
#include <sys/time.h>

// timers 
double start;
double end;

struct timeval total_time_start; // query time and index creation
struct timeval total_query_time_start; // total query time
struct timeval query_time_start; // query time

struct timeval current_time;


double total_time;
double total_query_time; // for all queries
double query_time; // for one query column (query set)


#define INIT_STATS() total_time = 0;\
                    total_query_time = 0;\
                    query_time = 0;\
            

#define COUNT_NEW_QUERY_TIME(qt) total_query_time += qt;

#define RESET_QUERY_TIME() query_time = 0;


#define COUNT_TOTAL_TIME_START gettimeofday(&total_time_start, NULL);
#define COUNT_QUERY_TIME_START gettimeofday(&query_time_start, NULL);

#define COUNT_TOTAL_TIME_END  gettimeofday(&current_time, NULL); \
                                start = total_time_start.tv_sec * 1000000 + (total_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                total_time += (end - start);

#define COUNT_QUERY_TIME_END  gettimeofday(&current_time, NULL); \
                                start = query_time_start.tv_sec * 1000000 + (query_time_start.tv_usec); \
                                end = current_time.tv_sec * 1000000  + (current_time.tv_usec); \
                                query_time += (end - start);
