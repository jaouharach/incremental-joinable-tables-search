//
//  ts.c
//  isaxlib
//
//  Created by Kostas Zoumpatianos on 3/7/12.
//  Copyright 2012 University of Trento. All rights reserved.
//
#include "../config.h"
#include "../globals.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ts.h"

/**
 This function converts a string of floats seperated by a delimeter into a ts 
 record of a size ts_size.
 @param char ts_str[]
 @param int ts_size
 @param const char * delims
 @return *ts
 */
enum response ts_parse_str(char ts_str[], ts_type * ts_out, int ts_size, const char * delims)
{
    int index=0;
    char *result = strtok( ts_str, delims );
	while( result != NULL ) {
   	        ts_out[index] =  atof(result);
		result = strtok( NULL, delims );
#ifdef SANITY_CHECK
        if (index >= ts_size)
        {
            fprintf(stderr, "Error in ts.c: Time series bigger than limit of %d.\n", ts_size);
            return FAILURE; 
        }
#endif
        index++;
	}
    free(result);
    return SUCCESS;
}

float ts_euclidean_distance(ts_type * t, ts_type * s, int size) {
    float distance = 0;
    while (size > 0) {
        size--;
        distance += (t[size] - s[size]) * (t[size] - s[size]);
    }
    //distance = sqrtf(distance);
    
    return distance;
}

// Caution: q is ordered
ts_type ts_warping_distance(ts_type * q, ts_type * t , int offset , int size, ts_type bsf, int * order, float warping)
{   
    // Window fraction = 5%
    float window_frac = warping;

    // unorder q into q_temp
    ts_type *q_temp = malloc(sizeof(ts_type) * size);
    for(int i = 0; i < size; ++i){
        q_temp[order[i]] = q[i];
    }

    // from https://github.com/ncanac/dtw/blob/master/dtw.c 
    
    /* DTW distance algorithm based on
    https://en.wikipedia.org/wiki/Dynamic_time_warping */
    const float LARGE_VALUE = __FLT_MAX__;
    const int min_window = 0;
    int i, j, minj, maxj, window;
    float dist, min;
    float **distances = malloc((size + 1) * sizeof(float *));
    for(i = 0; i < size + 1; ++i)
        distances[i] = malloc((size + 1) * sizeof(float));

    window = window_frac*size;

    for(i = 0; i <= size; ++i)
        for(j = 0; j <= size; ++j)
            distances[i][j] = LARGE_VALUE;
    distances[0][0] = 0.0;

    for(i = 0; i < size; ++i)
    {
        minj = i - window;
        if(minj < 0)
            minj = 0;
        maxj = i + window;
        if(maxj > size-1)
            maxj = size-1;
        // minimum cost in a row
        float rowmin = __FLT_MAX__;
        for(j = minj; j <= maxj; ++j)
        {
            dist = (q_temp[i] - t[j])*(q_temp[i] - t[j]);
            min = distances[i][j];
            if(min > distances[i][j+1])
                min = distances[i][j+1];
            if(min > distances[i+1][j])
                min = distances[i+1][j];
            distances[i+1][j+1] = dist + min;
            if(distances[i+1][j+1] < rowmin) 
                rowmin = distances[i+1][j+1];
        }
        // early termination
        if(rowmin > bsf){
            for(i = 0; i < size + 1; ++i)
                free(distances[i]);
            free(distances);
            free(q_temp);
            return rowmin;
        }
    }

    dist = distances[size][size];

    for(i = 0; i < size + 1; ++i)
        free(distances[i]);
    free(distances);
    free(q_temp);

    return dist;
}

ts_type ts_euclidean_distance_reordered(ts_type * q, ts_type * t , int j , int size ,ts_type bsf, int * order)
{
    int i;
    ts_type sum = 0;
    for ( i = 0 ; i < size && sum < bsf ; i++ )
    {
      ts_type x = t[order[i]];
      sum += (x-q[i])*(x-q[i]);      
    }
    return sum;
}

/** 
 This function prints a ts record of a size.
 @param ts *ts
 @param int size
*/
void ts_print(ts_type *ts, int size) 
{
    int i;
    for (i=0; i < size; i++) {
        printf("%lf", ts[i]);
    }
    printf("\n");
}
