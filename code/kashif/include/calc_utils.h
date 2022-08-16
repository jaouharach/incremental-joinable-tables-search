//
//  cals_utils.h
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2012 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#ifndef dstreelib_dstree_calc_utils_h
#define dstreelib_dstree_calc_utils_h

#include "../config.h"
#include "../globals.h"
#include "math.h"
#include "dstree_file_buffer_manager.h"

ts_type calc_mean (ts_type *, int start, int end);
ts_type calc_stdev (ts_type *, int start, int end);
int compare(const void *a, const void *b);
short compare_short (const void * a, const void * b);
//short compare_file_map_entries (const void * a, const void * b);
void get_current_time(char * time_buf);
short compare_file_buffer (const void * a, const void * b);
ts_type calculate_mean_std_dev_range(struct segment_sketch sketch, int len);
int get_segment_start(short * points, int idx);
int get_segment_end(short * points, int idx);
int get_segment_length(short * points, int i);
int znorm_comp(const void *a, const void* b);
ts_type get_max(ts_type * series, int start, int end);
void assign_max_segments(ts_type * u, ts_type * u_hat, short * points, int num_points);
ts_type get_min(ts_type * series, int start, int end);
void assign_min_segments(ts_type * l, ts_type * l_hat, short * points, int num_points);
void lower_upper_lemire(float *t, int len, int r, float *l, float *u);
// float get_dtw(float *q_temp, int *order1, float*buffer, int r, int lengthSeriesDataset, float bsf);

#endif
