//
//  calc_utils.c
//  ds-tree C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//



#include "../config.h"
#include "../globals.h"
#include "../include/calc_utils.h"
#include "math.h"


ts_type calc_mean (ts_type * series, int start, int end)
{

  ts_type mean = 0;

  if (start >= end)
  {
    int j = 0;
    j++;
    printf("error start > end \n");
  }
  else
  {
   for (int i=start; i < end; i++) 
   {
     mean += series[i];
   }
  
   mean /= (end-start); 
  }

  return mean;
  
}

/*
  Using the stdev computational formula.

*/

ts_type calc_stdev (ts_type * series, int start, int end)
{

  ts_type sum_x_squares=0, sum_x=0, stdev = 0; //sum of x's and sum of x squares
  int i, count_x;

  if (start >= end)
  {
    printf ("error in stdev start >= end\n");
  }
  else
  {
    count_x = end-start; //size of the series
  
    for (int i=start; i<end; i++) 
    {
     sum_x += series[i];
     sum_x_squares += pow(series[i],2);
    }
  
    sum_x_squares -= (pow(sum_x,2) / count_x);

    stdev = sqrt(sum_x_squares/count_x);
  }
  
  return stdev;
  
}

ts_type calc_mean_per_segment (ts_type * series, short * segments, ts_type *means, int size)
{

  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    means[i] = calc_mean (series, start, end);
    start = end;
  }
  
}

ts_type calc_stdev_per_segment (ts_type * series, short * segments, ts_type *stdevs, int size)
{

  int start=0, end;
  for (int i=0; i< size; i++) 
  {
    end = segments[i];
    stdevs[i] = calc_stdev (series, start, end);
    start = end;
  }
  
}



/* 
 This is the compare function used in the binary search code
 */

  

ts_type compare_ts_type (const void * a, const void * b)
{
  return ( *(ts_type*)a - *(ts_type*)b );
}

  

short compare_short (const void * a, const void * b)
{
  if (*(short*)a < *(short*)b )
    return -1;
  else if (*(short*)a == *(short*)b )
    return 0;
  else
    return 1;

}
/*
short compare_file_map_entries (const void * a, const void * b)
{
  char * entry_a = (char *) a;
  struct dstree_file_map *entry_b = (struct dstree_file_map*) b;

  return ( strcmp(entry_a, entry_b->filename));

}
*/
short compare_file_buffer (const void * a, const void * b)
{
  struct dstree_file_buffer * entry_a = *((struct dstree_file_buffer**) a);
  struct dstree_file_buffer * entry_b = *((struct dstree_file_buffer**) b);

  if (entry_a->buffered_list_size < entry_b->buffered_list_size )
    return 1;
  else if  (entry_a->buffered_list_size == entry_b->buffered_list_size )
    return 0;
  else
    return -1;
}

/*
  returns the current time in string format.
*/

void get_current_time(char * time_buf)
{
    time_t timer;
    
    struct tm* tm_info;

    time(&timer);
    tm_info = localtime(&timer);

    strftime(time_buf, 26, "%Y-%m-%d %H:%M:%S", tm_info);
}

ts_type calculate_mean_std_dev_range(struct segment_sketch sketch, int len){

  ts_type mean_width = sketch.indicators[0]-sketch.indicators[1];

  ts_type stdev_upper = sketch.indicators[2];
  ts_type stdev_lower = sketch.indicators[3];

  return (len * (mean_width * mean_width + stdev_upper * stdev_upper));
  
}


int get_segment_start(short * points, int idx)
{
  if (idx == 0 )
    return 0;
  else
    return points[idx-1];
}

int get_segment_end(short * points, int idx)
{
  return points[idx];
}

int get_segment_length(short * points, int i)
{

  if (i == 0)
    return points[i];
  else
    return points[i] - points[i-1];
}

ts_type get_max(ts_type * series, int start, int end)
{
    ts_type max = series[start];

    if (start >= end)
    {
      int j = 0;
      j++;
      printf("error start > end \n");
    }
    else
    {
      for (int i=start; i < end; i++) 
      {
        if(series[i] > max)
          max = series[i];
      }
    }
    return max;
}

void assign_max_segments(ts_type * u, ts_type * u_hat, short * points, int num_points)
{
    int start=0, end;
    for (int i=0; i< num_points; i++) 
    {
      end = points[i];
      u_hat[i] = get_max (u, start, end);
      start = end;
    }
}

ts_type get_min(ts_type * series, int start, int end)
{
    ts_type min = series[start];

    if (start >= end)
    {
      int j = 0;
      j++;
      printf("error start > end \n");
    }
    else
    {
      for (int i=start; i < end; i++) 
      {
        if(series[i] < min)
          min = series[i];
      }
    }
    return min;
}

void assign_min_segments(ts_type * l, ts_type * l_hat, short * points, int num_points)
{
    int start=0, end;
    for (int i=0; i< num_points; i++) 
    {
      end = points[i];
      l_hat[i] = get_min (l, start, end);
      start = end;
    }
}

/// This is from UCR SUITE code

  /// Data structure (circular array) for finding minimum and maximum for LB_Keogh envolop
  struct deque
  {   int *dq;
      int size,capacity;
      int f,r;
  };

  /// Initial the queue at the begining step of envelop calculation
  void init(struct deque *d, int capacity)
  {
      d->capacity = capacity;
      d->size = 0;
      d->dq = (int *) malloc(sizeof(int)*d->capacity);
      d->f = 0;
      d->r = d->capacity-1;
  }

  /// Destroy the queue
  void destroy(struct deque *d)
  {
      free(d->dq);
  }

  /// Insert to the queue at the back
  void push_back(struct deque *d, int v)
  {
      d->dq[d->r] = v;
      d->r--;
      if (d->r < 0)
          d->r = d->capacity-1;
      d->size++;
  }

  /// Delete the current (front) element from queue
  void pop_front(struct deque *d)
  {
      d->f--;
      if (d->f < 0)
          d->f = d->capacity-1;
      d->size--;
  }

  /// Delete the last element from queue
  void pop_back(struct deque *d)
  {
      d->r = (d->r+1)%d->capacity;
      d->size--;
  }

  /// Get the value at the current position of the circular queue
  int front(struct deque *d)
  {
      int aux = d->f - 1;

      if (aux < 0)
          aux = d->capacity-1;
      return d->dq[aux];
  }

  /// Get the value at the last position of the circular queueint back(struct deque *d)
  int back(struct deque *d)
  {
      int aux = (d->r+1)%d->capacity;
      return d->dq[aux];
  }

  /// Check whether or not the queue is empty
  int empty(struct deque *d)
  {
      return d->size == 0;
  }

  /// Finding the envelop of min and max value for LB_Keogh
  /// Implementation idea is intoruduced by Danial Lemire in his paper
  /// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
  void lower_upper_lemire(float *t, int len, int r, float *l, float *u)
  {
      struct deque du, dl;

      init(&du, 2*r+2);
      init(&dl, 2*r+2);

      push_back(&du, 0);
      push_back(&dl, 0);

      for (int i = 1; i < len; i++)
      {
          if (i > r)
          {
              u[i-r-1] = t[front(&du)];
              l[i-r-1] = t[front(&dl)];
          }
          if (t[i] > t[i-1])
          {
              pop_back(&du);
              while (!empty(&du) && t[i] > t[back(&du)])
                  pop_back(&du);
          }
          else
          {
              pop_back(&dl);
              while (!empty(&dl) && t[i] < t[back(&dl)])
                  pop_back(&dl);
          }
          push_back(&du, i);
          push_back(&dl, i);
          if (i == 2 * r + 1 + front(&du))
              pop_front(&du);
          else if (i == 2 * r + 1 + front(&dl))
              pop_front(&dl);
      }
      for (int i = len; i < len+r+1; i++)
      {
          u[i-r-1] = t[front(&du)];
          l[i-r-1] = t[front(&dl)];
          if (i-front(&du) >= 2 * r + 1)
              pop_front(&du);
          if (i-front(&dl) >= 2 * r + 1)
              pop_front(&dl);
      }
      destroy(&du);
      destroy(&dl);
  }