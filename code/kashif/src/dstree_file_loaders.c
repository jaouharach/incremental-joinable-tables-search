//
//  dstree_file_loaders.c
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
#include <stdio.h>
#include "../include/dstree_file_loaders.h"
#include "../include/kashif_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "../include/dstree_query_engine.h"
#include <sys/types.h>
#include <unistd.h>
#include <unistd.h>

//CORRECT STATS FOR ASCII
enum response dstree_query_ascii_file(const char *ifilename, int q_num, 
                  const char delimiter, struct dstree_index *index,
		  float minimum_distance, ts_type epsilon, ts_type r_delta)
{
    FILE * ifile;
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    ifile = fopen (ifilename,"r");
    COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n", ifilename);
        return FAILURE;
    }
    
    char *ts_str = NULL; 
    size_t linecap = 0;
    ssize_t linelen;
    unsigned int q_loaded = 0;
    ts_type * ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    if (ts == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: Querying..\
                         Could not allocate memory for time series!\n");
        return FAILURE;
    }

    
    COUNT_PARTIAL_INPUT_TIME_START
    while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 && q_loaded < q_num)
    {
        COUNT_PARTIAL_INPUT_TIME_END
	COUNT_PARTIAL_SEQ_INPUT
        if(ts_str == NULL)
        {
          fprintf(stderr,"Error in dstree_file_loaders.c: Querying..\
                         Could not get the time series from file %s.\n", ifilename);
          return FAILURE;	
        }
        if (!ts_parse_str(ts_str, ts, index->settings->timeseries_size, &delimiter))
        { 
           fprintf(stderr, "Error in dstree_file_loaders.c:  Querying..Could not parse the time series.\n");
           return FAILURE;              
        }
        
	  printf("\n\n");

        q_loaded++;
        COUNT_PARTIAL_INPUT_TIME_START    
    }



    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the query filename %s", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END
      
    free(ts);
    free(ts_str);    
    return SUCCESS;    
}

enum response dstree_query_binary_file(const char *ifilename, int q_num, struct dstree_index *index,
						   float minimum_distance, ts_type epsilon,
				       ts_type r_delta )
{
    RESET_PARTIAL_COUNTERS()
    COUNT_PARTIAL_TIME_START

     
    FILE * ifile;
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    ifile = fopen (ifilename,"rb");
    COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    fseek(ifile, 0L, SEEK_SET);
    COUNT_PARTIAL_INPUT_TIME_END        
    unsigned int ts_length = index->settings->timeseries_size;
    file_position_type total_records = sz/ts_length * sizeof(ts_type);
    unsigned int offset = 0;

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    unsigned int q_loaded = 0;
    ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
    ts_type * query_ts_reordered = calloc(1,sizeof(ts_type) * ts_length);
    int * query_order = calloc(1,sizeof(int) * ts_length);
        if( query_order == NULL)
           return FAILURE;

    while (q_loaded < q_num)
    {
  
        RESET_QUERY_COUNTERS ()
	  
        COUNT_PARTIAL_SEQ_INPUT      
        COUNT_PARTIAL_INPUT_TIME_START
	fread(query_ts, sizeof(ts_type), ts_length, ifile);
        COUNT_PARTIAL_INPUT_TIME_END


	reorder_query(query_ts,query_ts_reordered,query_order,ts_length);        

      
	struct query_result result = exact_search(query_ts, query_ts_reordered, query_order, offset, index, minimum_distance, epsilon, r_delta);	

        q_loaded++;
 
        get_query_stats(index,1);	
        print_query_stats(index,q_loaded, 1, ifilename);
	RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START
	  
    }
    COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
      
     free(query_ts);
     free(query_ts_reordered);
     free(query_order);

    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the query filename %s", ifilename);
        return FAILURE;
    }

     return SUCCESS;
}

enum response dstree_knn_query_binary_file(const char *ifilename,
					   int q_num,
					   struct dstree_index *index,
					   float minimum_distance,
					   ts_type epsilon,
					   ts_type r_delta,
					   unsigned int k,
					   boolean track_bsf,
					   boolean track_pruning,
					   boolean all_mindists,
					   boolean max_policy,
					   unsigned int nprobes,
					   unsigned char incremental,
             float warping)
{

  struct bsf_snapshot ** bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
  if(track_bsf)
  {
    max_bsf_snapshots = 10000;
    cur_bsf_snapshot = 0;

    bsf_snapshots = calloc(k, sizeof(struct bsf_snapshot*));      
    for (unsigned int i = 0; i < k; ++i)
    {
      bsf_snapshots[i] = calloc(max_bsf_snapshots, sizeof(struct bsf_snapshot));
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
      {
	bsf_snapshots[i][j].distance = FLT_MAX;
	bsf_snapshots[i][j].time = FLT_MAX;
	bsf_snapshots[i][j].series = NULL;	
	bsf_snapshots[i][j].checked_nodes = -1;
	bsf_snapshots[i][j].label = 0;		
      }      
    }
  }

  
    RESET_PARTIAL_COUNTERS()

    COUNT_PARTIAL_TIME_START

    FILE * ifile;
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    ifile = fopen (ifilename,"rb");
    COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    fseek(ifile, 0L, SEEK_SET);
    COUNT_PARTIAL_INPUT_TIME_END        
    unsigned int ts_length = index->settings->timeseries_size;
    file_position_type total_records = sz/ts_length * sizeof(ts_type);
    unsigned int offset = 0;

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    unsigned int q_loaded = 0;
    ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
    ts_type * query_ts_reordered = calloc(1,sizeof(ts_type) * ts_length);
    int * query_order = calloc(1,sizeof(int) * ts_length);
        if( query_order == NULL)
           return FAILURE;

 /*
     const char *filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
     filename = strcpy(filename, index->settings->root_directory);
     filename = strcat(filename, "raw_series.csv\0");

   printf ("series_file = %s\n", filename);
   printf ("dataset_file = %s\n", index->settings->dataset);

    FILE *series_file = fopen(filename, "a");
    FILE *dataset_file = fopen(index->settings->dataset, "rb");
*/
    FILE *series_file = NULL;
    FILE *dataset_file = NULL;
    while (q_loaded < q_num)
    {
  
        RESET_QUERY_COUNTERS ()
	  
        COUNT_PARTIAL_SEQ_INPUT      
        COUNT_PARTIAL_INPUT_TIME_START
	fread(query_ts, sizeof(ts_type), ts_length, ifile);
        COUNT_PARTIAL_INPUT_TIME_END


	reorder_query(query_ts,query_ts_reordered,query_order,ts_length);        

        q_loaded++;
	
        if (track_bsf)
	{
	  cur_bsf_snapshot = 0;
	  if (incremental)
	    {
	      exact_de_incr_progressive_knn_search (query_ts, query_ts_reordered,query_order, offset,
						    index, minimum_distance,epsilon, r_delta,
						    k, q_loaded, ifilename,
						    bsf_snapshots, &cur_bsf_snapshot, warping, dataset_file, series_file);
	    }
	  else
	    {
	      exact_de_progressive_knn_search (query_ts, query_ts_reordered,query_order, offset,
					       index, minimum_distance,epsilon, r_delta,
					       k, q_loaded, ifilename,
					       bsf_snapshots, &cur_bsf_snapshot);
	    }
	  for (unsigned int i = 0; i < k; ++i)
	  {
	    for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
	    {
	      bsf_snapshots[i][j].distance = FLT_MAX;
	      bsf_snapshots[i][j].time = FLT_MAX;
	      bsf_snapshots[i][j].series = NULL;
	      bsf_snapshots[i][j].checked_nodes = -1;
	      bsf_snapshots[i][j].label = 0;	      	      	      	      
	    }      
	  }
	}
	else 
	{
	  exact_de_knn_search(query_ts, query_ts_reordered, query_order, offset,
			   index, minimum_distance, epsilon, r_delta, k,
			   q_loaded, ifilename);		  
	}
  
	RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START
	  
    }
    COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
      
     free(query_ts);
     free(query_ts_reordered);
     free(query_order);

    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the query filename %s", ifilename);
        return FAILURE;
    }

    if (track_bsf)
    {
      for (unsigned int i = 0; i < k; ++i)
	{
	  for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
	    {
	      if(bsf_snapshots[i][j].series != NULL)
		free(bsf_snapshots[i][j].series);
	    }        
	  free(bsf_snapshots[i]);
	}
	free(bsf_snapshots);
    }
      //fclose(series_file);
      //fclose(dataset_file);

    return SUCCESS;
}
enum response dstree_knn_query_gt_binary_file(const char *ifilename,
					      int q_num,
					      struct dstree_index *index,
					      float minimum_distance,
					      ts_type epsilon,
					      ts_type r_delta,
					      unsigned int k,
					      boolean track_bsf,
					      boolean track_pruning,
					      boolean all_mindists,
					      boolean max_policy,
					      unsigned int nprobes,
					      unsigned char incremental)					   
{

  struct bsf_snapshot ** bsf_snapshots = NULL;
  unsigned int max_bsf_snapshots;
  unsigned int cur_bsf_snapshot;
  if(track_bsf)
  {
    max_bsf_snapshots = 10000;
    cur_bsf_snapshot = 0;

    bsf_snapshots = calloc(k, sizeof(struct bsf_snapshot*));      
    for (unsigned int i = 0; i < k; ++i)
    {
      bsf_snapshots[i] = calloc(max_bsf_snapshots, sizeof(struct bsf_snapshot));
      for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
      {
	bsf_snapshots[i][j].distance = FLT_MAX;
	bsf_snapshots[i][j].time = FLT_MAX;
	bsf_snapshots[i][j].series = NULL;
	bsf_snapshots[i][j].checked_nodes = -1;		
      }      
    }
  }
     /*
     const char *filename = malloc(sizeof(char) * (strlen(index->settings->root_directory) + 18));
     filename = strcpy(filename, index->settings->root_directory);
     filename = strcat(filename, "raw_series.csv\0");

    printf ("series_file = %s\n", filename);
    printf ("dataset_file = %s\n", index->settings->dataset);

    FILE *series_file = fopen(filename, "w");
    FILE *dataset_file = fopen(index->settings->dataset, "rb");
    */
    FILE * series_file = NULL;
    FILE * dataset_file = NULL;
  
    RESET_PARTIAL_COUNTERS()

    COUNT_PARTIAL_TIME_START

    FILE * ifile;
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    ifile = fopen (ifilename,"rb");
    COUNT_PARTIAL_INPUT_TIME_END    
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    fseek(ifile, 0L, SEEK_SET);
    COUNT_PARTIAL_INPUT_TIME_END        
    unsigned int ts_length = index->settings->timeseries_size;
    file_position_type total_records = sz/ts_length * sizeof(ts_type);
    unsigned int offset = 0;

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    unsigned int q_loaded = 0;
    ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
    ts_type * query_ts_reordered = calloc(1,sizeof(ts_type) * ts_length);
    int * query_order = calloc(1,sizeof(int) * ts_length);
        if( query_order == NULL)
           return FAILURE;

    while (q_loaded < q_num)
    {
  
        RESET_QUERY_COUNTERS ()
	  
        COUNT_PARTIAL_SEQ_INPUT      
        COUNT_PARTIAL_INPUT_TIME_START
	fread(query_ts, sizeof(ts_type), ts_length, ifile);
        COUNT_PARTIAL_INPUT_TIME_END


	reorder_query(query_ts,query_ts_reordered,query_order,ts_length);        

        q_loaded++;
	
        if (track_bsf)
	{
	  cur_bsf_snapshot = 0;
	  if (incremental)
	    {
	      exact_de_incr_progressive_knn_search (query_ts, query_ts_reordered,query_order, offset,
						    index, minimum_distance,epsilon, r_delta,
						    k, q_loaded, ifilename,
						    bsf_snapshots, &cur_bsf_snapshot,0, dataset_file, series_file);
	    }
	  else
	    {
	      exact_de_progressive_knn_search (query_ts, query_ts_reordered,query_order, offset,
					       index, minimum_distance,epsilon, r_delta,
					       k, q_loaded, ifilename,
					       bsf_snapshots, &cur_bsf_snapshot);
	    }
	  for (unsigned int i = 0; i < k; ++i)
	  {
	    for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
	    {
	      bsf_snapshots[i][j].distance = FLT_MAX;
	      bsf_snapshots[i][j].time = FLT_MAX;
	      bsf_snapshots[i][j].series = NULL;
	      bsf_snapshots[i][j].checked_nodes = -1;	      	      	      
	    }      
	  }
	}
	else 
	{
	  exact_de_knn_search(query_ts, query_ts_reordered, query_order, offset,
			   index, minimum_distance, epsilon, r_delta, k,
			   q_loaded, ifilename);		  
	}
  
	RESET_PARTIAL_COUNTERS()
	COUNT_PARTIAL_TIME_START
	  
    }
    COUNT_PARTIAL_TIME_END
    RESET_PARTIAL_COUNTERS()
      
     free(query_ts);
     free(query_ts_reordered);
     free(query_order);

    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the query filename %s", ifilename);
        return FAILURE;
    }

    if (track_bsf)
    {
      for (unsigned int i = 0; i < k; ++i)
	{
	  for (unsigned int j = 0; j < max_bsf_snapshots; ++j)
	    {
	      if(bsf_snapshots[i][j].series != NULL)
		free(bsf_snapshots[i][j].series);
	    }        
	  free(bsf_snapshots[i]);
	}
	free(bsf_snapshots);
    }
//fclose (dataset_file);
//fclose (series_file);

    return SUCCESS;
}

enum response dstree_tlb_binary_file(const char *ifilename, int q_num, struct dstree_index *index,
						   float minimum_distance)
{

    FILE * ifile;

    ifile = fopen (ifilename,"rb");

    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    fseek(ifile, 0L, SEEK_SET);
    unsigned int ts_length = index->settings->timeseries_size;
    file_position_type total_records = sz/ts_length * sizeof(ts_type);
    unsigned int offset = 0;

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    unsigned int q_loaded = 0;
    ts_type * query_ts = calloc(1,sizeof(ts_type) * ts_length);
      
    while (q_loaded < q_num)
    {
        total_tlb = 0;
	total_ts_count = 0;
	leaf_nodes_count = 0;

	fread(query_ts, sizeof(ts_type), ts_length, ifile);
      
	dstree_calc_tlb(query_ts, index, index->first_node);
	
        q_loaded++;
        print_tlb_stats(index,q_loaded, ifilename);

    }
      
     free(query_ts);

     if(fclose(ifile))
       {   
	 fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the query filename %s", ifilename);
	 return FAILURE;
       }
     
     return SUCCESS;
}

enum response dstree_index_ascii_file(const char *ifilename, file_position_type ts_num, 
                           const char delimiter,struct  dstree_index *index)
{
    double parse_time = 0;
    int ts_length = index->settings->timeseries_size;
    ts_type * ts = NULL;
    
    ts = NULL;
    ts = malloc (sizeof(ts_type) * ts_length); 

    if(ts == NULL)
    {
          fprintf(stderr,"Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
          return FAILURE;	
    }
    
    FILE * ifile;
    COUNT_INPUT_TIME_START
    ifile = fopen (ifilename,"r");
    COUNT_INPUT_TIME_END
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    
	
    char *ts_str = NULL; //= malloc(sizeof(char) * 2000);
    size_t linecap = 0;
    ssize_t linelen;
    
    file_position_type ts_loaded = 0;    
    

    int percentage = 100;
    if (percentage > 100)
    {
       percentage = (int) (ts_num / (file_position_type) 100);
    }


    COUNT_INPUT_TIME_START
    while ((linelen = getline(&ts_str, &linecap, ifile)) > 0 && ts_loaded<ts_num)
    {
        COUNT_INPUT_TIME_END
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
        printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
        COUNT_PARSE_TIME_START
        ts_parse_str(ts_str, ts, index->settings->timeseries_size, &delimiter);
        COUNT_PARSE_TIME_END


	if (!dstree_index_insert(index,ts))
	{
           fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
           return FAILURE;              
        }
    
        ts_loaded++;

	if(ts_loaded % percentage == 0)
            {
                float distance = 0;
            }

        COUNT_INPUT_TIME_START  
     }
    
    free(ts_str);
    free(ts);
    COUNT_INPUT_TIME_START
	fclose(ifile);
    COUNT_INPUT_TIME_END

      return SUCCESS;
}

enum response dstree_index_binary_file(const char *ifilename, file_position_type ts_num, struct dstree_index *index)
{
    double parse_time = 0;
    int ts_length = index->settings->timeseries_size;
    ts_type * ts = NULL;
    
    ts = NULL;
    ts = malloc (sizeof(ts_type) * ts_length); 

    if(ts == NULL)
    {
          fprintf(stderr,"Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
          return FAILURE;	
    }
    
    FILE * ifile;
      
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    ifile = fopen (ifilename,"rb");
    COUNT_PARTIAL_INPUT_TIME_END
    if (ifile == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    COUNT_PARTIAL_INPUT_TIME_END    
    file_position_type total_records = sz/index->settings->timeseries_size * sizeof(ts_type);

    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_SET);
    COUNT_PARTIAL_INPUT_TIME_END        
    if (total_records < ts_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        return FAILURE;
    }
    	
    file_position_type ts_loaded = 0;    
    
    int percentage = 100;
    if (percentage > 100)
    {
       percentage = (int) (ts_num / (file_position_type) 100);
    }

    while (ts_loaded<ts_num)
    {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
        printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
        COUNT_PARTIAL_SEQ_INPUT
	COUNT_PARTIAL_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_PARTIAL_INPUT_TIME_END

	if (!dstree_index_insert(index,ts))
	  {
           fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
           return FAILURE;              
        }
            COUNT_PARTIAL_TIME_END

	    index->stats->idx_building_total_time  += partial_time;	
	    index->stats->idx_building_input_time  += partial_input_time;
	    index->stats->idx_building_output_time += partial_output_time;

	    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_building_seq_output_count += partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

        ts_loaded++;

     
        if(ts_loaded % percentage == 0)
        {
           float distance = 0;
	   //PRINT_STATS(distance)
        }

     }

    free(ts);
    COUNT_PARTIAL_INPUT_TIME_START    
    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the filename %s", ifilename);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END    

            COUNT_PARTIAL_TIME_END
	    index->stats->idx_building_total_time  += partial_time;	
	    index->stats->idx_building_input_time  += partial_input_time;
	    index->stats->idx_building_output_time += partial_output_time;
	    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_building_seq_output_count += partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

      return SUCCESS;      
}

enum response dstree_index_classify_binary_file(const char *ifilename, file_position_type ts_num, struct dstree_index *index)
{
    double parse_time = 0;
    int ts_length = index->settings->timeseries_size;
    ts_type * ts = NULL;

    int filename_size = strlen(ifilename)+4;     
    const char *gt_filename = malloc(sizeof(char) * filename_size);
    gt_filename = strcpy(gt_filename, "\0");
    gt_filename = strcat(gt_filename, ifilename);
    gt_filename = strcat(gt_filename, ".gt\0");
    
    ts = NULL;
    ts = malloc (sizeof(ts_type) * ts_length); 

    //unsigned char ts_gt = 0;
    label_type ts_gt = 0;

    if(ts == NULL)
    {
          fprintf(stderr,"Error in dstree_file_loaders.c: Could not allocate memory for ts.\n");
          return FAILURE;	
    }
    
    FILE * ifile;
    FILE * gt_file;
    
      
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    ifile = fopen (ifilename,"rb");
    gt_file = fopen (gt_filename,"rb");

    
    COUNT_PARTIAL_INPUT_TIME_END
    if (ifile == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n", ifilename);
        return FAILURE;
    }
    if (gt_file == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n", gt_filename);
        return FAILURE;
    }
    free (gt_filename);
    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START    
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    COUNT_PARTIAL_INPUT_TIME_END    
    file_position_type total_records = sz/index->settings->timeseries_size * sizeof(ts_type);

    COUNT_PARTIAL_RAND_INPUT
    COUNT_PARTIAL_INPUT_TIME_START
    fseek(ifile, 0L, SEEK_SET);
    COUNT_PARTIAL_INPUT_TIME_END        
    if (total_records < ts_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        return FAILURE;
    }
    	
    //file_position_type ts_loaded = 0;
    unsigned int ts_loaded = 0;    
    
    int percentage = 100;
    if (percentage > 100)
    {
       percentage = (int) (ts_num / (file_position_type) 100);
    }

    while (ts_loaded<ts_num)
    {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
        printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
        COUNT_PARTIAL_SEQ_INPUT
	COUNT_PARTIAL_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        //fread(&ts_gt, sizeof(unsigned char), 1, gt_file);	
	fread(&ts_gt, sizeof(label_type), 1, gt_file);	
        COUNT_PARTIAL_INPUT_TIME_END

	  if (!index->settings->track_file_pos)
	    {
	      if (!dstree_index_classify_insert(index,ts, ts_gt,-1))
		{
		  fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
		  return FAILURE;              
		}	      
	    }
	  else
	    {
	      if (!dstree_index_classify_insert(index,ts, ts_gt, (ts_loaded+1) ))
		{
		  fprintf(stderr, "Error in dstree_file_loaders.c:  Could not \
                           add the time series to the index.\n");
		  return FAILURE;              
       	}	      	      
	    }	    

            COUNT_PARTIAL_TIME_END

	    index->stats->idx_building_total_time  += partial_time;	
	    index->stats->idx_building_input_time  += partial_input_time;
	    index->stats->idx_building_output_time += partial_output_time;

	    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_building_seq_output_count += partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

        ts_loaded++;

     
        if(ts_loaded % percentage == 0)
        {
           float distance = 0;
	   //PRINT_STATS(distance)
        }

     }

    free(ts);
    COUNT_PARTIAL_INPUT_TIME_START    
    if(fclose(ifile))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the filename %s", ifilename);
        return FAILURE;
    }

    if(fclose(gt_file))
    {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the filename %s", gt_file);
        return FAILURE;
    }
    COUNT_PARTIAL_INPUT_TIME_END    

            COUNT_PARTIAL_TIME_END
	    index->stats->idx_building_total_time  += partial_time;	
	    index->stats->idx_building_input_time  += partial_input_time;
	    index->stats->idx_building_output_time += partial_output_time;
	    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_building_seq_output_count += partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_building_rand_output_count += partial_rand_output_count;	    
 
	    RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START

      return SUCCESS;      
}


enum response reorder_query(ts_type * query_ts, ts_type * query_ts_reordered, int * query_order, int ts_length)
{
  
        q_index *q_tmp = malloc(sizeof(q_index) * ts_length);
        int i;
	
        if( q_tmp == NULL )
	  return FAILURE;

	for( i = 0 ; i < ts_length ; i++ )
        {
          q_tmp[i].value = query_ts[i];
          q_tmp[i].index = i;
        }
	
        qsort(q_tmp, ts_length, sizeof(q_index),znorm_comp);

        for( i=0; i<ts_length; i++)
        {
          
	  query_ts_reordered[i] = q_tmp[i].value;
          query_order[i] = q_tmp[i].index;
        }
        free(q_tmp);

	return SUCCESS;
}


int znorm_comp(const void *a, const void* b)
{
    q_index* x = (q_index*)a;
    q_index* y = (q_index*)b;


    if (fabsf(y->value) > fabsf(x->value) )
       return 1;
    else if (fabsf(y->value) == fabsf(x->value))
      return 0;
    else
      return -1;

}

/* start kashif changes */
enum response dstree_index_multiple_binary_files(const char * bin_files_directory, unsigned int total_data_files, struct dstree_index * index)
{
  int vector_length = index->settings->timeseries_size;
  int opened_files = 0; 

  // allocate memory for vector
	struct vector v; 
  v.values = (ts_type *) malloc(sizeof(ts_type) * vector_length);
  ts_type val; 
  unsigned int nvec = 0u;

	// open source dir
  struct dirent *dfile;
  DIR *dir = opendir(bin_files_directory);
  if (!dir)
  {
    printf("Error in dstree_file_loader.c: Unable to open directory stream!");
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0)
  {
    if(is_binaryfile(dfile->d_name))
    {
      total_data_files --;
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
        fprintf(stderr, "Error in dstree_file_loaders.c:  Vector length passed in argumentes (--timeseries-size %d) does not match vector length in file (%d) %s.\n", vector_length_in_filename, vector_length, bin_file_path);
        return FAILURE;
      }

      COUNT_PARTIAL_RAND_INPUT
      COUNT_PARTIAL_INPUT_TIME_START
      FILE * bin_file = fopen (bin_file_path, "rb");
      COUNT_PARTIAL_INPUT_TIME_END

      if (bin_file == NULL)
      {
        fprintf(stderr, "Error in dstree_file_loaders.c: File %s not found!\n", bin_file_path);
        return FAILURE;
      }

      /* Start processing file: read every vector in binary file*/
      int i = 0 , j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
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

            set_id += 1;
            
        }
        else if(i <= (unsigned int)nvec * vector_length)
        {
            // end of vector but still in current set
            if(j > (vector_length - 1)){
                j = 0; 
                /*Index v in dstree */
                if (!dstree_index_insert_vector(index,v.values, v.table_id, v.set_id))
                {
                fprintf(stderr, "Error in dstree_file_loaders.c:  Could not add the time series to the index.\n");
                return FAILURE;              
                }

                index->stats->idx_building_total_time  += partial_time;	
                index->stats->idx_building_input_time  += partial_input_time;
                index->stats->idx_building_output_time += partial_output_time;

                index->stats->idx_building_seq_input_count  += partial_seq_input_count;
                index->stats->idx_building_seq_output_count += partial_seq_output_count;
                index->stats->idx_building_rand_input_count  += partial_rand_input_count;
                index->stats->idx_building_rand_output_count += partial_rand_output_count;
            }
            COUNT_PARTIAL_SEQ_INPUT
	          COUNT_PARTIAL_INPUT_TIME_START 
            fread((void*)(&val), sizeof(val), 1, bin_file);
            COUNT_PARTIAL_INPUT_TIME_END 
            total_bytes--;
            v.values[j] = val;

            // end of last vector in current  set
            if(i == (unsigned int)nvec * vector_length)
            {   
                /*Index v in dstree */
                if (!dstree_index_insert_vector(index,v.values, v.table_id, v.set_id))
                {
                  fprintf(stderr, "Error in dstree_file_loaders.c:  Could not add the time series to the index.\n");
                  return FAILURE;              
                }
                index->stats->idx_building_total_time  += partial_time;	
                index->stats->idx_building_input_time  += partial_input_time;
                index->stats->idx_building_output_time += partial_output_time;

                index->stats->idx_building_seq_input_count  += partial_seq_input_count;
                index->stats->idx_building_seq_output_count += partial_seq_output_count;
                index->stats->idx_building_rand_input_count  += partial_rand_input_count;
                index->stats->idx_building_rand_output_count += partial_rand_output_count;

                i = 0; j = 0;
                nvec = 0u;
                continue;
            }
            i++;
            j++;
        }
      }
      
      COUNT_PARTIAL_INPUT_TIME_START
      if(fclose(bin_file))
      {   
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not close the filename %s", bin_file_path);
        return FAILURE;
      }
      COUNT_PARTIAL_INPUT_TIME_END
      COUNT_PARTIAL_TIME_END
      /* end read processing file */
      index->stats->idx_building_total_time  += partial_time;	
	    index->stats->idx_building_input_time  += partial_input_time;
	    index->stats->idx_building_output_time += partial_output_time;
	    index->stats->idx_building_seq_input_count  += partial_seq_input_count;
	    index->stats->idx_building_seq_output_count += partial_seq_output_count;
	    index->stats->idx_building_rand_input_count  += partial_rand_input_count;
	    index->stats->idx_building_rand_output_count += partial_rand_output_count;
      RESET_PARTIAL_COUNTERS()
	    COUNT_PARTIAL_TIME_START
    }
  }
  closedir(dir);
  free(v.values);
  if(opened_files == 0)
  {
    fprintf(stderr, "Error in dstree_file_loaders.c:  Could not find any binary file in directory %s.\n", bin_files_directory);
    return FAILURE;
  }   
}
/* end kashif changes */