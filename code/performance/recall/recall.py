# evaluate query time for k = {10, 100, 1.000, 10.000}

import os
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

BRUTE_FORCE_ALGO_NAME = 'bfed'
NUM_TOP = 10
def make_file(output_file):
    header = ['algorithm', 'total-files', 'data-gb-size', 'k','num-top', 'TQ:Q','recall']
    with open(output_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
    f.close()
    return True

def append_evaluation(algo, total_files, data_gb_size, k, top, tqq, recall, output_file):
    with open(output_file, "a") as f:
        writer = csv.writer(f)
        writer.writerow([algo, total_files, data_gb_size, k, top, tqq, recall])
    f.close()
    return True

def get_recall_evaluation(nqueries, source_dir, ground_truth_dir, output_file):
    k = 0
    k_count = 0
    queries = list()
    for subdir, dirs, files in os.walk(source_dir):
        if re.search(f'nn', os.path.basename(subdir)):
            k = int(re.findall(r"(\d+)nn", subdir)[0])
            # print(f"Current k = {k}")
            if k == 10:
                k_count = 1
            else:
                k_count += 4
            
        if re.search(f'_{nqueries}q_min', os.path.basename(subdir)) and k_count == 1:
            # get algorithm name
            _currentdir = os.path.basename(subdir)
            _ = _currentdir.split('_')
            algo = _[0]

            for file in files:
                if re.search(f'_runtime', file):
                    _ = file.split('_')
                    qtable_id = _[0].replace('TQ', '')
                    qset_id =  _[1].replace('Q', '')
                    qsize = _[2].replace('qsize','') # get number of candidate tables
                    total_files = _[3].replace('l','') # get number of candidate tables
                    data_gb_size = _[4].replace('dlsize','') # get data lake size in GB
                    # _len= _[5].replace('len','') # vector length
                    # _querytime = _[6].replace('runtime','') # get query runtime in sec
                    # _ndistcalc = _[7].replace('ndistcalc','') # get number of distance calculations
                    # _dataaccess = _[8].replace('dataaccess','') # total checked vector
                    # _dataaccess = _dataaccess.replace('.csv','') # remove extension '.csv'
                   
                    query = (qtable_id, qset_id, qsize)
                    queries.append(query)
                    # print("compute recall...")
                    dir = str(os.path.basename(subdir)).replace(algo, BRUTE_FORCE_ALGO_NAME)
                    compute_query_recall(query, source_dir+f'/{k}nn/{os.path.basename(subdir)}/'+file, ground_truth_dir+f'/{k}nn/{dir}/', algo, total_files, data_gb_size, k, NUM_TOP)
                    
    return queries, k_count

def compute_query_recall(query, query_results_csv_file, ground_truth_dir, algo, total_files, data_gb_size, k, num_top):
    # print(ground_truth_dir)
    # exit(1)
    for _, _, files in os.walk(ground_truth_dir, topdown=False):
        for gt_file in files:
            if re.search(f'_runtime', gt_file):
                _ = gt_file.split('_')
                qtable_id = _[0].replace('TQ', '')
                qset_id =  _[1].replace('Q', '')
                qsize = _[2].replace('qsize','') # get number of candidate tables

                gt_query = (qtable_id, qset_id, qsize)

                if gt_query == query:
                    # print(f"query = {gt_query} \t gt query = {gt_query}")
                    recall = compute_recall(query_results_csv_file, ground_truth_dir+gt_file)
                    append_evaluation(algo, total_files, data_gb_size, k, num_top, f"TQ{qtable_id}Q{qset_id}", recall, output_file)
                    # print(f"query {query} recall = {recall}.")
                # else:
                    # print(f"query = {gt_query} !=\t gt query = {gt_query}")


def compute_recall(csv_file, gt_csv_file):
    gt_results = pd.read_csv(gt_csv_file, delimiter=",") 
    results = pd.read_csv(csv_file, delimiter=",") 
    total_gt_rows = len(gt_results)
    total_rows = len(results)

    gt_results.columns = [c.replace(' ', '') for c in gt_results.columns]
    results.columns = [c.replace(' ', '') for c in results.columns]

    gt_results = gt_results.sort_values('TS:S')
    results = results.sort_values('TS:S')

    match_count = 0
    for gt_r in range(total_gt_rows):
        for r in range(total_rows):
            if results.loc[r, "TS:S"] == gt_results.loc[gt_r, "TS:S"] and results.loc[r, "qindex"] == gt_results.loc[gt_r, "qindex"]:
                if  results.loc[r, "sindex"] == gt_results.loc[gt_r, "sindex"] and results.loc[r, "d"] == gt_results.loc[gt_r, "d"]:
                    match_count += 1
                    break

    return match_count/total_gt_rows


source_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/kashif/results/100k-tables/"
ground_truth_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/bf/results/100k-tables"
output_file = "./csv/recall_eval.csv"
nqueries = 10

make_file(output_file)
queries, k_count = get_recall_evaluation(nqueries, source_dir, ground_truth_dir, output_file)

# print(queries)
# print(len(queries))
# compute_recall(csv_file, gt_csv_file)