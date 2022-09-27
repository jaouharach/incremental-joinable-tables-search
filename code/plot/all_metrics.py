# evaluate query time for k = {10, 100, 1.000, 10.000}

import os
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


# import modin.pandas as pd

BRUTE_FORCE_ALGO_NAME = 'bf'
NUM_TOP = 10
def make_file(output_file):
    header = ['algorithm', 'total-files', 'data-gb-size', 'k', 'TQ:Q','recall', 'precision', 'querytime']
    with open(output_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
    f.close()
    return True

def append_evaluation(output_file, algo, total_files, data_gb_size, k, tqq, recall, precision, querytime,):
    with open(output_file, "a") as f:
        writer = csv.writer(f)
        writer.writerow([algo, total_files, data_gb_size, k, tqq, recall, precision, querytime])
    f.close()
    return True

def get_recall_evaluation(nqueries, source_dir, ground_truth_dir, output_file):
    k_values = set()
    for subdir, dirs, files in os.walk(source_dir):
        if re.search(f'_{nqueries}q_min', os.path.basename(subdir)):
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
                    len = _[5].replace('len','')
                    k = _[6].replace('k','')
                    runtime = _[7].replace('runtime','')

                    if k not in k_values:
                        k_values.add(k)
                    
                    query = (qtable_id, qset_id, qsize)
                    query_results_csv_file = source_dir+f'/{os.path.basename(subdir)}/'+file
                    compute_all_query_metrics(ground_truth_dir, query_results_csv_file, output_file, query, runtime, algo, total_files, data_gb_size, int(k))

    return k_values

def compute_all_query_metrics(ground_truth_dir, query_results_csv_file, output_file, query, runtime, algo, total_files, data_gb_size, k):
    for gt_file in os.listdir(ground_truth_dir):
        if gt_file.startswith(f"TQ{query[0]}_Q{query[1]}"):
            print(f"[k = {k}] query {query}")
            print(gt_file)
            recall, precision = 0, 0
            recall, precision = compute_recall_and_precision(query_results_csv_file, ground_truth_dir+"/"+gt_file, k)
            append_evaluation(output_file, algo, total_files, data_gb_size, k, f"TQ{query[0]}Q{query[1]}", recall, precision, runtime)
            print(f"recall = {recall}, precision = {precision}.")

def compute_recall_and_precision(csv_file, gt_csv_file, k):
    gt_results = pd.read_csv(gt_csv_file, delimiter=",") 
    results = pd.read_csv(csv_file, delimiter=",") 

    # remove space in column names
    gt_results.columns = [c.replace(' ', '') for c in gt_results.columns]
    results.columns = [c.replace(' ', '') for c in results.columns]

    total_gt_rows = len(gt_results)
    total_rows = len(results)

    gt_results = gt_results.sort_values('q_pos')
    results = results.sort_values('q_pos')

    gt_results = gt_results.sort_values('TS:S')
    results = results.sort_values('TS:S')

    gt_results = gt_results.sort_values('s_pos')
    results = results.sort_values('s_pos')

    match_count = 0

    # coumpte total rows for precision, since it is possible that some vectors have no exact matches we should evalues kashif precision based to vectors for which bf found an exact match
    gt_result_count = gt_results.groupby(['q_pos']).size().reset_index(name='count')
    gt_result_count.loc[gt_result_count["count"] > k, "count"] = k
    exact_result_total_rows = gt_result_count['count'].sum()

    for r in range(total_rows):

        found = ((gt_results["TS:S"] == results.loc[r, "TS:S"]) 
        & (gt_results["q_pos"] == results.loc[r, "q_pos"])
        & (gt_results["s_pos"] == results.loc[r, "s_pos"])).any()

        if found:
            match_count += 1

    recall = match_count/total_gt_rows
    precision = match_count/exact_result_total_rows

    return recall, precision

def plot_results(csv_file, output_dir):
    # set text font
    plt.style.use('classic')
    # plt.rc('axes', axisbelow=True)
    # plt.rcParams.update({
    #                  "text.usetex": True,
    #                  "font.family": "serif",
    #                  "font.serif": "Computer Modern",
    #                  "savefig.dpi": 130})

    df = pd.read_csv(csv_file)

    # querytime
    querytime_df = df.groupby(
        ['k']
    ).agg(
        querytime = ('querytime','mean'),
    ).reset_index()

    # recall
    recall_df = df.groupby(
        ['k']
    ).agg(
        recall = ('recall','mean'),
    ).reset_index()

    # precision
    precision_df = df.groupby(
        ['k']
    ).agg(
        precision = ('precision','mean'),
    ).reset_index()
    precision_df = precision_df.round(2)

    print(querytime_df.head(13))

    # PLOT
    _, ax1 = plt.subplots()
    ax1.yaxis.grid(True, color="lightgray")
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, linestyle='dashed', color="lightgray")
    sns.barplot(data = recall_df,x="k", y="recall", ax=ax1, errorbar=None, color="cornflowerblue")
    sns.lineplot(data = querytime_df["querytime"], sort = False, ax=ax2, color="red")
    
    # put precision as labels on top of bar.
    labels = precision_df["precision"]
    ax1.bar_label(ax1.containers[0], labels=labels, color="cornflowerblue")
    ax1.margins(y=0.1)

    ax1.set_ylabel(r'$\mathrm{Mean \ recall}$')
    ax1.xaxis.label.set_size(16)
    ax1.yaxis.label.set_size(16)
    ax2.set_ylabel(r'$\mathrm{Mean \ query \ time \ (sec)}$')
    ax2.yaxis.label.set_size(16)
    ax1.set_ylim([0, 1.5])
    # ax2.set_yscale('linear')
    # ax2.set_ylim([75, 200])
    # ax2.legend(loc=0)
    ax2.set_xlabel(r'$\mathrm{k}$')
    

    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    # plt.show()
    plt.savefig(f"{output_dir}/kashif_querytime_recall.png")
    plt.close()


source_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/kashif/results/100k-res/"
ground_truth_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/bf/results/100k-exact-results/bf_l100000_10q_min50_max100/"
output_img_dir = "./img/"

csv_file = "./csv/all_metrics_eval.csv"
querytime_csv_file = "../query-time/csv/kashif_querytime_100k_new.csv"
nqueries = 10
k_count = 10
# make_file(csv_file)
# k_values = get_recall_evaluation(nqueries, source_dir, ground_truth_dir, csv_file)
plot_results(csv_file, output_img_dir)
