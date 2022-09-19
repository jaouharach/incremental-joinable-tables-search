# evaluate query time for k = {10, 100, 1.000, 10.000}

import os
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

#USE ONLY ONE OF THESE:
# os.environ["MODIN_ENGINE"] = "ray"  # Modin will use Ray
# os.environ["MODIN_ENGINE"] = "dask"  # Modin will use Dask

# import modin.pandas as pd

BRUTE_FORCE_ALGO_NAME = 'bf'
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
            k_count += 1
            
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

                    query = (qtable_id, qset_id, qsize)
                    queries.append(query)
                    dir = str(os.path.basename(subdir)).replace(algo, BRUTE_FORCE_ALGO_NAME)
                    print(f"----------- k = {k} ------------\n")
                    compute_query_recall(query, source_dir+f'/{k}nn/{os.path.basename(subdir)}/'+file, ground_truth_dir+f'/{dir}/', algo, total_files, data_gb_size, k, NUM_TOP, output_file)
                    
    return queries, k_count

def compute_query_recall(query, query_results_csv_file, ground_truth_dir, algo, total_files, data_gb_size, k, num_top, output_file):
    for _, _, files in os.walk(ground_truth_dir, topdown=False):
        for gt_file in files:
            if re.search(f'_runtime', gt_file):
                _ = gt_file.split('_')
                qtable_id = _[0].replace('TQ', '')
                qset_id =  _[1].replace('Q', '')
                qsize = _[2].replace('qsize','') # get number of candidate tables
            
                gt_query = (qtable_id, qset_id, qsize)
                if gt_query == query:
                    print(f"query = {gt_query} \t gt query = {gt_query}")
                    recall = compute_recall(query_results_csv_file, ground_truth_dir+gt_file)
                    append_evaluation(algo, total_files, data_gb_size, k, num_top, f"TQ{qtable_id}Q{qset_id}", recall, output_file)
                    print(f"query {query} recall = {recall}.")


def compute_recall(csv_file, gt_csv_file):
    gt_results = pd.read_csv(gt_csv_file, delimiter=",") 
    results = pd.read_csv(csv_file, delimiter=",") 

    # remove space in column names
    gt_results.columns = [c.replace(' ', '') for c in gt_results.columns]
    results.columns = [c.replace(' ', '') for c in results.columns]

    # correct column name (for result of older code)
    gt_results.rename(columns = {'qindex':'q_pos', 'sindex':'s_pos'}, inplace = True)
    results.rename(columns = {'qindex':'q_pos', 'sindex':'s_pos'}, inplace = True)

    total_gt_rows = len(gt_results)
    total_rows = len(results)

    gt_results = gt_results.sort_values('q_pos')
    results = results.sort_values('q_pos')

    gt_results = gt_results.sort_values('TS:S')
    results = results.sort_values('TS:S')

    gt_results = gt_results.sort_values('s_pos')
    results = results.sort_values('s_pos')

    match_count = 0
    for r in range(total_rows):
        # if gt_results.loc[r, "q_pos"] == results.loc[r, "q_pos"] and gt_results.loc[r, "s_pos"] == results.loc[r, "s_pos"] and gt_results.loc[r, "TS:S"] == results.loc[r, "TS:S"]:
        #     match_count += 1

        found = ((gt_results["TS:S"] == results.loc[r, "TS:S"]) 
        & (gt_results["q_pos"] == results.loc[r, "q_pos"])
        & (gt_results["s_pos"] == results.loc[r, "s_pos"])).any()

        if found:
            match_count += 1

    return match_count/total_gt_rows

def plot_results(csv_file, querytime_csv_file, output_dir, k_count):
    # set text font
    plt.rcParams.update({
                     "text.usetex": True,
                     "font.family": "serif",
                     "font.serif": "Computer Modern",
                     "savefig.dpi": 130})

    data1 = pd.read_csv(csv_file)
    data2 = pd.read_csv(querytime_csv_file)
    
    # using merge function by setting how='inner'
    df = pd.merge(data1, data2, 
                    on=['TQ:Q','k', 'total-files','data-gb-size', 'num-top'], 
                    how='inner')

    df = df.drop(['dataaccess'], axis = 1)
    df = df.drop(['ndistcalc'], axis = 1)
    df = df.drop(['num-top'], axis = 1)

    plt.rc('axes', axisbelow=True)

    print(df)
  
    recall_df = df.groupby(
        ['total-files', 'data-gb-size', 'k']
    ).agg(
        recall = ('recall','mean'),
    ).reset_index()

    querytime_df = df.groupby(
        ['total-files', 'data-gb-size', 'k']
    ).agg(
        querytime = ('querytime','mean'),
    ).reset_index()

    # BAR PLOT
    fig, ax1 = plt.subplots()
    ax1.yaxis.grid(True, color="lightgray")
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, linestyle='dashed', color="lightgray")
    sns.barplot(data = df,x="k", y="recall", ax=ax1, ci=None,)
    sns.lineplot(data = querytime_df['querytime'], sort = False, ax=ax2, color="blue")
    

    ax1.set_ylabel(r'$\mathrm{Mean \ recall}$')
    ax2.set_ylabel(r'$\mathrm{Mean \ query \ time \ (sec)}$')
    ax1.set_ylim([0, 1.5])
    ax2.set_yscale('linear')
    # ax2.legend(loc=0)
    ax2.set_xlabel(r'$\mathrm{k}$')

    # Set these based on your column counts
    columncounts = [25 for i in range(k_count)]

    # Maximum bar width is 1. Normalise counts to be in the interval 0-1. Need to supply a maximum possible count here as maxwidth
    def normaliseCounts(widths,maxwidth):
        widths = np.array(widths)/float(maxwidth)
        return widths

    widthbars = normaliseCounts(columncounts,100)

    # Loop over the bars, and adjust the width (and position, to keep the bar centred)
    for bar,newwidth in zip(ax2.patches,widthbars):
        x = bar.get_x()
        width = bar.get_width()
        centre = x+width/2.

        bar.set_x(centre-newwidth/2.)
        bar.set_width(newwidth)

    plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    # plt.legend(loc='upper left')
    # plt.title("Kashif: mean recall (10 query columns of size [50 - 100])")
    plt.savefig(f"{output_dir}/kashif_recall_querytime_on_vid.png")
    plt.close()


source_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/kashif/results/100k-results/"
ground_truth_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/bf/results/100k-results/"
output_img_dir = "./img/"

csv_file = "./csv/recall.csv"
querytime_csv_file = "../query-time/csv/kashif_querytime_100k_new.csv"
nqueries = 10
k_count = 10
make_file(csv_file)
queries, k_count = get_recall_evaluation(nqueries, source_dir, ground_truth_dir, csv_file)

plot_results(csv_file,querytime_csv_file, output_img_dir, k_count)