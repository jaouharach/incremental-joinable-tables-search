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

def make_file(output_file):
    header = ['algorithm', 'k', 'tq:q', 'recall', 'precision']
    with open(output_file, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
    f.close()
    return True

def append_evaluation(algo, k, tqq, recall, precision, output_file):
    with open(output_file, "a") as f:
        writer = csv.writer(f)
        writer.writerow([algo, k, tqq, recall, precision])
    f.close()
    return True

def get_recall(tqq, gt_file, algo_file):
    gt_results = pd.read_csv(gt_file, delimiter=",") 
    results = pd.read_csv(algo_file, delimiter=",") 

    total_gt_rows = len(gt_results)
    total_rows = len(results)

    gt_results = gt_results.sort_values('ts:s')
    results = results.sort_values('ts:s')

    match_count = 0
    tp_fn = 0
    tp_fp = len(results[results["tq:q"] == tqq])
    for r in range(total_gt_rows):
        if gt_results.loc[r, "tq:q"] == tqq:
            tp_fn += 1
            found = ((results["ts:s"] == gt_results.loc[r, "ts:s"])
            & (results["tq:q"] == tqq)).any()
            if found:
                match_count += 1

    if tp_fn != 0:
        # if tp_fp > tp_fn:
        #     tp_fp = tp_fn
        return match_count/tp_fn, match_count/tp_fp
    return 1.0, 1.0 


# change bars width
def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)

# plot recall for diff k values 
def plot_results(csv_file, output_dir, k_count):
    # set text font
    plt.rcParams.update({
                     "text.usetex": True,
                     "font.family": "serif",
                     "font.serif": "Computer Modern",
                     "savefig.dpi": 130})
    plt.rc('axes', axisbelow=True)

    df = pd.read_csv(csv_file)
    
  
    recall_df = df.groupby(['algorithm', 'k']
    ).agg(
        recall = ('recall','mean'),
        precision = ('precision','mean')
    ).reset_index()
    print(recall_df)

    recall_df = recall_df.drop(["algorithm"], axis=1)
    
    recall_df = recall_df.melt('k', var_name='metric', value_name='avg.val')
    print(recall_df)
    
    # BAR PLOT
    fig, ax1 = plt.subplots()
    ax1.yaxis.grid(True, color="lightgray")
    # sns.barplot(data = recall_df,x="k", y="recall", ax=ax1, ci=None, color="cornflowerblue")
    sns.factorplot(data=recall_df,  ax=ax1, x='k', y='avg.val', hue='metric', kind='bar',
     palette=["mediumaquamarine", "cornflowerblue"])

    ax1.set_ylabel(r'$\mathrm{avg. \ value}$')
    ax1.set_ylim([0, 1])
    ax1.set_xlabel(r'$\mathrm{k}$')

    change_width(ax1, .25)
    # plt.savefig(f"{output_dir}/kashif_avg_column_recall_pexeso_T1%.png")
    plt.show()
    plt.close()



k_settings = ["nok", "1", "10", "100", "1000"]
fout = "recall_diff_settings.csv"
gt_file = "pexeso_results_tau=2%_T=20%.csv"

make_file(fout)

for k in k_settings:
    results_file = f"kashif_results_k={k}.csv"
    df = pd.read_csv(results_file, delimiter=",") 
    tqqs = df['tq:q'].unique()
    print(f"settings : k = {k}")
    for tqq in tqqs:
        recall, precision = get_recall(tqq, gt_file, results_file)
        append_evaluation("kashif", k, tqq, recall, precision, fout)
        # print(f"tqq : {tqq}\t\trecall = {recall}, precision = {precision}\n")

plot_results(fout, "./img/", len(k_settings))