# evaluate query time for k = {10, 100, 1.000, 10.000}

import os
import csv
import re
import matplotlib.pyplot as plt
import pandas
import seaborn as sns
import numpy as np

def make_file(_output_file):
    _header = ['algorithm', 'total-files', 'data-gb-size', 'k','num-top', 'TQ:Q','querytime', 'ndistcalc', 'dataaccess']
    with open(_output_file, "w") as _f:
        writer = csv.writer(_f)
        writer.writerow(_header)
    _f.close()

def append_evaluation(_algo, _total_files, _dlsize, _k, _x, _tqq, _querytime, _ndistcalc, _dataaccess, _output_file):
    with open(_output_file, "a") as _f:
        writer = csv.writer(_f)
        writer.writerow([_algo, _total_files, _dlsize, _k, _x, _tqq, _querytime, _ndistcalc, _dataaccess])
    _f.close()
    return _output_file



def summarize_results_to_csv(_nqueries, _source_dir, _output_file, _algo_var):
    _k = 0
    _k_count = 0
    _k_values = set()
    for _subdir, _dirs, _files in os.walk(_source_dir):
        
        if re.search(f'_{_nqueries}q_min', os.path.basename(_subdir)):
            print(f"\nCurrent Directory: {_subdir}")
            _currentdir = os.path.basename(_subdir)

            # get algorithm name
            _ = _currentdir.split('_')
            _algo = _[0]

            for _file in _files:
                if re.search(f'_runtime', _file):
                    _ = _file.split('_')
                    _tqq = _[0] + _[1]
                    _qsize = _[2].replace('qsize','') # get number of candidate tables
                    _l = _[3].replace('l','') # get number of candidate tables
                    _dlsize = _[4].replace('dlsize','') # get data lake size in GB
                    _len = _[5].replace('len','') # vector length
                    _k = _[6].replace('k','') # vector length
                    _querytime = _[7].replace('runtime','') # get query runtime in sec
                    _ndistcalc = _[8].replace('ndistcalc','') # get number of distance calculations
                    _dataaccess = _[9].replace('dataaccess','') # total checked vector
                    _dataaccess = _dataaccess.replace('.csv','') # remove extension '.csv'

                    _k_values.add(_k)
                    print(f"_tqq {_tqq}, _qsize {_qsize}, _dlsize {_dlsize}, _k {_k}, _dataaccess {_dataaccess}\n")
                    # exit(1)
                    # save evaluation to csv file
                    append_evaluation(_algo_var, _l, _dlsize, _k, "-", _tqq, _querytime, _ndistcalc, _dataaccess, _output_file)
    return _k_values

def plot_results(csv_file, output_dir, k_values):
    k_count = len(k_values)
    # set text font
    plt.rcParams.update({
                     "text.usetex": True,
                     "font.family": "serif",
                     "font.serif": "Computer Modern",
                     "savefig.dpi": 130})


    sqa_df = pandas.read_csv(csv_file)
    sqa_df = sqa_df.drop(['TQ:Q'], axis = 1)
    sqa_df = sqa_df.drop(['num-top'], axis = 1)

    sqa_df_new = sqa_df.groupby(
        ['algorithm', 'total-files', 'data-gb-size', 'k']
    ).agg(
        querytime = ('querytime','mean'),
    ).reset_index()

    print(sqa_df_new)
    print(sqa_df_new)
    sqa_df_new[["k"]] = sqa_df_new[["k"]].astype(str) 

    # PLOT
    sns.color_palette("Set2")
    sns.set_style("whitegrid")
    ax = sns.lineplot(data = sqa_df_new, x="k", y="querytime", hue="algorithm", palette="Set2",style="algorithm", sort= False, markers=['o', 'v'])
    # ax = sns.barplot(data = sqa_df_new, x="k", y="querytime")

    plt.ylabel(r'$\mathrm{Mean \ query \ time \ (sec)}$', fontsize = 11)
    plt.xlabel(r'$\mathrm{k}$', fontsize = 11)
    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    # plt.legend(loc='upper left')
    plt.title(r'$\mathrm{5 \ query \ columns, \  query \ size \ \in \ [10 \ - \ 50]}$')
    plt.savefig(f"{output_dir}/kashif_querytime_sqa_vs_pqa.png")
    plt.close()



output_dir = "./img/"
# sequential query answers
sqa_dir = "/home/jaouhara/Documents/Projects/sqa-piqa/sqa/"
sqa_csv_file = "./csv/querytime_sqa_pqa.csv"

# parallel query answers
pqa_dir = "/home/jaouhara/Documents/Projects/sqa-piqa/piqa-1m-results/"
# pqa_csv_file = "./csv/querytime_pqa.csv"

nqueries = 5

make_file(sqa_csv_file)

k_values = summarize_results_to_csv(nqueries, sqa_dir, sqa_csv_file, "sqa")
k_values = summarize_results_to_csv(nqueries, pqa_dir, sqa_csv_file, "pqa")
print(sorted(k_values))
_df_new = plot_results(sqa_csv_file, output_dir, k_values)

