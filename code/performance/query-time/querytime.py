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



def summarize_results_to_csv(_nqueries, _source_dir, _output_file, _num_top=10):
    make_file(_output_file)
    _k = 0
    _k_count = 0
    for _subdir, _dirs, _files in os.walk(_source_dir):
        
        if re.search(f'nn', os.path.basename(_subdir)):
            _k = int(re.findall(r"(\d+)nn", _subdir)[0])
            print(f"Current k = {_k}")
            _k_count += 1
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
                    _len= _[5].replace('len','') # vector length
                    _querytime = _[6].replace('runtime','') # get query runtime in sec
                    _ndistcalc = _[7].replace('ndistcalc','') # get number of distance calculations
                    _dataaccess = _[8].replace('dataaccess','') # total checked vector
                    _dataaccess = _dataaccess.replace('.csv','') # remove extension '.csv'
                   
                    # save evaluation to csv file
                    append_evaluation(_algo, _l, _dlsize, _k, _num_top, _tqq, _querytime, _ndistcalc, _dataaccess, _output_file)
    return _k_count

def plot_results(_csv_file, _output_dir, _k_count):
    # set text font
    plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "cm"})
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

    _df = pandas.read_csv(_csv_file)
    _df = _df.drop(['TQ:Q'], axis = 1)
    _df = _df.drop(['num-top'], axis = 1)

    _df_new = _df.groupby(
        ['algorithm', 'total-files', 'data-gb-size', 'k']
    ).agg(
        querytime = ('querytime','mean'),
    ).reset_index()

    # BAR PLOT
    sns.set_style("whitegrid")
    ax = sns.barplot(data=_df_new, x="k", y="querytime")

    # Set these based on your column counts
    columncounts = [25 for i in range(_k_count)]

    # Maximum bar width is 1. Normalise counts to be in the interval 0-1. Need to supply a maximum possible count here as maxwidth
    def normaliseCounts(widths,maxwidth):
        widths = np.array(widths)/float(maxwidth)
        return widths

    widthbars = normaliseCounts(columncounts,100)

    # Loop over the bars, and adjust the width (and position, to keep the bar centred)
    for bar,newwidth in zip(ax.patches,widthbars):
        x = bar.get_x()
        width = bar.get_width()
        centre = x+width/2.

        bar.set_x(centre-newwidth/2.)
        bar.set_width(newwidth)
    
    plt.ylabel("Mean query time (sec)", fontsize = 11)
    plt.xlabel("k", fontsize = 11)
    plt.xticks(fontsize = 11)
    plt.yticks(fontsize = 11)
    # plt.legend(loc='upper left')
    plt.title("Kashif: mean query time (10 query columns of size [5 - 10])")
    plt.savefig(f"{_output_dir}/kashif_querytime.png")
    plt.close()


_source_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/kashif/results/100k-tables/"
_output_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/performance/query-time/img/"
_csv_file = "./csv/querytime.csv"
_nqueries = 10

_k_count = summarize_results_to_csv(_nqueries, _source_dir, _csv_file, _num_top=10)
_df_new = plot_results(_csv_file, _output_dir, _k_count)

