# evaluate query time for k = {10, 100, 1.000, 10.000}

import os
import csv
import re


def make_evaluation_file(_output_file):
    _header = ['algorithm', 'total-files', 'dlsize', 'k','num-top', 'TQ:Q','querytime', 'ndistcalc', 'dataaccess']
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



def save_evaluations(_nqueries, _source_dir, _output_file, _num_top=10):
    make_evaluation_file(_output_file)
    _k = 0
    for _subdir, _dirs, _files in os.walk(_source_dir):
        
        if re.search(f'nn', os.path.basename(_subdir)):
            _k = int(re.findall(r"(\d+)nn", _subdir)[0])
            print(f"Current k = {_k}")
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


_source_dir = "/home/jaouhara/Documents/Projects/iqa-demo/code/search-algorithms/bf/results/"
_output_file = "./csv/querytime.csv"
_nqueries = 10
save_evaluations(_nqueries, _source_dir, _output_file, _num_top=10)
