from pickle import NONE
from flask import Flask, render_template, request, redirect
from matplotlib.pyplot import table
import pandas as pd
import subprocess, re
import time, json, os


import sys
sys.path.append('./utils')
from query_to_bin import *

ALLOWED_EXTENSIONS = {'csv'}

RESULTS_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/query_results/"
UPLOAD_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/uploads/"
TMP_FOLDER = UPLOAD_FOLDER + "tmp/"
BIN_FOLDER = UPLOAD_FOLDER + "bins/"

EMBEDDING_DIM = 50
KASHIF_BIN = "../search-algorithms/kashif/bin/dstree"

KASHIF_IDX = "../../100-idx/" # storing 100 tables
RAW_DATA_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/raw-tables/"
CLEAN_DATA_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/metadata/"

EMBEDDING_MODEL = 'glove' # 'fasttext' or 'glove'
PATH_TO_MODEL = {"ar": "", "en": "/home/jaouhara/Documents/Projects/embedding_models/glove/glove.6B.50d.txt","fr": ""}

# EMBEDDING_MODEL = 'fasttext' 
# PATH_TO_MODEL = {"ar": "/home/jaouhara/Documents/Projects/embedding_models/fasttext/cc.ar.300.bin", "en": "/home/jaouhara/Documents/Projects/embedding_models/fasttext/cc.en.300.bin","fr": "/home/jaouhara/Documents/Projects/embedding_models/fasttext/cc.fr.300.bin"}


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = "super secret key"

# read raw file into pandas data frame
def read_raw_file(raw_filename, source_dir):
    raw_filename = os.path.basename(raw_filename)
    data = None
    filepath = source_dir + '/' + raw_filename

    try:
        f = open(filepath, 'r')
        data = f.read()
        table = json.loads(data)
        f.close()

        cols = table['relation']
        ncols = len(table['relation'])
        num_rows = [len(col) for col in cols]
        max_num_rows = max(num_rows)
        
        # df = pd.DataFrame(table['cols'])
        return ncols, cols, max_num_rows,  ""

    except OSError:
        return -1, -1, -1, f"Couldn't open raw data file {filepath}."

# check if table has header
def has_header(table):
    if table['hasHeader'] and table['headerPosition'] == 'FIRST_ROW':
        return True
    return False

# read raw file into pandas data frame
def read_metadata_file(binary_filename, clean_data_dir, raw_data_dir):
    # raw files are named t[table_id]_c[column_pos].json, mus extract table id and colun pos from binary raw filename
    table_id = int(re.search(r'_t(.*)c', binary_filename).group(1))
    column_pos = int(re.search(fr'_t{table_id}c(.*?)_', binary_filename).group(1))
    filename = f"t{table_id}_c{column_pos}.json"

    data = None
    filepath = clean_data_dir + '/' + filename
    raw_file_name = ""

    try:
        f = open(filepath, 'r')
        data = f.read()
        table = json.loads(data)
        raw_file_name = table["raw_file"]
        f.close()
    except OSError:
        return -1, f"Couldn't open clean data file {filepath}."
    
    raw_filepath = raw_data_dir + '/' + raw_file_name

    try:
        f = open(raw_filepath, 'r')
        raw_data = f.read()
        table = json.loads(raw_data)
        f.close()

        metadata = {"page_title" : "", "title": "", "header": [], "url": "",
        "text_before_table" : "", "text_after_table": "", "last_modified": "", 
        "raw_data_file": raw_file_name}
        
        if "pageTitle" in table:
            metadata["page_title"]  = table["pageTitle"]
        if "title" in table:
            metadata["title"]  = table["title"]

        if has_header(table):
            for col in table['relation']:
                metadata["header"].append(col[0])

        if "url" in table:
            metadata["url"]  = table["url"]
        if "textBeforeTable" in table:
            metadata["text_before_table"]  = table["textBeforeTable"]
        if "textAfterTable" in table:
            metadata["text_after_table"]  = table["textAfterTable"]
        if "lastModified" in table:
            metadata["last_modified"]  = table["lastModified"]
        
        print(metadata)
        return metadata, ""

    except OSError:
        return -1, f"Couldn't open raw data file {raw_filepath}."

# create temp file for the query column retrieved from a csv file
def csv_to_bin_file(query_file , column_idx, query_lang):
    df = pd.read_csv(query_file, engine='python')
    if column_idx < 0 or column_idx > df.shape[1] - 1:
        return -1, f"Wrong value for column idx. file only contains {df.shape[1]} columns."
    else:
        query_column = df.iloc[:, column_idx].values
        query_table_id = re.search(r't(.*).csv', query_file.filename)
        
        if query_table_id:
            query_table_id = int(query_table_id.group(1))
        else:
            query_table_id = 0
        query_size, msg = query_to_bin(query_column, TMP_FOLDER, BIN_FOLDER, EMBEDDING_MODEL, PATH_TO_MODEL[query_lang], EMBEDDING_DIM, query_table_id)
        
        return query_size, msg

def keywords_to_bin_file(keywords, query_lang):
    kws = list(keywords.split(" "))
    if not keywords:
        return -1, f"No keywwords were found."
    else:
        query_size, msg = query_to_bin(kws, TMP_FOLDER, BIN_FOLDER, EMBEDDING_MODEL, PATH_TO_MODEL[query_lang], EMBEDDING_DIM)
        return query_size, msg

def get_keyword_query_results(query_size, top, k, approx_error):
    kashif_output = run_kashif(KASHIF_BIN, KASHIF_IDX, BIN_FOLDER, query_size, RESULTS_FOLDER, BIN_FOLDER,  EMBEDDING_DIM, top, k, approx_error, keyword_query=1, join_query=0)
    print("prog output:")
    print(f"**{kashif_output}**")

    files = list(re.findall(r'@@(.*?)\$', kashif_output)) # get file names without duplicates
    table_ids = list(re.findall(r'table-(.*?)\-', kashif_output))
    min_distances = list(re.findall(r'min_distance=(.*?)\§', kashif_output))
    num_closest = list(re.findall(r'num_closest=(.*?)\#', kashif_output))
    total_matches = list(re.findall(r'total_matches=(.*?)\µ', kashif_output))
    query_time = "{:.2f}".format(float(re.search(r'query_time=(.*?)sec', kashif_output).group(1)))
    print(min_distances)
    print(f"Query time = {query_time} sec\n")
    temp = set(list(zip(files, table_ids, min_distances, num_closest, total_matches)))
    visited = set()
    results = list()
    for (f, table_id, d, nc, t) in temp:
        if (f, table_id) not in visited: 
            visited.add((f, table_id))
            results.append((f, int(table_id), float(d), int(nc), int(t)))
        
    # sort results by overlap 
    results.sort(key=lambda x: x[2])
    results.sort(key=lambda x: x[3], reverse=True)
    # results.sort(key=lambda x: x[4], reverse=True)

    metadata = list()
    for result in results:
        metad, msg = read_metadata_file(result[0], CLEAN_DATA_FOLDER, RAW_DATA_FOLDER)
        if metad == -1:
            print(msg)
            return render_template('index.html', error=msg)
        else:
            metadata.append(metad)
    
    results = list(zip(results, metadata))
    total_results = len(results)

    return total_results, results, metadata, query_time

def get_join_query_results(query_size, top, k, approx_error):
    kashif_output = run_kashif(KASHIF_BIN, KASHIF_IDX, BIN_FOLDER, query_size, RESULTS_FOLDER, BIN_FOLDER,  EMBEDDING_DIM, top, k, approx_error, keyword_query=0, join_query=1)
    print("prog output:")
    print(f"**{kashif_output}**")

    files = list(re.findall(r'@@(.*?)\$', kashif_output)) # get file names without duplicates
    table_ids = list(re.findall(r'table-(.*?)\-', kashif_output))
    column_pos = list(re.findall(r'column-(.*?)\-', kashif_output))
    overlaps = list(re.findall(r'overlap=(.*?)\§', kashif_output))
    query_time = "{:.2f}".format(float(re.search(r'query_time=(.*?)sec', kashif_output).group(1)))

    print(f"Query time = {query_time} sec\n")
    temp = set(list(zip(files, table_ids, column_pos, overlaps)))
    visited = set()
    results = list()
    for (f, tid, sid, o) in temp:
        if o == '0':
            continue
        if (f, tid, sid) not in visited: 
            visited.add((f, tid, sid))
            results.append((f,  int(tid), int(sid)+1, int(o)))
            print(f"{f}, {tid}, {sid}, {o}")
        
    # sort results by overlap 
    results.sort(key=lambda x:x[3], reverse=True)
    metadata = list()
    for result in results:
        metad, msg = read_metadata_file(result[0], CLEAN_DATA_FOLDER, RAW_DATA_FOLDER)
        if metad == -1:
            return render_template('index.html', error=msg)
        else:
            metadata.append(metad)
    
    results = list(zip(results, metadata))
    total_results = len(results)

    return total_results, results, metadata, query_time

def run_kashif(kashif_bin, kashif_idx, bin_folder, query_size, result_dir, dataset_folder, embedding_size, num_top, k, approx_error, keyword_query=0, join_query=1): 
    if(join_query == 1):
        query = subprocess.check_output([kashif_bin, '--index-path', kashif_idx, '--queries', bin_folder,
        '--nq', '1', '--queries-size',  str(query_size), 
        '--min-qset-size', str(query_size), '--max-qset-size', str(query_size+1),
        '--dataset', dataset_folder, '--total-data-files', '100000', 
        '--dataset-GB-size', '1', '--dataset-size', '5032000',
        '--result-dir', result_dir, '--k', str(k),
        '--top',  str(num_top), ' --delta', '1',
        '--epsilon',  str(approx_error), '--timeseries-size', str(embedding_size),
        '--track-bsf', '--incremental', '--leaf-size', '1000',
        '--buffer-size', '100',
        '--mode', '1',  '--warping', '0.0', '--ascii-input', '0',
        '--track-vector', '1'
        ])
    elif(keyword_query == 1):
        query = subprocess.check_output([kashif_bin, '--index-path', kashif_idx, '--queries', bin_folder,
        '--nq', '1', '--queries-size',  str(query_size), 
        '--min-qset-size', str(query_size), '--max-qset-size', str(query_size+1),
        '--dataset', dataset_folder, '--total-data-files', '100000', 
        '--dataset-GB-size', '1', '--dataset-size', '5032000',
        '--result-dir', result_dir, '--k', str(k),
        '--top',  str(num_top), ' --delta', '1',
        '--epsilon',  str(approx_error), '--timeseries-size', str(embedding_size),
        '--track-bsf', '--incremental', '--leaf-size', '1000',
        '--buffer-size', '100',
        '--mode', '1',  '--warping', '0.0', '--ascii-input', '0',
        '--track-vector', '1',
        '--keyword-search', '1',
        ])

    # for the index storing 500k tables
    # query = subprocess.check_output([kashif_bin, '--index-path', kashif_idx, '--queries', bin_folder,
    # '--nq', '1', '--queries-size',  str(query_size), 
    # '--min-qset-size', str(query_size), '--max-qset-size', str(query_size+1),
    # '--dataset', dataset_folder, '--total-data-files', '500000', 
    # '--dataset-GB-size', '3', '--dataset-size', '14574511',
    # '--result-dir', result_dir, '--k', str(k),
    # '--top',  str(num_top), ' --delta', '1',
    # '--epsilon',  str(approx_error), '--timeseries-size', str(embedding_size),
    # '--track-bsf', '--incremental', '--leaf-size', '100',
    # '--buffer-size', '500000',
    # '--mode', '1',  '--warping', '0.0', '--ascii-input', '0',
    # '--track-vector', '1'
    # ])

    # stdout, stderr = query.communicate()
    return query.decode('utf-8')

@app.route('/404')
def error():
   return render_template('404.html')


@app.route('/')
def home():
   return render_template('index.html')

@app.route('/query', methods = ['GET', 'POST'])
def process_query():
    if request.method == 'POST':
        keyword_query = 0
        join_query = 0
        # check if the post request has the file part
        if 'query_file' not in request.files and 'query_keywords' not in request.form:
            print('No query was submited')
            return redirect(request.url)
        
        # read form data
        input_file = request.files['query_file']
        column_idx = int(request.form['column_idx']) - 1
        top = int(request.form['top'])
        approx_error = float(request.form['approx_error'])
        k = int(request.form['k'])
        query_lang = str(request.form['query_lang'])

        if 'query_file' in request.files and input_file.filename != '':
            join_query = 1
        elif 'query_keywords' in request.form:
            keyword_query = 1
        
        total_results, results, metadata, query_time = None, None, None, None

        if(join_query):
            # upload query to server
            start = time.time()
            query_size, msg= csv_to_bin_file(request.files.get('query_file'), column_idx, query_lang)
            end = time.time()
            query_cleaning_time =  "{:.2f}".format(end - start)
            if msg:
                return render_template("index.html", error=msg)
            
            print("cleaned join query\n")
            total_results, results, metadata, query_time = get_join_query_results(query_size, top, k, approx_error)


        elif(keyword_query):
            # upload query to server
            start = time.time()
            query_size, msg= keywords_to_bin_file(str(request.form['query_keywords']), query_lang)
            end = time.time()
            query_cleaning_time =  "{:.2f}".format(end - start)
            if msg:
                return render_template("index.html", error=msg)
            print("cleaned keyword query\n")
            total_results, results, metadata, query_time = get_keyword_query_results(query_size, top, k, approx_error)
        
        else:
            return render_template("index.html", error="Couldn't know which type of query to run.")

        return render_template("index.html", total_results=total_results, results=results, metadata=metadata, query_cleaning_time=query_cleaning_time, query_time=query_time, keyword_query=keyword_query)

@app.route('/view-dataset', methods=['GET', 'POST'])
def view_dataset():
    if request.method == 'POST':
        filename = request.form['file_name']
        column_idx = request.form['column_idx']
        table_id = request.form['table_id']

        print(f"(+) table id = {table_id}")

        if not filename or not column_idx:
            return render_template('404.html')

        else:
            metadata, msg = read_metadata_file(filename, CLEAN_DATA_FOLDER, RAW_DATA_FOLDER)
            if metadata == -1:
                return render_template('view-dataset.html', error=msg)
            
            ncols, cols, max_num_rows, msg = read_raw_file(metadata["raw_data_file"], RAW_DATA_FOLDER)
            if ncols == -1:
                return render_template('view-dataset.html', error=msg)
            

            print(metadata)
            
            return render_template('view-dataset.html', metadata=metadata, file_name=request.form['file_name'], table_id=table_id, match_col=column_idx, cols=cols, ncols=ncols, max_num_rows=max_num_rows)


if __name__ == '__main__':
   app.run(debug = True)