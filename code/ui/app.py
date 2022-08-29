from flask import Flask, flash, render_template, request, redirect
import pandas as pd
import subprocess, re
import time


import sys
sys.path.append('./utils')
from query_to_bin import *

RESULTS_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/query_results/"
UPLOAD_FOLDER = "/home/jaouhara/Documents/Projects/iqa-demo/code/ui/data/uploads/"
TMP_FOLDER = UPLOAD_FOLDER + "tmp/"
BIN_FOLDER = UPLOAD_FOLDER + "bins/"
JSON_FOLDER = UPLOAD_FOLDER + "json/"
ALLOWED_EXTENSIONS = {'csv', 'json'}
GLOVE_PATH = "./glove/glove.6B.50d.txt"
EMBEDDING_DIM = 50
KASHIF_BIN = "../kashif/bin/dstree"
KASHIF_IDX = "../kashif_idx/"

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    
# create temp file for the query column
def create_query_bin_file(query_file , column_idx):
    df = pd.read_csv(query_file)
    if column_idx < 0 or column_idx > df.shape[1] - 1:
        return -1, f"Wrong value for column idx. file only contains {df.shape[1]} columns"
    else:
        query_column = df.iloc[:, column_idx].values
        query_size, msg = query_to_bin(query_column, TMP_FOLDER, BIN_FOLDER, GLOVE_PATH, EMBEDDING_DIM)
        
        return query_size, msg


def run_kashif(kashif_bin, kashif_idx, bin_folder, query_size, result_dir, dataset_folder, embedding_size, num_top, k, approx_error): 
    query = subprocess.check_output([kashif_bin, '--index-path', kashif_idx, '--queries', bin_folder,
    '--nq', '1', '--queries-size',  str(query_size), 
    '--min-qset-size', str(query_size), '--max-qset-size', str(query_size+1),
    '--dataset', dataset_folder, '--total-data-files', '100', 
    '--dataset-GB-size', '1', '--dataset-size', '120',
    '--result-dir', result_dir, '--k', str(k),
    '--top',  str(num_top), ' --delta', '1',
    '--epsilon',  str(approx_error), '--timeseries-size', str(embedding_size),
    '--track-bsf', '--incremental', '--leaf-size', '100',
    '--buffer-size', '100',
    '--mode', '1',  '--warping', '0.0', '--ascii-input', '0',
    '--track-vector', '1'
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
        # check if the post request has the file part
        if 'query_file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        
        # read form data
        input_file = request.files['query_file']
        column_idx = int(request.form['column_idx']) - 1
        top = int(request.form['top'])
        approx_error = float(request.form['approx_error'])
        k = int(request.form['k'])

        # file without a filename.
        if input_file.filename == '':
            flash('No selected file')
            return redirect(request.url)

        # upload query column to server
        start = time.time()
        query_size, msg= create_query_bin_file(request.files.get('query_file'), column_idx)
        end = time.time()
        query_cleaning_time =  "{:.2f}".format(end - start)

        if msg:
            return render_template("index.html", error=msg)
        
        kashif_output = run_kashif(KASHIF_BIN, KASHIF_IDX, BIN_FOLDER, query_size, RESULTS_FOLDER, BIN_FOLDER,  EMBEDDING_DIM, top, k, approx_error)
        print("prog output:")
        print(f"**{kashif_output}**")

        files = list(re.findall(r'@@(.*?)\$', kashif_output)) # get file names without duplicates
        column_pos = list(re.findall(r'column-(.*?)\-', kashif_output))
        overlaps = list(re.findall(r'overlap=(.*?)\§', kashif_output))
        query_time = "{:.2f}".format(float(re.search(r'query_time=(.*?)sec', kashif_output).group(1)))

        print(f"Query time = {query_time} sec\n")
        temp = set(list(zip(files, column_pos, overlaps)))
        visited = set()
        results = list()
        for (f, sid, o) in temp:
            if o == '0':
                continue
            if (f, sid) not in visited: 
                visited.add((f, sid))
                results.append((f, int(sid)+1, int(o)))
                print(f"{f}, {sid}, {o}")
            
        # sort results by overlap 
        results.sort(key=lambda x:x[2], reverse=True)
        return render_template("index.html", results=results, query_cleaning_time=query_cleaning_time, query_time=query_time)

@app.route('/view-dataset', methods=['GET'])
def view_dataset():
   return render_template('view-dataset.html')

if __name__ == '__main__':
   app.run(debug = True)