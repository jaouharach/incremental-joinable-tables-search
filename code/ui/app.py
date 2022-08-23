import os
from flask import Flask, flash, render_template, request, redirect
import pandas as pd
import subprocess
import shutil

import sys
sys.path.append('./utils/clean')
from query_to_json import *

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

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# convert query file into a binary file in the format [total_vectors_in column][vector_1][vector_2]...
# every [value] is stored in 4 bytes
def create_bin_query_file(tmp_query_file, glove_file):
    embedding_process = subprocess.Popen(['python3', './utils/embed/embed_2.py', TMP_FOLDER, JSON_FOLDER,  GLOVE_PATH, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = embedding_process.communicate()
    
    tobin_process = subprocess.Popen(['python3', './utils/tobin/tobin.py', JSON_FOLDER, BIN_FOLDER, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = tobin_process.communicate()
    
    return stdout
    
# create temp file for the query column
def create_temp_query_file(query_file , column_idx):
    df = pd.read_csv(query_file)
    print(f"shape = {df.shape}")
    if(column_idx > df.shape[1]):
        render_template("index.html", data = f"Wrong value for column idx. file only contains {df.shape[1]} columns")
    else:
        query_column = df.iloc[:, column_idx].values
        print(f"col = {query_column}")

        query_file, query_size = query_to_json(query_column, TMP_FOLDER)
        if(query_file == -1):
            render_template("index.html", data = "Column has no text value!")
            print("Column has no text value!")
        else:
            return query_file, query_size

def clear_folders(upload_folders):
    for folder in upload_folders:
        with os.scandir(folder) as entries:
            for entry in entries:
                if entry.is_dir() and not entry.is_symlink():
                    shutil.rmtree(entry.path)
                else:
                    os.remove(entry.path)
    return True

def run_kashif(kashif_bin, kashif_idx, bin_folder, query_size, result_dir, dataset_folder, embedding_size, top_k, approx_error):
    
    query_process = subprocess.call([kashif_bin, '--index-path', kashif_idx, '--queries', bin_folder,
    '--nq', '1', '--queries-size',  str(query_size), 
    ' --min-qset-size', str(query_size), '--max-qset-size', str(query_size),
    '--dataset', dataset_folder, '--total-data-files', '100',
    '--result-dir', result_dir, '--k', '100', '--top',  str(top_k), ' --delta', '1',
    '--epsilon',  str(approx_error), '--timeseries-size', str(embedding_size),
    '--track-bsf', '--incremental', '--leaf-size', '100',
    '--mode', '1',  '--warping', '0.0', '--ascii-input', '0'
    ])
    # query_process.communicate()

    
    
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
        column_idx = int(request.form['column_idx'])
        top = int(request.form['top'])
        approx_error = float(request.form['approx_error'])

        # file without a filename.
        if input_file.filename == '':
            flash('No selected file')
            return redirect(request.url)

        # empty upload files
        if(clear_folders([TMP_FOLDER, JSON_FOLDER, BIN_FOLDER]) != True):
            return render_template("index.html", data="Could't clear upload folders")

        # upload query column to server
        tmp_query_file, query_size = create_temp_query_file(request.files.get('query_file'), column_idx)

        # process the query and save it to a binary file
        create_bin_query_file(tmp_query_file, GLOVE_PATH)

        kashif_output = run_kashif(KASHIF_BIN, KASHIF_IDX, BIN_FOLDER, query_size, RESULTS_FOLDER, BIN_FOLDER,  EMBEDDING_DIM, top, approx_error)
        print(kashif_output)


        return render_template("index.html")
        
if __name__ == '__main__':
   app.run(debug = True)