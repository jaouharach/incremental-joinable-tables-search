from flask import Flask, flash, render_template, request, redirect, url_for
import pandas as pd
import subprocess, os, signal
import shutil, re


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
def create_bin_query_file():
    embedding_process = subprocess.Popen(['python3', './utils/embed/embed_2.py', TMP_FOLDER, JSON_FOLDER,  GLOVE_PATH, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = embedding_process.communicate()
    
    tobin_process = subprocess.Popen(['python3', './utils/tobin/tobin.py', JSON_FOLDER, BIN_FOLDER, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = tobin_process.communicate()
    
    bin_dir = os.listdir(BIN_FOLDER)
    json_dir = os.listdir(JSON_FOLDER)

    if len(bin_dir) == 0 or len(json_dir) == 0:
        return -1
    return 0
    
# create temp file for the query column
def create_temp_query_file(query_file , column_idx):
    df = pd.read_csv(query_file)
    print(f"col idx = {column_idx}")
    print(query_file)
    if column_idx < 0 or column_idx > df.shape[1] - 1:
        return -1, df.shape[1]
    else:
        query_column = df.iloc[:, column_idx].values
        print(f"col = {query_column}")

        query_file, query_size = query_to_json(query_column, TMP_FOLDER)
        if(query_file == -1):
            return -2, ""
        else:
            return query_file, query_size

# empty folder content
def clear_folders(upload_folders):
    for folder in upload_folders:
        with os.scandir(folder) as entries:
            for entry in entries:
                if entry.is_dir() and not entry.is_symlink():
                    shutil.rmtree(entry.path)
                else:
                    os.remove(entry.path)
    return True

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

        # empty upload files
        if(clear_folders([TMP_FOLDER, JSON_FOLDER, BIN_FOLDER]) != True):
            return render_template("index.html", error="Could't clear upload folders")

        # upload query column to server
        _, size = create_temp_query_file(request.files.get('query_file'), column_idx)
        if _ == -1:
            return render_template("index.html", error=f"Wrong value for column idx. file only contains {size} columns")
        elif _ == -2:
            return render_template("index.html", error="Please choose a column with text data.")
        
        # if upload query column to server succeeded
        query_size = size

        # process the query and save it to a binary file, if binary file not created then the chosen column had no text data
        if create_bin_query_file() == -1:
            return render_template("index.html", error="Please choose a column with text data.")


        kashif_output = run_kashif(KASHIF_BIN, KASHIF_IDX, BIN_FOLDER, query_size, RESULTS_FOLDER, BIN_FOLDER,  EMBEDDING_DIM, top, k, approx_error)
        # print("prog output:")
        print(f"**{kashif_output}**")

        files = list(re.findall(r'@@(.*?)\$', kashif_output)) # get file names without duplicates
        column_pos = list(re.findall(r'column-(.*?)\-', kashif_output))
        overlaps = list(re.findall(r'overlap=(.*?)\§', kashif_output))

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

        # print("files:")
        # print(files)
        # print("\n")
        # print("column postions:")
        # print(column_pos)
        # print("\n")
        # print("overlaps:")
        # print(overlaps)
        # print("\n")
        # print("results:")
        # print(results)
        # print("\n")


        return render_template("index.html", results=results)

@app.route('/view-dataset', methods=['GET'])
def view_dataset():
   return render_template('view-dataset.html')

if __name__ == '__main__':
   app.run(debug = True)