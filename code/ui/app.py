import os
from flask import Flask, flash, render_template, request, redirect
import pandas as pd
import subprocess

import sys
sys.path.append('./utils/clean')
from query_to_json import *

UPLOAD_FOLDER = "./uploads/"
TMP_FOLDER = UPLOAD_FOLDER + "tmp/"
BIN_FOLDER = UPLOAD_FOLDER + "bins/"
JSON_FOLDER = UPLOAD_FOLDER + "json/"
ALLOWED_EXTENSIONS = {'csv', 'json'}
GLOVE_PATH = "./glove/glove.6B.50d.txt"
EMBEDDING_DIM = 50


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def create_bin_query_file(tmp_query_file, glove_file):
    embedding_process = subprocess.Popen(['python3', './utils/embed/embed_2.py', TMP_FOLDER, JSON_FOLDER,  GLOVE_PATH, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = embedding_process.communicate()
    
    tobin_process = subprocess.Popen(['python3', './utils/tobin/tobin.py', JSON_FOLDER, BIN_FOLDER, str(EMBEDDING_DIM)],
                                stdout=subprocess.PIPE)
    stdout = tobin_process.communicate()
    
    return stdout
    


def create_temp_query_file(query_file , column_idx):
    df = pd.read_csv(query_file)
    print(f"shape = {df.shape}")
    if(column_idx > df.shape[1]):
        render_template("index.html", data = f"Wrong value for column idx. file only contains {df.shape[1]} columns")
    else:
        query_column = df.iloc[:, column_idx].values
        print(f"col = {query_column}")

        query_file = query_to_json(query_column, TMP_FOLDER)
        if(query_file == -1):
            render_template("index.html", data = "Column has no text value!")
            print("Column has no text value!")
        else:
            return query_file


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

        print(f"idx = {column_idx}")
        # empty file without a filename.
        if input_file.filename == '':
            flash('No selected file')
            return redirect(request.url)

        # upload query column to server
        tmp_query_file = create_temp_query_file(request.files.get('query_file'), column_idx)
        print(tmp_query_file)

        # process the query and save it to a binary file
        print("process:\n")
        bin_query_file = create_bin_query_file(tmp_query_file, GLOVE_PATH)
        print(bin_query_file)
        return render_template("index.html")
        
if __name__ == '__main__':
   app.run(debug = True)