# JSON table encoder, each vector is converted to 32 bits (base 2) the binary file stores vectors column wise, for every column we store [total vecotors][vector 1][vector 2]...[vector n]
# Binary conversion: little indian.
import struct
import numpy as np
import json
import time, os, sys, glob
import math 
# read table from json file
def get_json_table(file):
    data = None
    with open(file, 'r') as f:
        data = f.read()
    table = json.loads(data)
    return table

# save to .bin file
def save_to_file(target_dir, tableid, bin_table, totalvec, ncols, embedding_dim):
    file_name = f'{target_dir}/data_size{totalvec}_t{tableid}c{ncols}_len{embedding_dim}_noznorm.bin'
    f = open(file_name, "wb")
    f.write(bin_table)
    f.close()

def float_to_bytes(num):  # For testing.
    return  struct.pack('<f', num)

def int_to_bytes(num):
    return num.to_bytes(4,'little')

def vec_to_bin(vec):
    bin_vector = b''
    for _float in vec:
        if math.isnan(float(_float)):
            print(f"nan value in {_float}")
            exit(1)
        bin_vector = bin_vector + float_to_bytes(float(_float))
    return bin_vector

def col_to_bin(col):
    nvec = len(col)
    bin_col = int_to_bytes(int(nvec))
    for vec in col:
        bin_vec = vec_to_bin(vec)
        bin_col = bin_col + bin_vec 
    return bin_col, nvec

def table_to_bin(table):
    bin_table = b''
    totalvec = 0
    for i in range(table['ncols']):
        bin_col, nvec = col_to_bin(table['cols'][i])
        bin_table = bin_table + bin_col
        totalvec = totalvec + nvec

    return bin_table, table['ncols'], totalvec

def encode_all(source_dir, target_dir, embedding_dim):
    table_count = 0

    for input_file in glob.glob(os.path.join(source_dir, '*.json')):
        print(f"> Processing file: {input_file}")
        table = get_json_table(input_file)
        table_id = table['id']
        

        # convert table to binary format
        bin_table, ncols, totalvec = table_to_bin(table)
        print(f'total vectors = {totalvec}')
        # Save to file
        save_to_file(target_dir, table_id, bin_table, totalvec, ncols, embedding_dim)
        table_count = table_count + 1

        print(f"< File {input_file} has been processed.\n")

    return table_count

if len(sys.argv) >= 4:
    source_dir = sys.argv[1] # file where embeddings data was stored
    target_dir = sys.argv[2] # directory to store binary files
    embedding_dim = sys.argv[3] # dimension of embeddings

    # create target directory if doesn't exist
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print("Output directory " , target_dir ,  " created\n\n\n")
    else:    
        print("Output directory " , target_dir ,  " already exists\n\n\n")

    start = time.time()
    num_tables = encode_all(source_dir, target_dir, embedding_dim)
    runtime = time.time() - start

    print(f'Congrats! Generate binary files process is done. {num_tables} files are created')
    print(f"Execution time: {runtime} seconds")

else:
    print('Please provide source directory, target directory and embbeding dimension as args!')


