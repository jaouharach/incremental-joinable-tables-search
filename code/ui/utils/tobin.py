# JSON table encoder, each vector is converted to 32 bits (base 2) the binary file stores vectors column wise, for every column we store [total vecotors][vector 1][vector 2]...[vector n]
# Binary conversion: little indian.
import struct
import json

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
    return 1

def float_to_bytes(num):  # For testing.
    return  struct.pack('<f', num)

def int_to_bytes(num):
    return num.to_bytes(4,'little')

def vec_to_bin(vec):
    bin_vector = b''
    for _float in vec:
        bin_vector = bin_vector + float_to_bytes(_float)
    return bin_vector

def col_to_bin(col):
    nvec = len(col)
    bin_col = int_to_bytes(int(nvec))
    for vec in col:
        bin_col = bin_col +  vec_to_bin(vec)
    return bin_col, nvec

def table_to_bin(table):
    bin_table = b''
    totalvec = 0
    for i in range(table['ncols']):
        bin_col, nvec = col_to_bin(table['cols'][i])
        bin_table = bin_table + bin_col
        totalvec = totalvec + nvec

    return bin_table, table['ncols'], totalvec

def encode_query(query_embeddings, target_dir, embedding_dim):
    query_id = query_embeddings['id']
    # convert table to binary format
    bin_query, ncols, totalvec = table_to_bin(query_embeddings)
    # Save to file
    if not save_to_file(target_dir, query_id, bin_query, totalvec, ncols, embedding_dim):
        return False
    return True

 
