# CONVERT SETS OF TOKENS TO SETS OF VECTORS (EMBEDDINGS)
# Start process with command: python3 generate_embeddings.py [path_to_source_file] [path_to_target_file]
import json
import glob
import fasttext
import fasttext.util
import os
import time
import sys
import numpy as np

FASTTEXT_PATH = '/home/jchanchaf/local/src/fasttext-en/cc.en.300.bin'

# read table from json file
def get_json_table(file):
    data = None
    with open(file, 'r') as f:
        data = f.read()
    table = json.loads(data)
    return table

# save table to target dir
def save_table(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()
    print(f'Table {table["id"]} with {table["ncols"]} columns has been processed.')
    

def encode_col(col, ft, normalize_vectors):
    col_embeddings = list()
    if normalize_vectors:
        for token in col:
            vec = ft.get_word_vector(token).tolist()

            # normalize vector
            vec = np.array(vec)
            norm = np.linalg.norm(vec)
            norm_vec = list(vec/norm)
            
            if not np.isnan(norm_vec).any():
                col_embeddings.append(norm_vec)
    else:
        for token in col:
            vec = ft.get_word_vector(token).tolist()

            if not np.isnan(vec).any():
                col_embeddings.append(vec)

    return col_embeddings
    

def encode_tables(source_dir, target_dir, path_to_fasttext_model=FASTTEXT_PATH, embeddings_dim=100, normalize_vectors=True):
    # load Fast text model
    ft = fasttext.load_model(path_to_fasttext_model)
    fasttext.util.reduce_model(ft, embeddings_dim)

    # iterate over files in source directory
    for input_file in glob.glob(os.path.join(source_dir, '*.json')):
        print(f"Processing file: {input_file}")
        table = get_json_table(input_file)
        new_table = {'id': table['id'], 'ncols':  0, 'cols': []}

        for col in table['cols']:
            vectors_col = encode_col(col, ft, normalize_vectors)
            if(len(vectors_col) > 0):
                new_table['ncols'] += 1
                new_table['cols'].append(vectors_col)

        # save new table, each table is saved in a separate json file
        save_table(new_table, target_dir)


#command: python3 generate_embeddings.py '/mnt/beegfs/jchanchaf/src/wdc-experiment/extracted-data/tables' '/mnt/beegfs/jchanchaf/src/wdc-experiment/embedding-data/tables'

if len(sys.argv) >= 5:
    source_dir = sys.argv[1] #file where clean data was stored
    target_dir = sys.argv[2] #file to store data after converting each token to a vector
    path_to_fasttext_model = sys.argv[3]
    embedding_dim = int(sys.argv[4])
    normalize_vectors = False # normalize embedding vectors to be unit vectors

    # exit if source directory does not exist
    if not os.path.exists(source_dir):
        print("Error: Source directory " , source_dir ,  " doesn't exist\n")
        exit(1)

    # create target directory if doesn't exist
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print("Output directory " , target_dir ,  " does not exist. Output directory has been created\n")
    

    start = time.time()
    encode_tables(source_dir, target_dir, path_to_fasttext_model, embedding_dim, normalize_vectors)
    end = time.time()

    print(f'Success! All files have been procecessed. table embeddings are stored in {target_dir}\n')
    print(f'Embeddding process took {end - start} seconds.')

else:
    print('Please provide source directory, target directory, path to fasttext model (e.g. cc.en.300.bin), and embedding dimension as arguments!')
