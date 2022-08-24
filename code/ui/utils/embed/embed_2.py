# CONVERT SETS OF TOKENS TO SETS OF VECTORS (EMBEDDINGS) USING GLOVE MODEL
import numpy as np
import json
import os, sys, glob
import time
import math

# read table from json file
def get_json_table(file):
    data = None
    with open(file, 'r') as f:
        data = f.read()
    table = json.loads(data)
    return table

def load_embeddings_index(path_to_glove_file):
    # print("Load embeddings index...")
    embeddings_index = {}
    with open(path_to_glove_file) as f:
        for line in f:
            word, coefs = line.split(maxsplit=1)
            coefs = np.fromstring(coefs, "f", sep=" ")
            embeddings_index[word] = coefs
    # print("Found %s word vectors." % len(embeddings_index))
    # print("Done.")
    return embeddings_index

def encode_col(col, embeddings_index, embedding_dim=50):
    col_embeddings = list()
    # compute the average vector (embedding) for each token, don't count words for which there is no embedding in glove index
    for token in col:
        token_embedding = np.zeros(embedding_dim, dtype = float) # zeros array 
        words = token.split(' ')
        total_token_embeddings = 0
        for word in words:
            embedding = embeddings_index.get(word)
            if embedding is not None:
                embedding = np.array(embedding)

                if(np.isnan(embedding).any() == True): # if embedding is in the form [..., nan, ...]
                    continue
                else:
                    # print("original vector:")
                    # print(embedding)
                    # print("normalized vector:")
                    norm = np.linalg.norm(embedding)
                    normalized_embedding = embedding/norm
                    # print(normalized_embedding)
                    token_embedding = np.add(token_embedding, normalized_embedding)
                    total_token_embeddings += 1

         # append average vector to column
        token_embedding = np.true_divide(token_embedding, total_token_embeddings)
        # sanity check: add embedding if it does not contain  any nan value
        if(len([x for x in token_embedding if math.isnan(x)]) == 0):
            if(np.all(token_embedding == 0) == False):
                col_embeddings.append(list(token_embedding))
    
    return col_embeddings
    

def save_table(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()
    # print(f'Table {table["id"]} with {table["ncols"]} columns has been processed.')
    


def encode_tables(source_dir, target_dir, path_to_glove_file, embedding_dim):
    # check if input dir exists
    if not os.path.exists(f"{source_dir}"):
        exit(f"Coudn't find directory {source_dir}!")

    # load embeddings
    embeddings_index = load_embeddings_index(path_to_glove_file)
    embedding_dim = embedding_dim

    # iterate over files in source directory
    for input_file in glob.glob(os.path.join(source_dir, '*.json')):
        # print(f"Processing file: {input_file}")
        table = get_json_table(input_file)
        new_table = {'id': table['id'], 'ncols':  0, 'cols': []}

        for col in table['cols']:
            vectors_col = encode_col(col, embeddings_index, embedding_dim)
            if(len(vectors_col) > 0):
                new_table['ncols'] += 1
                new_table['cols'].append(vectors_col)

        # no columns in table
        if(new_table['ncols'] == 0):
            return -1
        else:
            # save new table, each table is saved in a separate json file
            save_table(new_table, target_dir)


if len(sys.argv) >= 5:
    source_dir = sys.argv[1]
    target_dir = sys.argv[2] # directory where all json files will be stored
    path_to_glove_file = sys.argv[3] # file of glove embeddings
    embedding_dim = int(sys.argv[4])

    # create target directory if doesn't exist
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)

    start = time.time()
    encode_tables(source_dir, target_dir, path_to_glove_file, embedding_dim)
    end = time.time()

    exit(0)

else:
    print('Please provide source directory, target directory, path to glove embeddings file, embedding dimension as arguments!')
        # ex: python3  tables-to-vectors.py ../data/clean-tables ../data/embeddings ./glove/glove.6B.50d.txt 50
