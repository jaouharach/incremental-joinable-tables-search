# CONVERT SETS OF TOKENS TO SETS OF VECTORS (EMBEDDINGS) USING GLOVE MODEL
import numpy as np
import math

def embed_query(query, path_to_glove_file, embedding_dim):
    # load embeddings
    embeddings_index = load_embeddings_index(path_to_glove_file)
    new_query = {'id': query['id'], 'ncols':  0, 'cols': []}

    for col in query['cols']:
        vectors_col = encode_col(col, embeddings_index, embedding_dim)
        if(len(vectors_col) > 0):
            new_query['ncols'] += 1
            new_query['cols'].append(vectors_col)

    # no columns in table
    if(new_query['ncols'] == 0):
        return False
    else:
        return new_query
        

def load_embeddings_index(path_to_glove_file):
    # print("Load embeddings index...")
    embeddings_index = {}
    with open(path_to_glove_file) as f:
        for line in f:
            word, coefs = line.split(maxsplit=1)
            coefs = np.fromstring(coefs, "f", sep=" ")
            embeddings_index[word] = coefs
    return embeddings_index

def encode_col(col, embeddings_index, embedding_dim=50):
    col_embeddings = list()
    # compute the average vector (embedding) for each token, don't count words for which there is no embedding in glove index
    for token in col:
        token_embedding = np.zeros(embedding_dim, dtype = float) # zeros array 
        words = token.split(' ')
        total_token_embeddings = 0
        for word in words:
            if not word:
                continue
            embedding = embeddings_index[word]
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
    