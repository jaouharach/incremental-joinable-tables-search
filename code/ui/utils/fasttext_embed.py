# CONVERT SETS OF TOKENS TO SETS OF VECTORS (EMBEDDINGS) USING GLOVE MODEL
import fasttext

FASTTEXT_PATH = '/home/jchanchaf/local/src/fasttext-en/cc.en.300.bin'

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
        return -1
    else:
        return new_query
        

def load_embeddings_index(path_to_fasttext_model=FASTTEXT_PATH, embeddings_dim=100):
    # load Fast text model
    ft = fasttext.load_model(path_to_fasttext_model)
    fasttext.util.reduce_model(ft, embeddings_dim)

def encode_col(col, ft):
    col_embeddings = [ft.get_word_vector(token).tolist() for token in col]
    return col_embeddings
    