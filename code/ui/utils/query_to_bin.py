# EXTRACT TABLES FROM WDC COMPLETE CORPUS 2015 DATASET
import json
import os, re
import shutil
import embed, tobin

def clear_folder(folder_path):
    with os.scandir(folder_path) as entries:
        for entry in entries:
            if entry.is_dir() and not entry.is_symlink():
                shutil.rmtree(entry.path)
            else:
                os.remove(entry.path)
    return True

def clean_array(array):
    # remove digit values
    no_digit = [re.sub('\W+',' ', str(x).lower()) for x in array if not (str(x) == "" or str(x) =="nan" or str(x).isdigit() or str(x)[0] == '-' and str(x[1:]).isdigit())]
    # remove duplicates
    distinct_values = list(set(no_digit))
    return distinct_values


def save_file(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()

    return output_file

def save_bin_file(target_dir, tableid, bin_table, totalvec, ncols, embedding_dim):
    file_name = f'{target_dir}/data_size{totalvec}_t{tableid}c{ncols}_len{embedding_dim}_noznorm.bin'
    f = open(file_name, "wb")
    f.write(bin_table)
    f.close()

def query_to_bin(query_array, text_target_dir, bin_target_dir, path_to_glove_file, embedding_dim):
    # check if target dir exists
    if not os.path.exists(text_target_dir) or not os.path.exists(bin_target_dir):
        return -1, "Target directory doesn't exist"

    # 1- clean data
    new_array = clean_array(query_array)
    query_size = 0
     # skip empty tables
    if not new_array or len(new_array) == 0:
        return -1, "Query column is a numerical column!"
    else:
        query_size = len(new_array)
        query_data = {'id': 0, 'ncols':  1, 'cols': [new_array]}
        
        # delete all temp files in temp folder
        if(clear_folder(text_target_dir) and clear_folder(bin_target_dir)):
            # save query (in text format) to temp folder
            file = save_file(query_data, text_target_dir)
            if not file:
                return -1, "Failed to save temp query file."
        else:
            return -1, "Failed to clear target folders."


    print(f"\nquery size = {query_size}")
    print(query_data)
    print("\n")
    
     # 2- embed
    query_embeddings = embed.embed_query(query_data, path_to_glove_file, embedding_dim)
    if query_embeddings == -1:
        return -1, "Faild to generate query embeddings."

    # 3- create binary file
    _ = tobin.encode_query(query_embeddings, bin_target_dir, embedding_dim)
    if not _:
        return -1, "Faild to create binary file."
    
    return query_size, ""



   
    



