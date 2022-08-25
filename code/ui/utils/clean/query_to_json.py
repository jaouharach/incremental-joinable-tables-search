# EXTRACT TABLES FROM WDC COMPLETE CORPUS 2015 DATASET
import json
import os, re
import shutil

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

def query_to_json(query_array, target_dir):
    # check if target dir exists
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print("Error in query_to_json.py target directory does not exist.")

    
    new_array = clean_array(query_array)
     # skip empty tables
    if not new_array or len(new_array) == 0:
        return -1, -1 # query column is a numerical column!
    else:
        query_size = len(new_array)
        data = {'id': 0, 'ncols':  1, 'cols': [new_array]}

        
        # delete all temp files in temp folder
        if(clear_folder(target_dir)):
            # save table to temp folder
            file = save_file(data, target_dir)
            print(f'Success! output file has been created, file name : {file}\n')
            return file, query_size
        else:
            print("Error query_to_json.py: couldn't clear target directory.")

   
    



