# EXTRACT TABLES FROM WDC COMPLETE CORPUS 2015 DATASET
import json
import re
import os, sys, glob
import time

from importlib_metadata import metadata

def get_json_table(file):
    data = None
    with open(file, 'r') as f:
        data = f.read()
    table = json.loads(data)
    return table

def remove_header(table):
    if has_header(table):
        for col in table['relation']:
            col.pop(0)
        # print('header has been removed.')

    return table['relation']

def get_metadata(table, raw_data_file):
    metadata = {"page_title" : "", "title": "", "header": [],"url": "",
     "text_before_table" : "", "text_after_table": "", "last_modified": "", "raw_data_file": raw_data_file}
     
    if "pageTitle" in table:
        metadata["page_title"]  = table["pageTitle"]
    if "title" in table:
        metadata["title"]  = table["title"]

    if has_header(table):
        for col in table['relation']:
            metadata["header"].append(col[0])

    if "url" in table:
        metadata["url"]  = table["url"]
    if "textBeforeTable" in table:
        metadata["text_before_table"]  = table["textBeforeTable"]
    if "textAfterTable" in table:
        metadata["text_after_table"]  = table["textAfterTable"]
    if "lastModified" in table:
        metadata["last_modified"]  = table["lastModified"]
    
    return metadata
    

def has_header(table):
    if table['hasHeader'] and table['headerPosition'] == 'FIRST_ROW':
        return True
    return False


def remove_digit_cols(rel):
    # remove empty and digit columns + sort columns
    col_to_delete = []

    for col in range(len(rel)):
        rel[col].sort()  # sort column
        c = re.sub('[^A-Za-z0-9 ]+', '', ''.join(rel[col]))
        if c.isdigit() or c.strip() == '':
            col_to_delete.append(rel[col])

    for col in col_to_delete:
        rel.remove(col)

    return rel

def remove_duplicates(rel):
    # remove special charcters and duplicates and convert all to lowercase
    for col in range(len(rel)):
        c = [re.sub('[^A-Za-z0-9 ]+', ' ', token).lower().rstrip('\r\n') for token in rel[col]]
        c = filter(lambda word: word.strip(), c) # remove blanck spaces and empty words
        rel[col] = c
        rel[col] = list(dict.fromkeys(rel[col]))  # remove duplicates
    
    return rel


def clean_relation(rel):
    rel = remove_digit_cols(rel)
    rel = remove_duplicates(rel)
    return rel

def save_raw_table(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()
    print(f'Raw table {table["id"]} with {table["ncols"]} columns has been processed.')

def save_table(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()
    print(f'Table {table["id"]} with {table["ncols"]} columns has been processed.')

def save_table_metadata(table, metadata, metadata_target_dir):
    output_file = f'{metadata_target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, ensure_ascii=False)
        f.close()
    print(f'Metadata file has been created.')


def extract_tables(source_dir, target_dir):
    tableid = 1
    # iterate over json files in source directory
    for input_file in glob.glob(f'{source_dir}/**/*.json', recursive = True):
        print(f"Processing file: {input_file}")
        table = get_json_table(input_file)
        metadata = get_metadata(table, input_file)
        
        # skip vertical tables
        if table['tableOrientation'] == 'VERTICAL':
            continue

        else:
            # clean data
            raw_relation = table['relation']
            # relation = remove_header(table)
            relation = clean_relation(raw_relation)
            # skip empty tables
            if len(relation) == 0:
                continue
                
            raw_file_name = os.path.basename(input_file)
            table = {'id': tableid, 'ncols':  len(relation), 'cols': relation, 'raw_file' : raw_file_name}

            # save table to disk
            save_table(table, target_dir)
            # save metadata to disk
            tableid = tableid + 1

if len(sys.argv) >= 3:
    source_dir = sys.argv[1]
    target_dir = sys.argv[2] # directory where all json files will be stored

    # create target directory if doesn't exist
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print("Output directory " , target_dir ,  " created\n\n\n")

    start = time.time()
    extract_tables(source_dir, target_dir)
    end = time.time()

    print(f'Success! All files have been procecessed.\nClean tables are stored in {target_dir}\n')
    print(f'Cleaning process took {end-start} seconds.')

else:
    print('Please provide source directory, target directory and raw data directory and metadata directory as args!')
    # ex: python3  clean-tables.py ../data/raw-files ../data/clean-tables ../data/meta

