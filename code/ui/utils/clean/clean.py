# EXTRACT TABLES FROM WDC COMPLETE CORPUS 2015 DATASET
import json
import re
import os, sys, glob
import time

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

def save_table(table, target_dir):
    output_file = f'{target_dir}/t{table["id"]}_c{table["ncols"]}.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(table, f, ensure_ascii=False)
        f.close()
    print(f'Table {table["id"]} with {table["ncols"]} columns has been processed.')

def extract_tables(source_dir, target_dir):
    tableid = 1
    # iterate over json files in source directory
    for input_file in glob.glob(f'{source_dir}/*.json', recursive = True):
        print(f"Processing file: {input_file}")
        table = get_json_table(input_file)

        # skip vertical tables
        if table['tableOrientation'] == 'VERTICAL':
            continue

        else:
            # clean data
            relation = remove_header(table)
            relation = clean_relation(relation)

            # skip empty tables
            if len(relation) == 0:
                continue

            table = {'id': tableid, 'ncols':  len(relation), 'cols': relation}

            # save table to disk
            save_table(table, target_dir)
            tableid = tableid + 1

if len(sys.argv) >= 3:
    source_dir = sys.argv[1]
    target_dir = sys.argv[2] #directory where all json files will be stored

    # create target directory if doesn't exist
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
        print("Output directory " , target_dir ,  " created\n\n\n")
    else:    
        print("Output directory " , target_dir ,  " already exists\n\n\n")

    start = time.time()
    extract_tables(source_dir, target_dir)
    end = time.time()

    print(f'Success! All files have been procecessed.\nClean tables are stored in {target_dir}\n')
    print(f'Cleaning process took {end-start} seconds.')

else:
    print('Please provide source directory and target directory as args!')
    # ex: python3  clean-tables.py ../data/raw-files ../data/clean-tables
