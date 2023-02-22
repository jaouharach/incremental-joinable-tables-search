
# Clean WDC webtables corpus

  

Scripts to extract and clean WDC web tables corpus, the following script will:

- Extract tables (horizontal tables only)
- Remove duplicates in the same column
- Remove numerical columns
- Remove empty tables
- Associate an incremental id to each table
- Store all tables (in JSON format) in one big file (you may want to change that and store each table in a separate JSON file especially if the number of tables to be cleaned is very large the output file may not fit into memory!

> **Requirments:** Python 3.x, libraries: os, sys, re, json, jsonlines, 

  

## Step 1: Download raw data

Download raw data, Example:
Create directory for raw data:

`mkdir $HOME/work-dir $HOME/work-dir/wdc2015-raw-data`
`cd "$HOME/work-dir/wdc2015-raw-data"`

Download and extract tables:
`wget http://data.dws.informatik.uni-mannheim.de/webtables/2015-07/englishCorpus/compressed/00.tar.gz && gunzip -k 00.tar.gz && tar -xf 00.tar`

## Step 2: Clean data
After extraction of raw data, go to the raw data directory, you'll notice that raw data files are grouped in a folder (or multiple folders) called 0, 1, 2 ...

Make a file called **files.txt** in the raw data folder that will contains names of all files in the set.
`(cd $HOME/work-dir/wdc2015-raw-data/`**`0`**`/ && ls > files.txt)` 

Check that the file contains raw files names
`(cd $HOME/work-dir/wdc2015-raw-data/`**`0`**`/ && head -5 files.txt)`

Run the bash script clean.sh using the name of the folder where the raw files are stored e.g. `$HOME/work-dir/wdc2015-raw-data/`**`0`**`/` 


Call the clean script using the following command: 
`bash clean.sh [raw_data_dir] [target_dir]`
  
>`[src_dir]`:  Directory where raw web table files are stored, if you down load multiple sets run this script for each raw web table files directory.

>`[target_dir]`: Directory where the clean data will be stored (you don't need to create the target directory beforehand)

Example: 
`bash clean.sh $HOME/work-dir/wdc2015-raw-data/0/ $HOME/work-dir/wdc2015-clean-data/`

After execution, go to `[target_dir]` you should see a file called **"clean_tables"** similar to the **"clean_tables_sample"** file, where each JSON object is a table, and each table is represented with a unique identifier "id", the total number of columns "ncols" and a list of columns (list of lists) "cols".