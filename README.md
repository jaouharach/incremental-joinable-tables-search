# Kashif: Incremental Joinable Tables Search
This is the README file for Kashif a joinable tables discovery framework based on dstree C implementation.

The dstree is a Data-adaptive and Dynamic Segmentation  Index for Whole Matching on Time Series.

URL link for the DSTREE paper (The original implementation was in Java):\
> https://www.cs.sfu.ca/~jpei/publications/time%20series%20index%20VLDB13.pdf 


This C version was developed by: Karima Echihabi on 18/12/2016 and Jaouhara Chanchaf on 02/22/2023


Copyright 2016 Paris Descartes University. All rights reserved.\
Copyright 2023 Mohammed 6 Polytechnic University. All rights reserved.\
This code is based on the isax code written by Zoumpatianos on 3/12/12.\
Copyright 2012 University of Trento. All rights reserved.

> To compile and run, go to code directory: 

#### DO ONCE TO SET UP ENVIRONMENT
1. Run 'aclocal' to generate an m4 environment for autotools to use:\
`aclocal`

2. Then run autoconf to turn our configure.ac into a configure script:\
`autoconf`

3. Then run automake to turn our Makefile.am into a Makefile.in:\
`automake --add-missing`

4. Compile\
`./configure`\
`make`

5. Run\
`bin/dstree --help`

#### DATA PREPARATION
1. Download and extract tables from :
> https://data.dws.informatik.uni-mannheim.de/webtables/2015-07/englishCorpus/compressed/

2. Clean data: \
To remove numerical data go to data-preparation/clean and run:\
`python3 clean.py [directory_of_raw_tables] [directory_for_clean_tables]`

3. Generate embeddings: \
    i. download fasttext binary file from:\
        > https://fasttext.cc/docs/en/crawl-vectors.html \
    ii. then go to data-preparation/embed/fasttext.py and update the following line with the correct path to fasttext file, if you want to normalize vectors set normalize_vectors to True:\
        `FASTTEXT_PATH = '/home/jchanchaf/local/src/fasttext-en/cc.en.300.bin`\
        `normalize_vectors = False` \
    iii. Finally, run:\
        `fasttext_embed.py [directory_of_clean_tables] [directory_to_store_embbedings] [path_to_fasttext_binary] [embeddings_length]`\

4. Generate Binary files:\
Go to data-preparation/bin and run: \
`tobin.py [directory_of_table_embeddings] [directory_to_store_binary_files] [embeddings_length]`

#### EXPERIMENTS
1. Build the index:\
`bin/dstree --dataset /data/real/jchanchaf/wdc-2015-en-full/subsets/100-nonorm/ --total-data-files 100 --buffer-size 600 --index-path /home/jaouhara.chanchaf/work-dir/indexes/100-nonorm-idx/ --ascii-input 0 --mode 0 --track-bsf   --incremental --delta 1 --epsilon 0 --k 100 --timeseries-size 50 --track-vector 1 --warping 0.00 --leaf-size 100`

2. Query the index:\
`bin/dstree --queries /data/real/jchanchaf/wdc-2015-en-full/10m/query/10q-size100/ --total-data-files 100000 --buffer-size 600 --index-path /home/jaouhara.chanchaf/work-dir/indexes/100k-nonorm-idx/ --ascii-input 0 --mode 1 --track-bsf   --incremental --delta 1 --epsilon 0 --timeseries-size 50 --track-vector 1 --warping 0.00 --leaf-size 100000  --top 10 --dataset /data/real/jchanchaf/wdc-2015-en-full/subsets/100-nonorm/ --result-dir /home/jaouhara.chanchaf/work-dir/exp-results/kashif-search 100tk/stop-mode-0/ --ground-truth-dir ./na/  --nq 1 --min-qset-size 100 --max-qset-size 100  --parallel --incremental --store-results-in-disk 1 --stop-when-nn-dist-changes 0 --k  10 --k-values "5,10" --knn-data-structure minmax-heap`

3. Parameters:\
* `--dataset [string]` : directory where data files are stored.
* `--queries [string]` : directory where query files are stored.
* `--mode [0/1]` : 0 for index building , 1 for querying.
* `--track-vector [0/1]` : keep track of nearest neigbor ids.
* `--top [x]` : nb of matching columns to be retrieved.
* `--timeseries-size [x]` : length of the query vectors.
* `--nq [x]` : number of query column to execute.
* `--max-qset-size [x]` : Maximum query column size.
* `--min-qset-size [x]` : Minumum query column size.
* `--num-threads [x]` : number of worker threads. 
* `--k [x]`: (kmax) number of nearest neigbors to retrieve for each query vector in the query column.
* `--k-values [string]` : intermediate k values for which we would like to measure kashif performance (query time, recall) (only works for incremental search).

* `--knn-data-structure [string]` : data structure used to store nearest neighbors. Options: sorted-arr, minmax-heap, ostree.
* `--stop-when-nn-dist-changes [0/1]` : only return nns of the same distance as the first nn.
* `--result-dir [string]` : directory where query results will be stored
* `--ground-truth-dir [string]`: directory where ground truth results are stored (necessary when measuring recall)
* `--store-results-in-disk [0/1]` : set it to 0 is you don't want query results to be stored in disk

- The parameter `--buffer-size [x]`is in MBs so in this example, the buffer is ~ 0.6 GB.


