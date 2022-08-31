> To compile code (from src folder): 
gcc utils.c bf.c -lm -w -g -o ../bin/bf

> To run code (from bf folder): 
bin/bf --dataset ../../data-samples/ --queries ../../data-samples/ --nq 1 --min-qset-size 2 --max-qset-size 10 --top 3 --result-dir ./results/ --total-data-files 50 --k 10 --vector-length 50

<!-- in server -->
bin/bf --dataset /home/jaouhara.chanchaf/work-dir/data/wdc-2015/clean-tables-beta2/bins/  --total-data-files 100000 --queries /home/jaouhara.chanchaf/work-dir/data/wdc-2015/clean-tables-beta2/query-columns[10]-size[5-10]/ --nq 10 --min-qset-size 5 --max-qset-size 10 --top 10 --result-dir /home/jaouhara.chanchaf/work-dir/exp-results/bf-search/ --k 10 --vector-length 50
