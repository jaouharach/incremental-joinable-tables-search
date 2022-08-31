> To compile code (from src folder): 
gcc utils.c bf.c -lm -w -g -o ../bin/bf

> To run code (from bf folder): 
bin/bf --dataset ../../data-samples/ --queries ../../data-samples/ --nq 1 --min-qset-size 2 --max-qset-size 10 --top 3 --result-dir ./results/ --total-data-files 50 --k 10 --vector-length 50