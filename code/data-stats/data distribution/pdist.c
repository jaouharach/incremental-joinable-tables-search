#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <pthread.h>
     
int num_threads = 32; 
char *dataset_dir = "/home/jaouhara.chanchaf/work-dir/data/wdc-2015/clean-tables-beta2/bins/";
unsigned int total_files = 100000;
unsigned long long total_vectors = 5027911;
unsigned int vector_length = 50;
char * output_file = "/data/real/jchanchaf/wdc-2015-en-full/pairwise_distances_100k.csv";


// Total files: 500000
// Total columns 2474940
// Total  vectors:	25279858

// Total files: 100000
// Total columns 494672
// Total  vectors:	5027911

// Total files: 10
// Total columns 37
// Total  vectors:	833
// Size in GB:	0

struct vector { 
  unsigned int table_id;
  unsigned int set_id;
  unsigned int pos;
  float * coords;
} vector;

struct pair_dist { 
  unsigned int src_table_id;
  unsigned int src_set_id;
  unsigned int src_pos;

  unsigned int dest_table_id;
  unsigned int dest_set_id;
  unsigned int dest_pos;

  float  dist;
} pair_dist;

struct vector *read_all_vectors(char *dataset_dir, unsigned int total_files, unsigned long long total_vectors, 
                                unsigned int vector_length, unsigned long long * total_read_vectors);
bool is_binaryfile(const char *filename);

void save_to_csv(struct vector *src, struct vector * dest, float dist, FILE* out);

float euclidean_distance(struct vector *a, struct vector *b, int len);

int main()
{
    clock_t start, start_full, end;
    double time_used, cpu_time_used;
     
    start = clock();
    start_full = clock();
    unsigned long long total_read_vectors = 0;
    struct vector * dataset_vectors = read_all_vectors(dataset_dir, total_files, total_vectors, 
                                                        vector_length, &total_read_vectors);

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Finished reading vectors, process took %.2f seconds.\n", cpu_time_used);

    unsigned long long num_distances = (unsigned long long) ((total_read_vectors * (total_read_vectors - 1)) / 2);
    unsigned long long total_dist_computations = 0;
    printf("Total vectors = %lld, computing %lld distances...\n", total_read_vectors, num_distances);

    // struct pair_dist * pairwise_distances = malloc(sizeof(struct pair_dist)* num_distances);
    double * pairwise_distances = malloc(sizeof(double)* num_distances);

    printf("1 pair requires %ld bytes\n", sizeof(struct pair_dist));

    if(pairwise_distances == NULL)
    {
        fprintf(stderr, "Error: couldn't allocate memory for dist array.\n");
        exit(1);
    }

    FILE *fp;
    fp = fopen(output_file, "wb");
    if (fp == NULL) 
    {
        fprintf(stderr, "Error: Could not open file %s!\n",
                output_file);
        exit(1);
    }
    // write header
    // fprintf(fp, "src, dest, d");

    start = clock();
    for(unsigned long long i = 0, s = 1, dist = 0; dist < num_distances; i++, s++)
    {
        struct vector * src = &dataset_vectors[i];
        for(unsigned long long j = s; j < total_read_vectors; j++)
        {
            struct vector * dest = &dataset_vectors[j];
            printf("[%lld/%lld] (%u, %u, %u) <-- d --> (%u, %u, %u)\n ", i+1, total_read_vectors,
                    src->table_id, src->set_id, src->pos, dest->table_id, dest->set_id, dest->pos);
            
            // pairwise_distances[dist].src_table_id = src->table_id;
            // pairwise_distances[dist].src_set_id = src->set_id;
            // pairwise_distances[dist].src_pos = src->pos;
            // pairwise_distances[dist].dest_table_id = dest->table_id;
            // pairwise_distances[dist].dest_set_id = dest->set_id;
            // pairwise_distances[dist].dest_pos = dest->pos;
            // pairwise_distances[dist].dist = euclidean_distance(src, dest, vector_length);
            
            pairwise_distances[dist] = (double) euclidean_distance(src, dest, vector_length);

            dist++;
            total_dist_computations++; 
        }       
    }
    fwrite(pairwise_distances, sizeof(*pairwise_distances), num_distances, fp);

    time_used = ((double) (clock() - start_full)) / CLOCKS_PER_SEC;
    
    printf("Total time: process took %.2f seconds.\n", time_used);
    printf("Finished computing distances, process performed %lld dist calc in %.2f seconds.\n", total_dist_computations, time_used);

    fclose(fp);
    free(dataset_vectors);
    return 0;
}

bool is_binaryfile(const char *filename) 
{
  // check if filename has bin extesion.
  char *ext = ".bin";
  size_t nl = strlen(filename), el = strlen(ext);
  return nl >= el && !strcmp(filename + nl - el, ext);
}

struct vector *read_all_vectors(char *dataset_dir, unsigned int total_files, unsigned long long total_vectors, 
                                unsigned int vector_length, unsigned long long *total_read_vectors)
{
    unsigned long long read_so_far = 0;
    struct vector *dataset_vectors = malloc(sizeof(struct vector) * total_vectors);
    if(dataset_vectors == NULL)
    {
        fprintf(stderr, "Error: Unable to allocate memory for all dataset vectors.");
        exit(1);
    }
    for(int i = 0; i < total_vectors; i++)
    {
        dataset_vectors[i].coords = malloc(sizeof(float) * vector_length);
        if(dataset_vectors[i].coords == NULL)
        {
            fprintf(stderr, "Error: Unable to allocate memory for all dataset vectors.");
            exit(1);
        }
    }
    struct vector temp; 

    // open source dir
    struct dirent *dfile;
    DIR *dir = opendir(dataset_dir);
    if (!dir) {
        fprintf(stderr, "Error: Unable to open directory stream! %s", dataset_dir);
        exit(1);
    }

    while ((dfile = readdir(dir)) != NULL && total_files > 0 && read_so_far < total_vectors)
    {
        if (is_binaryfile(dfile->d_name)) 
        {
            total_files--;
            // printf("file %d\n", total_files);
            // get full path of bin file
            char bin_file_path[PATH_MAX + 1] = "";
            strcat(bin_file_path, dataset_dir);
            strcat(bin_file_path, "/");
            strcat(bin_file_path, dfile->d_name);

            // get binary table info
            int datasize, table_id, nsets, vector_length_in_filename;
            sscanf(dfile->d_name, "data_size%d_t%dc%d_len%d_noznorm.bin", &datasize,
                    &table_id, &nsets, &vector_length_in_filename);

            // check if vector length in file name matches vector length passed as
            // argument
            if (vector_length_in_filename != vector_length) 
            {
                fprintf(stderr,
                        "Error: Vector length %d does not match vector length in file %s.\n", vector_length, bin_file_path);
                exit(1);
            }
            
            FILE *bin_file = fopen(bin_file_path, "rb");
            if (bin_file == NULL) {
                fprintf(stderr, "Error: Could not open file '%s' for this reason %s\n", bin_file_path, strerror(errno));
                exit(1);
            }

            /* Start processing file: read every vector in binary file*/
            unsigned int i = 0, j = 0, set_id = 0, total_bytes = (datasize * vector_length) + nsets;
            unsigned int nvec = 0u;
            float val;

            // counts 4 bytes as one because every 
            // vector coordinate is stored in 4 bytes
            while (total_bytes) 
            {
                // read a new column
                if (i == 0)
                {
                    fread(&nvec, sizeof(nvec), 1, bin_file);
                    // printf("|column| = %d\n", nvec);
                    total_bytes--;

                    // set new set id
                    dataset_vectors[read_so_far].table_id = table_id;
                    dataset_vectors[read_so_far].set_id = set_id;
                    dataset_vectors[read_so_far].pos = 0;

                    // printf("\n+ vec (%d, %d, %d) ...\n", table_id, set_id, 0);

                    set_id += 1;
                    i++;
                    j = 0;
                }
                // still read curr column
                else if (i <= (unsigned int)nvec * vector_length) 
                {
                    // end of vector but still in current column
                    if (j > (vector_length - 1)) {
                        // printf("(");
                        // for(int c = 0; c < vector_length; c++)
                        // {
                        //     printf("%.3f, ", dataset_vectors[read_so_far].coords[c]);
                        // }
                        // printf(")\n");

                        j = 0;
                        read_so_far++;

                        if(read_so_far >= total_vectors)
                        {
                            printf("Done.\n");
                            break; // we have read all vectors
                        }

                        dataset_vectors[read_so_far].table_id = dataset_vectors[read_so_far - 1].table_id;
                        dataset_vectors[read_so_far].set_id = dataset_vectors[read_so_far - 1].set_id;
                        dataset_vectors[read_so_far].pos = dataset_vectors[read_so_far - 1].pos + 1;

                        // printf("\n+ vec (%d, %d, %d) ...\n", dataset_vectors[read_so_far].table_id, dataset_vectors[read_so_far].set_id, dataset_vectors[read_so_far].pos);
                    }

                    fread((void *)(&val), sizeof(val), 1, bin_file);
                    total_bytes--;
                    dataset_vectors[read_so_far].coords[j] = val;

                    // end of last vector in current column
                    if (i == (unsigned int)nvec * vector_length)
                    {
                        // printf("(");
                        // for(int c = 0; c < vector_length; c++)
                        // {
                        //     printf("%.3f, ", dataset_vectors[read_so_far].coords[c]);
                        // }
                        // printf(")\n");

                        dataset_vectors[read_so_far].pos = 0;
                        i = 0;
                        j = 0;
                        nvec = 0u;
                        continue;
                    }
                    i++;
                    j++;
                }
            }
            fclose(bin_file);
        }
    }

    *total_read_vectors = read_so_far; 
    printf("Total read vector  = %lld/%lld\n", read_so_far, total_vectors);
    return dataset_vectors;
}

float euclidean_distance(struct vector *a, struct vector *b, int len)
{
    float d = 0;
    for(int i = 0; i < len; i++)
    {
        d += ((a->coords[i] - b->coords[i]) * (a->coords[i] - b->coords[i]));
    }
    // return sqrt(d);
    return d;

}