#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct point{
    char label;
    float coord [2];
} point;

struct pairwise_distance{
    struct point *src;
    struct point *dest;
    float d;
} pairwise_distance;

float euclidean_distance(struct point *a, struct point *b, int len)
{
    float d = 0;
    for(int i = 0; i < len; i++)
    {
        d += ((a->coord[i] - b->coord[i]) * (a->coord[i] - b->coord[i]));
    }
    return sqrt(d);

}

void save_to_csv(struct pairwise_distance * pds, int num_distances)
{
    FILE *fp;
    int i;
    char * filename = "pdist.csv";
    fp = fopen(filename, "w+");
    
    if (fp == NULL) {
        fprintf(stderr, "Error in dstree_file_loaders.c: Could not open file %s!\n",
                filename);
        exit(1);
    }
    // write header
    fprintf(fp, "src, dest, d");

    for(int i = 0; i < num_distances; i++)
    {
        fprintf(fp, "\n");
        fprintf(fp, "%c, %c, %f", pds[i].src->label, pds[i].dest->label, pds[i].d);
    }
    fclose(fp);
}
int main()
{
    int num_points = 11;
    int len = 2;

    struct point a = {'A', {3.2, 1.7}};
    struct point b = {'B', {0.1, 2.3}};
    struct point c = {'C', {1.3, 4.7}};
    struct point d = {'D', {0.01, 2.2}};
    struct point e = {'E', {4.2, 1.2}};
    struct point f = {'F', {0.2684, 0.1523}};
    struct point g = {'G', {0.7365, 0.502}};
    struct point h = {'H', {0.3942, 0.9667}};
    struct point i = {'I', {0.554, 0.25}};
    struct point j = {'J', {0.819, 0.618}};
    struct point k = {'K', {0.450, 0.23}};

    struct point dataset [11] = {a, b, c, d, e, f, g, h, i, j, k};

    for(int i = 0; i < num_points; i++)
    {
        printf("point %c  (%f, %f)\n", dataset[i].label, dataset[i].coord[0], dataset[i].coord[1]);
    }
    
    int num_distances = num_points * (num_points - 1) / 2;

    struct pairwise_distance * pds = malloc(sizeof(struct pairwise_distance) * num_distances);

    for(int i = 0, s = 1, dist = 0; dist < num_distances; i++, s++)
    {
        struct point * src = &dataset[i];
        for(int j = s; j < num_points; j++)
        {
            struct point * dest = &dataset[j];
            float d = euclidean_distance(src, dest, len);
            pds[dist].src = src;
            pds[dist].dest = dest;
            pds[dist].d = d;
            dist++;
        }       
    }

    printf("pairwise distance:\n");
    for(int dist = 0; dist < num_distances; dist++)
    {
        printf("src: %c, dest: %c , d = %f\n", pds[dist].src->label, pds[dist].dest->label, pds[dist].d);
    }
    save_to_csv(pds, num_distances);

    free(pds);
    return 0;
}


