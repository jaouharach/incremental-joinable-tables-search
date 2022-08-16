/*SCRIPT TO CONVERT ONE BINARY FILE TO TABLE, SET AND VECTOR STRUCTURES*/
#include <dirent.h>
#include <libgen.h>
#include <linux/limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#define NBITS 64
typedef float ts_type;

typedef float ts_type;

typedef struct vector {
  int table_id;
  int set_id;
  ts_type *values;
} vector;

/* convert IEEE double precision (64 bits) to float (/double) */
ts_type todecimal(char *s) {
  // printf("\ninput: %s\n", s);
  ts_type f;
  int sign, exp = 0;
  unsigned long long mant = 0;
  int i;

  sign = s[0] - '0';

  for (i = 1; i <= 11; i++)
    exp = exp * 2 + (s[i] - '0');

  if (exp > -1022) {
    mant = 1;
    exp -= 1023;
  } else {
    // Subnormal numbers
    mant = 0;
    exp = -1022;
  }
  for (i = 12; i < 64; i++) {
    mant = mant * 2 + (s[i] - '0');
  }
  /*

  printf("s = %d\n", sign);
  printf("exp = %d\n", exp);
  printf("m = %llu\n", mant);

  */

  f = pow(-1, sign) * pow(2, exp) * (mant / pow(2, 52));
  return f;
}

/* CONVERT BIN TO INT */
int toint(char *n) {
  int integer = 0;
  for (int i = 0; i < strlen(n); i++)
    integer = integer * 2 + (n[i] - '0');
  return integer;
}

bool is_binaryfile(const char *filename) {
  // check if filename has bin extesion.
  char *ext = ".bin";
  size_t nl = strlen(filename), el = strlen(ext);
  return nl >= el && !strcmp(filename + nl - el, ext);
}

// extract set bytes from bin file
char *get_set_bytes(FILE *f, int nvec, int start, int nbits, int vlen) {
  start += nbits;                  // skip bytes for num vectors
  int total = nvec * vlen * nbits; // total bits occupied by set values

  char *bytes = (char *)malloc(sizeof(char) * total);
  fseek(f, start, SEEK_SET);
  int i;
  for (i = 0; i < total; i++) {
    bytes[i] = fgetc(f);
  }
  bytes[i] = '\0';
  // //printf("\n%s\n", bytes);
  return bytes;
}

// get number of vectors in the current set (set that starts at 'start')
int get_num_vec(FILE *f, int start, int nbits) {
  char *bytes = (char *)malloc(sizeof(char) * nbits);
  int i;
  fseek(f, start, SEEK_SET);
  for (i = 0; i < nbits; i++) {
    bytes[i] = fgetc(f);
  }
  bytes[i] = '\0';
  // //printf("\n%s\n", bytes);
  return toint(bytes);
}

// read datasize from diffrent files and get total datasize of a datalake
unsigned int get_total_data_vectors(char *bindir, unsigned int total_data_files) {
  struct dirent *dfile;
  DIR *dir = opendir(bindir);
  unsigned int total_datasize = 0;
  unsigned int datasize;
  if (!dir) {
    printf("Unable to open directory stream!");
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && total_data_files > 0) {
    // // skip directories
    // if (dfile->d_type != DT_REG)
    //   continue;
    if (is_binaryfile(dfile->d_name)) {
      total_data_files--;
      // ** get binary table info
      sscanf(dfile->d_name, "data_size%d%*[^0123456789]", &datasize);
      total_datasize += datasize;
    }
  }
  closedir(dir);
  return total_datasize;
}

// get data lake size in GB
unsigned int get_dlsize(char *dl_dir, unsigned int l) {
  struct dirent *dfile;
  DIR *dir = opendir(dl_dir);
  float total_dlsize = 0.0;
  FILE *fp;

  if (!dir) {
    printf("Unable to open directory stream!");
    exit(1);
  }

  while ((dfile = readdir(dir)) != NULL && l > 0) {
    // // skip directories
    // if (dfile->d_type != DT_REG)
    //   continue;

    if (is_binaryfile(dfile->d_name))
    {
      char bin_file_path[PATH_MAX + 1] = "";
      strcat(bin_file_path, dl_dir);
      strcat(bin_file_path, "/");
      strcat(bin_file_path, dfile->d_name);
      l--;
      // ** get binary file size
      fp = fopen(bin_file_path, "r");
      if (fp == NULL) {
        printf("Could not get size of binary file %s. File not found!\n",
               dfile->d_name);
        exit(1);
      }
      fseek(fp, 0L, SEEK_END);
      total_dlsize += ftell(fp);
      fclose(fp);
    }
  }
  closedir(dir);
  return (unsigned int)round(total_dlsize / 1073741824);
}