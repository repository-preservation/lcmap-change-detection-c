#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
//#include "envi_header.h"
#include "const.h"

/* Input file type definition */
typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_MAX
} Input_type_t;

/* Structure for the metadata */
typedef struct {
    int lines;            /* number of lines in a scene */ 
    int samples;          /* number of samples in a scene */
    int data_type;        /* envi data type */
    int byte_order;       /* envi byte order */
    int utm_zone;         /* UTM zone; use a negative number if this is a
                             southern zone */
    int pixel_size;       /* pixel size */
    char interleave[MAX_STR_LEN];  /* envi save format */ 
    int  upper_left_x;    /* upper left x coordinates */                         
    int  upper_left_y;    /* upper left y coordinates */ 
} Input_meta_t;

/* Structure for the 'input' data type */
typedef struct {
  Input_type_t file_type;  /* Type of the input image files */
  Input_meta_t meta;       /* Input metadata */
  FILE *fp_bin[TOTAL_BANDS][MAX_SCENE_LIST];           
} Input_t;

/* Prototypes */
FILE *open_raw_binary
(
    char *infile,        /* I: name of the input file to be opened */
    char *access_type    /* I: string for the access type for reading the
                               input file */
);

void close_raw_binary
(
    FILE *fptr      /* I: pointer to raw binary file to be closed */
);

int write_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to write to the file */
    int nsamps,         /* I: number of samples to write to the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    void *img_array     /* I: array of nlines * nsamps * size to be written
                              to the raw binary file */
);

int read_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to read from the file */
    int nsamps,         /* I: number of samples to read from the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    void *img_array     /* O: array of nlines * nsamps * size to be read from
                              the raw binary file (sufficient space should
                              already have been allocated) */
);

int read_envi_header
(
    char *filename,     /* I: header filename with full path */
    Input_meta_t *meta  /* O: saved header file info */
);


void usage ();

#endif
