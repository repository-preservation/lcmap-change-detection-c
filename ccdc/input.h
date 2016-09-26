#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>
#include "const.h"

/* All defines in one file to avoid conflicts */
#include "defines.h"

/* Input file type definition */
typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_MAX
} Input_type_t;

/* Structure for the metadata */
typedef struct {
    int lines;          /* number of lines in a scene                         */ 
    int samples;        /* number of samples in a scene                       */
    int data_type;      /* envi data type                                     */
    int byte_order;     /* envi byte order                                    */
    int utm_zone;       /* UTM zone; use a negative number if this is a
                             southern zone                                    */
    int pixel_size;     /* pixel size                                         */
    char interleave[MAX_STR_LEN];  /* envi save format                        */ 
    int  upper_left_x;  /* upper left x coordinates                           */                         
    int  upper_left_y;  /* upper left y coordinates                           */ 
} Input_meta_t;

/* Structure for the 'input' data type */
typedef struct {
  Input_type_t file_type;  /* Type of the input image files */
  Input_meta_t meta;       /* Input metadata */
  FILE *fp_bin[TOTAL_BANDS][MAX_SCENE_LIST]; /* file pointer for image files  */
} Input_t;

/* Prototypes */
FILE *open_raw_binary
(
    char *infile,       /* I: name of the input file to be opened             */
    char *access_type   /* I: string for the access type for reading the
                               input file                                     */
);

void close_raw_binary
(
    FILE *fptr          /* I: pointer to raw binary file to be closed         */
);

int write_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file                  */
    int nlines,         /* I: number of lines to write to the file            */
    int nsamps,         /* I: number of samples to write to the file          */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8))   */
    void *img_array     /* I: array of nlines * nsamps * size to be written
                              to the raw binary file                          */
);

int read_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file                  */
    int nlines,         /* I: number of lines to read from the file           */
    int nsamps,         /* I: number of samples to read from the file         */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8))   */
    void *img_array     /* O: array of nlines * nsamps * size to be read from
                              the raw binary file (sufficient space should
                              already have been allocated)                    */
);

int read_envi_header
(
    char *data_type,    /* I: input data type                                 */
    char *scene_name,   /* I: scene name                                      */
    Input_meta_t *meta  /* O: saved header file info                          */
);

int read_cfmask
(
    int  curr_scene_num, /* I:   current num. in list of scenes to read       */
    char *data_type,     /* I:   type of flies, tifs or single BIP            */
    char **scene_list,   /* I:   current scene name in list of sceneIDs       */
    int  row,            /* I:   the row (Y) location within img/grid         */
    int  col,            /* I:   the col (X) location within img/grid         */
    int  nrows,          /* I:   number of rows to read in a tile of some data*/
    int  ncols,          /* I:   number of cols to read in a line of bip data */
    int  num_samples,    /* I:   number of image samples (X width)            */
    FILE ***fp_tifs,     /* I/O: file ptr array for tif band file names       */
    FILE **fp_bip,       /* I/O: file pointer array for BIP file names        */
    unsigned char *fmask_buf,/* O:   pointer to cfmask band values            */
                         /* I/O: Worldwide Reference System path and row for  */
                         /* I/O: the current swath, this group of variables   */
                         /* I/O: is for filtering out swath overlap, and      */
    int *prev_wrs_path,  /* I/O: using the first of two scenes in a swath,    */
    int *prev_wrs_row,   /* I/O: , because it is recommended to use the meta  */
    int *prev_year,      /* I/O: data from the first for things like sun      */
    int *prev_jday,      /* I/O: angle, etc. However, always removing a       */
    unsigned char *prev_fmask_buf,/* I/O: redundant x/y location specified    */
    int *valid_scene_count,/* I/O: x/y is not always valid for gridded data,  */
    int *swath_overlap_count,/* I/O: it may/may not be in swap overlap area.  */
    char **valid_scene_list,/* I/O: 2-D array for list of filtered            */
    int *clear_sum,      /* I/O: Total number of clear cfmask pixels          */
    int *water_sum,      /* I/O: counter for cfmask water pixels.             */
    int *shadow_sum,     /* I/O: counter for cfmask shadow pixels.            */
    int *sn_sum,         /* I/O: Total number of snow cfmask pixels           */
    int *cloud_sum,      /* I/O: counter for cfmask cloud pixels.             */
    int *fill_sum,       /* I/O: counter for cfmask fill pixels.              */
    int *all_sum,        /* I/O: Total of all cfmask pixels                   */
    unsigned char *updated_fmask_buf, /* I/O: new entry in valid fmask values */
    int *updated_sdate_array, /* I/O: new buf of valid date values            */
    int *sdate,          /* I:   Original array of julian date values         */
    int *valid_num_scenes/* I/O: number of valid scenes after reading cfmask  */
);


int read_stdin
(
    int           *updated_sdate_array, /* O:   pointer to date values buffer */
    int           *buf,                 /* O:   pointer to image bands buffer */
    unsigned char *updated_cfmask_buf,  /* O:   pointer to cfmask pixel buffer*/
    int           num_bands,            /* I:   total number of bands         */
    int           *clear_sum,           /* O:   accumulator for clear  pixels */
    int           *water_sum,           /* O:   accumulator for water  pixels */
    int           *shadow_sum,          /* O:   accumulator for shadow pixels */
    int           *snow_sum,            /* O:   accumulator for snow   pixels */
    int           *cloud_sum,           /* O:   accumulator for cloud  pixels */
    int           *fill_sum,            /* O:   accumulator for fill   pixels */
    int           *all_sum,             /* O:   accumulator for all    pixels */
    int           *valid_num_scenes,    /* O:   total scenes read             */
    bool          debug                 /* I:   flag for printing mesgs/info  */
);


int assign_cfmask_values
(
    unsigned char cfmask_value,         /* I: current cfmask pixel value.     */
    int           *clear_sum,           /* O: accumulator for clear  pixels   */
    int           *water_sum,           /* O: accumulator for water  pixels   */
    int           *shadow_sum,          /* O: accumulator for shadow pixels   */
    int           *snow_sum,            /* O: accumulator for snow   pixels   */
    int           *cloud_sum,           /* O: accumulator for cloud  pixels   */
    int           *fill_sum,            /* O: accumulator for full   pixels   */
    int           *all_sum              /* O: accumulator for all    pixels   */
);


int read_tifs
(
    char *sceneID_name,  /* I:   current file name in list of sceneIDs        */
    FILE ***fp_tifs,     /* I/O: file pointer array for band file names       */
    int  curr_scene_num, /* I:   current num. in list of scenes to read       */
    int  row,            /* I:   the row (Y) location within img/grid         */
    int  col,            /* I:   the col (X) location within img/grid         */
    int  nrows,          /* I:   number of rows to read in tile of data */
    int  ncols,          /* I:   number of colums to read (1 for bip)   */
    int  num_samples,    /* I:   number of image samples (X width)            */
    bool debug,          /* I:   flag for printing debug messages             */
    int  *image_buf      /* O:   pointer to 2-D image band values array       */
);


int read_bip
(
    char *current_scene_name, /* I:   current file name in list of sceneIDs  */
    FILE **fp_bip,            /* I/O: file pointer array for BIP  file names */
    int  curr_scene_num,      /* I:   current num. in list of scenes to read */
    int  row,                 /* I:   the row (Y) location within img/grid   */
    int  col,                 /* I:   the col (X) location within img/grid   */
    int  nrows,               /* I:   number of rows to read in tile of data */
    int  ncols,               /* I:   number of colums to read (1 for bip)   */
    int  num_samples,         /* I:   number of image samples (X width)      */
    int  *image_buf           /* O:   pointer to 2-D image band values array */
);


void usage ();

#endif
