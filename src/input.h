#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "date.h"
//#include "envi_header.h"
#include "input.h"

#define INPUT_FILL (0)
#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)
#define NBAND_THM_MAX 2
#define CCDC_VERSION "1.0.0"
#define MINSIGMA 1e-5
#define MAX_STR_LEN 510

typedef signed short int16;
typedef unsigned char uint8;

typedef enum
{
    HDF_FILE = 0,
    BINARY_FILE
}File_type;

/* Input file type definition */
typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_MAX
} Input_type_t;

/* Satellite type definition */
typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_1 = 0, 
  SAT_LANDSAT_2, 
  SAT_LANDSAT_3, 
  SAT_LANDSAT_4, 
  SAT_LANDSAT_5, 
  SAT_LANDSAT_7, 
  SAT_LANDSAT_8, 
  SAT_MAX
} Sat_t;

/* Instrument type definition */
typedef enum {
  INST_NULL = -1,
  INST_MSS = 0, 
  INST_TM,
  INST_ETM, 
  INST_OLI_TIRS, 
  INST_MAX
} Inst_t;

/* World Reference System (WRS) type definition */
typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

/* Band gain settings (ETM+ only) */
typedef enum {
  GAIN_NULL = -1,
  GAIN_HIGH = 0, 
  GAIN_LOW, 
  GAIN_MAX
} Gain_t;

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
                               input file; use the raw_binary_format array
                               at the top of this file */
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
    char *scene_name,      /* I: scene name*/
    Input_metad_t *meta    /* O: saved header file info */
)

void split_filename 
(
    const char *filename,       /* I: Name of file to split */
    char *directory,            /* O: Directory portion of file name */
    char *scene_name,           /* O: Scene name portion of the file name */
    char *extension             /* O: Extension portion of the file name */
);

char *sub_string
(
    const char *source,
    size_t start,
    size_t length
); 

int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    int row,                  /* I: row number for the pixel */
    int col,                  /* I: col number for the pixel */
    float *min_rmse,          /* I: minimum rmse threshold value */
    float *t_cg,              /* I: chi-square inversed T_cg */
    float *t_max_cg,          /* I: chi-square inversed T_max_cg for 
                                    last step noise removal */
    int * conse,              /* I: number of points used for change detection */ 
    bool *verbose             /* O: verbose flag */
);

int create_scene_list
(
    const char *item,         /* I: string of file items be found */
    char **scene_list,        /* O: scene_list used for ccdc processing */ 
    bool *verbose             /* O: verbose flag */
);

int convert_year_doy_to_jday_from_0000
(
    int year,      /* I: year */
    int doy,       /* I: day of the year */
    int *jday      /* O: julian date since year 0000 */
);

int convert_jday_from_0000_to_year_doy
(
    int jday,      /* I: julian date since year 0000 */
    int *year,     /* O: year */
    int *doy       /* O: day of the year */
);

int sort_scene_based_on_year_doy
(
    char **scene_list,      /* I/O: scene_list, sorted as output */
    int num_scenes,         /* I: number of scenes in the scene list */
    int *sdate              /* O: year plus date since 0000 */
);

void update_cft
(
    int i_span,
    int n_times,
    int min_num_c,
    int mid_num_c,
    int max_num_c,
    int num_c,
    int *update_number_c
);

void median_filter
(
    int16 *array,      /* I: input array */
    int array_len,     /* I: number of elements in input array */
    int n,             /* I: output order N, here is an odd number */
    int *output_array  /* O: output array */
);

void array_intersection
(
    int *array1,       /* I: input array 1 */
    int array_len1,    /* I: number of elements in input array1 */
    int *array2,       /* I: input array 2 */
    int array_len2,    /* I: number of elements in input array2 */
    int *output_array, /* O: output array */
    int *output_len    /* O: output array length */
);

void matlab_norm
(
    float *array,        /* I: input array */
    int array_len,       /* I: number of elements in input array */
    float  *output_norm  /* O: output norm value */
);

void matlab_2d_norm
(
    float *array,        /* I: input array */
    int array_dim1,      /* I: number of input elements in 1st dim */
    int array_dim2,      /* I: number of input elements in 2nd dim */
    float  *output_norm  /* O: output norm value */
);

void square_root_mean
(
    float **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements */
    float **fit_ctf,     /* I: */
    float  *rmse         /* O: output rmse value */
);

void matlab_2d_array_mean
(
    float **array,       /* I: input array */
    int dim1_number,     /* I: first dimension number used */   
    int array_len2,      /* I: number of input elements in 1st dim */
    float  *output_mean  /* O: output norm value */
);

void usage ();

#endif
