#ifndef INPUT_H
#define INPUT_H

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "date.h"
//#include "envi_header.h"
#include "input.h"
#include "const.h"

#define INPUT_FILL (0)
#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)
#define NBAND_THM_MAX 2
#define CCDC_VERSION "1.0.0"
#define MINSIGMA 1e-5

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
    char *scene_name,      /* I: scene name*/
    Input_meta_t *meta    /* O: saved header file info */
);

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
    int *row,                 /* O: row number for the pixel */
    int *col,                 /* O: col number for the pixel */
    float *t_cg,              /* O: chi-square inversed T_cg */
    bool *verbose             /* O: verbose flag */
);

int create_scene_list
(
    const char *item,         /* I: string of file items be found */
    int num_scenes,           /* I/O: number of scenes */
    char **scene_list         /* O: scene_list used for ccdc processing */ 
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

int median_variogram
(
    int16 **array,      /* I: input array */
    int dim1_start,     /* I: dimension 1 start index */
    int dim1_end,       /* I: dimension 1 end index */
    int dim2_len,       /* I: dimension 2 length in input array */
    float *output_array /* O: output array */
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

void square_root_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements */
    float **fit_ctf,     /* I: */
    float  *rmse         /* O: output rmse value */
);

void partial_square_root_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
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

void matlab_2d_int_array_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float  *output_mean  /* O: output norm value */
);

void matlab_2d_partial_mean
(
    float **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
);

void matlab_int_2d_partial_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
);

void matlab_2d_partial_square_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
);

void matlab_2d_array_norm
(
    float **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float  *output_norm  /* O: output norm value */
);

void get_ids_length
(
    int *id_array,        /* I: input array */
    int start,            /* I: array start index */
    int end,              /* I: array end index */
    int *id_len           /* O: number of non-zero number in the array */
);

void quick_sort_index (float arr[], int idx[], int left, int right);

int auto_mask
(
    int *clrx,
    int16 **clry,
    int start,
    int end,
    float years,
    int t_b1,
    int t_b2,
    int n_t,
    int *bl_ids
);

int auto_ts_fit
(
    int *clrx,
    int16 **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse
);

int auto_ts_fit_full
(
    int *clrx,
    int16 **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse,
    float **v_dif
);

void auto_ts_predict
(
    int *clrx,
    float **coefs,
    int band_index,
    int start,
    int end,
    float *pred_y
);

void usage ();

#endif
