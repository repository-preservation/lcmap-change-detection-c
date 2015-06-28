#ifndef CCDC_H
#define CCDC_H


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>


#include "input.h"


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

void rmse_from_square_root_mean
(
    int16 **array,      /* I: input array */
    float fit_cft,      /* I: input fit_cft value */
    int dim1_len,       /* I: dimension 1 length */
    int dim2_index,     /* I: dimension 2 index in input array */
    float *rmse         /* O: output rmse */
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

void matlab_2d_array_median
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float *output_median /* O: output norm value */
);

void matlab_2d_int_array_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int din1_len,        /* I: number of input elements in 1st dim */
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

void rmse_from_square_root_mean
(
    int16 **array,      /* I: input array */
    float fit_cft,      /* I: input fit_cft value */
    int dim1_len,       /* I: dimension 1 length */
    int dim2_index,     /* I: dimension 2 index in input array */
    float *rmse         /* O: output rmse */
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

void matlab_2d_array_median
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float *output_median /* O: output norm value */
);

void matlab_2d_int_array_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int din1_len,        /* I: number of input elements in 1st dim */
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

#endif /* CCDC_H */
