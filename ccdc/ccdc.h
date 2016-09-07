#ifndef CCDC_H
#define CCDC_H


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#include "input.h"

int get_args
(
    int argc,              /* I: number of cmd-line args                    */
    char *argv[],          /* I: string of cmd-line args                    */
    int *row,              /* O: row number for the pixel                   */
    int *col,              /* O: col number for the pixel                   */
    int *nrows,            /* O: number of Y pixels for the tile            */
    int *ncols,            /* O: number of X pixels for the tile            */
    char *in_path,         /* O: directory location of input data           */
    char *out_path,        /* O: direcotry location of output files         */
    char *data_type,       /* O: data type: tifs, bip, stdin, bip_lines.    */
    char *scene_list_file, /* O: optional file name of list of sceneIDs     */
    bool *verbose          /* O: verbose flag                               */
);

void get_scenename
(
    const char *filename, /* I: Name of file to split               */
    char *directory,      /* O: Directory portion of file name      */
    char *scene_name,     /* O: Scene name portion of the file name */
    char *appendix        /* O: Appendix portion of the file name   */
);

int create_scene_list
(
    const char *item,         /* I: string of file items be found          */
    int *num_scenes,          /* I/O: number of scenes                     */
    char *sceneListFileName   /* I: file name containing list of scene IDs */
);

int convert_year_doy_to_jday_from_0000
(
    int year,      /* I: year                        */
    int doy,       /* I: day of the year             */
    int *jday      /* O: julian date since year 0000 */
);

int sort_scene_based_on_year_doy_row
(
    char **scene_list,      /* I/O: scene_list, sorted as output     */
    int num_scenes,         /* I: number of scenes in the scene list */
    int *sdate              /* O: year plus date since 0000          */
);

void quick_sort_2d_float
(
    float arr[],
    float *brr[],
    int left,
    int right
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
    float **array,      /* I: input array                       */
    int dim1_len,       /* I: dimension 1 length in input array */
    int dim2_start,     /* I: dimension 2 start index           */
    int dim2_end,       /* I: dimension 2 end index             */
    float *output_array /* O: output array                      */
);

void split_directory_scenename
(
    const char *filename,       /* I: Name of scene with path to split    */
    char *directory,            /* O: Directory portion of file name      */
    char *scene_name            /* O: Scene name portion of the file name */
);

void rmse_from_square_root_mean
(
    float **array,      /* I: input array                      */
    float fit_cft,      /* I: input fit_cft value              */
    int dim1_index,     /* I: dimension 1 index in input array */
    int dim2_len,       /* I: dimension 2 length               */
    float *rmse         /* O: output rmse                      */
);

void partial_square_root_mean
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim   */
    float **fit_ctf,     /* I:                                     */
    float  *rmse         /* O: output rmse value                   */
);

void matlab_2d_array_mean
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */   
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_mean   /* O: output norm value                   */
);

void matlab_2d_float_median
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */   
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_median /* O: output norm value                   */
);

void matlab_2d_partial_mean
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */
    int start,           /* I: number of start elements in 2nd dim */
    int end,             /* I: number of end elements in 2nd dim   */
    float *output_mean   /* O: output norm value                   */
);

void matlab_float_2d_partial_median
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */
    int start,           /* I: number of start elements in 2nd dim */
    int end,             /* I: number of end elements in 2nd dim   */
    float *output_median /* O: output norm value                   */
);

void matlab_2d_partial_square_mean
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */   
    int start,           /* I: number of start elements in 2nd dim */
    int end,             /* I: number of end elements in 2nd dim   */
    float  *output_mean  /* O: output norm value                   */
);

void matlab_2d_array_norm
(
    float **array,       /* I: input array                         */
    int dim1_index,      /* I: 1st dimension index                 */
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float  *output_norm  /* O: output norm value                   */
);

void get_ids_length
(
    int *id_array,        /* I: input array */
    int start,            /* I: array start index */
    int end,              /* I: array end index */
    int *id_len           /* O: number of non-zero number in the array */
);

void matlab_unique
(
    int *clrx,
    float **clry,
    int nums,
    int *new_nums
);

int auto_mask
(
    int *clrx,
    float **clry,
    int start,
    int end,
    float years,
    float t_b1,
    float t_b2,
    float n_t,
    int *bl_ids
);

int auto_ts_fit
(
    int *clrx,
    float **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse,
    float **v_dif
);

int auto_ts_predict
(
    int *clrx,
    float **coefs,
    int df,
    int band_index,
    int start,
    int end,
    float *pred_y
);

extern void elnet_(
    
// input:

    int *ka,		//   ka = algorithm flag
			//      ka=1 => covariance updating algorithm
			//      ka=2 => naive algorithm
    double *parm,	//   parm = penalty member index (0 <= parm <= 1)
			//        = 0.0 => ridge
			//        = 1.0 => lasso
    int *no,		//   no = number of observations
    int *ni,		//   ni = number of predictor variables
    double *x,		//   x[ni][no] = predictor data matrix flat file (overwritten)
    double *y,		//   y[no] = response vector (overwritten)
    double *w,		//   w[no]= observation weights (overwritten)
    int *jd,		//   jd(jd(1)+1) = predictor variable deletion flag
			//      jd(1) = 0  => use all variables
			//      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
    double *vp,		//   vp(ni) = relative penalties for each predictor variable
			//      vp(j) = 0 => jth variable unpenalized
    double cl[][2],	//   cl(2,ni) = interval constraints on coefficient values (overwritten)
			//      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
			//      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
    int *ne,		//   ne = maximum number of variables allowed to enter largest model
			//        (stopping criterion)
    int *nx,		//   nx = maximum number of variables allowed to enter all models
			//        along path (memory allocation, nx > ne).
    int *nlam,		//   nlam = (maximum) number of lamda values
    double *flmin,	//   flmin = user control of lamda values (>=0)
			//      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
			//      flmin >= 1.0 => use supplied lamda values (see below)
    double *ulam,	//   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
    double *thr,	//   thr = convergence threshold for each lamda solution.
			//      iterations stop when the maximum reduction in the criterion value
			//      as a result of each parameter update over a single pass
			//      is less than thr times the null criterion value.
			//      (suggested value, thr=1.0e-5)
    int *isd,		//   isd = predictor variable standarization flag:
			//      isd = 0 => regression on original predictor variables
			//      isd = 1 => regression on standardized predictor variables
			//      Note: output solutions always reference original
			//            variables locations and scales.
    int *intr,		//   intr = intercept flag
			//      intr = 0/1 => don't/do include intercept in model
    int *maxit,		//   maxit = maximum allowed number of passes over the data for all lambda
			//      values (suggested values, maxit = 100000)

// output:

    int *lmu,		//   lmu = actual number of lamda values (solutions)
    double *a0,		//   a0(lmu) = intercept values for each solution
    double *ca,		//   ca(nx,lmu) = compressed coefficient values for each solution
    int *ia,		//   ia(nx) = pointers to compressed coefficients
    int *nin,		//   nin(lmu) = number of compressed coefficients for each solution
    double *rsq,	//   rsq(lmu) = R**2 values for each solution
    double *alm,	//   alm(lmu) = lamda values corresponding to each solution
    int *nlp,		//   nlp = actual number of passes over the data for all lamda values
    int *jerr		//   jerr = error flag:
			//      jerr  = 0 => no error
			//      jerr > 0 => fatal error - no output returned
			//         jerr < 7777 => memory allocation error
			//         jerr = 7777 => all used predictors have zero variance
			//         jerr = 10000 => maxval(vp) <= 0.0
			//      jerr < 0 => non fatal error - partial output:
			//         Solutions for larger lamdas (1:(k-1)) returned.
			//         jerr = -k => convergence for kth lamda value not reached
			//            after maxit (see above) iterations.
			//         jerr = -10000-k => number of non zero coefficients along path
			//            exceeds nx (see above) at kth lamda value.
);

extern void spelnet_(
    
// input:

    int *ka,		//   ka = algorithm flag
			//      ka=1 => covariance updating algorithm
			//      ka=2 => naive algorithm
    double *parm,	//   parm = penalty member index (0 <= parm <= 1)
			//        = 0.0 => ridge
			//        = 1.0 => lasso
    int *no,		//   no = number of observations
    int *ni,		//   ni = number of predictor variables
    double *x,		//   x[ni][no] = predictor data matrix flat file (overwritten)
    double *y,		//   y[no] = response vector (overwritten)
    double *w,		//   w[no]= observation weights (overwritten)
    int *jd,		//   jd(jd(1)+1) = predictor variable deletion flag
			//      jd(1) = 0  => use all variables
			//      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
    double *vp,		//   vp(ni) = relative penalties for each predictor variable
			//      vp(j) = 0 => jth variable unpenalized
    double cl[][2],	//   cl(2,ni) = interval constraints on coefficient values (overwritten)
			//      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
			//      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
    int *ne,		//   ne = maximum number of variables allowed to enter largest model
			//        (stopping criterion)
    int *nx,		//   nx = maximum number of variables allowed to enter all models
			//        along path (memory allocation, nx > ne).
    int *nlam,		//   nlam = (maximum) number of lamda values
    double *flmin,	//   flmin = user control of lamda values (>=0)
			//      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
			//      flmin >= 1.0 => use supplied lamda values (see below)
    double *ulam,	//   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
    double *thr,	//   thr = convergence threshold for each lamda solution.
			//      iterations stop when the maximum reduction in the criterion value
			//      as a result of each parameter update over a single pass
			//      is less than thr times the null criterion value.
			//      (suggested value, thr=1.0e-5)
    int *isd,		//   isd = predictor variable standarization flag:
			//      isd = 0 => regression on original predictor variables
			//      isd = 1 => regression on standardized predictor variables
			//      Note: output solutions always reference original
			//            variables locations and scales.
    int *intr,		//   intr = intercept flag
			//      intr = 0/1 => don't/do include intercept in model
    int *maxit,		//   maxit = maximum allowed number of passes over the data for all lambda
			//      values (suggested values, maxit = 100000)

// output:

    int *lmu,		//   lmu = actual number of lamda values (solutions)
    double *a0,		//   a0(lmu) = intercept values for each solution
    double *ca,		//   ca(nx,lmu) = compressed coefficient values for each solution
    int *ia,		//   ia(nx) = pointers to compressed coefficients
    int *nin,		//   nin(lmu) = number of compressed coefficients for each solution
    double *rsq,	//   rsq(lmu) = R**2 values for each solution
    double *alm,	//   alm(lmu) = lamda values corresponding to each solution
    int *nlp,		//   nlp = actual number of passes over the data for all lamda values
    int *jerr		//   jerr = error flag:
			//      jerr  = 0 => no error
			//      jerr > 0 => fatal error - no output returned
			//         jerr < 7777 => memory allocation error
			//         jerr = 7777 => all used predictors have zero variance
			//         jerr = 10000 => maxval(vp) <= 0.0
			//      jerr < 0 => non fatal error - partial output:
			//         Solutions for larger lamdas (1:(k-1)) returned.
			//         jerr = -k => convergence for kth lamda value not reached
			//            after maxit (see above) iterations.
			//         jerr = -10000-k => number of non zero coefficients along path
			//            exceeds nx (see above) at kth lamda value.
);

/*--------------------------------------------------------------------
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx = input to elnet
c    lmu,ca,ia,nin = output from elnet
c
c output:
c
c    b(ni,lmu) = all elnet returned solutions in uncompressed format
----------------------------------------------------------------------*/
extern int solns_(
    int *ni,            //   ni = number of predictor variables
    int *nx,            //   nx = maximum number of variables allowed to enter all models
    int *lmu,           //   lmu = actual number of lamda values (solutions)
    double *ca,         //   ca(nx,lmu) = compressed coefficient values for each solution
    int *ia,            //   ia(nx) = pointers to compressed coefficients
    int *nin,           //   nin(lmu) = number of compressed coefficients for each solution
    double *b           //   b(ni,lmu) = compressed coefficient values for each solution
);

extern int c_glmnet(
    int no,             // number of observations (no)
    int ni,             // number of predictor variables (ni)
    double *x,          // input matrix, x[ni][no]
    double *y,          // response vaiable, of dimentions (no)
    int nlam,           // number of lambda values
    double *ulam,       // value of lambda values, of dimentions (nlam)
    double parm,        // the alpha variable

    int *lmu,           // lmu = actual number of lamda values (solutions)
    double cfs[nlam][ni+1] // results = cfs[lmu][ni + 1]
);

#endif /* CCDC_H */
