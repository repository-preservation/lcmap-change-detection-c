/* this is an effort to consolidate all defines.  Previously, they    */
/* scattered throughout .c and/or .h files, and not always included   */
/* and/or avilable everywhere or where needed.  Also, some reduncancy */
/* and conflicts existed.                                             */

/* from ccdc.c */
#define NUM_LASSO_BANDS 5 /* Number of bands for Least Absolute Shrinkage */
                          /* and Selection Operator LASSO regressions */
#define TOTAL_IMAGE_BANDS 7 /* Number of image bands, including thermal*/
#define TOTAL_BANDS 8     /* Total image plus mask bands, for loops.  */
#define MIN_NUM_C 4       /* Minimum number of coefficients           */
#define MID_NUM_C 6       /* Mid-point number of coefficients         */
#define MAX_NUM_C 8       /* Maximum number of coefficients           */
#define CONSE 6           /* No. of CONSEquential pixels 4 bldg. model*/
#define N_TIMES 3         /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25  /* average number of days per year          */
#define MAX_NUM_FC 30000  /* Values change with number of pixels run */
//#define NUM_FC 10         /* Values change with number of pixels run  */
#define T_CONST 4.89      /* Threshold for cloud, shadow, and snow detection */
#define MIN_YEARS 1       /* minimum year for model intialization     */
#define T_SN 0.75         /* no change detection for permanent snow pixels */ 
#define T_CLR 0.25        /* Fmask fails threshold                    */
#define T_CG 15.0863      /* chi-square inversed T_cg (0.99) for noise removal */
#define T_MAX_CG 35.8882  /* chi-square inversed T_max_cg (1e-6) for 
                             last step noise removal                  */


/* from 2darray.c */
/* Define a unique (i.e. random) value that can be used to verify a pointer
   points to an LSRD_2D_ARRAY. This is used to verify the operation succeeds to
   get an LSRD_2D_ARRAY pointer from a row pointer. */
#define SIGNATURE 0x326589ab

/* Given an address returned by the allocate routine, get a pointer to the
   entire structure. */
#define GET_ARRAY_STRUCTURE_FROM_PTR(ptr) \
    ((LSRD_2D_ARRAY *)((char *)(ptr) - offsetof(LSRD_2D_ARRAY, memory_block)))


/* from input.c */
//#define TOTAL_IMAGE_BANDS 7

/* from misc.c */
/* 12-31-1972 is 720624 in julian day since year 0000 */
#define JULIAN_DATE_LAST_DAY_1972 720624 
#define LANDSAT_START_YEAR 1973
#define LEAP_YEAR_DAYS 366
#define NON_LEAP_YEAR_DAYS 365
#define AVE_DAYS_IN_A_YEAR 365.25
#define ROBUST_COEFFS 5
#define LASSO_COEFFS 8
//#define TOTAL_IMAGE_BANDS 7

/* from input.h */
/* possible cfmask values */
#define CFMASK_CLEAR   0
#define CFMASK_WATER   1
#define CFMASK_SHADOW  2
#define CFMASK_SNOW    3
#define CFMASK_CLOUD   4
#define CFMASK_FILL  255
#define IMAGE_FILL -9999
#define CFMASK_BAND    7
#define THERMAL_BAND   6

/* from output.h */
#define FILL_VALUE 255
#define NUM_COEFFS 8
#define NUM_BANDS 7
