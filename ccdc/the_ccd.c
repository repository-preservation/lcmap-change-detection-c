/******************************************************************************/
/* File: the_ccd.c                                                            */
/******************************************************************************/

#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/timeb.h>

#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "ccdc.h"
//#include "matio.h" // for future matlab output.......
#include "defines.h"


/*******************************************************************************
MODULE: the_ccd

PURPOSE: receives the pixel pointer to all the pixels, and does the actual
         processing, then returns the coefficients and other values to main.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160923     Brian Davis      Original development

NOTES:
*******************************************************************************/

int the_ccd
(
    int valid_num_scenes,            /* Number of scenes after fmask filtering*/
    int row,                         /* The input row of the data frame.      */
    int col,                         /* The input column of the data frame.   */
    int nrows,                       /* The number of rows in "tile" of data  */
    int ncols,                       /* Number of columns in "tile" of data   */
    int num_samples,                 /* Number of samples in input grid       */
    int *updated_sdate_array,        /* sdate array after cfmask filtering    */
    unsigned char *updated_fmask_buf,/*sub-set of fmask buf, valid pixels only*/
    int *buf,                        /* This is the image bands buffer.       */
    int clr_sum,                     /* Total number of clear cfmask pixels   */
    int sn_sum,                      /* Total number of snow  cfmask pixels   */
    int water_sum,                   /* Total of cfmask water pixels.         */
    int all_sum,                     /* Total of all cfmask pixels            */
    bool verbose,                    /* Verbose flag for printing messages    */
    Output_t *rec_cg,                /* Output structure and metadata         */
    int *num_fc                      /* Number of fitted curves, cumulative   */
                                     /* for all rows/cols and models          */
)

{
    int status;                      /* Return value from function call       */
    char FUNC_NAME[] = "the_ccd";    /* For printing error messages           */
    char msg_str[MAX_STR_LEN];       /* Input data scene name                 */

    int i, k, m, b, k_new;           /* Loop counters                         */

    int num_c = 8;                   /* Max number of coefficients for model  */
    int rec_fc;                      /* Record num. of functional curves      */
    float v_start[NUM_LASSO_BANDS];  /* Vector for start of observation(s)    */
    float v_end[NUM_LASSO_BANDS];    /* Vector for end of observastion(s)     */
    float v_slope[NUM_LASSO_BANDS];  /* Vector for anormalized slope values   */
    float v_dif[NUM_LASSO_BANDS];    /* Vector for difference values          */
    float **v_diff;

    float clr_pct;                   /* Percent clear cfmask pixels           */
    float sn_pct;                    /* Percent snow cfmask pixels            */
    int n_sn = 0;                    /* Number of snow cfmask pixels          */
    int n_clr = 0;                   /* Number of clear cfmask pixels         */
    int *clrx;                       /* clear pixel curve in X direction ?    */
    float **clry;                    /* clear pixel curve in Y direction ?    */
    int *cpx;                        /* nunber of clear pixels X ?            */
    float **cpy;                     /* nunber of clear pixels Y ?            */
    int i_start;                     /* The first observation for TSFit       */
    int end;                         /* The end of clear observations of total*/
    int end2;                        /* place holder until final end is found */
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float *rmse;                     /* Root Mean Squared Error array.        */
    int i_span;                      /* index for span of consecutive obs. ?  */
    int update_num_c=MIN_NUM_C;      /* Number of coefficients to update      */
    int bl_train;                    /* Flag for which way to train the model.*/
    float time_span;                 /* Span of time in no. of years.         */
    int *bl_ids;
    int *id_range;
    int *ids;
    int *ids_old;
    int *rm_ids;
    int rm_ids_len;
    int ids_len;                     /*number of ids, incremented continuously*/
    int ids_old_len;
    int i_rec;                       /* start of model before noise removal   */
    float v_dif_norm;                /* vector for normalized difference vals */
    int i_count;                     /* Count difference of i each iteration  */
    float **v_dif_mag;               /* vector for magnitude of differences.  */
    int i_conse, i_b;
    float *vec_mag;                  /* permanent storage for magnitude vector*/
    float *vec_magg;                 /* temp storage for magnitude vecor      */
    float v_dif_mean;
    float vec_magg_min;
    float **rec_v_dif;               /* record of vector differences          */
    float **rec_v_dif_copy;          /* copy of record of vector differences  */
    float **temp_v_dif;              /* for the thermal band.......           */
    float adj_rmse[TOTAL_IMAGE_BANDS];/* Adjusted RMSE for all bands          */
    float mini_rmse;                 /* Mimimum RMSE                          */
    int bl_tmask;
    int n_rmse;                      /* number of RMSE values                 */
    int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5}; /* LASSO band indexes */
    float tmpcg_rmse[NUM_LASSO_BANDS]; /* to temporarily change RMSE          */
    int d_rt;
    float *d_yr;
    int id_last;                     /* The last stable id.                   */
    float ts_pred_temp;
    int i_break;                     /* for recording break points, i is index*/
    int i_ini;                       /* for recording begin of time,i is index*/
    int ini_conse;                   /* Initial CONSE.                        */
    float break_mag;

    bool debug = 1;                  /* This replaces "ifdef 0" convention.   */

    int row_inx, col_inx;            /* for looping through 2D array of pixels*/
    int row_count, col_count;        /*keep track of row/col relative to start*/
    int offset;                      /* to calc buf ptr location of pixels    */
    int scene_inx;                   /* for looping through scenes            */
    int row_size;                    /* size of a row in bytes.               */
    int col_size;                    /* size of a column in bytes.            */
    int scn_size;                    /* size of a scene in bytes.             */

    time_t now;                      /* For logging the start, stop, and some */
    time (&now);                     /*     intermediate times.               */


    /**************************************************************************/
    /*                                                                        */
    /* Allocate memmory for the variables need inside the algorithm.          */
    /*                                                                        */
    /**************************************************************************/

    clrx = malloc(valid_num_scenes * sizeof(int));
    if (clrx == NULL)
    {
        RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);
    }
    
    id_range = (int *)calloc(valid_num_scenes, sizeof(int));
    if (id_range == NULL)
    {
        RETURN_ERROR("ERROR allocating id_range memory", FUNC_NAME, FAILURE);
    }
    
    ids = (int *)calloc(valid_num_scenes, sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }
    
    ids_old = (int *)calloc(valid_num_scenes, sizeof(int));
    if (ids_old == NULL)
    {
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);
    }
    
    bl_ids = (int *)calloc(valid_num_scenes, sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }
    
    rm_ids = (int *)calloc(valid_num_scenes, sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }
    
    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));
    if (clry == NULL)
    {
        RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
    }
    
    fit_cft = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, MAX_NUM_C, 
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }
    
    rmse = (float *)calloc(TOTAL_IMAGE_BANDS, sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }
    
    vec_mag = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }
    
    v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, CONSE,
                sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory", FUNC_NAME, FAILURE);
    }
    
    rec_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                            sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }
    rec_v_dif_copy = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                                 sizeof (float));
    if (rec_v_dif_copy == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif_copy memory",FUNC_NAME, FAILURE);
    }
    
    temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                             sizeof (float));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }

    for (row_inx = (row - 1), row_count = 0; row_count < nrows; row_count++, row_inx++)
    {
      for (col_inx = (col - 1), col_count = 0; col_count < ncols; col_count++, col_inx++)
      {

        if (debug)
            printf("row_inx %d col_inx %d\n", row_inx, col_inx);

        /******************************************************************/
        /*                                                                */
        /* Percent of clear pixels: clear (0) or water (1).               */    
        /*                                                                */
        /******************************************************************/

        clr_pct = (float) (clr_sum + water_sum) / (float) all_sum;

        /******************************************************************/
        /*                                                                */
        /* percent of snow observations (3).                              */
        /*                                                                */
        /******************************************************************/

        if ((clr_sum + water_sum + sn_sum) != 0)
            sn_pct =  (float) sn_sum / (float)(clr_sum + water_sum + sn_sum); 
        else
            sn_pct = (float)sn_sum;

        //if (verbose)
        //{
        //    printf("  Number inputs specified      = %d\n", all_sum);
        //    printf("  Number of non-overlap pixels = %d\n", (all_sum - swath_overlap_count));
        //    printf("  Number of fill (255)  pixels = %d\n", fill_sum);
        //    printf("  Number of non-fill    pixels = %d\n", all_sum);
        //    printf("  Number of clear  (0)  pixels = %d\n", clr_sum);
        //    printf("  Number of water  (1)  pixels = %d\n", water_sum);
        //    printf("  Number of shadow (2)  pixels = %d\n", shadow_sum);
        //    printf("  Number of snow   (3)  pixels = %d\n", sn_sum);
        //    printf("  Number of cloud  (4)  pixels = %d\n", cloud_sum);
        //    printf("  Number of clear+water pixels = %d\n", clr_sum + water_sum);
        //    printf("  Percent of clear pixels      = %f (of non-fill pixels)\n", clr_pct);
        //    printf("  Percent of clear pixels      = %f (of non-fill, non-cloud, non-shadow pixels)\n",
        //           (float) clr_sum / (float) (clr_sum + water_sum + sn_sum) * 100);
        //    printf("  Percent of snow  pixels      = %f (of non-fill, non-cloud, non-shadow pixels)\n", (sn_pct * 100));
        //}

        row_size = ncols * valid_num_scenes * TOTAL_IMAGE_BANDS;
        col_size =         valid_num_scenes * TOTAL_IMAGE_BANDS;
        scn_size =                            TOTAL_IMAGE_BANDS;
        for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
        { 
            /**************************************************************/
            /*                                                            */
            /* pixel value ranges should follow physical rules.i          */
            /* Convert Kelvin to Celsius (for new espa data).             */
            /*                                                            */
            /**************************************************************/

            offset = (row_count * row_size) + (col_count * col_size) + (scene_inx * scn_size);

            if (buf[offset + THERMAL_BAND] != -9999)
                buf[offset + THERMAL_BAND] = (int16)(buf[offset + THERMAL_BAND] * 10 - 27315);

            if ((buf[offset + 0] > 0)     && (buf[offset + 0] < 10000) &&
                (buf[offset + 1] > 0)     && (buf[offset + 1] < 10000) &&
                (buf[offset + 2] > 0)     && (buf[offset + 2] < 10000) &&
                (buf[offset + 3] > 0)     && (buf[offset + 3] < 10000) &&
                (buf[offset + 4] > 0)     && (buf[offset + 4] < 10000) &&
                (buf[offset + 5] > 0)     && (buf[offset + 5] < 10000) &&
                (buf[offset + 6] > -9320) && (buf[offset + 6] < 7070))
            {
                id_range[scene_inx] = 1;
            }
            else
            {
                id_range[scene_inx] = 0;
            }
        }

        /**************************************************************/
        /*                                                            */
        /* Initilialize rec_cg                                        */
        /*                                                            */
        /**************************************************************/

        rec_cg[0].t_break = 0;

        /******************************************************************/
        /*                                                                */
        /* Fit permanent snow observations.                               */
        /*                                                                */
        /******************************************************************/

        if (clr_pct < T_CLR)
        {
            if (sn_pct > T_SN)
            {
                n_sn = 0;

                /**********************************************************/
                /*                                                        */
                /* Snow observations are "good" now.                      */
                /*                                                        */
                /**********************************************************/

                for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
                { 
    	            if (((updated_fmask_buf[scene_inx] == CFMASK_SNOW) || 
                         (updated_fmask_buf[scene_inx] < 2))           &&
                          id_range[scene_inx] == 1)
                    {
                        clrx[n_sn] = updated_sdate_array[scene_inx];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                        {
                            offset = (row_count * row_size) + (col_count * col_size) + (scene_inx * scn_size) + k;
                 	    clry[k][n_sn] = (float)buf[offset];
                        }
                        n_sn++;
                    }
                }  
                end = n_sn;

                if (n_sn < N_TIMES * MIN_NUM_C) // not enough snow pixels
                {
                    RETURN_ERROR ("Not enough good snow observations\n", 
                         FUNC_NAME, FAILURE);
                }

                /**********************************************************/
                /*                                                        */
                /* Remove repeated ids.                                   */
                /*                                                        */
                /**********************************************************/

                matlab_unique(clrx, clry, n_sn, &end);

                /**********************************************************/
                /*                                                        */
                /* Start model fit for snow persistent pixels.            */
                /*                                                        */
                /**********************************************************/

                if (verbose)
                    printf ("Fit permanent snow observations, now pixel = %f\n", 
                            100.0 * sn_pct); 

                i_start = 1; /* the first observation for TSFit */

                /**********************************************************/
                /*                                                        */
                /* Treat saturated and unsaturated pixels differently.    */
                /*                                                        */
                /**********************************************************/

                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
                    i_span = 0;
                    for (i = 0; i < end; i++)
                    {
                        if (k != TOTAL_IMAGE_BANDS - 1) // for optical bands
                        {
                            if (clry[i][k] > 0.0 && clry[i][k] < 10000.0)
                            {
                                clrx[i_span] = clrx[i];
                                clry[i_span][k] = clry[i][k];
                                i_span++;
                            }

                            if (i_span < MIN_NUM_C * N_TIMES)
                                fit_cft[i][k] = 10000; // fixed value for saturated pixels
                            else
                            {
                                status = auto_ts_fit(clrx, clry, k, 0, i_span-1, MIN_NUM_C, 
                                         fit_cft, &rmse[k], temp_v_dif); 
                                if (status != SUCCESS)  
                                    RETURN_ERROR ("Calling auto_ts_fit1\n", 
                                           FUNC_NAME, EXIT_FAILURE);

                            } 
                        }
                        else // for thermal band
                        {
                            if (clry[i][k] > -9300.0 && clry[i][k] < 7070.0)
                            {
                                clrx[i_span] = clrx[i];
                                clry[i_span][k] = clry[i][k];
                                i_span++;
                            }
                            
                            status = auto_ts_fit(clrx, clry, k, 0, i_span-1, MIN_NUM_C, fit_cft, 
                                     &rmse[k], temp_v_dif); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit2\n", 
                                      FUNC_NAME, EXIT_FAILURE);
                        }
                    }
                }

                /******************************************************/
                /*                                                    */
                /* NUM of Fitted Curves (num_fc)                      */
                /*                                                    */
                /******************************************************/

                (*num_fc)++;
                if (*num_fc >= MAX_NUM_FC)
                {
                    /**************************************************/
                    /*                                                */
                    /* Reallocate memory for rec_cg                   */
                    /*                                                */
                    /**************************************************/
                    rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
                    {
                        RETURN_ERROR("ERROR allocating rec_cg memory",
                                     FUNC_NAME, FAILURE);
                    }
                }

                /**********************************************************/
                /*                                                        */
                /* Update information at each iteration                   */
                /* Record time of curve start, end.                       */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].t_start = clrx[i_start-1]; 
                rec_cg[*num_fc].t_end = clrx[end-1]; 

                /**********************************************************/
                /*                                                        */
                /* No break at the moment.                                */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].t_break = 0; 

                /**********************************************************/
                /*                                                        */
                /* Record postion of the pixel.                           */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].pos = (row) * num_samples + col_inx + 1; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Record fitted coefficients.                    */
                        /*                                                */
                        /**************************************************/

                        rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record rmse of the pixel.                          */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].rmse[i_b] = rmse[i_b]; 
                }

                /**********************************************************/
                /*                                                        */
                /* Record change probability, number of observations.     */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].change_prob = 0.0; 
                rec_cg[*num_fc].num_obs = n_sn; 
                rec_cg[*num_fc].category = 50 + MIN_NUM_C; /* snow pixel */

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */ 
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].magnitude[i_b] = 0.0; 
                }

            }  // if sn_pct > T_SN

            else

            {

                /**********************************************************/
                /*                                                        */
                /* No change detection for clear observations, backup     */
                /* algorithm.                                             */
                /*                                                        */
                /**********************************************************/

                n_clr = 0;

                for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
                { 
                    if (id_range[scene_inx] == 1)
                    {
                        clrx[n_clr] = updated_sdate_array[scene_inx];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                            {
                            offset = (row_count * row_size) + (col_count * col_size) + (scene_inx * scn_size) + k;
                            clry[k][n_clr] = (float)buf[offset];
                            }
                        n_clr++;
                    }   
                }
                end = n_clr;

                /**********************************************************/
                /*                                                        */
                /* Remove repeated ids.                                   */
                /*                                                        */
                /**********************************************************/

                matlab_unique(clrx, clry, n_clr, &end);

                /**********************************************************/
                /*                                                        */
                /* Start model fit for clear persistent pixels.           */
                /*                                                        */
                /**********************************************************/

                if (verbose)
                    printf ("Fmask failed, clear pixel = %f\n", 100.0 * clr_pct); 

                n_clr = 0;
                float band2_median; // probably not good practice to declare here....
                quick_sort_float(clry[1], 0, end - 1);
                matlab_2d_float_median(clry, 1, end, &band2_median);

                for (i = 0; i < end; i++)
                {
                    if (clry[1][i] < (band2_median + 400.0))
                    {
                        clrx[n_clr] = clrx[i];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                        { 
                            clry[k][n_clr] = clry[k+1][i]; 
                        }
                        n_clr++;
                    }
                }
                end = n_clr;

                /**********************************************************/
                /*                                                        */
                /* The first observation for TSFit.                       */
                /*                                                        */
                /**********************************************************/

                i_start = 1; /* the first observation for TSFit */

                if (n_clr < N_TIMES * MIN_NUM_C)
                {
                    continue;
                }
                else
                {
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        status = auto_ts_fit(clrx, clry, i_b, 0, end-1, MIN_NUM_C, 
                                             fit_cft, &rmse[k], temp_v_dif); 
                        if (status != SUCCESS)
                        {  
                            RETURN_ERROR ("Calling auto_ts_fit for clear persistent pixels\n", 
                                          FUNC_NAME, FAILURE);
                        }
                    }
                }

                /******************************************************/
                /*                                                    */
                /* NUM of Fitted Curves (num_fc)                      */
                /*                                                    */
                /******************************************************/

                (*num_fc)++;
                if (*num_fc >= MAX_NUM_FC)
                {
                    /**************************************************/
                    /*                                                */
                    /* Reallocate memory for rec_cg                   */
                    /*                                                */
                    /**************************************************/
                    rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
                    {
                        RETURN_ERROR("ERROR allocating rec_cg memory",
                                     FUNC_NAME, FAILURE);
                    }
                }

                /**********************************************************/
                /*                                                        */
                /* Update information at each iteration.                  */
                /* Record time of curve start, time of curve end.         */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].t_start = clrx[i_start-1]; 
                rec_cg[*num_fc].t_end = clrx[end-1]; 

                /**********************************************************/
                /*                                                        */
                /* No break at the moment.                                */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].t_break = 0; 

                /**********************************************************/
                /*                                                        */
                /* Record postion of the pixel.                           */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].pos = (row-1) * num_samples + col_inx + 1; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Record fitted coefficients.                    */
                        /*                                                */
                        /**************************************************/

                        rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record rmse of the pixel.                          */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].rmse[i_b] = rmse[i_b]; 
                }

                /**********************************************************/
                /*                                                        */
                /* Record change probability, number of observations,     */
                /* fit category.                                          */
                /*                                                        */
                /**********************************************************/
                rec_cg[*num_fc].change_prob = 0.0; 
                rec_cg[*num_fc].num_obs = n_clr; 
                rec_cg[*num_fc].category = 40 + MIN_NUM_C; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */ 
                    /*                                                    */
                    /******************************************************/
                    rec_cg[*num_fc].magnitude[i_b] = 0.0; 
                }
            }
        }
        else /* normal CCDC procedure */
        {
            if (verbose)
            {
                printf("Seasonal Snow (Snow < %f)\n", 100.0 * sn_pct);
                printf("Fmask works, clear pixels (land/water) = %f\n", 100.0 * clr_pct);
            }

            n_clr = 0;
            for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
            { 
                if ((updated_fmask_buf[scene_inx] < 2) && (id_range[scene_inx] == 1)) // offsets
                {
                    clrx[n_clr] = updated_sdate_array[scene_inx]; // offsets
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                    {
                        offset = (row_count * row_size) + (col_count * col_size) + (scene_inx * scn_size) + k;
                        clry[k][n_clr] = (float)buf[offset];
                    }
                    n_clr++;
                }   
            }
            end = n_clr;
            if (debug)
            {
                printf("end_clr=%d\n",end);
            }

            /**************************************************************/
            /*                                                            */
            /* Remove repeated ids.                                       */
            /*                                                            */
            /**************************************************************/

            matlab_unique(clrx, clry, n_clr, &end);

            /**************************************************************/
            /*                                                            */
            /* Calculate median variogram.                                */
            /*                                                            */
            /**************************************************************/

            for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
            {
                adj_rmse[k] = 0.0;
            } 
            status = median_variogram(clry, TOTAL_IMAGE_BANDS, 0, end-1, adj_rmse);
            if (status != SUCCESS)
            {
                RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME, 
                             FAILURE);
            }

            if (debug)
            {
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                    printf("k,adj_rmse[k]=%d,%f\n",k,adj_rmse[k]);
            }

            /**************************************************************/
            /*                                                            */
            /* Start with mininum requirement of clear obs.               */
            /*                                                            */
            /**************************************************************/

            i = N_TIMES * MIN_NUM_C;

            /**************************************************************/
            /*                                                            */
            /* The first observation for TSFit.                           */
            /*                                                            */
            /**************************************************************/

            i_start = 1; 

            /**************************************************************/
            /*                                                            */
            /* Record the start of the model initialization               */
            /*     (0=>initial;1=>done)                                   */
            /*                                                            */
            /**************************************************************/

            bl_train = 0;

            /**********************************************************/
            /*                                                        */
            /* NUM of Fitted Curves (num_fc)                          */
            /*                                                        */
            /**********************************************************/

            (*num_fc)++;
            if (*num_fc >= MAX_NUM_FC)
            {
                /* Reallocate memory for rec_cg */
                rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
                {
                    RETURN_ERROR("ERROR allocating rec_cg memory",
                                 FUNC_NAME, FAILURE);
                }
            }

            /**************************************************************/
            /*                                                            */
            /* Record the num_fc at the beginning of each pixel.          */
            /*                                                            */
            /**************************************************************/
            rec_fc = *num_fc;

            /**************************************************************/
            /*                                                            */
            /* Record the start of Tmask (0=>initial;1=>done)             */
            /*                                                            */
            /**************************************************************/

            bl_tmask = 0;

            /**************************************************************/
            /*                                                            */
            /* If verbose, record the start time of just the CDCD         */
            /*     algorithm.  Up until here, it has all just been        */
            /*     setting it up......                                    */
            /*                                                            */
            /**************************************************************/

            if (verbose)
            {
                snprintf (msg_str, sizeof(msg_str), "CCDC init_time=%s\n", ctime (&now));
                LOG_MESSAGE (msg_str, FUNC_NAME);
            }

            /**************************************************************/
            /*                                                            */
            /* While loop - process til the last clear observation - CONSE*/
            /*                                                            */
            /**************************************************************/

            while (i <= end - CONSE)
            {
                /**********************************************************/
                /*                                                        */
                /* span of "i"                                            */
                /*                                                        */
                /**********************************************************/

                i_span = i - i_start + 1;

                /**********************************************************/
                /*                                                        */
                /* span of time (num of years)                            */
                /*                                                        */
                /**********************************************************/
                time_span = (float)(clrx[i-1] - clrx[i_start-1]) / NUM_YEARS;

                /**********************************************************/
                /*                                                        */
                /* basic requrirements: 1) enough observations;           */
                /*                      2) enough time                    */
                /*                                                        */
                /**********************************************************/

                if ((i_span >= N_TIMES * MIN_NUM_C) && (time_span >= (float)MIN_YEARS))
                {
                    /******************************************************/
                    /*                                                    */
                    /* Initializing model.                                */
                    /*                                                    */
                    /******************************************************/

                    if (bl_train == 0)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Step 1: noise removal.                         */ 
                        /*                                                */
                        /**************************************************/

                        status = auto_mask(clrx, clry, i_start-1, i+CONSE-1,
                                           (float)(clrx[i+CONSE-1]-clrx[i_start-1]) / NUM_YEARS, 
                                           adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR("ERROR calling auto_mask during model initilization", 
                                          FUNC_NAME, FAILURE);
                        }

                        /**************************************************/
                        /*                                                */
                        /* Clear the IDs buffers.                         */
                        /*                                                */
                        /**************************************************/

                        for (k = 0; k < valid_num_scenes; k++)
                            ids[k] = 0;

                        /**************************************************/
                        /*                                                */
                        /* IDs to be removed.                             */
                        /*                                                */
                        /**************************************************/

                        for (k = i_start-1; k < i+CONSE; k++)
                        {
                            ids[k-i_start+1] = k;
                        }
                        m = 0;
                        i_span = 0;
                        for (k = 0; k < i-i_start+1; k++)
                        {
                            if (bl_ids[k] == 1) 
                            {
                                rm_ids[m] = ids[k];
                                m++;
                            }
                            else
                                i_span++;  /* update i_span after noise removal */
                        }

                        rm_ids_len = m;

                        /**************************************************/
                        /*                                                */
                        /* Check if there are enough observation.         */
                        /*                                                */
                        /**************************************************/

                        if (i_span < (N_TIMES * MIN_NUM_C))
                        {
                            /**********************************************/
                            /*                                            */
                            /* Move forward to the i+1th clear observation*/
                            /*                                            */
                            /**********************************************/

                            i++;

                            /**********************************************/
                            /*                                            */
                            /* Not enough clear observations.             */
                            /*                                            */
                            /**********************************************/

                            continue;
                        }
                        else
                        {
    
                            if (end == 0)
                                RETURN_ERROR("No available data point", FUNC_NAME, FAILURE);

                            /**************************************************/
                            /*                                                */
                            /* Allocate memory for cpx, cpy.                  */
                            /*                                                */
                            /**************************************************/

                            cpx = malloc(end * sizeof(int));
                            if (cpx == NULL)
                                RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                            cpy = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, end,
                                             sizeof (float));
                            if (cpy == NULL)
                            {
                                RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                            }

                            /**************************************************/
                            /*                                                */
                            /* Remove noise pixels between i_start & i.       */
                            /*                                                */
                            /**************************************************/

                            m = 0;
                            for (k = 0, k_new=0; k < end; k++)
                            {
                                if (m < rm_ids_len && k == rm_ids[m])
                                {
                                    m++;
                                    continue;
                                }
                                cpx[k_new] = clrx[k];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                {
                                    cpy[b][k_new] = clry[b][k];
                                }
                                k_new++;
                            }
                            end2 = k_new;
                            end = end2;
    
                            /**************************************************/
                            /*                                                */
                            /* Record i before noise removal.                 */ 
                            /* This is very important, ie model is not yet    */
                            /* initialized.   The multitemporal masking shall */
                            /* be done again instead of removing outliers  In */
                            /* every masking.                                 */
                            /*                                                */
                            /**************************************************/

                            i_rec = i;

                            /**************************************************/
                            /*                                                */
                            /* Update i afer noise removal.                   */
                            /*     (i_start stays the same).                  */
                            /*                                                */
                            /**************************************************/

                            i = i_start + i_span - 1;

                            /**************************************************/
                            /*                                                */
                            /* Update span of time (num of years).            */
                            /*                                                */
                            /**************************************************/

                            time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                            /**************************************************/
                            /*                                                */
                            /* Check if there is enough time.                 */
                            /*                                                */
                            /**************************************************/

                            if (time_span < MIN_YEARS)
                            {
                                i = i_rec;   /* keep the original i */

                                /**********************************************/
                                /*                                            */
                                /* Move forward to the i+1th clear observation*/
                                /*                                            */
                                /**********************************************/

                                i++;        
                                free(cpx);
                                status = free_2d_array ((void **) cpy);
                                if (status != SUCCESS)
                                {
                                      RETURN_ERROR ("Freeing memory: cpy\n", 
                                            FUNC_NAME, FAILURE);
                                }
                                continue;    /* not enough time */
                            }
                            else
                            {

                                /**********************************************/
                                /*                                            */
                                /* Remove noise in original arrays.           */
                                /*                                            */
                                /**********************************************/

                                for (k = 0; k < end; k++)
                                {
                                    clrx[k] = cpx[k];
                                    for (m = 0; m < TOTAL_IMAGE_BANDS; m++)
                                    {
                                        clry[m][k] = cpy[m][k];
                                    }
                                }

                                free(cpx);
                                status = free_2d_array ((void **) cpy);
                                if (status != SUCCESS)
                                {
                                    RETURN_ERROR ("Freeing memory: cpy\n", 
                                         FUNC_NAME, FAILURE);
                                }

                                /**********************************************/
                                /*                                            */
                                /* Step 2) model fitting: initialize model    */
                                /*         testing variables defining         */
                                /*         computed variables.                */
                                /*                                            */
                                /**********************************************/

                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                {

                                    /******************************************/
                                    /*                                        */
                                    /* Initial model fit.                     */
                                    /*                                        */
                                    /******************************************/

                                    status = auto_ts_fit(clrx, clry, b, i_start-1, i-1, 
                                             MIN_NUM_C, fit_cft, &rmse[b], rec_v_dif); 
                                    if (status != SUCCESS)  
                                    {
                                        RETURN_ERROR ("Calling auto_ts_fit during model initilization\n", 
                                             FUNC_NAME, FAILURE);
                                    }
                                }

                                v_dif_norm = 0.0;
                                for(b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                    /******************************************/
                                    /*                                        */
                                    /* Calculate min. rmse.                   */
                                    /*                                        */
                                    /******************************************/

                                    mini_rmse = max(adj_rmse[lasso_blist[b]], rmse[lasso_blist[b]]);

                                    /******************************************/
                                    /*                                        */
                                    /* Compare the first observation.         */
                                    /*                                        */
                                    /******************************************/

                                    v_start[b] = rec_v_dif[lasso_blist[b]][0] 
                                            / mini_rmse;

                                    /******************************************/
                                    /*                                        */
                                    /* Compare the last clear observation.    */
                                    /*                                        */
                                    /******************************************/

                                    v_end[b] = rec_v_dif[lasso_blist[b]][i-i_start]
                                                            / mini_rmse;

                                    /******************************************/
                                    /*                                        */
                                    /* Anormalized slope values.              */
                                    /*                                        */
                                    /******************************************/

                                    v_slope[b] = fit_cft[lasso_blist[b]][1] *
                                                    (clrx[i-1]-clrx[i_start-1])/mini_rmse;

                                    /******************************************/
                                    /*                                        */
                                    /* Difference in model intialization.     */
                                    /*                                        */
                                    /******************************************/

                                    v_dif[b] = fabs(v_slope[b]) + fabs(v_start[b]) + fabs(v_end[b]); 
                                    v_dif_norm += v_dif[b] * v_dif[b];               
                                }
    
                                /**********************************************/
                                /*                                            */
                                /* Find stable start for each curve.          */
                                /*                                            */
                                /**********************************************/

                                if (v_dif_norm > T_CG)
                                {
                                    /******************************************/
                                    /*                                        */
                                    /* Start from next clear observation.     */
                                    /*                                        */
                                    /******************************************/

                                    i_start++;

                                    /******************************************/
                                    /*                                        */
                                    /* Move forward to the i+1th clear        */
                                    /* observation.                           */
                                    /*                                        */
                                    /******************************************/

                                    i++;

                                    /******************************************/
                                    /*                                        */
                                    /* Keep all data and move to the next obs.*/
                                    /*                                        */
                                    /******************************************/

                                    continue;
                                }
                                else
                                {

                                    /******************************************/
                                    /*                                        */
                                    /* Model is ready.                        */
                                    /*                                        */
                                    /******************************************/

                                    bl_train = 1;

                                    /******************************************/
                                    /*                                        */
                                    /* Count difference of i for each         */
                                    /* iteration.                             */
                                    /*                                        */
                                    /******************************************/

                                    i_count = 0;

                                    /******************************************/
                                    /*                                        */
                                    /* Find the previous break point.         */
                                    /*                                        */
                                    /******************************************/

                                    if (*num_fc == rec_fc)
                                    {
                                        i_break = 1; /* first curve */
                                    }
                                    else
                                    {
                                        /**************************************/
                                        /*                                    */
                                        /* After the first curve, compare     */
                                        /* rmse to determine which curve to   */
                                        /* determine t_break.                 */
                                        /*                                    */
                                        /**************************************/

                                        for (k = 0; k < end; k++) 
                                        {
                                            if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                                            {
                                                i_break = k + 1;
                                                break;
                                            }
                                        }
                                    }

                                    if (i_start > i_break)
                                    {
                                        /**************************************/
                                        /*                                    */
                                        /* Model fit at the beginning of the  */
                                        /* time series.                       */
                                        /*                                    */
                                        /**************************************/

                                        for(i_ini = i_start-2; i_ini >= i_break-1; i_ini--)
                                        {
                                            if ((i_start - i_break) < CONSE)
                                            {
                                                ini_conse = i_start - i_break;
                                            }
                                            else
                                            {
                                                ini_conse = CONSE;
                                            }


                                            /**********************************/
                                            /*                                */
                                            /* Allocate memory for            */ 
                                            /* model_v_dif, v_diff, vec_magg  */ 
                                            /* for the non-stdin branch here. */
                                            /*                                */
                                            /**********************************/

                                            v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                                        ini_conse, sizeof (float));
                                            if (v_diff == NULL)
                                            {
                                                RETURN_ERROR ("Allocating v_diff memory", 
                                                               FUNC_NAME, FAILURE);
                                            }
 
                                            vec_magg = (float *) malloc(ini_conse * sizeof (float));
                                            if (vec_magg == NULL)
                                            {
                                                RETURN_ERROR ("Allocating vec_magg memory", 
                                                              FUNC_NAME, FAILURE);
                                            }

                                            /**********************************/
                                            /*                                */
                                            /* Detect change.                 */
                                            /* value of difference for CONSE  */
                                            /* obsservations                  */
                                            /* Record the magnitude of change */
                                            /*                                */
                                            /**********************************/

                                            vec_magg_min = 9999.0;
                                            for (i_conse = 0; i_conse < ini_conse; i_conse++)
                                            {
                                                v_dif_norm = 0.0;
                                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                                {
                                                    /**************************/
                                                    /*                        */
                                                    /* Absolute differences.  */
                                                    /*                        */
                                                    /**************************/

                                                    auto_ts_predict(clrx, fit_cft, MIN_NUM_C, i_b, i_ini-i_conse,
                                                                    i_ini-i_conse, &ts_pred_temp);
                                                    v_dif_mag[i_b][i_conse] = (float)clry[i_b][i_ini-i_conse] - 
                                                                       ts_pred_temp;

                                                    /**************************/
                                                    /*                        */
                                                    /* Normalize to z-score.  */
                                                    /*                        */
                                                    /**************************/

                                                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                                                    {
                                                        if (i_b == lasso_blist[b])
                                                        {
                                                            /******************/
                                                            /*                */
                                                            /* Minimum rmse.  */ 
                                                            /*                */
                                                            /******************/

                                                            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);

                                                            /******************/
                                                            /*                */
                                                            /* z-scores.      */
                                                            /*                */
                                                            /******************/

                                                            v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] 
                                                                                          / mini_rmse;
                                                            v_dif_norm += v_diff[b][i_conse] * v_diff[b][i_conse];
                                                        }
                                                    }
                                                }
                                                vec_magg[i_conse] = v_dif_norm; 

                                                if (vec_magg_min > vec_magg[i_conse])
                                                {
                                                    vec_magg_min =  vec_magg[i_conse];
                                                }
                                            }

                                            /**********************************/
                                            /*                                */
                                            /* Change detection.              */
                                            /*                                */
                                            /**********************************/

                                            if (vec_magg_min > T_CG) /* change detected */
                                            {
                                                break;
                                            }
                                            else if (vec_magg[0] > T_MAX_CG) /* false change */
                                            {
                                                for (k = i_ini; k < end - 1; k++)
                                                {
                                                    clrx[k] = clrx[k+1];
                                                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                                    {
                                                        clry[b][k] = clry[b][k+1];
                                                    }
                                                }
                                                i--;
                                                end--;

                                                /******************************/
                                                /*                            */
                                                /* Update i_start if i_ini is */
                                                /* not a confirmed break.     */
                                                /*                            */
                                                /******************************/

                                                i_start = i_ini;
                                            }

                                            /**********************************/
                                            /*                                */
                                            /* Free the temporary memory.     */
                                            /*                                */
                                            /**********************************/

                                            free(vec_magg);
                                            status = free_2d_array ((void **) v_diff);
                                            if (status != SUCCESS)
                                            {
                                                RETURN_ERROR ("Freeing memory: v_diff\n", 
                                                              FUNC_NAME, FAILURE);
                                            }
                                        }
                                    }

                                    /******************************************/
                                    /*                                        */
                                    /* Enough to fit simple model and confirm */
                                    /* a break.                               */
                                    /*                                        */
                                    /******************************************/

                                    if ((*num_fc == rec_fc) && ((i_start - i_break) >= CONSE))
                                    {
                                        /**************************************/
                                        /*                                    */
                                        /* Defining computed variables.       */
                                        /*                                    */
                                        /**************************************/

                                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                        {
                                            status = auto_ts_fit(clrx, clry, i_b, i_break-1, i_start-2, 
                                                     MIN_NUM_C, fit_cft, &rmse[i_b], temp_v_dif); 
                                            if (status != SUCCESS)
                                            {  
                                                  RETURN_ERROR ("Calling auto_ts_fit with enough observations\n", 
                                                             FUNC_NAME, FAILURE);
                                            }
                                        }

                                        /**************************************/
                                        /*                                    */
                                        /* Record time of curve end,          */
                                        /* postion of the pixels.             */
                                        /*                                    */
                                        /**************************************/

                                        rec_cg[*num_fc].t_end = clrx[i_start-2]; 
                                        rec_cg[*num_fc].pos = (row-1) * num_samples + col_inx + 1; 

                                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                        {
                                            for (k = 0; k < MIN_NUM_C; k++)
                                            {
                                                /******************************/
                                                /*                            */
                                                /* Record fitted coefficients.*/
                                                /*                            */
                                                /******************************/

                                                rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                             }

                                            /**********************************/
                                            /*                                */
                                            /* Record rmse of the pixel.      */                         
                                            /*                                */
                                            /**********************************/

                                            rec_cg[*num_fc].rmse[i_b] = rmse[i_b]; 
                                        }

                                        /**************************************/
                                        /*                                    */
                                        /* Record break time, fit category,   */
                                        /* change probability, time of curve  */
                                        /* start, number of observations,     */
                                        /* change magnitude.                  */
                                        /*                                    */
                                        /**************************************/

                                        rec_cg[*num_fc].t_break = clrx[i_start -1];
                                        rec_cg[*num_fc].category = 10 + MIN_NUM_C;
                                        rec_cg[*num_fc].change_prob = 1.0;
                                        rec_cg[*num_fc].t_start = clrx[0];
                                        rec_cg[*num_fc].num_obs = i_start - i_break;

                                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                        {
                                            quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                                            matlab_2d_float_median(v_dif_mag, i_b, ini_conse, 
                                                                  &v_dif_mean);
                                            rec_cg[*num_fc].magnitude[i_b] = -v_dif_mean; 
                                        }

                                        /**************************************/
                                        /*                                    */
                                        /* Identified and move on for the     */
                                        /* nex tfunctional curve.             */
                                        /*                                    */
                                        /**************************************/

                                        (*num_fc)++;  
                                        if (*num_fc >= MAX_NUM_FC)
                                        {
                                            /**********************************/
                                            /*                                */
                                            /* Reallocate memory for rec_cg.  */ 
                                            /*                                */
                                            /**********************************/

                                            rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                                            if (rec_cg == NULL)
                                            {
                                                RETURN_ERROR("ERROR allocating rec_cg memory", 
                                                             FUNC_NAME, FAILURE);
                                            }
                                        }           
                                    }
                                }
                            }
                        }
                    } /* end of initializing model */
 
                    /******************************************************/
                    /*                                                    */
                    /* Allocate memory for v_diff for the non-stdin branch*/ 
                    /*                                                    */
                    /******************************************************/

                    v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                      CONSE, sizeof (float));
                    if (v_diff == NULL)
                    {
                        RETURN_ERROR ("Allocating v_diff memory", 
                                      FUNC_NAME, FAILURE);
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Continuous monitoring started!!!                   */
                    /*                                                    */
                    /******************************************************/

                    if (bl_train == 1)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Clears the IDs buffers.                        */
                        /*                                                */
                        /**************************************************/

                        for (k = 0; k < valid_num_scenes; k++)
                        {
                            ids[k] = 0;
                        }

                        /**************************************************/
                        /*                                                */
                        /* All IDs.                                       */
                        /*                                                */
                        /**************************************************/

                        ids_len = 0;
                        for (k = i_start-1; k < i; k++)
                        {
                            ids[k-i_start+1] = k;
                            ids_len++;
                        }
                        i_span = i - i_start +1;

                        /**************************************************/
                        /*                                                */
                        /* Determine the time series model.               */
                        /*                                                */
                        /**************************************************/

                        update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C, 
                                   num_c, &update_num_c);

                        /******************************************************/
                        /*                                                    */
                        /* initial model fit when there are not many observations.*/
                        /* if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C))*/
                        /*                                                    */
                        /******************************************************/

                        if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                        {
                            /**********************************************/
                            /*                                            */
                            /* update i_count at each iteration.          */
                            /*                                            */
                            /**********************************************/

                            i_count = clrx[i-1] - clrx[i_start-1];

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                status = auto_ts_fit(clrx, clry, i_b, i_start-1, i-1, update_num_c, 
                                                     fit_cft, &rmse[i_b], rec_v_dif); 
                                if (status != SUCCESS) 
                                { 
                                    RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n", 
                                                  FUNC_NAME, FAILURE);
                                }
                            }

                            /**********************************************/
                            /*                                            */
                            /* Updating information for the first         */
                            /* iteration.  Record time of curve start and */
                            /* time of curve end.                         */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].t_start = clrx[i_start-1]; 
                            rec_cg[*num_fc].t_end = clrx[i-1]; 

                            /**********************************************/
                            /*                                            */
                            /* No break at the moment.                    */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].t_break = 0; 

                            /**********************************************/
                            /*                                            */
                            /* Record postion of the pixel.               */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].pos = (row-1) * num_samples + col_inx + 1; // offsets

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                for (k = 0; k < MAX_NUM_C; k++)
                                {
                                    /**************************************/
                                    /*                                    */
                                    /* Record fitted coefficients.        */
                                    /*                                    */
                                    /**************************************/

                                    rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                }
                                /******************************************/
                                /*                                        */
                                /* Record rmse of the pixel.              */
                                /*                                        */
                                /******************************************/

                                rec_cg[*num_fc].rmse[i_b] = rmse[i_b]; 
                            }
                            /**********************************************/
                            /*                                            */
                            /* Record change probability, number of       */
                            /* observations, fit category.                */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].change_prob = 0.0; 
                            rec_cg[*num_fc].num_obs = i-i_start+1; 
                            rec_cg[*num_fc].category = 0 + update_num_c; 

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                /******************************************/
                                /*                                        */
                                /* Record change magnitude.               */
                                /*                                        */
                                /******************************************/
                                rec_cg[*num_fc].magnitude[i_b] = 0.0; 
                            }

                            /**********************************************/
                            /*                                            */
                            /* Detect change, value of difference for     */
                            /* CONSE observations.                        */
                            /*                                            */
                            /**********************************************/

                            for (i_conse = 0; i_conse < CONSE; i_conse++)
                            {
                                v_dif_norm = 0.0;
                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    /**************************************/
                                    /*                                    */
                                    /* Absolute differences.              */
                                    /*                                    */
                                    /**************************************/

                                    auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i+i_conse, i+i_conse, 
                                                    &ts_pred_temp);
                                    v_dif_mag[i_b][i_conse] = (float)clry[i_b][i+i_conse] - ts_pred_temp; 

                                    /**************************************/
                                    /*                                    */
                                    /* Normalize to z-score.              */
                                    /*                                    */
                                    /**************************************/

                                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                                    {
                                        if (i_b == lasso_blist[b])
                                        {
                                            /******************************/
                                            /*                            */
                                            /* Minimum rmse, z-scores.    */ 
                                            /*                            */
                                            /******************************/

                                            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
                                            v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] / mini_rmse;
                                            v_dif_norm += v_diff[b][i_conse] * v_diff[b][i_conse];
                                        }
                                    }
                                }
                                vec_mag[i_conse] = v_dif_norm;
                            }

                            /**********************************************/
                            /*                                            */
                            /* Clears the IDs_old buffers.                */
                            /*                                            */
                            /**********************************************/

                            for (k = 0; k < ids_len; k++)
                            {
                                ids_old[k] = 0;
                            }

                            /**********************************************/
                            /*                                            */
                            /* IDs that have not been updated.            */
                            /*                                            */
                            /**********************************************/

                            for (k = 0; k < ids_len; k++)
                            {
                                ids_old[k] = ids[k];
                            }
                            ids_old_len = ids_len;

                        }
                        else
                        {
                            if ((float)(clrx[i-1] - clrx[i_start-1]) >= (1.33*(float)i_count))
                            {
                                /******************************************/
                                /*                                        */
                                /* Update i_count at each iteration year. */
                                /*                                        */
                                /******************************************/

                                i_count = clrx[i-1] - clrx[i_start-1];

                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    status = auto_ts_fit(clrx, clry, i_b, i_start-1, i-1, update_num_c, 
                                                         fit_cft, &rmse[i_b], rec_v_dif); 
                                    if (status != SUCCESS)  
                                    {
                                        RETURN_ERROR ("Calling auto_ts_fit for change detection with "
                                             "enough observations\n", FUNC_NAME, FAILURE);
                                    }
                                }

                                /******************************************/
                                /*                                        */
                                /* Record fitted coefficients.            */
                                /*                                        */
                                /******************************************/
                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    for (k = 0; k < MAX_NUM_C; k++)
                                    {
                                        /**********************************/
                                        /*                                */
                                        /* Record fitted coefficients.    */
                                        /*                                */
                                        /**********************************/

                                        rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                                    } 
                                    /**************************************/
                                    /*                                    */
                                    /* Record rmse of the pixel.          */
                                    /*                                    */
                                    /**************************************/

                                    rec_cg[*num_fc].rmse[i_b] = rmse[i_b]; 
                                }
                                /******************************************/
                                /*                                        */
                                /* Record number of observations, fit     */
                                /* category.                              */
                                /*                                        */
                                /******************************************/

                                rec_cg[*num_fc].num_obs = i-i_start+1; 
                                rec_cg[*num_fc].category = 0 + update_num_c; 

                                /******************************************/
                                /*                                        */
                                /* Clears the IDs_Old buffers.            */
                                /*                                        */
                                /******************************************/

                                for (k = 0; k < ids_len; k++)
                                {
                                    ids_old[k] = 0;
                                }

                                /******************************************/
                                /*                                        */
                                /* IDs that have not been updated.        */
                                /*                                        */
                                /******************************************/

                                for (k = 0; k < ids_len; k++)
                                {
                                    ids_old[k] = ids[k];
                                }
                                ids_old_len = ids_len;

                            }

                            /**********************************************/
                            /*                                            */
                            /* Record time of curve end.                  */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].t_end = clrx[i-1];

                            /**********************************************/
                            /*                                            */
                            /* Use fixed number for RMSE computing.       */
                            /*                                            */
                            /**********************************************/

                            n_rmse = N_TIMES * rec_cg[*num_fc].category;

                            /**********************************************/
                            /*                                            */
                            /* Better days counting for RMSE calculating  */
                            /* relative days distance.                    */
                            /*                                            */
                            /**********************************************/

                            if (ids_old_len == 0)
                            {
                                RETURN_ERROR ("No data points for RMSE calculating", 
                                             FUNC_NAME, FAILURE);
                            }

                            d_yr = malloc(ids_old_len * sizeof(float));
                            if (d_yr == NULL)
                            {
                                RETURN_ERROR ("Allocating d_yr memory", 
                                             FUNC_NAME, FAILURE);
                            }
 
                            for(m = 0; m < ids_old_len; m++)
                            {
                                d_rt = clrx[ids_old[m]] - clrx[i+CONSE-1]; 
                                d_yr[m] = fabs(round((float)d_rt / NUM_YEARS) * NUM_YEARS - (float)d_rt);
                            }

                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            {
                                for (m = 0; m < ids_old_len; m++)
                                    rec_v_dif_copy[b][m] = rec_v_dif[b][m];
                            }

                            /**********************************************/
                            /*                                            */
                            /* Sort the rec_v_dif based on d_yr.          */
                            /*                                            */
                            /**********************************************/

                            quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1);
                            for(b = 0; b < NUM_LASSO_BANDS; b++)
                                tmpcg_rmse[b] = 0.0;

                            /**********************************************/
                            /*                                            */
                            /* Temporarily changing RMSE.                 */
                            /*                                            */
                            /**********************************************/

                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                matlab_2d_array_norm(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                                 &tmpcg_rmse[b]);
                                tmpcg_rmse[b] /= sqrt(n_rmse - rec_cg[*num_fc].category);
                            }

                            /**********************************************/
                            /*                                            */
                            /* Free allocated memories.                   */
                            /*                                            */
                            /**********************************************/
                            free(d_yr);

                            /**********************************************/
                            /*                                            */
                            /* Move the ith col to i-1th col.             */
                            /*                                            */
                            /**********************************************/

                            for (m = 0; m < CONSE-1; m++)
                            {
                                vec_mag[m] = vec_mag[m+1];
                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                    v_diff[b][m] = v_diff[b][m+1];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                    v_dif_mag[b][m] = v_dif_mag[b][m+1];
                            }

                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                                v_diff[b][CONSE-1] = 0.0;
                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                v_dif_mag[b][CONSE-1] = 0.0;
                            vec_mag[CONSE-1] = 0.0;

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                /******************************************/
                                /*                                        */
                                /* Absolute difference for all bands.     */
                                /*                                        */
                                /******************************************/

                                auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i+CONSE-1, 
                                                i+CONSE-1, &ts_pred_temp);
                                v_dif_mag[i_b][CONSE-1] = clry[i_b][i+CONSE-1] - ts_pred_temp;

                                /******************************************/
                                /*                                        */
                                /* Normalized to z-scores.                */
                                /*                                        */
                                /******************************************/
                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                    if (i_b == lasso_blist[b])
                                    {
                                        /**********************************/
                                        /*                                */
                                        /* Minimum rmse.                  */
                                        /*                                */
                                        /**********************************/
                                        mini_rmse = max(adj_rmse[i_b], tmpcg_rmse[b]);

                                        /**********************************/
                                        /*                                */
                                        /* Z-score.                       */
                                        /*                                */
                                        /**********************************/

                                        v_diff[b][CONSE-1] = v_dif_mag[i_b][CONSE-1] / mini_rmse;
                                        vec_mag[CONSE-1] += v_diff[b][CONSE-1] * v_diff[b][CONSE-1]; 
                                    }         
                                }
                            }
                        }

                        break_mag = 9999.0;
                        for (m = 0; m < CONSE; m++)
                        {
                            if (break_mag > vec_mag[m])
                            {
                                break_mag = vec_mag[m];
                            }
                        }

                        if (break_mag > T_CG)
                        {

                            if (debug)
                            {
                                printf("Change Magnitude = %.2f\n", break_mag - T_CG);
                            }

                            /**********************************************/
                            /*                                            */
                            /* Record break time.                        */
                            /*                                            */
                            /**********************************************/

                            rec_cg[*num_fc].t_break = clrx[i];
                            rec_cg[*num_fc].change_prob = 1.0;

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                quick_sort_float(v_dif_mag[i_b], 0, CONSE-1);
                                matlab_2d_float_median(v_dif_mag, i_b, CONSE,
                                                       &rec_cg[*num_fc].magnitude[i_b]);
                            }
                            /**********************************************/
                            /*                                            */
                            /* Identified and move on for the next        */
                            /* functional curve.                          */
                            /*                                            */
                            /**********************************************/

                            (*num_fc)++;

                            if (*num_fc >= MAX_NUM_FC)
                            {
                                /******************************************/
                                /*                                        */
                                /* Reallocate memory for rec_cg.          */ 
                                /*                                        */
                                /******************************************/

                                rec_cg = realloc(rec_cg, (*num_fc + 1) * sizeof(Output_t));
                                if (rec_cg == NULL)
                                {
                                    RETURN_ERROR("ERROR allocating rec_cg memory", 
                                                         FUNC_NAME, FAILURE);
                                }
                            }

                            /**********************************************/
                            /*                                            */
                            /* Start from i+1 for the next functional     */
                            /* curve.                                     */
                            /*                                            */
                            /**********************************************/

                            i_start = i + 1;

                            /**********************************************/
                            /*                                            */
                            /* Start training again.                      */
                            /*                                            */
                            /**********************************************/

                            bl_train = 0;
                        }
                        else if (vec_mag[0] > T_MAX_CG)
                        {
                            /**********************************************/
                            /*                                            */
                            /* Remove noise.                              */
                            /*                                            */
                            /**********************************************/

                            for (m = i; m < end -1; m++)
                            {
                                clrx[m] = clrx[m+1];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                    clry[b][m] = clry[b][m+1];
                            }
                            end--; /* check if this is needed */

                            i--;   /* stay & check again after noise removal */
                        }
                        status = free_2d_array ((void **) v_diff);
                        if (status != SUCCESS)
                        {
                            RETURN_ERROR ("Freeing memory: v_diff\n", 
                                          FUNC_NAME, FAILURE);
                        }
                    } /* end of continuous monitoring */ 
                }  /* end of checking basic requrirements */ 

                /**********************************************************/
                /*                                                        */
                /* Move forward to the i+1th clear observation.           */
                /*                                                        */
                /**********************************************************/

                i++;

            } /* end of "while (i <= end - CONSE) */

            /**************************************************************/
            /*                                                            */
            /* Two ways for processing the end of the time series.        */ 
            /*                                                            */
            /**************************************************************/

            if (bl_train == 1)
            {

                /**********************************************************/
                /*                                                        */
                /* If no break, find at the end of the time series,       */
                /* define probability of change based on CONSE.           */
                /*                                                        */
                /**********************************************************/

                for (i_conse = CONSE - 1; i_conse >= 0; i_conse--)
                {
                    if (vec_mag[i_conse] <= T_CG)
                    {
                        /**************************************************/
                        /*                                                */
                        /* The last stable ID.                            */
                        /*                                                */
                        /**************************************************/

                        id_last = i_conse + 1;
                        break;
                    }
                } 

                /**********************************************************/
                /*                                                        */
                /* Update change probability, end time of the curve.      */
                /*                                                        */
                /**********************************************************/

                rec_cg[*num_fc].change_prob = (CONSE - id_last) / CONSE; 
                rec_cg[*num_fc].t_end = clrx[end - CONSE + id_last];

                /**********************************************************/
                /*                                                        */
                /* Mean value fit for the rest of the pixels < CONSE & > 1*/
                /*                                                        */
                /**********************************************************/

                if (CONSE > id_last)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Update time of the probable change.                */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].t_break = clrx[end-CONSE+id_last+1];

                    /******************************************************/
                    /*                                                    */
                    /* Update magnitude of change.                        */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], id_last, CONSE-2);
                        matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, CONSE-1,
                                                       &rec_cg[*num_fc].magnitude[i_b]);
                    }
                }
            }

            else if (bl_train == 0)

            {
                /**********************************************************/
                /*                                                        */
                /* If break found close to the end of the time series,    */ 
                /* use [CONSE,MIN_NUM_C*N_TIMES+CONSE) to fit curve.      */
                /*                                                        */
                /* Update i_start.                                        */
                /*                                                        */
                /**********************************************************/
                if (*num_fc == rec_fc)
                {
                    /******************************************************/
                    /*                                                    */
                    /* First curve.                                       */
                    /*                                                    */
                    /******************************************************/

                    i_start = 1;
                }
                else
                {
                    for (k = 0; k < valid_num_scenes; k++) 
                    {
                        if (clrx[k] >= rec_cg[*num_fc-1].t_break)
                        {
                            i_start = k + 1;
			    break;
                        }
                    }
                }

                for (m = 0; m < valid_num_scenes; m++)
                {
                    bl_ids[m] = 0;
                }

                if ((end - i_start + 1) > CONSE)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Multitemporal cloud mask.                          */
                    /*                                                    */
                    /******************************************************/

                    status = auto_mask(clrx, clry, i_start-1, end-1,
                                       (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS, 
                                       adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
                    if (status != SUCCESS)
                        RETURN_ERROR("ERROR calling auto_mask at the end of time series", 
                                      FUNC_NAME, FAILURE);

                    /******************************************************/
                    /*                                                    */
                    /* Clears the IDs buffers.                            */
                    /*                                                    */
                    /******************************************************/

                    for (m = 0; m < valid_num_scenes-1; m++)
                    {
                        ids[m] = 0;
                    }

                    /******************************************************/
                    /*                                                    */
                    /* IDs to be removed.                                 */
                    /*                                                    */
                    /******************************************************/

                    for (k = i_start-1; k < end; k++)
                    {
                        ids[k-i_start+1] = k;
                    }
                    m= 0;
                    i_span = 0;
                    for (k = 0; k < end-i_start+1; k++)
                    {
                        if (bl_ids[k] == 1) 
                        {
                            rm_ids[m] = ids[k];
                            m++;
                        }
                        else
                            i_span++;  /* update i_span after noise removal */
                    }
                    rm_ids_len = m;

                    /******************************************************/
                    /*                                                    */
                    /* Remove noise pixels between i_start & i.           */
                    /*                                                    */
                    /******************************************************/

                    m = 0;
                    for (k = 0, k_new=0; k < end-i_start+1; k++)
                    {
                        if (m < rm_ids_len && k == rm_ids[m])
                        {
                            m++;
                            end--;
                            continue;
                        }
                        clrx[k_new] = clrx[k];
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                             clry[i_b][k_new] = clry[i_b][k];
                        k_new++;
                    }
                }

                if ((end - i_start + 1) >= CONSE)
                {
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        status = auto_ts_fit(clrx, clry, i_b, i_start-1, end-1, MIN_NUM_C, 
                                             fit_cft, &rmse[i_b], temp_v_dif); 
                        if (status != SUCCESS)  
                        {
                             RETURN_ERROR ("Calling auto_ts_fit at the end of time series\n", 
                                    FUNC_NAME, FAILURE);
                        }
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record time of curve start, time of curve end,     */
                    /* break time, postion of the pixel.                  */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].t_start = clrx[i_start-1];
                    rec_cg[*num_fc].t_end = clrx[end-1];
                    rec_cg[*num_fc].t_break = 0;
                    rec_cg[*num_fc].pos = (row-1) * num_samples + col_inx + 1; // offsets

                    /******************************************************/
                    /*                                                    */
                    /* Record fitted coefficients.                        */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        for (k = 0; k < MAX_NUM_C; k++)
                        {
                            rec_cg[*num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                        }
                        rec_cg[*num_fc].rmse[i_b] = rmse[i_b];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record change probability, number of observations, */
                    /* fit category.                                      */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[*num_fc].change_prob = 0.0;
                    rec_cg[*num_fc].num_obs = i_span;
                    rec_cg[*num_fc].category = 20 + MIN_NUM_C; /* simple model fit at the end */

                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        rec_cg[*num_fc].magnitude[i_b] = 0.0; 
                    }
                }
            }
        }
      } // for ncols
    } // for nrows

    // add the free calls here......

    return (SUCCESS);
}
