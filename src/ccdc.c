#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "ccdc.h"

#define LASSO_BANDS 5
#define MAX_NUM_FC 10 /* Values change with number of pixels run */
int lasso_blist[LASSO_BANDS] = {1, 2, 3, 4, 6}; /* This is band index */
int cmpfunc (const void * a, const void * b)
{
   return ( *(float*)a - *(float*)b );
}

/******************************************************************************
METHOD:  ccdc

PURPOSE:  the main routine for CCDC (Continuous Change Detection and 
          Classification) in C

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the ccdc
SUCCESS         Processing was successful

PROJECT:  Land Change Monitoring, Assessment and Projection (LCMAP) Project

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/15/2013   Song Guo         Original Development

NOTES: type ./ccdc --help for information to run the code
******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";
    char msg_str[MAX_STR_LEN];  /* input data scene name */
    //    char filename[MAX_STR_LEN];         /* input binary filenames */
    float t_cg;
    int conse;
    int status;                 /* return value from function call */
    Output_t *rec_cg = NULL;    /* output structure and metadata */
    bool verbose;               /* verbose flag for printing messages */
    int i, k, m, b;
    //    char **scene_list = NULL;
    FILE *fd;
    int num_scenes = MAX_SCENE_LIST;
    int min_num_c = 4;
    int mid_num_c = 6;
    int max_num_c = 8;
    int num_c = 8;            /* max number of coefficients for the model */
    int n_times = 3;          /* number of clear observations/coefficients*/
    int num_fc = 0;           /* intialize NUM of Functional Curves */
    int rec_fc;
    float num_yrs = 365.25;   /* number of days per year */
    int t_const = 3;          /* Threshold for cloud, shadow, and snow detection */
    int mini_yrs = 1;         /* minimum year for model intialization */
    float v_start[LASSO_BANDS];
    float v_end[LASSO_BANDS];
    float v_slope[LASSO_BANDS];
    float v_dif[LASSO_BANDS];
    float t_sn = 0.75;        /* no change detection for permanent snow pixels */ 
    float t_clr = 0.25;       /* Fmask fails threshold */
    int *sdate;
    Input_meta_t *meta;
    int row, col;
    //    int landsat_number;
    int fmask_sum = 0;
    int clr_sum = 0;
    int sn_sum = 0;
    int all_sum = 0;
    float sn_pct;
    float clr_pct;
    int n_ws = 0;
    int n_sn = 0;
    int n_clr = 0;
    int *id_clr;
    int *id_all;
    int *id_sn;
    int *clrx;
    int16 **clry;
    int *cpx;
    int16 **cpy;
    int i_start;
    int end;
    float **fit_cft;
    float *rmse;
    int i_span;
    int update_num_c;
    int bl_train;
    float time_span;
    int *bl_ids;
    int *id_range;
    int *ids;
    int *ids_old;
    int *rm_ids;
    int rm_ids_len;
    int i_rec;
    float **v_diff;
    float v_dif_norm;
    int i_count;
    float **v_dif_mag;
    int i_conse, i_b;
    float *vec_mag;
    float *vec_magg;
    float v_dif_mean;
    float **rec_v_dif;
    float **rec_v_dif_temp;
    float adj_rmse[TOTAL_BANDS-1];
    float mini_rmse;
    int bl_tmask;
    int n_rmse;
    float tmpcg_rmse[TOTAL_BANDS-1];
    int d_rt;
    float *d_yr;
    float break_mag;
    float max_v_dif;
    int id_last;
    float ts_pred_temp;
    FILE *fp_bin_out;
    int ids_old_len;
    int i_break;
    int i_ini;
    int ini_conse;
    float vec_magg_min;
    float t_max_cg = 35.8882;    /* chi-square inversed T_max_cg for 
                                    last step noise removal */
    int ids_len;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    row = 5099;
    col = 3191;
    t_cg = 15.0863;
    conse = 6;

#if 0
    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    status = get_args (argc, argv, &row, &col, &t_cg, &conse, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }
#endif
#if 0
    /* allocate memory for scene_list */
    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                         sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR ("Allocating scene_list memory", FUNC_NAME, FAILURE);
    }

    /* check if scene_list.txt file exists, if not, create the scene_list
       from existing files in the current data working directory */
    if (access("scene_list.txt", F_OK) != -1) /* File exists */
    {
        num_scenes = MAX_SCENE_LIST;
    }
    else /* File not exists */
    {
        status = create_scene_list("L*_sr_band1.hdr", num_scenes, scene_list);
        if(status != SUCCESS)
            RETURN_ERROR("Running create_scene_list file", FUNC_NAME, FAILURE); 
    }

    fd = fopen("scene_list.txt", "r");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_scenes; i++)
    {
        if (fscanf(fd, "%s", scene_list[i]) == EOF)
            break;
    }
    num_scenes = i;
#endif

    num_scenes = 455;
    /* Allocate memory */
    sdate = malloc(num_scenes * sizeof(int));
    if (sdate == NULL)
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);

    clrx = malloc(num_scenes * sizeof(int));
    if (clrx == NULL)
        RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);

    id_range = (int *)calloc(num_scenes, sizeof(int));
    if (id_range == NULL)
        RETURN_ERROR("ERROR allocating id_range memory", FUNC_NAME, FAILURE);

    id_clr = (int *)calloc(num_scenes, sizeof(int));
    if (id_clr == NULL)
        RETURN_ERROR("ERROR allocating id_clr memory", FUNC_NAME, FAILURE);

    id_all = (int *)calloc(num_scenes, sizeof(int));
    if (id_all == NULL)
        RETURN_ERROR("ERROR allocating id_all memory", FUNC_NAME, FAILURE);

    id_sn = (int *)calloc(num_scenes, sizeof(int));
    if (id_sn == NULL)
        RETURN_ERROR("ERROR allocating id_sn memory", FUNC_NAME, FAILURE);

    ids = (int *)calloc(num_scenes, sizeof(int));
    if (ids == NULL)
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);

    ids_old = (int *)calloc(num_scenes, sizeof(int));
    if (ids_old == NULL)
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);

    bl_ids = (int *)calloc(num_scenes, sizeof(int));
    if (bl_ids == NULL)
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);

    rm_ids = (int *)calloc(num_scenes, sizeof(int));
    if (rm_ids == NULL)
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);

    clry = (int16 **) allocate_2d_array (num_scenes, TOTAL_BANDS - 1, 
                                         sizeof (int16));
    if (clry == NULL)
    {
        RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
    }

    fit_cft = (float **) allocate_2d_array (max_num_c, TOTAL_BANDS - 1, 
                                         sizeof (float));
    if (fit_cft == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    rmse = malloc((TOTAL_BANDS - 1) * sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    vec_mag = malloc(conse * sizeof(float));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }

    /* allocate memory for v_dif_mag */ 
    v_dif_mag = (float **) allocate_2d_array(conse,
                         TOTAL_BANDS - 1, sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory", 
                                 FUNC_NAME, FAILURE);
    }

#if 0
    /* sort scene_list based on year & julian_day */
    status = sort_scene_based_on_year_doy(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, EXIT_FAILURE);
    }
#endif
    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    //    status = read_envi_header(scene_list[0], meta);
    status = read_envi_header("LC80460272013120LGN01_MTLstack", meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header", 
                      FUNC_NAME, EXIT_FAILURE);
    }

    if (verbose)
    {
        /* Print some info to show how the input metadata works */
        printf ("DEBUG: Number of input lines: %d\n", meta->lines);
        printf ("DEBUG: Number of input samples: %d\n", meta->samples);
        printf ("DEBUG: UL_MAP_CORNER: %d, %d\n", meta->upper_left_x,
                meta->upper_left_y);
        printf ("DEBUG: ENVI data type: %d\n", meta->data_type);
        printf ("DEBUG: ENVI byte order: %d\n", meta->byte_order);
        printf ("DEBUG: UTM zone number: %d\n", meta->utm_zone);
        printf ("DEBUG: Pixel size: %d\n", meta->pixel_size);
        printf ("DEBUG: Envi save format: %s\n", meta->interleave);
    }
#if 0
    /* Open input files */
    FILE *fp_bin[num_scenes][TOTAL_BANDS];
    short int buf[num_scenes][TOTAL_BANDS-1];
    unsigned char fmask_buf[num_scenes];
    /* Open input files */
    for (i = 0; i < num_scenes; i++)
    {
        for (k = 0; k < TOTAL_BANDS; k++)
        {
            landsat_number = atoi(sub_string(scene_list[i],2,1));
            if (landsat_number != 8)
            {
                if (k == 5)
                    sprintf(filename, "%s_toa_band6.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+1);
            }
            else
            {
                if (k == 5)
                    sprintf(filename, "%s_toa_band10.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else if (k == 6)
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+1);
                else 
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+2);
            }

            fp_bin[i][k] = open_raw_binary(filename,"rb");
            if (fp_bin[i][k] == NULL)
                printf("error open %d scene, %d bands files\n",i, k+1);
            if (k != TOTAL_BANDS-1)
            {
                fseek(fp_bin[i][k], (row * meta->samples + col)*sizeof(short int), SEEK_SET);
                if (read_raw_binary(fp_bin[i][k], 1, 1,
                    sizeof(short int), &buf[i][k]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
            }
            else
            {
                fseek(fp_bin[i][k], (row * meta->samples + col)*sizeof(unsigned char), 
                    SEEK_SET);
                if (read_raw_binary(fp_bin[i][k], 1, 1,
                    sizeof(unsigned char), &fmask_buf[i]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
            }
            close_raw_binary(fp_bin[i][k]);
        }
    }
#endif
    /* Temporary code for testing purpose */
    unsigned char fmask_buf[num_scenes];
    int buf[num_scenes][TOTAL_BANDS+1];
    fd = fopen("pixel_inputs.txt","r");
    if (fd == NULL)
        RETURN_ERROR ("Open input file", FUNC_NAME, EXIT_FAILURE);

    for (i = 0; i < num_scenes; i++)
    {
        for (i_b = 0; i_b < TOTAL_BANDS+1; i_b++)
        {
            if (i_b == 0)
            {
                fscanf(fd, "%d", &buf[i][i_b]);
                clrx[i] = buf[i][i_b];
                sdate[i] = clrx[i];  //Changed
            }
            else if (i_b == 8)
            {
                fscanf(fd, "%d\n", &buf[i][i_b]);
                fmask_buf[i] = (unsigned char)buf[i][i_b];
            }
            else if (i_b == 6)
            {
                fscanf(fd,"%d", &buf[i][i_b+1]);
                clry[i][i_b] = (int16)buf[i][i_b+1];
            }
            else if (i_b == 7)
            {
                fscanf(fd,"%d", &buf[i][i_b-1]);
                clry[i][i_b-2] = (int16)buf[i][i_b-1];
            }
            else
            {
                fscanf(fd,"%d", &buf[i][i_b]);
                clry[i][i_b-1] = (int16)buf[i][i_b];
            }
        }
    }

    /* Only run CCDC for places where more than 50% of images has data */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] < 255)
            fmask_sum++;
    }
    if (fmask_sum < (int) 0.5 * num_scenes)
    {
        RETURN_ERROR ("Not enough clear-sky pisels", FUNC_NAME, EXIT_FAILURE);
    }
    else
    {
        printf("Clear-sky pixel percentage = %f\n", 
               (float)fmask_sum / (float)num_scenes);
    }

    /* CHANGE: need change back to 0-6 from 1-7 if original 
       inputs are used */

    /* pixel value ranges should follow physical rules */
    for (i = 0; i < num_scenes; i++)
    { 
        if ((buf[i][1] > 0) && (buf[i][1] < 10000) &&
            (buf[i][2] > 0) && (buf[i][2] < 10000) &&
            (buf[i][3] > 0) && (buf[i][3] < 10000) &&
            (buf[i][4] > 0) && (buf[i][4] < 10000) &&
            (buf[i][5] > 0) && (buf[i][5] < 10000) &&
            (buf[i][6] > -9320) && (buf[i][6] < 7070) &&
            (buf[i][7] > 0) && (buf[i][7] < 10000))
            id_range[i] = 1;
        else
            id_range[i] = 0;
    }

    /* Get each mask pixel totals */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] < 2)
            clr_sum++;
        if (fmask_buf[i] < 255)
            all_sum++;
        if (fmask_buf[i] == 3)
            sn_sum++;
    }


    /* percent of clear pixels */
    clr_pct = (float) clr_sum / (float) all_sum;

    /* percent of snow observations */
    sn_pct =  (float) sn_sum / (float)(clr_sum + sn_sum); 

    printf("clr_sum,all_sum,clr_pct=%d,%d,%f\n",clr_sum,all_sum,clr_pct);

    /* Allocate memory for rec_cg */ 
    rec_cg = malloc(MAX_NUM_FC * sizeof(Output_t));
    if (rec_cg == NULL)
        RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);

    /* fit permanent snow observations */
    if (clr_pct < t_clr)
    {
        if (sn_pct > t_sn)
        {
            if (n_sn < n_times * min_num_c ) /* not enough snow pixels */
            {
                RETURN_ERROR ("Not enough good snow observations\n", 
                     FUNC_NAME, EXIT_FAILURE);
            }
            else
            {
                /* start model fit for snow persistent pixels */
                printf ("Fit permanent snow observations, now pixel = %f\n", 
                       100.0 * sn_pct); 

                /* snow observations are "good" now */
                for (i = 0; i < num_scenes; i++)
                { 
                        if (fmask_buf[i] == 3)
                    {
                        clrx[n_sn] = sdate[i];
                        for (k = 0; k < TOTAL_BANDS - 1; k++)
                            clry[n_sn][k] = buf[i][k];
                        n_sn++;
                    }
                }   
            }
            end = n_sn;

            /* the first observation for TSFit */
            i_start = 1; /* the first observation for TSFit */
#if 0
            /* identified and move on for the next curve */
            num_fc++; /* not needed, as array index starts with zero */
#endif
            /* treat saturated and unsaturated pixels differently */
            for (k = 0; k < TOTAL_BANDS -1; k++)
            {
                i_span = 0;
                if (k != TOTAL_BANDS - 3) /* for optical bands */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > 0 && clry[i][k] < 10000)
                        {
                            clrx[i_span] = sdate[i];
                            for (k = 0; k < TOTAL_BANDS - 1; k++)
                                 clry[i_span][k] = buf[i][k];
                                i_span++;
                        }
                    }

                    if (i_span < min_num_c * n_times)
                        fit_cft[i][k] = 10000; /* fixed value for saturated pixels */
                    else
                    {
                        status = auto_ts_fit(clrx, clry, k, 0, i_span-1, min_num_c, 
                                 fit_cft, &rmse[k]); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit1\n", 
                                   FUNC_NAME, EXIT_FAILURE);

                    } 
                }
                else /* for thermal band */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > -9300 && clry[i][k] < 7070)
                        {
                            clrx[i_span] = sdate[i];
                            for (k = 0; k < TOTAL_BANDS - 1; k++)
                                 clry[i_span][k] = buf[i][k];
                            i_span++;
                        }
                            
                        status = auto_ts_fit(clrx, clry, k, 0, i_span-1, min_num_c, fit_cft, 
                                 &rmse[k]); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit2\n", 
                                 FUNC_NAME, EXIT_FAILURE);
                    }
                }
            }

            /* update information at each iteration */
            /* record time of curve start */
            rec_cg[num_fc].t_start = clrx[i_start-1]; 
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end-1]; 
            /* no break at the moment */
            rec_cg[num_fc].t_break = 0; 
            /* record postion of the pixel */
            rec_cg[num_fc].pos.row = row; 
            /* record postion of the pixel */
            rec_cg[num_fc].pos.col = col; 
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                for (k = 0; k < min_num_c; k++)
                    /* record fitted coefficients */
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                /* record rmse of the pixel */
                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0; 
            /* record number of observations */
            rec_cg[num_fc].num_obs = n_ws; 
            /* record fit category */
            rec_cg[num_fc].category = 50 + min_num_c; /* snow pixel */
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                /* record change magnitude */ 
                rec_cg[num_fc].magnitude[i_b] = 0.0; 
            /* NUM of Fitted Curves (num_fc) */
            num_fc++;   
        }
        else
        {
            /* no change detection for clear observations */
            if (clr_sum < (n_times * min_num_c)) /* not enough snow pixels */ 
            {
                RETURN_ERROR("Not enough good clear observations\n", 
                            FUNC_NAME, EXIT_FAILURE);
            }
            else
            {
                /* start model fit for snow persistent pixels */
                printf ("Fmask failed, clear pixel = %f\n", 
                       100.0 * clr_pct); 

                /* snow observations are "good" now */
                for (i = 0; i < num_scenes; i++)
                { 
                    if ((fmask_buf[i] == 0 || fmask_buf[i] == 1) && id_range[i] == 1)
                    {
                        clrx[n_sn] = sdate[i];
                        for (k = 0; k < TOTAL_BANDS - 1; k++)
                             clry[n_sn][k] = buf[i][k];
                        n_clr++;
                    }   
                }
                end = n_clr;

                /* the first observation for TSFit */
                i_start = 1; /* the first observation for TSFit */
#if 0
                /* identified and move on for the next curve */
                num_fc++; /* not needed, as array index starts with zero */
#endif
                /* Do lasso regression */
                for (k = 0; k < TOTAL_BANDS -1; k++)
                {
                    status = auto_ts_fit(clrx, clry, k, 0, num_scenes-1, min_num_c, 
                             fit_cft, &rmse[k]); 
                    if (status != SUCCESS)  
                        RETURN_ERROR ("Calling auto_ts_fit3\n", 
                               FUNC_NAME, EXIT_FAILURE);
                }

                /* update information at each iteration */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[i_start-1]; 
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end-1]; 
                /* no break at the moment */
                rec_cg[num_fc].t_break = 0; 
                /* record postion of the pixel */
                rec_cg[num_fc].pos.row = row; 
                /* record postion of the pixel */
                rec_cg[num_fc].pos.col = col; 
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                {
                    for (k = 0; k < min_num_c; k++)
                        /* record fitted coefficients */
                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    /* record rmse of the pixel */
                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                }
                /* record change probability */
                rec_cg[num_fc].change_prob = 0.0; 
                /* record number of observations */
                rec_cg[num_fc].num_obs = n_ws; 
                /* record fit category */
                rec_cg[num_fc].category = 40 + min_num_c; /* snow pixel */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    /* record change magnitude */ 
                    rec_cg[num_fc].magnitude[i_b] = 0.0; 
                /* NUM of Fitted Curves (num_fc) */
                num_fc++;   
            }
        }
    }
    else /* clear land or water pixels */
    {
        printf("Seasonal Snow (Snow < %f)\n", 100.0 * sn_pct);
        printf("Fmask works, clear pixels (land/water) = %f\n", 100.0 * clr_pct);

        for (i = 0; i < num_scenes; i++)
        { 
            if ((fmask_buf[i] < 2) && id_range[i] == 1)
            {
                clrx[n_clr] = sdate[i];
                for (k = 0; k < TOTAL_BANDS - 1; k++)
                     clry[n_clr][k] = clry[i][k];
                n_clr++;
            }   
        }
        end = n_clr;

        printf("end_clr=%d\n",end);
        /* calculate median variogram */
        status = median_variogram(clry, 0, end-1, TOTAL_BANDS-1, adj_rmse);
        if (status != SUCCESS)
            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME, 
                         FAILURE);
#if 0
        for (k = 0; k < TOTAL_BANDS-1; k++)
         printf("k,adj_rmse[k]=%d,%f\n",k,adj_rmse[k]);
#endif
        /* start with mininum requirement of clear obs */
        i = n_times * min_num_c;

        /* the first observation for TSFit */
        i_start = 1; 

        /* record the start of the model initialization (0=>initial;1=>done) */
        bl_train = 0;

        /* record the num_fc at the beginning of each pixel */
        rec_fc = num_fc;

        /* record the start of Tmask (0=>initial;1=>done) */
        bl_tmask = 0;

        /* while loop - process till the last clear observation - conse */
        while (i <= end - conse)
        {
            /* Allocate memory for rec_v_dif */
            rec_v_dif = (float **)allocate_2d_array(i-i_start+1, TOTAL_BANDS - 1,
                                                            sizeof (float));
            if (rec_v_dif == NULL)
                RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);

            /* span of "i" */
            i_span = i - i_start + 1;

            printf("end,conse=%d,%d\n",end,conse);
            printf("i_start,i,i_span=%d,%d,%d\n",i_start,i,i_span);

            /* span of time (num of years) */
            time_span = (float)(clrx[i-1] - clrx[i_start-1]) / num_yrs;

            /* basic requrirements: 1) enough observations; 2) enough time */
            if (i_span >= n_times * min_num_c && time_span >= (float)mini_yrs)

            /* initializing model */
            if (bl_train == 0)
            {
                /* step 1: noise removal */ 
                status = auto_mask(clrx, clry, i_start-1, i+conse-1,
                                   (float)(clrx[i+conse-1]-clrx[i_start-1])/num_yrs, 
                                   adj_rmse[1], adj_rmse[4], t_const, bl_ids);
                if (status != SUCCESS)
                    RETURN_ERROR("ERROR calling auto_mask routine", 
                                  FUNC_NAME, FAILURE);
#if 0
                for(k = i_start-1; k < i+conse; k++)
                     printf("k,bl_ids[k]=%d,%d\n",k,bl_ids[k]);
#endif
                /* Clears the IDs buffers */
                for (k = 0; k < num_scenes; k++)
                    ids[k] = 0;

                /* IDs to be removed */
                for (k = i_start-1; k < i+conse; k++)
                    ids[k-i_start+1] = k;
                m= 0;
                i_span = 0;
                for (k = 0; k < i-i_start+1; k++)
                {
                    if (bl_ids[k] == 1) 
                    {
                        rm_ids[m] = ids[k];
                        printf("m,rm_ids[m]=%d,%d\n",m,rm_ids[m]);
                        m++;
                    }
                    else
                        i_span++;  /* update i_span after noise removal */
                }

                rm_ids_len = m;
                /* update i_span after noise removal */
                // i_span = i - i_start +1 - rm_ids_len;

                printf("rm_ids_len, update i_span = %d,%d\n", rm_ids_len, i_span);
                /* check if there is enough observation */
                if (i_span < (n_times * min_num_c))
                {
                    /* move forward to the i+1th clear observation */
                    i++;
                    printf("i++=%d\n",i);
                    /* not enough clear observations */
                    continue;
                }
                else
                {
                    /* allocate memory for cpx, cpy */
                    cpx = malloc(end * sizeof(int));
                    if (cpx == NULL)
                        RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                    cpy = (int16 **) allocate_2d_array (end, TOTAL_BANDS - 1, 
                                         sizeof (int16));
                    if (cpy == NULL)
                    {
                        RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                    }

                    /* copy clrx and clry to cpx and cpy */
                    for (k = 0; k < end; k++)
                    {
                        cpx[k] = clrx[k];
                        for (m = 0; m < TOTAL_BANDS - 1; m++)
                            cpy[k][m] = clry[k][m];
                    }

                    printf("end1,rm_ids_len=%d,%d\n",end,rm_ids_len);
                    /* remove noise pixels between i_start & i */
                    for (m = 0; m < rm_ids_len; m++)
                    {
                        printf("rm_ids[m]=%d\n",rm_ids[m]);
                        if (m != 0)
                            rm_ids[m]--;
                        for (k = rm_ids[m]; k < end-1; k++)
                        {
                            cpx[k] = cpx[k+1];
                            for (b = 0; b < TOTAL_BANDS - 1; b++)
                                cpy[k][b] = cpy[k+1][b];
                        }
                        end--;
                    }

                    printf("end2=%d\n",end);
                    /* record i before noise removal 
                       This is very important as if model is not initialized 
                       the multitemporal masking shall be done again instead 
                       of removing outliers in every masking */
                    i_rec=i;

                    /* update i afer noise removal (i_start stays the same) */
                    i=i_start + i_span - 1;
                    printf("update i = %d\n",i);

                    /* update span of time (num of years) */
                    time_span=(cpx[i-1] - cpx[i_start-1]) / num_yrs;

                    /* check if there is enough time */
                    if (time_span < mini_yrs)
                    {
                        i = i_rec;   /* keep the original i */
                        /* move forward to the i+1th clear observation */
                        i++;        
                        printf("i++2=%d\n",i);
                        free(cpx);
                        status = free_2d_array ((void **) cpy);
                        if (status != SUCCESS)
                              RETURN_ERROR ("Freeing memory: cpy\n", 
                                        FUNC_NAME, EXIT_FAILURE);
                        continue;    /* not enough time */
                    }
                    else
                    {
                        /* remove noise */
                        for (k = 0; k < end; k++)
                        {
                            clrx[k] = cpx[k];
                            for (m = 0; m < TOTAL_BANDS - 1; m++)
                                clry[k][m] = cpy[k][m];
                        }

                        free(cpx);
                        status = free_2d_array ((void **) cpy);
                        if (status != SUCCESS)
                              RETURN_ERROR ("Freeing memory: cpy\n", 
                                        FUNC_NAME, EXIT_FAILURE);
#if 0
                        /* record the start of Tmask (0=>initial;1=>done) */
                        bl_tmask = 0;
#endif
                        /* Step 2: model fitting: initialize model testing variables
                           defining computed variables */
                        v_dif_norm = 0.0;
                        for (b = 0; b < TOTAL_BANDS-1; b++)
                        {
                            /* Initial model fit */
                            status = auto_ts_fit_full(clrx, clry, b, i_start-1, i-1, 
                                     min_num_c, fit_cft, &rmse[b], rec_v_dif); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit4\n", 
                                     FUNC_NAME, EXIT_FAILURE);
#if 0
                            for(k = 0; k < NUM_COEFFS; k++)
                             printf("lasso_blist[b],k,fit_cft[k][lasso_blist[b]],rmse[lasso_blist[b-1]]=%d,%d,%f,%f\n",lasso_blist[b],k,fit_cft[k][lasso_blist[b]],rmse[lasso_blist[b]]);
#endif
                        }
#if 0
printf("conse=%d\n",conse);
#endif
                        for(b = 0; b < LASSO_BANDS; b++)
                        {
                            /* calculate mini rmse */
                            mini_rmse = max(adj_rmse[lasso_blist[b]], rmse[lasso_blist[b]]);

                            /* compare the first clear obs */
                            auto_ts_predict(clrx, fit_cft, lasso_blist[b], i_start-1, i_start-1, 
                                            &ts_pred_temp);  
#if 0
printf("conse1.5=%d,%d\n",lasso_blist[b],conse);
 printf("clry[i_start-1][lasso_blist[b]],ts_pred_temp,mini_rmse=%d,%f,%f\n",
        clry[i_start-1][lasso_blist[b]],ts_pred_temp,mini_rmse);
#endif
                            v_start[lasso_blist[b]] = (clry[i_start-1][lasso_blist[b]] -
                               ts_pred_temp)/mini_rmse;
#if 0
printf("v_start[lasso_blist[b]]=%f\n",v_start[lasso_blist[b]]);
printf("conse2=%d,%d\n",lasso_blist[b],conse);
#endif
                            /* compare the last clear observation */
                            auto_ts_predict(clrx, fit_cft, lasso_blist[b], i-1, i-1, &ts_pred_temp);
                            v_end[lasso_blist[b]] = (clry[i-1][lasso_blist[b]]-
                                                     ts_pred_temp)/mini_rmse;

                            /* anormalized slope values */
                            v_slope[lasso_blist[b]] = fit_cft[1][lasso_blist[b]] *
                                    (clrx[i-1]-clrx[i_start-1])/mini_rmse;
                            
                            /* differece in model intialization */
                            v_dif[lasso_blist[b]] = fabs(v_slope[lasso_blist[b]]) + 
                                                    fabs(v_start[lasso_blist[b]]) + 
                                                    fabs(v_end[lasso_blist[b]]);  
                            v_dif_norm += v_dif[lasso_blist[b]] * v_dif[lasso_blist[b]];               
                        }
                        printf("v_dif_norm=%f\n",v_dif_norm);
#if 0
printf("conse3=%d\n",conse);
#endif
                        /* find stable start for each curve */
                        if (v_dif_norm > t_cg)
                        {
                            /* start from next clear obs */
                            i_start++;
                            printf("i_start1=%d\n",i_start);
                            /* move forward to the i+1th clear observation */
                            i++;
                            printf("i_start++,i++3=%d,%d\n",i_start,i);

                            /* keep all data and move to the next obs */
                            continue;
                        }
                        else
                        {
                            /* model ready */
                            bl_train = 1;

                            /* count difference of i for each iteration */
                            i_count = 0;

                            /* find the previous break point */
                            if (num_fc == rec_fc)
                                i_break = 1; /* first curve */
                            else
                            {
                                for (k = 0; k < num_scenes; k++) 
                                {
                                    if (clrx[k] > rec_cg[num_fc-1].t_end)
                                     i_break = k + 1;
                                }
                            printf("i_start,i_break2=%d,%d\n",i_start,i_break);
                            }

                            if (i_start > i_break)
                            {
                                /* model fit at the beginning of the time series */
                                for(i_ini = i_start-2; i_ini >= i_break-1; i_ini--)
                                {
                                    if ((i_start - i_break) < conse)
                                        ini_conse = i_start - i_break;
                                    else
                                        ini_conse = conse;
                            
                                    /* allocate memory for model_v_dif */ 
                                    v_dif_mag = (float **) allocate_2d_array(ini_conse,
                                              TOTAL_BANDS-1, sizeof (float));
                                    if (v_dif_mag == NULL)
                                    {
                                        RETURN_ERROR ("Allocating v_dif_mag memory", 
                                                      FUNC_NAME, FAILURE);
                                    }

                                    /* allocate memory for v_diff */ 
                                    v_diff = (float **) allocate_2d_array(ini_conse,
                                              LASSO_BANDS, sizeof (float));
                                    if (v_diff == NULL)
                                    {
                                        RETURN_ERROR ("Allocating v_diff memory", 
                                                      FUNC_NAME, FAILURE);
                                    }
 
                                    /* allocate memory for vec_magg */ 
                                    vec_magg = (float *) malloc(ini_conse * sizeof (float));
                                    if (vec_magg == NULL)
                                    {
                                        RETURN_ERROR ("Allocating vec_magg memory", 
                                                      FUNC_NAME, FAILURE);
                                    }
   
                                    /* detect change. 
                                       value of difference for conse obs
                                       record the magnitude of change */
                                    printf("ini_conse,i=%d,%d\n",ini_conse,i);
                                    vec_magg_min = 9999.0;
                                    for (i_conse = 0; i_conse < ini_conse-1; i_conse++)
                                    {
                                        v_dif_norm = 0.0;
                                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                                        {
                                            /* absolute differences */
                                            auto_ts_predict(clrx, fit_cft, i_b, i_ini-i_conse,
                                                            i_ini-i_conse, &ts_pred_temp);
                                            v_dif_mag[i_conse][i_b] = (float)clry[i_ini-i_conse][i_b] - 
                                                    ts_pred_temp;

                                            /* normalize to z-score */
                                            for (b = 0; b < LASSO_BANDS; b++)
                                            {
                                                if (i_b == lasso_blist[b])
                                                {
                                                    /* minimum rmse */ 
                                                    mini_rmse = max(adj_rmse[i_b], rmse[i_b]);

                                                    /* z-scores */
                                                    v_diff[i_conse][i_b] = v_dif_mag[i_conse][i_b] 
                                                                       / mini_rmse;
                                                    v_dif_norm += v_diff[i_conse][i_b] * v_diff[i_conse][i_b];
                                                }
                                            }
                                        }
                                        vec_magg[i_conse] = v_dif_norm; 

                                        printf("i_conse,vec_magg[i_conse]1 = %d,%f\n",i_conse,vec_magg[i_conse]);

                                        if (vec_magg_min < vec_magg[i_conse])
                                            vec_magg_min =  vec_magg[i_conse];
                                    }

                                    /* change detection */
                                    if (vec_magg_min > t_cg) /* change detected */
                                        break;
                                    else if (vec_magg[0] > t_max_cg) /*flase change */
                                    {
                                        for (k = i_ini; k < end; k++)
                                        {
                                            clrx[k] = clrx[k+1];
                                            for (b = 0; b < TOTAL_BANDS-1; b++)
                                                clry[k][b] = clry[k+1][b];
                                        }
                                        i--;
                                    }

                                    /* free the memory */
                                    free(vec_magg);
                                    status = free_2d_array ((void **) v_diff);
                                    if (status != SUCCESS)
                                         RETURN_ERROR ("Freeing memory: v_diff\n", 
                                             FUNC_NAME, EXIT_FAILURE);

                                    /* update i_start if i_ini is not a confirmed break */
                                    i_start = i_ini;
                                    printf("i_start=%d\n",i_start);

                                    /* update curves */
                                    for (i_b = 0; i_b < TOTAL_BANDS-1; i_b++)
                                    {
                                        status = auto_ts_fit(clrx, clry, b, i_start-1, i-1, 
                                                 min_num_c, fit_cft, &rmse[b]); 
                                        if (status != SUCCESS)  
                                            RETURN_ERROR ("Calling auto_ts_fit5\n", 
                                                FUNC_NAME, EXIT_FAILURE);
                                    }
                                }
                            }

printf("i_start,i_break,conse3=%d,%d,%d\n",i_start,i_break,conse);
                            /* enough to fit simple model and confirm a break */
                            if ((i_start - i_break) >= conse)
                            {
printf("i_start,i_break4=%d,%d\n",i_start,i_break);
                                /* defining computed variables */
                                for (i_b = 0; i_b < TOTAL_BANDS -1; i_b++)
                                {
                                    status = auto_ts_fit(clrx, clry, i_b, i_break-1, i_start-2, 
                                             min_num_c, fit_cft, &rmse[i_b]); 
                                    if (status != SUCCESS)  
                                        RETURN_ERROR ("Calling auto_ts_fit6\n", 
                                                 FUNC_NAME, EXIT_FAILURE);
                                }

                                /* record time of curve start */
                                rec_cg[num_fc].t_start = clrx[i_break-1]; 
                                /* record time of curve end */
                                rec_cg[num_fc].t_end = clrx[i_start-2]; 
                                /* record postion of the pixel */
                                rec_cg[num_fc].pos.row = row; 
                                /* record postion of the pixel */
                                rec_cg[num_fc].pos.col = col; 
                                for (i_b = 0; i_b < TOTAL_BANDS -1; i_b++)
                                {
                                    for (k = 0; k < min_num_c; k++)
                                    /* record fitted coefficients */
                                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                    /* record rmse of the pixel */                         
                                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                }

                                /* treat first curve different */
                                if (num_fc == rec_fc)                           
                                {
                                    /* record fit category */
                                    rec_cg[num_fc].category = 10 + min_num_c;
                                    /* record break time */
                                    rec_cg[num_fc].t_break = clrx[i_start-1];
                                    /* record change probability */
                                    rec_cg[num_fc].change_prob = 1.0;
                                }
                                else
                                {
                                    /* record fit category */
                                    rec_cg[num_fc].category = 30 + min_num_c;
                                    /* record break time */
                                    rec_cg[num_fc].t_break = 0;
                                    /* record change probability */
                                    rec_cg[num_fc].change_prob = 0.0;
                                } 
                                /* record number of observations */
                                rec_cg[num_fc].num_obs = i_start - i_break; 
                                /* record change magnitude */
                                for (i_b = 0; i_b < TOTAL_BANDS -1; i_b++)
                                {
                                    matlab_2d_array_mean(v_dif_mag, i_b, ini_conse, 
                                                         &v_dif_mean);  
                                    rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                                }
                                /* identified and move on for the next functional curve */
                                num_fc++;                                      
                            }
                            else if (i_start > i_break)
                            {
                                /* median value fit for the rest of the pixels < conse */
                                for (i_b = 0; i_b< TOTAL_BANDS - 1; i_b++)
                                {
                                    matlab_int_2d_partial_mean(clry, i_b, i_break-1, i_start-1, 
                                              &fit_cft[0][i_b]);
                                    partial_square_root_mean(clry, i_b, i_break-1, i_start-1, 
                                              fit_cft, &rmse[i_b]);
                                }
                                /* record time of curve start */
                                rec_cg[num_fc].t_start = clrx[i_break-1]; 
                                /* record time of curve end */  
                                rec_cg[num_fc].t_end = clrx[i_start-2]; 
                                /* record fitted coefficients */
                                rec_cg[num_fc].pos.row = row; 
                                /* record fitted coefficients */
                                rec_cg[num_fc].pos.col = col; 
                                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                                {
                                    for (k = 0; k < max_num_c; k++)
                                        /* record fitted coefficients */
                                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                    /* record rmse of the pixel */
                                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                }
                                /* treat first curve different */
                                if (num_fc == rec_fc)                           
                                {
                                    /* record fit category */
                                    rec_cg[num_fc].category = 10 + 1;
                                    /* record break time */
                                    rec_cg[num_fc].t_break = clrx[i_start-1];
                                    /* record change probability */
                                    rec_cg[num_fc].change_prob = (float)(i_start-i_break)/(float)conse;
                                }
                                else
                                {
                                    /* record fit category */
                                    rec_cg[num_fc].category = 30 + 1;
                                    /* record break time */
                                    rec_cg[num_fc].t_break = 0;
                                    /* record change probability */
                                    rec_cg[num_fc].change_prob = 0.0;
                                } 
                                /* record number of observations */
                                rec_cg[num_fc].num_obs = i_start-i_break; 
                                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                                {
                                    matlab_2d_array_mean(v_dif_mag, i_b, ini_conse, &v_dif_mean);
                                    /* record change magnitude */ 
                                    rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                                }
                                    /* NUM of Fitted Curves (num_fc) */
                                    num_fc++;
                            }
                        }   
                    }  /* end  of "if (time_span < mini_yrs)" */           
                } /* end of "if (i_span < n_times*min_num_c)" */  
            } /* end of initializing model */

            /* allocate memory for v_diff */ 
            v_diff = (float **) allocate_2d_array(conse,
                              TOTAL_BANDS - 1, sizeof (float));
            if (v_diff == NULL)
            {
                RETURN_ERROR ("Allocating v_diff memory", 
                              FUNC_NAME, FAILURE);
            }

            /* continuous monitoring started!!! */
            if (bl_train == 1)
            {
                /* Clears the IDs buffers */
                for (k = 0; k < num_scenes; k++)
                    ids[k] = 0;

                /* all IDs */
                for (k = i_start-1; k < i; k++)
                {
                    ids[k-i_start+1] = k;
                }
                i_span = i - i_start +1;

                printf("i_start, i, up-update i_span=%d,%d,%d\n",
                        i_start, i, i_span);

                /* determine the time series model */
                update_cft(i_span, n_times, min_num_c, mid_num_c, max_num_c, 
                           num_c, &update_num_c);

                /* dynamic model fit when there are not many obs */
                get_ids_length(ids_old, 0, num_scenes-1, &ids_old_len);
                printf("i_count,i_span,ids_old_len=%d,%d,%d\n",i_count,i_span,ids_old_len);
                if (i_count == 0 || i_span < (n_times * max_num_c))
                {
                    /* update i_count at each interation */
                    i_count = clrx[i-1] - clrx[i_start-1];

                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        status = auto_ts_fit_full(clrx, clry, i_b, i_start-1, i-1, update_num_c, 
                                             fit_cft, &rmse[i_b], rec_v_dif); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit7\n", 
                                     FUNC_NAME, EXIT_FAILURE);
#if 0
                            for(k = i_start-1; k < i-1; k++)
                             printf("i_b,rec_v_dif[k][i_b]1=%d,%f\n",i_b,rec_v_dif[k][i_b]);
#endif
                    }

                    /* updating information for the first iteration */
                    /* record time of curve start */
                    rec_cg[num_fc].t_start = clrx[i_start-1]; 
                    /* record time of curve end */
                    rec_cg[num_fc].t_end = clrx[i-1]; 
                    /* no break at the moment */
                    rec_cg[num_fc].t_break = 0; 
                    /* record postion of the pixel */
                    rec_cg[num_fc].pos.row = row; 
                    /* record postion of the pixel */
                    rec_cg[num_fc].pos.col = col; 
                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        for (k = 0; k < update_num_c; k++)
                            /* record fitted coefficients */
                            rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                            /* record rmse of the pixel */
                            rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                    }
                    /* record change probability */
                    rec_cg[num_fc].change_prob = 0.0; 
                    /* record number of observations */
                    rec_cg[num_fc].num_obs = i-i_start+1; 
                    /* record fit category */
                    rec_cg[num_fc].category = 0 + update_num_c; 
                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        /* record change magnitude */
                        rec_cg[num_fc].magnitude[i_b] = 0.0; 
                    }
                           
                    /* detect change, value of difference for conse obs */
                    for (i_conse = 0; i_conse < conse; i_conse++)
                    {
                        v_dif_norm = 0.0;
                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        {
                            /* absolute differences */
                            auto_ts_predict(clrx, fit_cft, i_b, i+i_conse-1, i+i_conse-1, 
                                &ts_pred_temp);
                            v_dif_mag[i_conse][i_b] = (float)clry[i+i_conse-1][i_b] - ts_pred_temp; 
       
                           /* normalize to z-score */
                            for (b = 0; b < LASSO_BANDS; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /* minimum rmse */ 
                                    mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
 
                                    /* z-scores */
                                    v_diff[i_conse][i_b] = v_dif_mag[i_conse][i_b] / mini_rmse;
                                    v_dif_norm += v_diff[i_conse][i_b] * v_diff[i_conse][i_b];
#if 0
                                    printf("i_conse,i_b,clry[i+i_conse][i_b],ts_pred_temp,mini_rmse,v_diff[i_conse][i_b]=%d,%d,%d,%f,%f,%f\n",i_conse,i_b,clry[i+i_conse][i_b],ts_pred_temp,mini_rmse,v_diff[i_conse][i_b]);
#endif
                                }
                            }
                        }
                        vec_mag[i_conse] = v_dif_norm;

                        printf("i_conse,vec_mag[i_conse]2=%d,%f\n",i_conse,vec_mag[i_conse]);

                    }
                    get_ids_length(ids_old, 0, num_scenes-1, &ids_old_len);
                    get_ids_length(ids, 0, num_scenes-1, &ids_len);
                    printf("length(ids)1,length(ids_old)1,i-i_start+1=%d,%d,%d\n",
                    ids_len,ids_old_len,i-i_start+1);

                    /* Clears the IDsOld buffers */
                    for (k = 0; k < num_scenes; k++)
                        ids_old[k] = 0;

                    /* IDs that haven't been updated */
                    for (k = 0; k < num_scenes-1; k++)
                        ids_old[k] = ids[k];
                }
                else
                {
                    if ((clrx[i-1] - clrx[i_start-1]) >= (i_count + i_count / 3))
                    {
                        /* update i_count at each interation year */
                        i_count = clrx[i-1] - clrx[i_start-1];

                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        {
                            status = auto_ts_fit_full(clrx, clry, i_b, i_start-1, i-1, update_num_c, 
                                             fit_cft, &rmse[i_b], rec_v_dif); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit8\n", 
                                     FUNC_NAME, EXIT_FAILURE);
#if 0
                            for(k = i_start-1; k < i-1; k++)
                             printf("i_b,rec_v_dif[k][i_b]2=%d,%f\n",i_b,rec_v_dif[k][i_b]);
#endif
                        }

                        /* record fitted coefficients */
                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        {
                            for (k = 0; k < update_num_c; k++)
                                /* record fitted coefficients */
                                rec_cg[num_fc].coefs[i][k] = fit_cft[i_b][k]; 
                                /* record rmse of the pixel */
                                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                        }
                        /* record number of observations */
                        rec_cg[num_fc].num_obs = i-i_start+1; 
                        /* record fit category */
                        rec_cg[num_fc].category = 0 + update_num_c; 

                        get_ids_length(ids_old, 0, num_scenes-1, &ids_old_len);
                        get_ids_length(ids, 0, num_scenes-1, &ids_len);
                        printf("length(ids)2,length(ids_old)2,i-i_start+1=%d,%d,%d\n",
                        ids_len,ids_old_len,i-i_start+1);

                        /* Clears the IDsOld buffers */
                        for (k = 0; k < num_scenes; k++)
                            ids_old[k] = 0;

                        /* IDs that haven't been updated */
                        for (k = 0; k < num_scenes; k++)
                            ids_old[k] = ids[k];
                    }

                    /* record time of curve end */
                    rec_cg[num_fc].t_end = clrx[i-1]; /* record time of curve end */

                    printf("i_start, i, k=%d,%d,%d\n",i_start, i,k);
                    get_ids_length(ids_old, 0, num_scenes-1, &ids_old_len);
                    printf("ids_old_len3, i- i_start +1=%d,%d\n",ids_old_len,i - i_start +1);
#if 0
                    /* use temporally-adjusted RMSE */
                    if (ids_old_len <= n_times * max_num_c)
                    {
                        /* number of observations for calculating RMSE */
                        n_rmse = ids_old_len;
                        for (b = 0; b < TOTAL_BANDS-1; b++)
                            tmpcg_rmse[b] = rmse[b];
                    }
                    else
                    {
#endif
                    /* use fixed number for RMSE computing */
                    n_rmse = n_times * max_num_c;

                    /* better days counting for RMSE calculating */
                    /* relative days distance */
                    d_yr = malloc(ids_old_len * sizeof(float));
                    if (d_yr == NULL)
                        RETURN_ERROR ("Allocating d_yr memory", 
                                      FUNC_NAME, FAILURE);
                    rec_v_dif_temp = (float **)allocate_2d_array(ids_old_len,  
                                      LASSO_BANDS, sizeof(float));
                    if (rec_v_dif_temp == NULL)
                        RETURN_ERROR ("Allocating rec_v_dif_temp memory", 
                                      FUNC_NAME, FAILURE);
 
                    for(m = 0; m < ids_old_len; m++)
                    {
#if 0
                     printf("m,ids_old[m],i+conse-1=%d,%d,%d\n",m,ids_old[m],i+conse-1);
#endif
                        d_rt = clrx[ids_old[m]] - clrx[i+conse-1]; 
                        d_yr[m] = fabs(round((float)d_rt/num_yrs)*num_yrs - (float)d_rt);
#if 0
                        printf("m,d_rt,d_yr[m]=%d,%d,%f\n",m,d_rt,d_yr[m]);
#endif
                    }

                    get_ids_length(ids, 0, num_scenes-1, &ids_len);
                    get_ids_length(ids_old, 0, num_scenes-1, &ids_old_len);
                    printf("ids_len,ids_old_len=====%d,%d\n",ids_len,ids_old_len);

                    /* sort the d_yr */
                    qsort(d_yr, ids_old_len, sizeof(float), cmpfunc);

                    for(b = 0; b < TOTAL_BANDS-1; b++)
                        tmpcg_rmse[b] = 0.0;

                    for(m = 0; m < n_rmse; m++)
                    {
                        for (b = 0; b < LASSO_BANDS; b++)
                        {
                            rec_v_dif_temp[m][lasso_blist[b]] = rec_v_dif[ids_old[m] - 
                                     ids_old[0]][lasso_blist[b]];
                            tmpcg_rmse[lasso_blist[b]] += rec_v_dif_temp[m][lasso_blist[b]] *
                                     rec_v_dif_temp[m][lasso_blist[b]]; 
                        }
                    }

                    /* temporarily changing RMSE */
                    for (b = 0; b < LASSO_BANDS; b++)
                    {
                        tmpcg_rmse[lasso_blist[b]] = sqrt(tmpcg_rmse[lasso_blist[b]]) / 
                                sqrt(n_rmse - update_num_c);  
                    }

                    /* free allocated memories */
                    free(d_yr);
                    status = free_2d_array ((void **) rec_v_dif_temp);
                    if (status != SUCCESS)
                        RETURN_ERROR ("Freeing memory: rec_v_dif_temp\n", 
                                  FUNC_NAME, EXIT_FAILURE);

                    /* move the ith col to i-1th col */
                    for (b = 0; b < TOTAL_BANDS-1; b++)
                    {
                        for (m = 0; m < conse-1; m++)
                        {
                            v_diff[m][b] = v_diff[m+1][b];
                            v_dif_mag[m][b] = v_dif_mag[m+1][b];
                            vec_mag[m] = vec_mag[m+1];
                        }
                        v_diff[conse-1][b] = 0.0;
                        v_dif_mag[conse-1][b] = 0.0;
                    }
                    vec_mag[conse-1] = 0.0;

                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        /* absolute difference for all bands */
                        auto_ts_predict(clrx, fit_cft, i_b, i+conse-1, 
                             i+conse-1, &ts_pred_temp);
                        v_dif_mag[conse-1][i_b] = clry[i+conse-1][i_b] - ts_pred_temp;

                        /* normalized to z-scores */
                        for (b = 0; b < LASSO_BANDS; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                /* mini rmse */
                                mini_rmse = max(adj_rmse[i_b], tmpcg_rmse[i_b]);

                                /* z-score */
                                v_diff[conse-1][i_b] = v_dif_mag[conse-1][i_b] / mini_rmse;
                                vec_mag[conse-1] += v_diff[conse-1][i_b] * v_diff[conse-1][i_b]; 
                            }
                        }

                    }
                }
                break_mag = 9999.0;
                for (m = 0; m < conse; m++)
                {
                 printf("m,vec_mag[m]3=%d,%f\n",m,vec_mag[m]);
                    if (break_mag > vec_mag[m])
                        break_mag = vec_mag[m];
                }

                printf("break_mag,vec_mag[0],t_max_cg=%f,%f,%f\n",break_mag,vec_mag[0],t_max_cg);
                if (break_mag > t_cg)
                {
                    printf("Change Magnitude = %.2f\n",break_mag);

                    /* record break time */
                    rec_cg[num_fc].t_break = clrx[i];
                    rec_cg[num_fc].change_prob = 1.0;
                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        matlab_2d_array_mean(v_dif_mag, 0, conse, 
                            &rec_cg[num_fc].magnitude[i_b]);

                    /* identified and move on for the next functional curve */
                    num_fc++;
                    /* start from i+1 for the next functional curve */
                    i_start = i + 1;
                    printf("i_start2=%d\n",i_start);
                    /* start training again */
                    bl_train = 0;
                }
                else if (vec_mag[0] > t_max_cg)
                {
                    /* remove noise */
                    for (m = i; m < num_scenes - 1; m++)
                    {
                        clrx[m] = clrx[m+1];
                        for (b = 0; b < TOTAL_BANDS - 1; b++)
                            clry[m][b] = clry[m+1][b];
                    }
                    end--; /* check if this is needed */

                    /* stay & check again after noise removal */
                    i--;
                    printf("i--=%d\n",i);
                }
            } /* end of continuous monitoring */ 
            /* move forward to the i+1th clear observation */
            i++;
            printf("i++4=%d\n",i);
            status = free_2d_array ((void **) rec_v_dif);
            if (status != SUCCESS)
                RETURN_ERROR ("Freeing memory: rec_v_dif\n", 
                           FUNC_NAME, EXIT_FAILURE);
        } /* end of "while (i <= end - conse) */

        /* Two ways for processing the end of the time series */ 
        if (bl_train == 1)
        {
            end = num_scenes - 1;
            /* if no break find at the end of the time series,
               define probability of change based on conse */
            for (i_conse = conse - 1; i_conse >= 0; i_conse--)
            {
                max_v_dif = 0.0;
                for (i_b = 0; i_b < LASSO_BANDS; i_b++)
                {
                    if (max_v_dif < v_diff[i_conse][i_b])
                        max_v_dif = v_diff[i_conse][i_b];
                    if (max_v_dif <= t_cg)
                    {
                        /* the last stable id */
                        id_last = i_conse + 1;
                        break;
                    }
                }
            } 

            printf("id_last = %d\n",id_last);
            /* update change probability */
            rec_cg[num_fc].change_prob = (conse - id_last) / conse; 
            /* update end time of the curve */
            rec_cg[num_fc].t_end = clrx[end - conse + id_last];

            /* mean value fit for the rest of the pixels < conse & > 1 */
            if (conse > id_last)
            {
                /* update time of the probable change */
                rec_cg[num_fc].t_break = clrx[end-conse+id_last+1];
                /* update magnitude of change */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    matlab_2d_partial_mean(v_dif_mag, i_b, id_last, conse-1, 
                                         &rec_cg[num_fc].magnitude[i_b]);

                /* median of the last < conse pixels */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                {
                    matlab_int_2d_partial_mean(clry, i_b, end-conse+id_last, end - 1, 
                          &fit_cft[0][i_b]);
                    matlab_2d_partial_square_mean(clry, i_b, end-conse+id_last, end - 1,  
                                                  &rmse[i_b]);
                }
                /* identified and move on for the next functional curve */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[end-conse+id_last];
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end-1];
                /* record break time */
                rec_cg[num_fc].t_break = 0;
                /* record postion of the pixel */
                rec_cg[num_fc].pos.row = row;
                rec_cg[num_fc].pos.col = col;
                /* record fitted coefficients */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                {
                    for (k = 0; k < update_num_c; k++)
                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    rec_cg[num_fc].rmse[i_b] = rmse[i_b];
                }
                /* record change probability */
                rec_cg[num_fc].change_prob = 0.0;
                /* record number of observations */
                rec_cg[num_fc].num_obs = conse - id_last;
                /* record fit category */
                rec_cg[num_fc].category = 20 + 1; /* mean value fit at the end */
                /* record change magnitude */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                rec_cg[num_fc].magnitude[i_b] = 0.0; /* record change magnitude */
                num_fc++;   
            }
        }
        else if (bl_train == 0)
        {
            /* if break found close to the end of the time series 
               Use [conse,min_num_c*n_times+conse) to fit curve */
            /* update i_start */
            if (num_fc == rec_fc)
            {
                /* first curve */
                i_start = 1;
                printf("update i_start2=%d\n",i_start);
            }
            else
            {
                for (k = 0; k < num_scenes; k++) 
                {
                     if (clrx[k] > rec_cg[num_fc-1].t_end)
                         i_start  = k + 1;
                }
                printf("update i_start3=%d\n",i_start);
            }

            get_ids_length(clrx, 0, num_scenes-1, &end);
            printf("i_start,end=%d,%d\n",i_start,end);
            /* multitemporal cloud mask */
            status = auto_mask(clrx, clry, i_start-1, end-1,
                               (float)(clrx[end-1]-clrx[i_start-1])/num_yrs, 
                               adj_rmse[1], adj_rmse[4], t_const, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask routine", 
                                  FUNC_NAME, FAILURE);

            /* update i_span after noise removal */
            i_span = 0;
            for (m = i_start-1; m < end; m++)
            {
                if (bl_ids[m] == 0)
                    i_span++;
                printf("i_span update =%d\n", i_span);
            }

            for (m = 0; m < num_scenes-1; m++)
                ids[m] = 0;

            for (m = i_start-1; m <= end - 1; m++)
            {
                ids[m-i_start+1] = m;
            }
            m= 0;
            for (k = 0; k < end-conse; k++)
            {
                if (bl_ids[k] == 1) 
                {
                    rm_ids[m] = ids[k];
                    m++;
                }
            }
            rm_ids_len = m;

            /* remove noise pixels between i_start & i */
            for (m = 0; m < rm_ids_len; m++)
            {
                for (k = rm_ids[m]; k < end - 1; k++)
                {
                    clrx[k] = clrx[k+1];
                    for (b = 0; b < TOTAL_BANDS - 1; b++)
                        clry[k][b] = clry[k+1][b];
                }
                end--;
            }

            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                status = auto_ts_fit(clrx, clry, i_b, i_start-1, end-1, min_num_c, 
                                     fit_cft, &rmse[i_b]); 
                if (status != SUCCESS)  
                     RETURN_ERROR ("Calling auto_ts_fit9\n", FUNC_NAME, EXIT_FAILURE);
            }

            /* record time of curve start */
            rec_cg[num_fc].t_start = clrx[i_start-1];
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end-1];
            /* record break time */
            rec_cg[num_fc].t_break = 0;
            /* record postion of the pixel */
            rec_cg[num_fc].pos.row = row;
            rec_cg[num_fc].pos.col = col;
            /* record fitted coefficients */
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                rec_cg[num_fc].rmse[i_b] = rmse[i_b];
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0;
            /* record number of observations */
            rec_cg[num_fc].num_obs = i_span + conse;
            /* record fit category */
            rec_cg[num_fc].category = 20 + min_num_c; /* simple model fit at the end */
            /* record change magnitude */
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                rec_cg[num_fc].magnitude[i_b] = 0.0; /* record change magnitude */
            num_fc++;
        }
    }

    /* Free memory allocation */
    free(sdate);
    free(clrx);
    free(rmse);
    free(id_range);
    free(id_clr);
    free(id_all);
    free(id_sn);
    free(ids);
    free(ids_old);
    free(bl_ids);
    free(rm_ids);
    free(vec_mag);
    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
        RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                   FUNC_NAME, EXIT_FAILURE);
    status = free_2d_array ((void **) v_diff);
        if (status != SUCCESS)
            RETURN_ERROR ("Freeing memory: v_diff\n", 
                FUNC_NAME, EXIT_FAILURE);
#if 0
    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      EXIT_FAILURE);
    }
#endif
    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Output rec_cg structure to the output file 
       note: can use fread to read out the structure from the output file */
    fp_bin_out = fopen("output.bin", "wb");
    if (fp_bin_out == NULL)
        RETURN_ERROR ("Opening output.bin file\n", FUNC_NAME,
                      EXIT_FAILURE);
    fwrite(rec_cg, sizeof(Output_t), 1, fp_bin_out);
    fclose(fp_bin_out);   

    if (verbose)
    {
        for (i = 0; i < num_fc; i++)
        {
            printf("i=%d\n",i);
            printf("rec_cg[i].t_start=%d\n",rec_cg[i].t_start);
            printf("rec_cg[i].t_end=%d\n",rec_cg[i].t_end);
            printf("rec_cg[i].t_break=%d\n",rec_cg[i].t_break);
            printf("rec_cg[i].pos.row=%d\n",rec_cg[i].pos.row);
            printf("rec_cg[i].pos.col=%d\n",rec_cg[i].pos.col);
            printf("rec_cg[i].num_obs=%d\n",rec_cg[i].num_obs);
            printf("rec_cg[i].category=%d\n",rec_cg[i].category);
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
                    printf("i_b,k,rec_cg[i].coefs[i_b][k] = %d,%d,%f\n", 
                         i_b,k,rec_cg[i].coefs[i_b][k]); 
                printf("rec_cg[i].rmse[i_b] = %f\n",rec_cg[i].rmse[i_b]);
                printf("rec_cg[i].magnitude[i_b]=%f\n",rec_cg[i].magnitude[i_b]); 
            }
        }
    }

    /* Free rec_cg memory*/
    free(rec_cg);

    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC end_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    return EXIT_SUCCESS;
    }


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/15/2013   Song Guo         Original Development
8/15/2013   Song Guo         Modified to use TOA reflectance file 
                             as input instead of metadata file
2/19/2014   Gail Schmidt     Modified to utilize the ESPA internal raw binary
                             file format

******************************************************************************/
void
usage ()
{
    printf ("Continuous Change Detection and Classification\n");
    printf ("\n");
    printf ("usage: ccdc"
            " --row=input row number"
            " --col=input col number"
            " --t_cg=chi-square inversed T_cg value" 
            " --conse = number of points used for change detection" 
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    -row: input row number\n");
    printf ("    --col=input col number\n");
    printf ("    --t_cg=chi-square inversed T_cg value\n");
    printf ("    --conse = number of points used for change detection\n");
    printf ("\n");
    printf ("where the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n");
    printf ("Example: ccdc"
            " --row=5099"
            " --col=3191"
            " --t_cg=15.0863"
            " --conse=6"
            " --verbose\n");
    printf ("Note: The ccdc must run from the directory"
            " where the input data are located.\n\n");
}
