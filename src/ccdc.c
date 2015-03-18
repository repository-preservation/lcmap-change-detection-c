#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>

#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "ccdc.h"

#define NUM_LASSO_BANDS 5
#define LASSO_BANDS 5
#define MAX_NUM_FC 10 /* Values change with number of pixels run */
int lasso_blist[LASSO_BANDS] = {2, 3, 4, 5, 7};

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
int
main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";
    char msg_str[MAX_STR_LEN];  /* input data scene name */
    char filename[MAX_STR_LEN];         /* input binary filenames */
    float t_cg;
    float t_max_cg;
    int conse;
    int status;                 /* return value from function call */
    Output_t *rec_cg = NULL;    /* output structure and metadata */
    bool verbose;               /* verbose flag for printing messages */
    int i, k, m, b;
    char **scene_list = NULL;
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
    int t_const = 400;        /* Threshold for cloud, shadow, and snow detection */
    int mini_yrs = 1;         /* minimum year for model intialization */
    int num_detect = LASSO_BANDS;/* number of bands for change detection */
    float p_min = 0.1;        /* percent of ref for mini_rmse */
    float t_ws = 0.95;        /* no change detection for permanent water pixels */
    float t_sn = 0.6;         /* no change detection for permanent snow pixels */ 
    float t_cs = 0.6;         /* Fmask fails threshold */
    int *sdate;
    Input_meta_t *meta;
    int row, col;
    int landsat_number;
    int fmask_sum = 0;
    int clear_water_sum = 0;
    int clear_land_sum = 0;
    int snow_sum = 0;
    int cloud_and_shadow_sum = 0;
    float ws_pct;
    float sn_pct;
    float cs_pct;
    int n_ws = 0;
    int n_sn = 0;
    int n_clr = 0;
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
    int v_qa;
    int bl_train;
    float time_span;
    int *bl_ids;
    int *id_range;
    int *ids;
    int *ids_old;
    int *clr_ids;
    int *rm_ids;
    int rm_ids_len;
    int i_rec;
    float *v_start;
    float *v_end;
    float *v_slope;
    float *v_dif;
    float **v_diff;
    float mean_v;
    float v_dif_norm;
    int i_count;
    int new_i_start;
    float **v_dif_mag;
    int i_conse, i_b;
    float *vec_mag;
    float v_dif_mean;
    float **rec_v_dif;
    float **rec_v_dif_temp;
    float rec_adj_mini[TOTAL_BANDS-1];
    float mini[TOTAL_BANDS-1];
    float mini_rmse;
    int bl_tmask;
    int n_rmse;
    float tmpcg_rmse[TOTAL_BANDS-1];
    int d_rt;
    float *d_yr;
    float break_mag;
    float max_v_dif;
    int id_last;
    float v_dif_temp;
    float ts_pred_temp;
    FILE *fp_bin_out;
    int *id_csc;
    float min_rmse;
    int *ids_old_len;
    int *d_yr_idx;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    status = get_args (argc, argv, &row, &col, &min_rmse, &t_cg, &t_max_cg, 
                       &conse, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /* minimum rmse */
    for (i = 0; i < TOTAL_BANDS - 2; i++) 
        mini[i] = 100.0 * min_rmse;     
    mini[TOTAL_BANDS-1] = 400.0 * min_rmse;

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

    ids_old = (int *)calloc(num_scenes, sizeof(int));
    if (ids_old == NULL)
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);

    bl_ids = (int *)calloc(num_scenes, sizeof(int));
    if (bl_ids == NULL)
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);

    clr_ids = (int *)calloc(num_scenes, sizeof(int));
    if (clr_ids == NULL)
        RETURN_ERROR("ERROR allocating clr_ids memory", FUNC_NAME, FAILURE);

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

    /* sort scene_list based on year & julian_day */
    status = sort_scene_based_on_year_doy(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    status = read_envi_header(scene_list[0], meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
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
        printf("Clear-sky pixel percentage = %f7.2\n", 
            (float)(fmask_sum / num_scenes));
    }

    /* pixel value ranges should follow physical rules */
    for (i = 0; i < num_scenes; i++)
    { 
        if ((buf[i][0] > 0) && (buf[i][0] < 10000) &&
            (buf[i][1] > 0) && (buf[i][1] < 10000) &&
            (buf[i][2] > 0) && (buf[i][2] < 10000) &&
            (buf[i][3] > 0) && (buf[i][3] < 10000) &&
            (buf[i][4] > 0) && (buf[i][4] < 10000) &&
            (buf[i][5] > -9320) && (buf[i][5] < 7070) &&
            (buf[i][6] > 0) && (buf[i][6] < 10000))
            id_range[i] = 1;
        else
            id_range[i] = 0;
    }

    /* Get each mask pixel totals */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] == 0)
            clear_land_sum++;
        if (fmask_buf[i] == 1)
            clear_water_sum++;
        if (fmask_buf[i] == 2 || fmask_buf[i] ==4)
            cloud_and_shadow_sum++;
        if (fmask_buf[i] == 3)
            snow_sum++;
    }

    /* percent of water pixels */
    ws_pct = (float) clear_water_sum / (float) (clear_land_sum + clear_water_sum);

    /* percent of snow observations */
    sn_pct =  (float) snow_sum / (float)(clear_land_sum + 
               clear_water_sum + snow_sum); 

    /* percent of cloud and cloud shadow pixels */
    cs_pct = (float)cloud_and_shadow_sum / (float)fmask_sum;  


    /* Allocate memory for rec_cg */ 
    rec_cg = malloc(MAX_NUM_FC * sizeof(Output_t));
    if (rec_cg == NULL)
        RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);


    /* fit only water pixels (permanet water) */
    if (ws_pct > t_ws)
    {
        for (i = 0; i < num_scenes; i++)
        { 
            if (fmask_buf[i] == 1 && id_range[i] == 1)
            {
                clrx[n_ws] = sdate[i];
                for (k = 0; k < TOTAL_BANDS - 1; k++)
                     clry[n_ws][k] = buf[i][k];
                n_ws++;
            }   
        }
        end = n_ws;

        printf("Permanent water pixel = %f\n", 100.0 * ws_pct);

        if (n_ws < n_times * min_num_c ) /* not enough water pixels */
        {
            RETURN_ERROR ("Not enough good water observations\n", 
                 FUNC_NAME, EXIT_FAILURE);
        }
        else
        {
            /* start model fit for snow persistent pixels */
            printf ("Fit permanent water observations\n"); 

            /* the first observation for TSFit */
            i_start = 0; /* the first observation for TSFit */

            /* use all water and snow values to fit */
            i_span = n_ws;
            update_cft(i_span, n_times, min_num_c, mid_num_c, max_num_c, num_c,
                       &update_num_c);

            for (k = 0; k < TOTAL_BANDS - 1; k++)
            {
                 status = auto_ts_fit(clrx, clry, k, 0, end-1, update_num_c, fit_cft, 
                                  &rmse[k], &v_dif_temp); 
                 if (status != SUCCESS)  
                     RETURN_ERROR ("Calling auto_ts_fit1\n", 
                          FUNC_NAME, EXIT_FAILURE);
            }

            /* update information at each iteration */
            /* record time of curve start */
            rec_cg[num_fc].t_start = clrx[i_start];
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end];
            /* no break at the moment */
            rec_cg[num_fc].t_break = 0;
            /* record postion of the pixel */
            rec_cg[num_fc].pos.row = row;
            /* record postion of the pixel */
            rec_cg[num_fc].pos.col = col;
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
                    /* record postion of the pixel */
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                /* record rmse of the pixel */
                rec_cg[num_fc].rmse[i_b] = rmse[i]; 
            }
            /* record change probability */            
            rec_cg[num_fc].change_prob = 0.0;
            /* record number of observations */
            rec_cg[num_fc].num_obs = n_ws;
            /* record fit category */
            rec_cg[num_fc].category = 20 + update_num_c;
            for (i = 0; i < TOTAL_BANDS - 1; i++)
                /* record change magnitude */
                rec_cg[num_fc].magnitude[i] = 0.0;
            /* NUM of Fitted Curves (num_fc) */
            num_fc++;   
        }
    }
    else if (sn_pct > t_sn) /* fit only snow pixels (permanet snow) */
    {
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
        end = n_sn;

        printf("Permanent snow pixel = %f\n", 100.0 * sn_pct);

        if (n_sn < n_times * min_num_c ) /* not enough snow pixels */
        {
            RETURN_ERROR ("Not enough good snow observations\n", 
                 FUNC_NAME, EXIT_FAILURE);
        }
        else
        {
            /* start model fit for snow persistent pixels */
            printf ("Fit permanent snow observations\n"); 

            /* the first observation for TSFit */
            i_start = 0; /* the first observation for TSFit */

            /* treat saturated and unsaturated pixels differently */
            for (k = 0; k < TOTAL_BANDS -1; k++)
            {
                i_span = 0;
                if (k != TOTAL_BANDS - 3) /* for optical bands */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > 0 && clry[i][k] < 10000)
                            i_span++;
                        if (i_span < min_num_c * n_times)
                            fit_cft[i][k] = 10000; /* fixed value for saturated
                                                       pixels */
                        else
                        {
                            update_cft(i_span, n_times, min_num_c, mid_num_c, 
                                       max_num_c, num_c, &update_num_c);

                            for (k = 0; k < TOTAL_BANDS - 1; k++)
                            {
                                status = auto_ts_fit(clrx, clry, k, 0, i_span-1, update_num_c, 
                                         fit_cft, &rmse[k], &v_dif_temp); 
                                if (status != SUCCESS)  
                                    RETURN_ERROR ("Calling auto_ts_fit2\n", 
                                           FUNC_NAME, EXIT_FAILURE);
                            }
                        }
                    } 
                }
                else /* for thermal band */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > -9300 && clry[i][k] < 7070)
                            i_span++;
                        update_cft(i_span, n_times, min_num_c, mid_num_c, 
                                   max_num_c, num_c, &update_num_c);

                        status = auto_ts_fit(clrx, clry, k, 0, i_span-1, update_num_c, fit_cft, 
                                 &rmse[k], &v_dif_temp); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit3\n", 
                                 FUNC_NAME, EXIT_FAILURE);
                    }
                }
            }

            /* update information at each iteration */
            /* record time of curve start */
            rec_cg[num_fc].t_start = clrx[i_start]; 
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end]; 
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
                rec_cg[num_fc].rmse[i] = rmse[i]; 
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0; 
            /* record number of observations */
            rec_cg[num_fc].num_obs = n_ws; 
            /* record fit category */
            rec_cg[num_fc].category = 20 + update_num_c; 
            for (i = 0; i < TOTAL_BANDS - 1; i++)
                /* record change magnitude */ 
                rec_cg[num_fc].magnitude[i] = 0.0; 
            /* NUM of Fitted Curves (num_fc) */
            num_fc++;   
        }
    }
    else /* clear land or water pixels */
    {
        printf("Land pixel (water/snow) = %f/%f\n", 100.0 * ws_pct, 100.0 * sn_pct);

        if (cs_pct < t_cs)  /* Fmask works */
        {
            for (i = 0; i < num_scenes; i++)
            { 
                if ((fmask_buf[i] == 0 || fmask_buf[i] == 1) && id_range[i] == 1)
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_BANDS - 1; k++)
                         clry[n_clr][k] = buf[i][k];
                    n_clr++;
                }   
            }
            end = n_clr;
            v_qa = 40; /* QA var for normal procedure */
        }
        else /* Fmask fails in persistent commission */
        {
            /* allocate temporary memory */
            int16 **cscy;
            cscy = (int16 **) allocate_2d_array (num_scenes, TOTAL_BANDS - 1, 
                                     sizeof (int16));
            if (cscy == NULL)
            {
                RETURN_ERROR ("Allocating cscy memory", FUNC_NAME, FAILURE);
            }
            int16 *cscy_fil;
            cscy_fil = (int16 *)malloc(num_scenes * sizeof (int16));
            if (cscy_fil == NULL)
            {
                RETURN_ERROR ("Allocating cscy_fil memory", FUNC_NAME, FAILURE);
            }

            id_csc = (int *)malloc(num_scenes * sizeof (int));
            if (id_csc == NULL)
            {
                RETURN_ERROR ("Allocating id_csc memory", FUNC_NAME, FAILURE);
            }

            /* catch back all valid wthin range observations */
            for (i = 0; i < num_scenes; i++)
            { 
                if ((fmask_buf[i] == 0 || fmask_buf[i] == 1 ||
                     fmask_buf[i] == 2 || fmask_buf[i] == 4) && (id_range[i] == 1))
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_BANDS - 1; k++)
                         cscy[n_clr][k] = buf[i][k];
                    id_csc[i] = 1;
                    n_clr++;
                }   
                else
                    id_csc[i] = 0;
            }

            /* more clear observation added */
            median_filter(cscy, n_clr, 2 * conse +1, cscy_fil);
            int n_clr_final = 0;
            for (i = 0; i < n_clr; i++)
            { 
                for (k = 0; k < TOTAL_BANDS - 1; k++)
                {
                    if ((cscy[i][k] < cscy_fil[i] + 400) || 
                        (id_csc[i] == 1 && fmask_buf[id_csc[i]] <= 1))
                    {
                        clry[i][k] = cscy[i][k];
                        id_range[i] = 1;
                        n_clr_final++;
                    }
                    else
                        id_range[i] = 0;
                }
            }
            end = n_clr_final;
            v_qa = 30; /* QA var for normal procedure */

            /* Free the temporary memory */
            status = free_2d_array ((void **) cscy);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: cscy\n", FUNC_NAME,
                              EXIT_FAILURE);
            }
            free(cscy_fil);
            free(id_csc);
        }

        /* start with mininum requirement of clear obs */
        i = n_times * min_num_c;

        /* the first observation for TSFit */
        i_start = 0; 

        /* record the start of the model initialization (0=>initial;1=>done) */
        bl_train = 0;

        /* record the num_fc at the beginning of each pixel */
        rec_fc = num_fc;

        /* record the start of Tmask (0=>initial;1=>done) */
        bl_tmask = 0;

        /* while loop - process till the last clear observation - conse */
        while (i < end - conse)
        {
            /* span of "i" */
            i_span = i - i_start + 1;

            /* span of time (num of years) */
            time_span = (float)(clrx[i] - clrx[i_start]) / num_yrs;

            /* basic requrirements: 1) enough observations; 2) enough time */
            if (i_span >= n_times * min_num_c && time_span >= (float)mini_yrs)

            /* initializing model */
            if (bl_train == 0)
            {
                /* do for the first time */
                if (bl_tmask == 0)
                {
                    for (k = i_start; k < end; k++)
                        clr_ids[k] = 1;
                    /* update BL_tmask value */
                    bl_tmask = 1;
                }

                /* Allocate memory for ids, rm_ids */
                ids = (int *)calloc(i, sizeof(int));
                if (ids == NULL)
                    RETURN_ERROR("ERROR allocating ids memory", 
                                  FUNC_NAME, FAILURE);

                rm_ids = (int *)calloc(i-i_start+1, sizeof(int));
                if (rm_ids == NULL)
                    RETURN_ERROR("ERROR allocating rm_ids memory", 
                                  FUNC_NAME, FAILURE);

                /* step 1: noise removal */ 
                status = auto_mask(clrx, clry, i_start, i+conse,
                                   (clrx[i+conse]-clrx[i_start])/num_yrs, t_const, bl_ids);
                if (status != SUCCESS)
                    RETURN_ERROR("ERROR calling auto_mask routine", 
                                  FUNC_NAME, FAILURE);
  
                for (k = i_start; k <= i; k++)
                    ids[k] = k;
                m= 0;
                for (k = 0; k < end-conse; k++)
                {
                    if (bl_ids[k] == 0) 
                    {            
                        clr_ids[k] = 1;
                    }
                    else
                    {
                        clr_ids[k] = 0;
                        rm_ids[m] = ids[k];
                        m++;
                    }
                }

                rm_ids_len = m;
                /* update i_span after noise removal */
                i_span = i - i_start + 1 - rm_ids_len;

                /* check if there is enough observation */
                if (i_span < n_times * min_num_c)
                {
                    /* move forward to the i+1th clear observation */
                    i++;
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

                    /* remove noise pixels between i_start & i */
                    for (m = 0; m < rm_ids_len; m++)
                    {
                        for (k = rm_ids[m]; k < end - 1; k++)
                        {
                            cpx[k] = cpx[k+1];
                            for (b = 0; b < TOTAL_BANDS - 1; b++)
                                cpy[k][b] = cpy[k][b];
                        }
                        end--;
                    }

                    /* record i before noise removal 
                       This is very important as if model is not initialized 
                       the multitemporal masking shall be done again instead 
                       of removing outliers in every masking */
                    i_rec=i;

                    /* update i afer noise removal (i_start stays the same) */
                    i=i_start + i_span - 1;

                    /* update span of time (num of years) */
                    time_span=(cpx[i] - cpx[i_start]) / num_yrs;

                    /* check if there is enough time */
                    if (time_span < mini_yrs)
                    {
                        i = i_rec;   /* keep the original i */
                        /* move forward to the i+1th clear observation */
                        i++;        
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

                        /* Step 2: model fitting: initialize model testing variables
                           defining computed variables */
                        v_start = malloc(num_detect * sizeof(float));
                        if (v_start == NULL)
                            RETURN_ERROR ("Allocating v_start memory", FUNC_NAME, FAILURE);
                        v_end = malloc(num_detect * sizeof(float));
                        if (v_end == NULL)
                            RETURN_ERROR ("Allocating v_end memory", FUNC_NAME, FAILURE);
                        v_slope = malloc(num_detect * sizeof(float));
                        if (v_slope == NULL)
                            RETURN_ERROR ("Allocating v_slope memory", FUNC_NAME, FAILURE);
                        v_dif = malloc(num_detect * sizeof(float));
                        if (v_dif== NULL)
                            RETURN_ERROR ("Allocating v_dif memory", FUNC_NAME, FAILURE);
 
                        for (b = 0; b < num_detect; b++)
                        { 
                            /* Initial model fit */
                            status = auto_ts_fit(clrx, clry, k, i_start, i, min_num_c, fit_cft, 
                                     &rmse[k], &v_dif[k]); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit4\n", 
                                     FUNC_NAME, EXIT_FAILURE);

                            /* calculate mini rmse with mean values & mini */
                            mean_v = fit_cft[0][lasso_blist[b]] + 
                                 fit_cft[1][lasso_blist[b]] * 
                             (float)(clrx[i_start]+clrx[i]) / 2.0;
                            mini_rmse = max(mean_v * p_min, mini[lasso_blist[b]]);
                            /* minimum rmse */
                            mini_rmse = max(mini_rmse, rmse[lasso_blist[b]]);

                            /* compare the first clear obs */
                            auto_ts_predict(clrx[i_start], fit_cft, b, &ts_pred_temp);  
                            v_start[lasso_blist[b]] = (clry[i_start][lasso_blist[b]] -
                               ts_pred_temp)/mini_rmse;

                            /* compare the last clear observation */
                            auto_ts_predict(clrx[i], fit_cft, b, &ts_pred_temp);
                            v_end[lasso_blist[b]] = (clry[i][lasso_blist[b]]-
                                                     ts_pred_temp)/mini_rmse;

                            /* anormalized slope values */
                            v_slope[lasso_blist[b]] = fit_cft[1][lasso_blist[b]] *
                                    (clrx[i]-clrx[i_start])/mini_rmse;
                            
                            /* differece in model intialization */
                            v_dif[lasso_blist[b]] = fabs(v_slope[lasso_blist[b]]) + 
                                                    fabs(v_start[lasso_blist[b]]) + 
                                                    fabs(v_end[lasso_blist[b]]);
                        }
                        matlab_norm(v_dif, LASSO_BANDS, &v_dif_norm);
                        v_dif_norm *= v_dif_norm; 

                        /* free allocated memories */
                        free(ids);
                        free(rm_ids);
                        free(v_start); 
                        free(v_end);
                        free(v_slope); 
                        free(v_dif);

                        /* find stable start for each curve */
                        if (v_dif_norm > t_cg)
                        {
                            /* start from next clear obs */
                            i_start++;

                            /* move forward to the i+1th clear observation */
                            i++;

                            /* keep all data and move to the next obs */
                            continue;
                        }
                        else
                        {
                            /* model ready */
                            bl_train = 1;

                            /* count difference of i for each iteration */
                            i_count = 0;

                            /* give new i_start for the new curve */
                            new_i_start = i_start;

                            /* model fit at the beginning of time series */
                            if ((num_fc == rec_fc) && (i_start > 1))
                            {
                                /* allocate memory for model_v_dif */ 
                                v_dif_mag = (float **) allocate_2d_array(i_start-1,
                                          num_detect, sizeof (float));
                                if (v_dif_mag == NULL)
                                {
                                    RETURN_ERROR ("Allocating v_dif_mag memory", 
                                                  FUNC_NAME, FAILURE);
                                }

                                /* allocate memory for v_diff */ 
                                v_diff = (float **) allocate_2d_array(i_start-1,
                                          num_detect, sizeof (float));
                                if (v_diff == NULL)
                                {
                                    RETURN_ERROR ("Allocating v_diff memory", 
                                                  FUNC_NAME, FAILURE);
                                }

                                /* detect change. 
                                   value of difference for conse obs
                                   record the magnitude of change */
                                for (i_conse = 0; i_conse < i_start-1; i++)
                                {
                                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                                    {
                                        /* absolute differences */
                                        auto_ts_predict(clrx[i_conse], fit_cft, i_b, 
                                            &ts_pred_temp);
                                        v_dif_mag[i_conse][i_b] = clry[i_conse][i_b] - 
                                            ts_pred_temp;

                                        /* normalize to z-score */
                                        for (b = 0; b < LASSO_BANDS - 1; b++)
                                        {
                                            if (i_b == lasso_blist[b])
                                            {
                                                /* calculate mini rmse with mean values & mini */
                                                mean_v = fit_cft[0][i_b] + fit_cft[1][i_b] * 
                                                    (float)(clrx[i_start] + clrx[i]) / 2.0;
                                                mini_rmse = max(mean_v * p_min, mini[i_b]);

                                                /* minimum rmse */ 
                                                mini_rmse = max(mini_rmse, rmse[i_b]);
 
                                                /* z-scores */
                                                v_diff[i_conse][i_b] = abs(v_dif_mag[i_conse][i_b]) 
                                                                       / mini_rmse;
                                            }
                                        }
                                    }
                                    matlab_2d_norm(v_diff, i_start-1, LASSO_BANDS, &v_dif_norm);
                                    vec_mag[i_conse] *= v_dif_norm; 

                                    if (vec_mag[i_conse] <= t_cg)
                                    {
                                        new_i_start = i_conse;
                                        break;
                                    }
                                }

                                if (new_i_start > conse)
                                {
                                    /* defining computed variables */
                                    for (i_b = 0; i_b < TOTAL_BANDS -1; i_b++)
                                    {
                                        status = auto_ts_fit(clrx, clry, i_b, 0, new_i_start-1, 
                                                 min_num_c, fit_cft, &rmse[i_b], &v_dif_temp); 
                                        if (status != SUCCESS)  
                                            RETURN_ERROR ("Calling auto_ts_fit5\n", 
                                                FUNC_NAME, EXIT_FAILURE);
                                        /* record time of curve start */
                                        rec_cg[num_fc].t_start = clrx[0]; 
                                        /* record time of curve end */
                                        rec_cg[num_fc].t_end = clrx[new_i_start-1]; 
                                        /* no break at the moment */
                                        rec_cg[num_fc].t_break = new_i_start; 
                                        /* record postion of the pixel */
                                        rec_cg[num_fc].pos.row = row; 
                                        /* record postion of the pixel */
                                        rec_cg[num_fc].pos.col = col; 
                                        for (k = 0; k < update_num_c; k++)
                                            /* record fitted coefficients */
                                            rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                        /* record rmse of the pixel */                         
                                        rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                        
                                        /* record change probability */
                                        rec_cg[num_fc].change_prob = 1.0; 
                                        /* record number of observations */
                                        rec_cg[num_fc].num_obs = new_i_start-1; 
                                        /* record fit category */
                                        rec_cg[num_fc].category = v_qa + min_num_c; 
                                        /* record change magnitude */
                                        matlab_2d_array_mean(v_dif_mag, i_b, conse, 
                                                             &v_dif_mean);  
                                        rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 

                                        /* identified and move on for the next functional curve */
                                        num_fc++;  
                                    }
                                    status = free_2d_array ((void **) v_dif_mag);
                                    if (status != SUCCESS)
                                        RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                                             FUNC_NAME, EXIT_FAILURE);
                                    status = free_2d_array ((void **) v_diff);
                                    if (status != SUCCESS)
                                         RETURN_ERROR ("Freeing memory: v_diff\n", 
                                             FUNC_NAME, EXIT_FAILURE);

                                }
                                else if (new_i_start > 1)
                                {
                                    /* median value fit for the rest of the pixels < conse */
                                    for (i_b = 0; i< TOTAL_BANDS - 1; i++)
                                    {
                                        matlab_2d_int_array_mean(clry, i_b, new_i_start-1, &fit_cft[0][i_b]);
                                        square_root_mean(clry, i_b, new_i_start-1, fit_cft, &rmse[i_b]);
                                    }
                                    /* record time of curve start */
                                    rec_cg[num_fc].t_start = clrx[0]; 
                                    /* record time of curve end */  
                                    rec_cg[num_fc].t_end = clrx[new_i_start-1]; 
                                    /* no break at the moment */
                                    rec_cg[num_fc].t_break = clrx[new_i_start]; 
                                    /* record fitted coefficients */
                                    rec_cg[num_fc].pos.row = row; 
                                    /* record fitted coefficients */
                                    rec_cg[num_fc].pos.col = col; 
                                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i++)
                                    {
                                        for (k = 0; k < update_num_c; k++)
                                            /* record fitted coefficients */
                                            rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                        /* record rmse of the pixel */
                                        rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                    }
                                    /* record change probability */
                                    rec_cg[num_fc].change_prob = -(new_i_start-1)/conse; 
                                    /* record number of observations */
                                    rec_cg[num_fc].num_obs = new_i_start-1; 
                                    /* record fit category */
                                    rec_cg[num_fc].category = v_qa + 1; /* record fit category */
                                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                                    {
                                        matlab_2d_array_mean(v_dif_mag, i_b, conse, &v_dif_mean);
                                        /* record change magnitude */ 
                                        rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                                    }
                                    /* NUM of Fitted Curves (num_fc) */
                                    num_fc++;
                                }    
                            }
                        }
                    }              
                } /* end of "if (v_dif > t_cg)" */  
            } /* end of initializing model */

            /* continuous monitoring started!!! */
            if (bl_train == 1)
            {
                /* all IDs */
                for (k = i_start; k < i; k++)
                {
                    ids[k] = k;
                }
                i_span = i - i_start + 1;

                /* determine the time series model */
                update_cft(i_span, n_times, min_num_c, mid_num_c, max_num_c, 
                           num_c, &update_num_c);

                /* defining computed variables */
                rec_v_dif = (float **)allocate_2d_array(k, TOTAL_BANDS - 1,
                                                            sizeof (float));
                if (rec_v_dif == NULL)
                    RETURN_ERROR ("Allocating rec_v_dif memory", 
                                  FUNC_NAME, FAILURE);

                /* dynamic model fit when there are not many obs */
                if (i_count == 0)
                {
                    /* update i_count at each interation */
                    i_count = clrx[i] - clrx[i_start];

                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        status = auto_ts_fit(clrx, clry, i_b, i_start, i, min_num_c, 
                                             fit_cft, &rmse[i_b], &v_dif_temp); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit6\n", 
                                     FUNC_NAME, EXIT_FAILURE);
                        for (b = 0; b < LASSO_BANDS - 1; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                /* calculate mini rmse with mean values & mini */
                                mean_v = fit_cft[0][i_b] + fit_cft[1][i_b] * 
                                    (float)(clrx[i_start] + clrx[i]) / 2.0;
                                rec_adj_mini[i_b] = max(mean_v * p_min, mini[i_b]);
                            }
                        }
                    }

                    /* updating information for the first iteration */
                    /* record time of curve start */
                    rec_cg[num_fc].t_start = clrx[new_i_start]; 
                    /* record time of curve end */
                    rec_cg[num_fc].t_end = clrx[i]; 
                    /* no break at the moment */
                    rec_cg[num_fc].t_break = 0; 
                    /* record postion of the pixel */
                    rec_cg[num_fc].pos.row = row; 
                    /* record postion of the pixel */
                    rec_cg[num_fc].pos.col = col; 
                    for (i = 0; i < TOTAL_BANDS - 1; i++)
                    {
                        for (k = 0; k < update_num_c; k++)
                            /* record fitted coefficients */
                            rec_cg[num_fc].coefs[i][k] = fit_cft[i][k]; 
                            /* record rmse of the pixel */
                            rec_cg[num_fc].rmse[i] = rmse[i]; 
                    }
                    /* record change probability */
                    rec_cg[num_fc].change_prob = 0.0; 
                    /* record number of observations */
                    rec_cg[num_fc].num_obs = i-i_start+1; 
                    /* record fit category */
                    rec_cg[num_fc].category = v_qa + update_num_c; 
                    for (i = 0; i < TOTAL_BANDS - 1; i++)
                    {
                        /* record change magnitude */
                        rec_cg[num_fc].magnitude[i] = 0.0; 
                    }

                    /* detect change, value of difference for conse obs */
                           
                    /* allocate memory for v_dif_mag */ 
                    v_dif_mag = (float **) allocate_2d_array(conse,
                                         TOTAL_BANDS - 1, sizeof (float));
                    if (v_dif_mag == NULL)
                    {
                        RETURN_ERROR ("Allocating v_dif_mag memory", 
                                             FUNC_NAME, FAILURE);
                    }

                    /* allocate memory for v_diff */ 
                    v_diff = (float **) allocate_2d_array(conse,
                                    TOTAL_BANDS - 1, sizeof (float));
                    if (v_diff == NULL)
                    {
                        RETURN_ERROR ("Allocating v_diff memory", 
                                       FUNC_NAME, FAILURE);
                    }

                    for (i_conse = 0; i_conse < i_start-1; i++)
                    {
                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        {
                            /* absolute differences */
                            auto_ts_predict(clrx[i_conse], fit_cft, i_b, &ts_pred_temp);
                            v_dif_mag[i_conse][i_b] = clry[i_conse][i_b] - ts_pred_temp; 
       
                           /* normalize to z-score */
                            for (b = 0; b < LASSO_BANDS - 1; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /* minimum rmse */ 
                                    mini_rmse = max(rec_adj_mini[i_b], rmse[i_b]);
 
                                    /* z-scores */
                                    v_diff[i_conse][i_b] = abs(v_dif_mag[i_conse][i_b]) / mini_rmse;
                                }
                            }
                        }
                        matlab_2d_norm(v_diff, conse, LASSO_BANDS, &v_dif_norm);
                        vec_mag[i_conse] = v_dif_norm * v_dif_norm; 
                    }

                    /* IDs that haven't been updated */
                    for (k = i_start; k < i; k++)
                        ids_old[k] = ids[k];

                    status = free_2d_array ((void **) v_dif_mag);
                    if (status != SUCCESS)
                        RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                               FUNC_NAME, EXIT_FAILURE);
                    status = free_2d_array ((void **) v_diff);
                    if (status != SUCCESS)
                        RETURN_ERROR ("Freeing memory: v_diff\n", 
                               FUNC_NAME, EXIT_FAILURE);
                }
                else
                {
                    if ((clrx[i] - clrx[i_start]) >= (i_count + i_count / 3))
                    {
                        /* update i_count at each interation year */
                        i_count = clrx[i] - clrx[i_start];

                        /* defining computed variables */
                        rec_v_dif = (float **)allocate_2d_array(b, TOTAL_BANDS - 1,
                                                            sizeof (float));
                        if (rec_v_dif == NULL)
                        {
                            RETURN_ERROR ("Allocating rec_v_dif memory", 
                                          FUNC_NAME, FAILURE);
                        }


                        for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                        {
                            status = auto_ts_fit(clrx, clry, i_b, i_start, i, min_num_c, 
                                             fit_cft, &rmse[i_b], &v_dif_temp); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit7\n", 
                                     FUNC_NAME, EXIT_FAILURE);

                            for (b = 0; b < LASSO_BANDS - 1; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /* calculate mini rmse with mean values & mini */
                                    mean_v = fit_cft[0][lasso_blist[b]] + 
                                             fit_cft[1][lasso_blist[b]] * 
                                             (float)(clrx[i_start]+clrx[i]) / 2.0;
                                    rec_adj_mini[i_b] = max(mean_v * p_min, mini[i_b]);
                                }
                            }
                        }

                        /* record fitted coefficients */
                        for (i = 0; i < TOTAL_BANDS - 1; i++)
                        {
                            for (k = 0; k < update_num_c; k++)
                                /* record fitted coefficients */
                                rec_cg[num_fc].coefs[i][k] = fit_cft[i][k]; 
                                /* record rmse of the pixel */
                                rec_cg[num_fc].rmse[i] = rmse[i]; 
                        }
                        /* record number of observations */
                        rec_cg[num_fc].num_obs = i-i_start+1; 
                        /* record fit category */
                        rec_cg[num_fc].category = v_qa + update_num_c; 

                        /* IDs that haven't been updated */
                        for (k = i_start; k < i; k++)
                            ids_old[k] = ids[k];
                    }

                    /* record time of curve end */
                    rec_cg[num_fc].t_end = clrx[i]; /* record time of curve end */

                    ids_old_len = get_ids_length(ids_old, k);
                    /* use temporally-adjusted RMSE */
                    if (*ids_old_len <= n_times * max_num_c)
                    {
                        /* number of observations for calculating RMSE */
                        n_rmse = *ids_old_len;
                        for (b = 0; b < TOTAL_BANDS-1; b++)
                            tmpcg_rmse[b] = rmse[b];
                    }
                    else
                    {
                        /* use fixed number for RMSE computing */
                        n_rmse = n_times * max_num_c;

                        /* better days counting for RMSE calculating */
                        /* relative days distance */
                        d_yr = malloc(*ids_old_len * sizeof(float));
                        if (d_yr == NULL)
                            RETURN_ERROR ("Allocating d_yr memory", 
                                          FUNC_NAME, FAILURE);
                        d_yr_idx = malloc(*ids_old_len * sizeof(int));
                        if (d_yr == NULL)
                            RETURN_ERROR ("Allocating d_yr_idx memory", 
                                          FUNC_NAME, FAILURE);
                        rec_v_dif_temp = (float **)allocate_2d_array(*ids_old_len,  
                                          TOTAL_BANDS - 1, sizeof(float));
                        if (rec_v_dif_temp == NULL)
                            RETURN_ERROR ("Allocating rec_v_dif_temp memory", 
                                          FUNC_NAME, FAILURE);
 
                        for(m = 0; m < *ids_old_len; m++)
                        {
                            d_rt = clrx[m] - clrx[i+conse]; 
                            d_yr[m] = fabs(round((float)(d_rt/num_yrs)*num_yrs - d_rt));
                            d_yr_idx[m] = m;
                        }

                        /* sort the d_yr */
                        quick_sort_index(d_yr, d_yr_idx, 0, *ids_old_len-1);

                        for(m = 0; m < *ids_old_len; m++)
                        {
                            for (b = 0; b < TOTAL_BANDS - 1; b++)
                            {
                                rec_v_dif_temp[m][b] = rec_v_dif[ids_old[d_yr_idx[m]] -
                                    ids_old[0]][b];
                            }
                        }
                        /* temporarily changing RMSE */
                        for (b = 0; b < TOTAL_BANDS-1; b++)
                        {
                            matlab_2d_array_norm(rec_v_dif_temp, b, *ids_old_len, 
                                                 &tmpcg_rmse[b]);
                            tmpcg_rmse[b] /= sqrt(n_rmse - update_num_c);  
                        }

                        /* free allocated memories */
                        free(d_yr);
                        free(d_yr_idx);
                        status = free_2d_array ((void **) rec_v_dif_temp);
                        if (status != SUCCESS)
                            RETURN_ERROR ("Freeing memory: rec_v_dif_temp\n", 
                                  FUNC_NAME, EXIT_FAILURE);
                    }

                    /* move the ith col to i-1th col */
                    for (b = 0; b < TOTAL_BANDS-1; b++)
                    {
                        for (m = 1; m < conse-1; m++)
                        {
                            v_diff[m][b] = v_diff[m+1][b];
                            v_dif_mag[m][b] = v_dif_mag[m+1][b];
                            vec_mag[m] = vec_mag[m+1];
                        }
                        v_diff[conse][b] = 0.0;
                        v_dif_mag[conse][b] = 0.0;
                    }
                    vec_mag[conse] = 0.0;

                    for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                    {
                        /* absolute difference for all bands */
                        auto_ts_predict(clrx[i+conse], fit_cft, i_b, &ts_pred_temp);
                        v_dif_mag[conse][i_b] = clry[i+conse][i_b] - ts_pred_temp;

                        /* normalized to z-scores */
                        for (b = 0; b < LASSO_BANDS - 1; b++)
                        {
                            if (i_b == lasso_blist[b])
                            {
                                /* mini rmse */
                                mini_rmse = max(rec_adj_mini[i_b], tmpcg_rmse[i_b]);

                                /* z-score */
                                v_diff[end][i_b] = (v_dif_mag[end][i_b]) / mini_rmse;
                            }
                        }

                    }
                    matlab_norm(v_diff[end], LASSO_BANDS, &vec_mag[conse]);
                    vec_mag[conse] *= vec_mag[conse];
                }
                break_mag = 0.0;
                for (m = 0; m < conse; m++)
                {
                    if (break_mag > vec_mag[m])
                        break_mag = vec_mag[m];
                }

                if (break_mag > t_cg)
                {
                    printf("Change Magnitude = %.2f\n",break_mag);

                    /* record break time */
                    rec_cg[num_fc].t_break = clrx[i+1];
                    rec_cg[num_fc].change_prob = 1.0;
                    for (i = 0; i < TOTAL_BANDS - 1; i++)
                        matlab_2d_array_mean(v_dif_mag, 0, conse, 
                            &rec_cg[num_fc].magnitude[i]);

                    /* identified and move on for the next functional curve */
                    num_fc++;
                    /* start from i+1 for the next functional curve */
                    i_start = i + 1;
                    /* start training again */
                    bl_train = 0;
                }
                else if (vec_mag[0] > t_max_cg)
                {
                    /* remove noise */
                    for (m = i; m < num_scenes; m++)
                    {
                        clrx[m] = clrx[m+1];
                        for (b = 0; b < TOTAL_BANDS - 1; b++)
                            clry[m][b] = clry[m+1][b];
                    }
                    end--; /* check if this is needed */

                    /* stay & check again after noise removal */
                    i--;

                }
                status = free_2d_array ((void **) v_dif_mag);
                if (status != SUCCESS)
                    RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                           FUNC_NAME, EXIT_FAILURE);
            }

            /* move forward to the i+1th clear observation */
            i++;
        }

        end = num_scenes;
        /* Two ways for processing the end of the time series */ 
        if (bl_train == 1)
        {
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
                        id_last = i_conse;
                        break;
                    }
                }
            } 

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
                                         &rec_cg[num_fc].magnitude[i]);

                /* median of the last < conse pixels */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                {
                    matlab_2d_int_partial_mean(clry, i_b, end-conse+id_last, end - 1, 
                          &fit_cft[0][i_b]);
                    matlab_2d_partial_square_mean(clry, i_b, end-conse+id_last, end - 1,  
                                                  &rmse[i_b]);
                }
                /* identified and move on for the next functional curve */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[end-conse+id_last+1];
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end];
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
                rec_cg[num_fc].category = v_qa + 1; /* mean value fit at the end */
                /* record change magnitude */
                for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
                rec_cg[num_fc].magnitude[i_b] = 0.0; /* record change magnitude */
                num_fc++;   
            }
        }
        else if (bl_train == 0)
        {
            /* if break find close to the end of the time series 
               Use [conse,min_num_c*n_times+conse) to fit curve */
            /* multitemporal cloud mask */
            status = auto_mask(clrx, clry, i_start, end,
                               (clrx[end]-clrx[i_start])/num_yrs, t_const, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask routine", 
                                  FUNC_NAME, FAILURE);

            /* update i_span after noise removal */
            i_span = 0;
            for (m = 0; m < end-conse; m++)
            {
                if (bl_ids[m] == 0)
                    i_span++;
            }

            /* Allocate memory for ids, rm_ids */
            ids = (int *)calloc(i+conse-1, sizeof(int));
            if (ids == NULL)
                RETURN_ERROR("ERROR allocating ids memory", 
                              FUNC_NAME, FAILURE);

            rm_ids = (int *)calloc(i+conse-1-i_start, sizeof(int));
            if (rm_ids == NULL)
                RETURN_ERROR("ERROR allocating rm_ids memory", 
                              FUNC_NAME, FAILURE);

            for (m = i_start; m <= i + conse - 1; m++)
            {
                ids[m] = m;
            }
            m= 0;
            for (k = 0; k < end-conse; k++)
            {
                if (bl_ids[k] == 0) 
                {            
                    clr_ids[k] = 1;
                }
                else
                {
                    clr_ids[k] = 0;
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
                    cpx[k] = cpx[k+1];
                    for (b = 0; b < TOTAL_BANDS - 1; b++)
                        cpy[k][b] = cpy[k][b];
                }
                end--;
            }

            
            for (i_b = 0; i_b < TOTAL_BANDS - 1; i_b++)
            {
                status = auto_ts_fit(clrx, clry, i_b, i_start, end, min_num_c, 
                                     fit_cft, &rmse[i_b], &v_dif_temp); 
                if (status != SUCCESS)  
                     RETURN_ERROR ("Calling auto_ts_fit8\n", FUNC_NAME, EXIT_FAILURE);
            }

            /* record time of curve start */
            rec_cg[num_fc].t_start = clrx[i_start];
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end];
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
            rec_cg[num_fc].category = v_qa + min_num_c; /* simple model fit at the end */
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
    free(bl_ids);
    free(clr_ids);
    free(vec_mag);
    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

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
            " --min_rmse=input minimum rmse value"
            " --t_cg=chi-square inversed T_cg value" 
            " --t_max_cg=chi-square inversed T_max_cg for " 
                         "last step noise removal"
            " --conse = number of points used for change detection" 
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    -row: input row number\n");
    printf ("    --col=input col number\n");
    printf ("    --min_rmse=input minimum rmse value\n");
    printf ("    --t_cg=chi-square inversed T_cg value\n");
    printf ("    --t_max_cg=chi-square inversed T_max_cg for " 
                         "last step noise removal\n");
    printf ("    --conse = number of points used for change detection\n");
    printf ("\n");
    printf ("where the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n");
    printf ("Example: ccdc"
            " --row=3845"
            " --col=2918"
            " --min_rmse=0.5"
            " --t_cg=15.0863"
            " --t_max_cg=35.8882"
            " --conse=6"
            " --verbose\n");
    printf ("Note: The ccdc must run from the directory"
            " where the input data are located.\n\n");
}
