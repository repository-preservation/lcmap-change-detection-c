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

#define NUM_LASSO_BANDS 5
#define TOTAL_IMAGE_BANDS 7
#define TOTAL_BANDS 8
#define MIN_NUM_C 4
#define MID_NUM_C 6
#define MAX_NUM_C 8
#define CONSE 6
#define N_TIMES 3 /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25 /* average number of days per year */
#define NUM_FC 10 /* Values change with number of pixels run */
#define T_CONST 3 /* Threshold for cloud, shadow, and snow detection */
#define MIN_YEARS 1 /* minimum year for model intialization */
#define T_SN 0.75        /* no change detection for permanent snow pixels */ 
#define T_CLR 0.25       /* Fmask fails threshold */
#define T_CG 15.0863     /* chi-square inversed T_cg for noise removal */
#define T_MAX_CG 35.8882 /* chi-square inversed T_max_cg for 
                            last step noise removal */
#define CFMASK_CLEAR 0
#define CFMASK_WATER 1
#define CFMASK_CLOUD 2
#define CFMASK_SNOW 3
#define CFMASK_SHADOW 4 

const char scene_list_name[] = {"scene_list.txt"};
int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 6}; /* This is band index */
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

The current category in output structure: 
first digit:          
0: normal model (no change)     
1: change at the beginning of time series model                  
2: change at the end of time series model                    
3: disturbance change in the middle               
4: fmask fail scenario
5: permanent snow scenario
second digit:
1: model has only constant term
4: model has 3 coefs + 1 const 
6: model has 5 coefs + 1 const
8: model has 7 coefs + 1 const
******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";
    char msg_str[MAX_STR_LEN];  /* input data scene name */
    //    char filename[MAX_STR_LEN];         /* input binary filenames */
    int status;                 /* return value from function call */
    Output_t *rec_cg = NULL;    /* output structure and metadata */
    bool verbose;               /* verbose flag for printing messages */
    int i, k, m, b, k_new;
    //    char **scene_list = NULL;
    FILE *fd;
    int num_scenes = MAX_SCENE_LIST;
    int num_c = 8;            /* max number of coefficients for the model */
    int num_fc = 0;           /* intialize NUM of Functional Curves */
    int rec_fc;
    float v_start[NUM_LASSO_BANDS];
    float v_end[NUM_LASSO_BANDS];
    float v_slope[NUM_LASSO_BANDS];
    float v_dif[NUM_LASSO_BANDS];
    float v_difs[CONSE][NUM_LASSO_BANDS];
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
    float **v_diff = NULL;
    float v_dif_norm;
    int i_count;
    float **v_dif_mag;
    int i_conse, i_b;
    float *vec_mag;
    float *vec_magg;
    float v_dif_mean;
    float **rec_v_dif;
    float **temp_v_dif;
    float adj_rmse[TOTAL_IMAGE_BANDS];
    float mini_rmse;
    int bl_tmask;
    int n_rmse;
    float tmpcg_rmse[TOTAL_IMAGE_BANDS];
    int d_rt;
    float *d_yr;
    float break_mag;
    int id_last;
    float ts_pred_temp;
    FILE *fp_bin_out;
    int ids_old_len;
    int i_break;
    int i_ini;
    int ini_conse;
    float vec_magg_min;
    int ids_len;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Read the command-line arguments */
    status = get_args (argc, argv, &row, &col, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    printf("row,col,verbose=%d,%d,%d\n",row,col,verbose);

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
    if (access(scene_list_name, F_OK) != -1) /* File exists */
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
    {
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
    }

    clrx = malloc(num_scenes * sizeof(int));
    if (clrx == NULL)
    {
        RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);
    }

    id_range = (int *)calloc(num_scenes, sizeof(int));
    if (id_range == NULL)
    {
        RETURN_ERROR("ERROR allocating id_range memory", FUNC_NAME, FAILURE);
    }

    id_clr = (int *)calloc(num_scenes, sizeof(int));
    if (id_clr == NULL)
    {
        RETURN_ERROR("ERROR allocating id_clr memory", FUNC_NAME, FAILURE);
    }

    id_all = (int *)calloc(num_scenes, sizeof(int));
    if (id_all == NULL)
    {
        RETURN_ERROR("ERROR allocating id_all memory", FUNC_NAME, FAILURE);
    }

    id_sn = (int *)calloc(num_scenes, sizeof(int));
    if (id_sn == NULL)
    {
        RETURN_ERROR("ERROR allocating id_sn memory", FUNC_NAME, FAILURE);
    }

    ids = (int *)calloc(num_scenes, sizeof(int));
    if (ids == NULL)
    {
        RETURN_ERROR("ERROR allocating ids memory", FUNC_NAME, FAILURE);
    }

    ids_old = (int *)calloc(num_scenes, sizeof(int));
    if (ids_old == NULL)
    {
        RETURN_ERROR("ERROR allocating ids_old memory", FUNC_NAME, FAILURE);
    }

    bl_ids = (int *)calloc(num_scenes, sizeof(int));
    if (bl_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating bl_ids memory", FUNC_NAME, FAILURE);
    }

    rm_ids = (int *)calloc(num_scenes, sizeof(int));
    if (rm_ids == NULL)
    {
        RETURN_ERROR("ERROR allocating rm_ids memory", FUNC_NAME, FAILURE);
    }

    clry = (int16 **) allocate_2d_array (num_scenes, TOTAL_IMAGE_BANDS, 
                                         sizeof (int16));
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

    rmse = malloc((TOTAL_IMAGE_BANDS) * sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    vec_mag = malloc(CONSE * sizeof(float));
    if (vec_mag == NULL)
    {
        RETURN_ERROR ("Allocating vec_mag memory", FUNC_NAME, FAILURE);
    }

    /* allocate memory for v_dif_mag */ 
    v_dif_mag = (float **) allocate_2d_array(CONSE,
                         TOTAL_IMAGE_BANDS, sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory", 
                                 FUNC_NAME, FAILURE);
    }

    /* Allocate memory for rec_v_dif */
    rec_v_dif = (float **)allocate_2d_array(num_scenes, TOTAL_IMAGE_BANDS,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    /* Allocate memory for temp_v_dif */
    temp_v_dif = (float **)allocate_2d_array(num_scenes, TOTAL_IMAGE_BANDS,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }
#if 0
    /* sort scene_list based on year & julian_day */
    status = sort_scene_based_on_year_doy(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, FAILURE);
    }
#endif
    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
    {
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
    }

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    //    status = read_envi_header(scene_list[0], meta);
    status = read_envi_header("LC80460272013120LGN01_MTLstack", meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling read_envi_header", 
                      FUNC_NAME, FAILURE);
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
    short int buf[num_scenes][TOTAL_IMAGE_BANDS];
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
            if (k != TOTAL_IMAGE_BANDS)
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
    int buf[num_scenes][TOTAL_BANDS + 1];
    fd = fopen("pixel_inputs.txt","r");
    if (fd == NULL)
        RETURN_ERROR ("Open input file", FUNC_NAME, FAILURE);

    for (i = 0; i < num_scenes; i++)
    {
        for (i_b = 0; i_b < TOTAL_BANDS + 1; i_b++)
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
	{
            fmask_sum++;
	}
    }
    if (fmask_sum < (int) 0.5 * num_scenes)
    {
        RETURN_ERROR ("Not enough clear-sky pixels", FUNC_NAME, FAILURE);
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
	{
            id_range[i] = 1;
	}
        else
	{
            id_range[i] = 0;
	}
    }

    /* Get each mask pixel totals */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] < CFMASK_CLOUD)
	{
            clr_sum++;
	}
        if (fmask_buf[i] < FILL_VALUE)
	{
            all_sum++;
	}
        if (fmask_buf[i] == CFMASK_SNOW)
	{
            sn_sum++;
	}
    }


    /* percent of clear pixels */    
    clr_pct = (float) clr_sum / (float) all_sum;

    /* percent of snow observations */
    sn_pct =  (float) sn_sum / (float)(clr_sum + sn_sum); 

    printf("clr_sum,all_sum,clr_pct=%d,%d,%f\n",clr_sum,all_sum,clr_pct);

    /* Allocate memory for rec_cg */ 
    rec_cg = malloc(NUM_FC * sizeof(Output_t));
    if (rec_cg == NULL)
    {
        RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);
    }

    /* fit permanent snow observations */
    if (clr_pct < T_CLR)
    {
        if (sn_pct > T_SN)
        {
            /* snow observations are "good" now */
            for (i = 0; i < num_scenes; i++)
            { 
                if (fmask_buf[i] == CFMASK_SNOW && id_range[i] == 1)
                {
                    clrx[n_sn] = sdate[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		    {
                        clry[n_sn][k] = buf[i][k];
		    }
                    n_sn++;
                }
            }  
            end = n_sn;

            if (n_sn < 1 ) /* not enough snow pixels */
            {
                RETURN_ERROR ("Not enough good snow observations\n", 
                     FUNC_NAME, FAILURE);
            }

            /* start model fit for snow persistent pixels */
            printf ("Fit permanent snow observations, now pixel = %f\n", 
                   100.0 * sn_pct); 


            /* the first observation for TSFit */
            i_start = 1; /* the first observation for TSFit */
#if 0
            /* identified and move on for the next curve */
            num_fc++; /* not needed, as array index starts with zero */
#endif
            if (n_sn < N_TIMES * MIN_NUM_C)
            {
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    matlab_2d_array_median(clry, i_b, end, &fit_cft[0][i_b]);
		    rmse_from_square_root_mean(clry, fit_cft[0][i_b], end, i_b, &rmse[i_b]); 
		}
                rec_cg[num_fc].category = 50 + 1; /* snow pixel */
	    }
            else
            {
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    status = auto_ts_fit(clrx, clry, i_b, 0, end-1, MIN_NUM_C, 
                                         fit_cft, &rmse[k], temp_v_dif); 
                    if (status != SUCCESS) 
		    { 
                        RETURN_ERROR ("Calling auto_ts_fit for permanent snow\n", 
                                      FUNC_NAME, FAILURE);
		    }
		}
                rec_cg[num_fc].category = 50 + MIN_NUM_C; /* snow pixel */
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
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < MIN_NUM_C; k++)
		{
                    /* record fitted coefficients */
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
		}
                /* record rmse of the pixel */
                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0; 
            /* record number of observations */
            rec_cg[num_fc].num_obs = n_sn; 
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                /* record change magnitude */ 
                rec_cg[num_fc].magnitude[i_b] = 0.0; 
            /* NUM of Fitted Curves (num_fc) */
            num_fc++;
	    if (num_fc >= 10)
	    {
                /* Reallocate memory for rec_cg */ 
		rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
		{
                    RETURN_ERROR("ERROR allocating rec_cg memory", 
                                 FUNC_NAME, FAILURE);
		}
            }   
        }
        else
        {
            /* snow observations are "good" now */
            for (i = 0; i < num_scenes; i++)
            { 
                if ((fmask_buf[i] == CFMASK_CLEAR || fmask_buf[i] == CFMASK_WATER)
                      && id_range[i] == 1)
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                         clry[n_clr][k] = buf[i][k];
                    n_clr++;
                }   
            }
            end = n_clr;

            /* no change detection for clear observations */
            if (n_clr < 1) /* not enough clear pixels */ 
            {
                RETURN_ERROR("Not enough good clear observations\n", 
                            FUNC_NAME, FAILURE);
            }

            /* start model fit for clear persistent pixels */
            printf ("Fmask failed, clear pixel = %f\n", 
                   100.0 * clr_pct); 

            /* the first observation for TSFit */
            i_start = 1; /* the first observation for TSFit */
#if 0
            /* identified and move on for the next curve */
            num_fc++; /* not needed, as array index starts with zero */
#endif
            if (n_clr < N_TIMES * MIN_NUM_C)
            {
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    matlab_2d_array_median(clry, i_b, end, &fit_cft[0][i_b]);
       	            rmse_from_square_root_mean(clry, fit_cft[0][i_b], end, i_b, &rmse[i_b]); 
		}
                rec_cg[num_fc].category = 40 + 1; /* clear pixel */
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
                rec_cg[num_fc].category = 40 + MIN_NUM_C; /* clear pixel */
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
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < MIN_NUM_C; k++)
		{
                    /* record fitted coefficients */
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
		}
                /* record rmse of the pixel */
                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0; 
            /* record number of observations */
            rec_cg[num_fc].num_obs = n_clr; 
            /* record fit category */
            rec_cg[num_fc].category = 40 + MIN_NUM_C; /* snow pixel */
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
	    {
                /* record change magnitude */ 
                rec_cg[num_fc].magnitude[i_b] = 0.0; 
	    }
            /* NUM of Fitted Curves (num_fc) */
            num_fc++;   
       	    if (num_fc >= 10)
	    {
                /* Reallocate memory for rec_cg */ 
                rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
		{
                    RETURN_ERROR("ERROR allocating rec_cg memory", 
                                 FUNC_NAME, FAILURE);
		}
            }
        }
    }
    else /* clear land or water pixels */
    {
        printf("Seasonal Snow (Snow < %f)\n", 100.0 * sn_pct);
        printf("Fmask works, clear pixels (land/water) = %f\n", 100.0 * clr_pct);

        for (i = 0; i < num_scenes; i++)
        { 
            if ((fmask_buf[i] == CFMASK_CLEAR || fmask_buf[i] == CFMASK_WATER) 
                && id_range[i] == 1)
            {
                clrx[n_clr] = sdate[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		{
                     clry[n_clr][k] = clry[i][k];
		}
                n_clr++;
            }   
        }
        end = n_clr;

        printf("end_clr=%d\n",end);
        /* calculate median variogram */
        status = median_variogram(clry, 0, end-1, TOTAL_IMAGE_BANDS, adj_rmse);
        if (status != SUCCESS)
	{
            RETURN_ERROR("ERROR calling median_variogram routine", FUNC_NAME, 
                         FAILURE);
	}
#if 0
        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
         printf("k,adj_rmse[k]=%d,%f\n",k,adj_rmse[k]);
#endif
        /* start with mininum requirement of clear obs */
        i = N_TIMES * MIN_NUM_C;

        /* the first observation for TSFit */
        i_start = 1; 

        /* record the start of the model initialization (0=>initial;1=>done) */
        bl_train = 0;

        /* record the num_fc at the beginning of each pixel */
        rec_fc = num_fc;

        /* record the start of Tmask (0=>initial;1=>done) */
        bl_tmask = 0;

        /* while loop - process till the last clear observation - CONSE */
        while (i <= end - CONSE)
        {
            /* span of "i" */
            i_span = i - i_start + 1;

            printf("end,CONSE=%d,%d\n",end,CONSE);
            printf("i_start,i,i_span=%d,%d,%d\n",i_start,i,i_span);

            /* span of time (num of years) */
            time_span = (float)(clrx[i-1] - clrx[i_start-1]) / NUM_YEARS;

            /* basic requrirements: 1) enough observations; 2) enough time */
            if ((i_span >= N_TIMES * MIN_NUM_C) && (time_span >= (float)MIN_YEARS))
            {
                /* initializing model */
                if (bl_train == 0)
                {
                    /* step 1: noise removal */ 
                    status = auto_mask(clrx, clry, i_start-1, i+CONSE-1,
                                   (float)(clrx[i+CONSE-1]-clrx[i_start-1]) / NUM_YEARS, 
                                   adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
                    if (status != SUCCESS)
		    {
                        RETURN_ERROR("ERROR calling auto_mask during model initilization", 
                                      FUNC_NAME, FAILURE);
		    }
#if 0
                    for(k = i_start-1; k < i+CONSE; k++)
                        printf("k,bl_ids[k]=%d,%d\n",k,bl_ids[k]);
#endif
                    /* Clears the IDs buffers */
                    for (k = 0; k < num_scenes; k++)
                        ids[k] = 0;

                    /* IDs to be removed */
                    for (k = i_start-1; k < i+CONSE; k++)
		    {
                        ids[k-i_start+1] = k;
		    }
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
                    if (i_span < (N_TIMES * MIN_NUM_C))
                    {
                        /* move forward to the i+1th clear observation */
                        i++;
                        printf("i++=%d\n",i);
                        /* not enough clear observations */
                        continue;
                    }

                    /* allocate memory for cpx, cpy */
                    cpx = malloc(end * sizeof(int));
                    if (cpx == NULL)
                        RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                    cpy = (int16 **) allocate_2d_array (end, TOTAL_IMAGE_BANDS, 
                                     sizeof (int16));
                    if (cpy == NULL)
                    {
                        RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                    }

                    /* copy clrx and clry to cpx and cpy */
                    for (k = 0; k < end; k++)
                    {
                        cpx[k] = clrx[k];
                        for (m = 0; m < TOTAL_IMAGE_BANDS; m++)
			{
                            cpy[k][m] = clry[k][m];
			}
                    }

                    printf("end1,rm_ids_len=%d,%d\n",end,rm_ids_len);
                    /* remove noise pixels between i_start & i */
                    m = 0;
                    for (k = 0, k_new=0; k < end; k++)
                    {
                        if (m < rm_ids_len && k == rm_ids[m])
                        {
                            m++;
                            continue;
                        }
                        cpx[k_new] = clrx[k];
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
			{
                            cpy[k_new][i_b] = clry[k][i_b];
			}
                        k_new++;
                    }
                    end = k_new;
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
                    time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                    /* check if there is enough time */
                    if (time_span < MIN_YEARS)
                    {
                        i = i_rec;   /* keep the original i */
                        /* move forward to the i+1th clear observation */
                        i++;        
                        printf("i++2=%d\n",i);
                        free(cpx);
                        status = free_2d_array ((void **) cpy);
                        if (status != SUCCESS)
			{
                              RETURN_ERROR ("Freeing memory: cpy\n", 
                                    FUNC_NAME, FAILURE);
			}
                        continue;    /* not enough time */
                    }

                    /* remove noise */
                    for (k = 0; k < end; k++)
                    {
                        clrx[k] = cpx[k];
                        for (m = 0; m < TOTAL_IMAGE_BANDS; m++)
			{
                            clry[k][m] = cpy[k][m];
			}
                    }

                    free(cpx);
                    status = free_2d_array ((void **) cpy);
                    if (status != SUCCESS)
		    {
                        RETURN_ERROR ("Freeing memory: cpy\n", 
                             FUNC_NAME, FAILURE);
		    }
#if 0
                    /* record the start of Tmask (0=>initial;1=>done) */
                    bl_tmask = 0;
#endif

                    /* Step 2: model fitting: initialize model testing variables
                               defining computed variables */
                    v_dif_norm = 0.0;
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
                        /* Initial model fit */
                        status = auto_ts_fit(clrx, clry, b, i_start-1, i-1, 
                                 MIN_NUM_C, fit_cft, &rmse[b], rec_v_dif); 
                        if (status != SUCCESS)  
			{
                            RETURN_ERROR ("Calling auto_ts_fit during model initilization\n", 
                                 FUNC_NAME, FAILURE);
			}
#if 0
                        for(k = 0; k < NUM_COEFFS; k++)
                        printf("lasso_blist[b],k,fit_cft[k][lasso_blist[b]],rmse[lasso_blist[b-1]]=%d,%d,%f,%f\n",lasso_blist[b],k,fit_cft[k][lasso_blist[b]],rmse[lasso_blist[b]]);
#endif
                    }

                    for(b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        /* calculate mini rmse */
                        mini_rmse = max(adj_rmse[lasso_blist[b]], rmse[lasso_blist[b]]);

                        /* compare the first observation */
                        v_start[lasso_blist[b]] = rec_v_dif[0][lasso_blist[b]] 
                                / mini_rmse;

                        /* compare the last clear observation */
                        v_end[lasso_blist[b]] = rec_v_dif[i-i_start][lasso_blist[b]]
                                                / mini_rmse;

                        /* anormalized slope values */
                        v_slope[lasso_blist[b]] = fit_cft[1][lasso_blist[b]] *
                                        (clrx[i-1]-clrx[i_start-1])/mini_rmse;
                            
                        /* difference in model intialization */
                        v_dif[lasso_blist[b]] = fabs(v_slope[lasso_blist[b]]) + 
                                                fabs(v_start[lasso_blist[b]]) + 
                                                fabs(v_end[lasso_blist[b]]);  
                        v_dif_norm += v_dif[lasso_blist[b]] * v_dif[lasso_blist[b]];               
                    }
                    printf("v_dif_norm=%f\n",v_dif_norm);

                    /* find stable start for each curve */
                    if (v_dif_norm > T_CG)
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

                    /* model ready */
                    bl_train = 1;

                    /* count difference of i for each iteration */
                    i_count = 0;

                    /* find the previous break point */
                    if (num_fc == rec_fc)
		    {
                        i_break = 1; /* first curve */
		    }
                    else
                    {
                        for (k = 0; k < num_scenes; k++) 
                        {
                            if (clrx[k] > rec_cg[num_fc-1].t_break)
                            {
                                i_break = k + 1;
                                break;
                            }
                        }
                        printf("i_start,i_break2=%d,%d\n",i_start,i_break);
                    }

                    if (i_start > i_break)
                    {
                        /* model fit at the beginning of the time series */
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
                            
                            /* allocate memory for model_v_dif */ 
                            v_dif_mag = (float **) allocate_2d_array(ini_conse,
                                              TOTAL_IMAGE_BANDS, sizeof (float));
                            if (v_dif_mag == NULL)
                            {
                                RETURN_ERROR ("Allocating v_dif_mag memory", 
                                              FUNC_NAME, FAILURE);
                            }

                            /* allocate memory for v_diff */ 
                            v_diff = (float **) allocate_2d_array(ini_conse,
                                        NUM_LASSO_BANDS, sizeof (float));
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
                               value of difference for CONSE obs
                               record the magnitude of change */
                            printf("ini_conse,i_ini,i=%d,%d,%d\n",ini_conse,i_ini,i);
                            vec_magg_min = 9999.0;
                            for (i_conse = 0; i_conse < ini_conse; i_conse++)
                            {
                                v_dif_norm = 0.0;
                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    /* absolute differences */
                                    auto_ts_predict(clrx, fit_cft, i_b, i_ini-i_conse,
                                                    i_ini-i_conse, &ts_pred_temp);
                                    v_dif_mag[i_conse][i_b] = (float)clry[i_ini-i_conse][i_b] - 
                                                       ts_pred_temp;

                                    /* normalize to z-score */
                                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                                    {
                                        if (i_b == lasso_blist[b])
                                        {
                                            /* minimum rmse */ 
                                            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);

                                            /* z-scores */
                                            v_diff[i_conse][i_b] = v_dif_mag[i_conse][i_b] 
                                                                          / mini_rmse;
                                            v_dif_norm += v_diff[i_conse][i_b] * v_diff[i_conse][i_b];
#if 0
                                            printf("=%d,%d,%d,%d,%f,%f,%f,%f\n",i_conse,i_b,i_ini,clry[i_ini-i_conse][i_b],ts_pred_temp,v_dif_mag[i_conse][i_b],mini_rmse,v_diff[i_conse][i_b]);
#endif
					    break;
                                        }
                                    }
                                }
                                vec_magg[i_conse] = v_dif_norm; 

                                printf("i_conse,vec_magg[i_conse]1 = %d,%f\n",i_conse,vec_magg[i_conse]);

                                if (vec_magg_min > vec_magg[i_conse])
				{
                                    vec_magg_min =  vec_magg[i_conse];
				}
                            }

                            /* change detection */
                            if (vec_magg_min > T_CG) /* change detected */
			    {
                                break;
			    }
                            else if (vec_magg[0] > T_MAX_CG) /*false change */
                            {
                                for (k = i_ini; k < end - 1; k++)
                                {
                                    clrx[k] = clrx[k+1];
                                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
				    {
                                        clry[k][b] = clry[k+1][b];
				    }
                                }
                                i--;
#if 0
                                /* remove noise pixels */
                                m = 0;
                                for (k = 0, k_new=0; k < end; k++)
                                {
                                    if (m < i_ini && k == i_ini)
                                    {
                                        m++;
                                        continue;
                                    }
                                    cpx[k_new] = clrx[k];
                                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			            {
                                        cpy[k_new][b] = clry[k][b];
			            }
                                    k_new++;
                                }
                                end = k_new;
                                i--;
#endif
				printf("end,i=%d,%d\n",end,i);
                            }

                            /* free the memory */
                            free(vec_magg);
                            status = free_2d_array ((void **) v_diff);
                            if (status != SUCCESS)
			    {
                                RETURN_ERROR ("Freeing memory: v_diff\n", 
                                              FUNC_NAME, FAILURE);
			    }

                            /* update i_start if i_ini is not a confirmed break */
                            i_start = i_ini;
                            printf("i_start=%d\n",i_start);

                            /* update curves */
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                status = auto_ts_fit(clrx, clry, i_b, i_start-1, i-1, 
						   MIN_NUM_C, fit_cft, &rmse[b], temp_v_dif); 
                                if (status != SUCCESS)  
				{
                                    RETURN_ERROR ("Calling auto_ts_fit at the beginning of "
                                        "the time series\n", FUNC_NAME, FAILURE);
				}
                            }
                        }
                    }

printf("i_start,i_break,conse3=%d,%d,%d\n",i_start,i_break,CONSE);
                    /* enough to fit simple model and confirm a break */
                    if ((i_start - i_break) >= CONSE)
                    {
                        /* defining computed variables */
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

                        /* record time of curve end */
                        rec_cg[num_fc].t_end = clrx[i_start-2]; 
                        /* record postion of the pixel */
                        rec_cg[num_fc].pos.row = row; 
                        /* record postion of the pixel */
                        rec_cg[num_fc].pos.col = col; 
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            for (k = 0; k < MIN_NUM_C; k++)
			    {
                                /* record fitted coefficients */
                                rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
			    }
                            /* record rmse of the pixel */                         
                            rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                        }
			/* recored break time */
			rec_cg[num_fc].t_break = clrx[i_start -1];
			/* record change probability */
			rec_cg[num_fc].change_prob = 1.0;

                        /* treat first curve different */
                        if (num_fc == rec_fc)                           
                        {
         		    /* record time of curve start */ 
         		    rec_cg[num_fc].t_start = clrx[0];
                            /* record fit category */
                            rec_cg[num_fc].category = 10 + MIN_NUM_C;
                            /* record change probability */
                            rec_cg[num_fc].change_prob = 1.0;
                        }
                        else
                        {
         		    /* record time of curve start */ 
         		    rec_cg[num_fc].t_start = rec_cg[num_fc-1].t_break;
                            /* record fit category */
                            rec_cg[num_fc].category = 30 + MIN_NUM_C;
                            /* record change probability */
                            rec_cg[num_fc].change_prob = 0.0;
                        } 
                        /* record number of observations */
                        rec_cg[num_fc].num_obs = i_start - i_break; 
                        /* record change magnitude */
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            matlab_2d_array_mean(v_dif_mag, i_b, ini_conse, 
                                                 &v_dif_mean);  
                            rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                        }
                        /* identified and move on for the next functional curve */
                        num_fc++;                                      
                        if (num_fc >= 10)
		        {
                            /* Reallocate memory for rec_cg */ 
		            rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                            if (rec_cg == NULL)
			    {
                                RETURN_ERROR("ERROR allocating rec_cg memory", 
                                             FUNC_NAME, FAILURE);
			    }
		        }
                    }
                    else if (i_start > i_break)
                    {
                        /* median value fit for the rest of the pixels < CONSE */
                        for (i_b = 0; i_b< TOTAL_IMAGE_BANDS; i_b++)
                        {
                            matlab_int_2d_partial_mean(clry, i_b, i_break-1, i_start-1, 
                                  &fit_cft[0][i_b]);
                            partial_square_root_mean(clry, i_b, i_break-1, i_start-1, 
                                  fit_cft, &rmse[i_b]);
                        }

                        /* record time of curve end */  
                        rec_cg[num_fc].t_end = clrx[i_start-2]; 
                        /* record fitted coefficients */
                        rec_cg[num_fc].pos.row = row; 
                        /* record fitted coefficients */
                        rec_cg[num_fc].pos.col = col; 
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            for (k = 0; k < MAX_NUM_C; k++)
			    {
                                /* record fitted coefficients */
                                rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
			    } 
                            /* record rmse of the pixel */
                            rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                        }
			/* record break time */
			rec_cg[num_fc].t_break = clrx[i_start - 1];
			/* record change probability */
			rec_cg[num_fc].change_prob = (float)(i_start-i_break) / (float)CONSE;
                        /* treat first curve different */
                        if (num_fc == rec_fc)                           
                        {
                            /* record time of curve start */
               		    rec_cg[num_fc].t_start = clrx[0];
                            /* record fit category */
                            rec_cg[num_fc].category = 10 + 1;
                            /* record break time */
                            rec_cg[num_fc].t_break = clrx[i_start-1];
                            /* record change probability */
                            rec_cg[num_fc].change_prob = (float)(i_start-i_break)/(float)CONSE;
                        }
                        else
                        {
                            /* record time of curve start */
               		    rec_cg[num_fc].t_start = rec_cg[num_fc-1].t_break;
                            /* record fit category */
                            rec_cg[num_fc].category = 30 + 1;
                            /* record change probability */
                            rec_cg[num_fc].change_prob = 0.0;
                        } 
                        /* record number of observations */
                        rec_cg[num_fc].num_obs = i_start-i_break; 
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            matlab_2d_array_mean(v_dif_mag, i_b, ini_conse, &v_dif_mean);
                            /* record change magnitude */ 
                            rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                        }

                        /* NUM of Fitted Curves (num_fc) */
                        num_fc++;
                        if (num_fc >= 10)
		        {
                            /* Reallocate memory for rec_cg */ 
		            rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                            if (rec_cg == NULL)
			    {
                                RETURN_ERROR("ERROR allocating rec_cg memory", 
                                             FUNC_NAME, FAILURE);
			    }
                        }           
                    }
                } /* end of initializing model */
 
                /* continuous monitoring started!!! */
                if (bl_train == 1)
                {
                    /* Clears the IDs buffers */
                    for (k = 0; k < num_scenes; k++)
		    {
                        ids[k] = 0;
		    }

                    /* all IDs */
                    ids_len = 0;
                    for (k = i_start-1; k < i; k++)
                    {
                        ids[k-i_start+1] = k;
                        ids_len++;
                    }
                    i_span = i - i_start +1;

                printf("i_start, i, up-update i_span=%d,%d,%d\n",
                       i_start, i, i_span);

                    /* determine the time series model */
                    update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C, 
                              num_c, &update_num_c);

                printf("update_cft=%d\n",update_num_c);
                    /* dynamic model fit when there are not many obs */
                    //                if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C))
                    if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                    {
                       /* update i_count at each iteration */
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
#if 0
                            for(k = 0; k < i-i_start+1; k++)
                             printf("i_b,k,rec_v_dif[k][i_b]1=%d,%d,%f\n",i_b,k,rec_v_dif[k][i_b]);
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
                       for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                       {
                           for (k = 0; k < update_num_c; k++)
			   {
                               /* record fitted coefficients */
                               rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
			   }
                           /* record rmse of the pixel */
                           rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                       }
                       /* record change probability */
                       rec_cg[num_fc].change_prob = 0.0; 
                       /* record number of observations */
                       rec_cg[num_fc].num_obs = i-i_start+1; 
                       /* record fit category */
                       rec_cg[num_fc].category = 0 + update_num_c; 
                       for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                       {
                           /* record change magnitude */
                           rec_cg[num_fc].magnitude[i_b] = 0.0; 
                       }
                           
                       /* detect change, value of difference for CONSE obs */
                       for (i_conse = 0; i_conse < CONSE; i_conse++)
                       {
                           v_dif_norm = 0.0;
                           for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                           {
                               /* absolute differences */
                               auto_ts_predict(clrx, fit_cft, i_b, i+i_conse, i+i_conse, 
                                    &ts_pred_temp);
                               v_dif_mag[i_conse][i_b] = (float)clry[i+i_conse][i_b] - ts_pred_temp; 
       
                               /* normalize to z-score */
                               for (b = 0; b < NUM_LASSO_BANDS; b++)
                               {
                                   if (i_b == lasso_blist[b])
                                   {
                                       /* minimum rmse */ 
                                       mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
 
                                       /* z-scores */
                                       v_difs[i_conse][i_b] = v_dif_mag[i_conse][i_b] / mini_rmse;
                                       v_dif_norm += v_difs[i_conse][i_b] * v_difs[i_conse][i_b];
#if 0
printf("=%d,%d,%d,%f,%f,%f,%f\n",i_conse,i_b,clry[i+i_conse][i_b],ts_pred_temp,
       v_dif_mag[i_conse][i_b],mini_rmse,v_difs[i_conse][i_b]);
#endif
                                        break;
                                   }
                               }
                            }
                            vec_mag[i_conse] = v_dif_norm;

                            printf("i_conse,vec_mag[i_conse]2=%d,%f\n",i_conse,vec_mag[i_conse]);

                        }

                        /* Clears the IDsOld buffers */
                        for (k = 0; k < ids_len; k++)
			{
                            ids_old[k] = 0;
			}

                        /* IDs that haven't been updated */
                        for (k = 0; k < ids_len; k++)
			{
                            ids_old[k] = ids[k];
			}
                        ids_old_len = ids_len;

                    printf("length(ids)1,length(ids_old)1,i-i_start+1=%d,%d,%d\n",
                    ids_len,ids_old_len,i-i_start+1);
                    }
                    else
                    {
                        if ((clrx[i-1] - clrx[i_start-1]) >= (i_count + i_count / 3))
                        {
                            /* update i_count at each iteration year */
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
#if 0
                                for(k = i_start-1; k < i-1; k++)
                                 printf("i_b,rec_v_dif[k][i_b]2=%d,%f\n",i_b,rec_v_dif[k][i_b]);
#endif
                            }

                            /* record fitted coefficients */
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                for (k = 0; k < update_num_c; k++)
				{
                                    /* record fitted coefficients */
                                    rec_cg[num_fc].coefs[i][k] = fit_cft[i_b][k];
				} 
                                /* record rmse of the pixel */
                                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                            }
                            /* record number of observations */
                            rec_cg[num_fc].num_obs = i-i_start+1; 
                            /* record fit category */
                            rec_cg[num_fc].category = 0 + update_num_c; 

                            /* Clears the IDsOld buffers */
                            for (k = 0; k < ids_len; k++)
			    {
                                ids_old[k] = 0;
			    }

                            /* IDs that haven't been updated */
                            for (k = 0; k < ids_len; k++)
			    {
                                ids_old[k] = ids[k];
			    }
                            ids_old_len = ids_len;

                            printf("length(ids)2,length(ids_old)2,i-i_start+1=%d,%d,%d\n",
                               ids_len,ids_old_len,i-i_start+1);
                        }

                        /* record time of curve end */
                        rec_cg[num_fc].t_end = clrx[i-1]; /* record time of curve end */

                        printf("length(ids_len)3, length(ids_old_len)3, i- i_start +1=%d,%d,%d\n",
                               ids_len, ids_old_len, i - i_start +1);
#if 0
                        /* use temporally-adjusted RMSE */
                        if (ids_old_len <= N_TIMES * MAX_NUM_C)
                        {
                            /* number of observations for calculating RMSE */
                            n_rmse = ids_old_len;
                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			    {
                                tmpcg_rmse[b] = rmse[b];
			    }
                        }
                        else
                        {
#endif
                        /* use fixed number for RMSE computing */
                        n_rmse = N_TIMES * MAX_NUM_C;

                        /* better days counting for RMSE calculating */
                        /* relative days distance */
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
#if 0
                            printf("m,d_rt,d_yr[m]=%d,%d,%f\n",m,d_rt,d_yr[m]);
#endif
                        }

                        printf("length(IDs),length(IDsOld)=%d,%d\n",  ids_len, ids_old_len); 
                        /* sort the d_yr */
                        qsort(d_yr, ids_old_len, sizeof(float), cmpfunc);
#if 0
                       for(m = 0; m < ids_old_len; m++)
                         printf("m,ids_old[m],d_yr[m]=%d,%d,%f\n",m,ids_old[m],d_yr[m]);
#endif

                       for(b = 0; b < TOTAL_IMAGE_BANDS; b++)
                           tmpcg_rmse[b] = 0.0;

                        /* temporarily changing RMSE */
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            matlab_2d_array_norm(rec_v_dif, lasso_blist[b], n_rmse,
                                             &tmpcg_rmse[lasso_blist[b]]);
                            tmpcg_rmse[lasso_blist[b]] /= sqrt(n_rmse - update_num_c);
#if 0
                            for (m = 0; m < n_rmse; m++)
                         printf("m,b,rec_v_dif[m]===%d,%d,%f\n",m,b,rec_v_dif[m][lasso_blist[b]]);
#endif
                            printf("====%d,%f\n",b,tmpcg_rmse[lasso_blist[b]]); 

                        }

                        /* free allocated memories */
                        free(d_yr);

                        /* move the ith col to i-1th col */
                        for (m = 0; m < CONSE-1; m++)
                        {
                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                            {
                                v_difs[m][b] = v_difs[m+1][b];
                                v_dif_mag[m][b] = v_dif_mag[m+1][b];
                            }
                            vec_mag[m] = vec_mag[m+1];
#if 0
                     printf("m,vec_mag[m]3=%d,%f\n",m,vec_mag[m]);
#endif
                            v_difs[CONSE-1][b] = 0.0;
                            v_dif_mag[CONSE-1][b] = 0.0;
                        }
                        vec_mag[CONSE-1] = 0.0;

                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            /* absolute difference for all bands */
                            auto_ts_predict(clrx, fit_cft, i_b, i+CONSE-1, 
                                 i+CONSE-1, &ts_pred_temp);
                            v_dif_mag[CONSE-1][i_b] = (float)clry[i+CONSE-1][i_b] - ts_pred_temp;
#if 0
                            printf("===%d,%d,%f,%f\n",i_b,clry[i+CONSE-1][i_b],ts_pred_temp,
                              v_dif_mag[CONSE-1][i_b]);
#endif
                            /* normalized to z-scores */
                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /* mini rmse */
                                    mini_rmse = max(adj_rmse[i_b], tmpcg_rmse[i_b]);

                                    /* z-score */
                                    v_difs[CONSE-1][i_b] = v_dif_mag[CONSE-1][i_b] / mini_rmse;
                                    vec_mag[CONSE-1] += v_difs[CONSE-1][i_b] * v_difs[CONSE-1][i_b]; 
				    break;
                                }
                            }

                        }
printf("m,vec_mag[conse-1]3=%d,%f\n",CONSE-1,vec_mag[CONSE-1]);
                    }
                    break_mag = 9999.0;
                    for (m = 0; m < CONSE; m++)
                    {
                        printf("m,vec_mag[m]4=%d,%f\n",m,vec_mag[m]);
                        if (break_mag > vec_mag[m])
			{
                            break_mag = vec_mag[m];
			}
                    }

                    printf("break_mag,vec_mag[0],T_MAX_CG=%f,%f,%f\n",break_mag,vec_mag[0],T_MAX_CG);
                    if (break_mag > T_CG)
                    {
                        printf("Change Magnitude = %.2f\n",break_mag);

                        /* record break time */
                        rec_cg[num_fc].t_break = clrx[i];
                        rec_cg[num_fc].change_prob = 1.0;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            matlab_2d_array_mean(v_dif_mag, 0, CONSE, 
                                &rec_cg[num_fc].magnitude[i_b]);

                        /* identified and move on for the next functional curve */
                        num_fc++;
                        if (num_fc >= 10)
		        {
                            /* Reallocate memory for rec_cg */ 
		            rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                            if (rec_cg == NULL)
			    {
                                RETURN_ERROR("ERROR allocating rec_cg memory", 
                                                     FUNC_NAME, FAILURE);
			    }
		        }

                        /* start from i+1 for the next functional curve */
                        i_start = i + 1;
                        printf("i_start2=%d\n",i_start);
                        /* start training again */
                        bl_train = 0;
                    }
                    else if (vec_mag[0] > T_MAX_CG)
                    {
#if 0
                        /* remove noise */
                        m = 0;
                        for (k = 0, k_new=0; k < end; k++)
                        {
                            if (m < i && k == i)
                            {
                                m++;
                                continue;
                            }
                            cpx[k_new] = clrx[k];
                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			    {
                                cpy[k_new][b] = clry[k][b];
			    }
                            k_new++;
                        }
                        end = k_new;
#endif
                        for (m = i; m < num_scenes - 1; m++)
                        {
                            clrx[m] = clrx[m+1];
                            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			    {
                                clry[m][b] = clry[m+1][b];
			    }
                        }
                        end--; /* check if this is needed */

                        /* stay & check again after noise removal */
                        i--;
                        printf("i--=%d\n",i);
                    }
                }    
            } /* end of continuous monitoring */ 
            /* move forward to the i+1th clear observation */
            i++;
            printf("i++4=%d\n",i);
        } /* end of "while (i <= end - CONSE) */

        /* Two ways for processing the end of the time series */ 
        if (bl_train == 1)
        {
            end = num_scenes - 1;
            /* if no break find at the end of the time series,
               define probability of change based on CONSE */
            for (i_conse = CONSE - 1; i_conse >= 0; i_conse--)
            {
                if (vec_mag[i_conse] <= T_CG)
                {
                    /* the last stable id */
                    id_last = i_conse + 1;
                    break;
                }
            } 

            printf("id_last = %d\n",id_last);
            /* update change probability */
            rec_cg[num_fc].change_prob = (CONSE - id_last) / CONSE; 
            /* update end time of the curve */
            rec_cg[num_fc].t_end = clrx[end - CONSE + id_last];

            /* mean value fit for the rest of the pixels < CONSE & > 1 */
            if (CONSE > id_last)
            {
                /* update time of the probable change */
                rec_cg[num_fc].t_break = clrx[end-CONSE+id_last+1];
                /* update magnitude of change */
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
		{
                    matlab_2d_partial_mean(v_dif_mag, i_b, id_last, CONSE-1, 
                                         &rec_cg[num_fc].magnitude[i_b]);
		}

                /* median of the last < CONSE pixels */
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    matlab_int_2d_partial_mean(clry, i_b, end-CONSE+id_last, end - 1, 
                          &fit_cft[0][i_b]);
                    matlab_2d_partial_square_mean(clry, i_b, end-CONSE+id_last, end - 1,  
                                                  &rmse[i_b]);
                }
                /* identified and move on for the next functional curve */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[end-CONSE+id_last];
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end-1];
                /* record break time */
                rec_cg[num_fc].t_break = 0;
                /* record postion of the pixel */
                rec_cg[num_fc].pos.row = row;
                rec_cg[num_fc].pos.col = col;
                /* record fitted coefficients */
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < update_num_c; k++)
		    {
                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
		    }
                    rec_cg[num_fc].rmse[i_b] = rmse[i_b];
                }
                /* record change probability */
                rec_cg[num_fc].change_prob = 0.0;
                /* record number of observations */
                rec_cg[num_fc].num_obs = CONSE - id_last;
                /* record fit category */
                rec_cg[num_fc].category = 20 + 1; /* mean value fit at the end */
                /* record change magnitude */
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
		{
                    rec_cg[num_fc].magnitude[i_b] = 0.0; /* record change magnitude */
		}
                num_fc++;   
                if (num_fc >= 10)
		{
                    /* Reallocate memory for rec_cg */ 
		    rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
		    {
                        RETURN_ERROR("ERROR allocating rec_cg memory", 
                                     FUNC_NAME, FAILURE);
		    }
		}
            }
        }
        else if (bl_train == 0)
        {
            /* if break found close to the end of the time series 
               Use [CONSE,MIN_NUM_C*N_TIMES+CONSE) to fit curve */
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
                     if (clrx[k] > rec_cg[num_fc-1].t_break)
		     {
                         i_start  = k + 1;
		     }
                }
                printf("update i_start3=%d\n",i_start);
            }

            for (m = 0; m < num_scenes; m++)
	    {
                bl_ids[m] = 0;
	    }

	    //            get_ids_length(clrx, 0, num_scenes-1, &end);
            printf("i_start,end=%d,%d\n",i_start,end);
            /* multitemporal cloud mask */
            status = auto_mask(clrx, clry, i_start-1, end-1,
                               (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS, 
                               adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
            if (status != SUCCESS)
                RETURN_ERROR("ERROR calling auto_mask at the end of time series", 
                                  FUNC_NAME, FAILURE);

            /* update i_span after noise removal */
            i_span = 0;
            for (m = i_start-1; m < end; m++)
            {
                if (bl_ids[m] == 0)
		{
                    i_span++;
		}
            }

            for (m = 0; m < num_scenes-1; m++)
	    {
                ids[m] = 0;
	    }

            ids_len = 0;
            for (m = i_start-1; m <= end - 1; m++)
            {
                ids[m-i_start+1] = m;
                ids_len++;
            }
            m= 0;
            for (k = 0; k < end-CONSE; k++)
            {
                if (bl_ids[k] == 1) 
                {
                    rm_ids[m] = ids[k];
                    printf("m,rm_ids[m]=%d,%d\n",m,rm_ids[m]);
                    m++;
                }
            }
            rm_ids_len = m;

            printf("rm_ids_len, update i_span2 =%d,%d\n", rm_ids_len, i_span);
            /* remove noise pixels between i_start & i */
            m = 0;
            for (k = 0, k_new=0; k < end; k++)
            {
                if (m < rm_ids_len && k == rm_ids[m])
                {
                    m++;
                    continue;
                }
                clrx[k_new] = clrx[k];
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
		{
                    clry[k_new][i_b] = clry[k][i_b];
		}
                k_new++;
            }
            end = k_new;

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

	    /* record time of curve start */
	    if (num_fc == rec_fc)
	    {
                rec_cg[num_fc].t_start = clrx[0];
	    }
	    else
	    {
                rec_cg[num_fc].t_start = rec_cg[num_fc-1].t_break;
	    }
            /* record time of curve end */
            rec_cg[num_fc].t_end = clrx[end-1];
            /* record break time */
            rec_cg[num_fc].t_break = 0;
            /* record postion of the pixel */
            rec_cg[num_fc].pos.row = row;
            rec_cg[num_fc].pos.col = col;
            /* record fitted coefficients */
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
		{
                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
		}
                rec_cg[num_fc].rmse[i_b] = rmse[i_b];
            }
            /* record change probability */
            rec_cg[num_fc].change_prob = 0.0;
            /* record number of observations */
            rec_cg[num_fc].num_obs = i_span;
            /* record fit category */
            rec_cg[num_fc].category = 20 + MIN_NUM_C; /* simple model fit at the end */
            /* record change magnitude */
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
	    {
                rec_cg[num_fc].magnitude[i_b] = 0.0; /* record change magnitude */
	    }
            num_fc++;
            if (num_fc >= 10)
	    {
                /* Reallocate memory for rec_cg */ 
		rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
		{
                    RETURN_ERROR("ERROR allocating rec_cg memory", 
                                 FUNC_NAME, FAILURE);
		}
	    }
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
    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n", 
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n", 
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                   FUNC_NAME, FAILURE);
    }
#if 0
    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      FAILURE);
    }
#endif
    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      FAILURE);
    }

    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      FAILURE);
    }

    /* Output rec_cg structure to the output file 
       note: can use fread to read out the structure from the output file */
    fp_bin_out = fopen("output.bin", "wb");
    if (fp_bin_out == NULL)
    {
        RETURN_ERROR ("Opening output.bin file\n", FUNC_NAME,
                      FAILURE);
    }
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
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
		{
                    printf("i_b,k,rec_cg[i].coefs[i_b][k] = %d,%d,%f\n", 
                         i_b,k,rec_cg[i].coefs[i_b][k]); 
		}
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

    return SUCCESS;
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/15/2015   Song Guo         Original Development

******************************************************************************/
void
usage ()
{
    printf ("Continuous Change Detection and Classification\n");
    printf ("\n");
    printf ("usage:\n");
    printf ("ccdc"
            " --row=<input row number>"
            " --col=<input col number>"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --row=: input row number\n");
    printf ("    --col=: input col number\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("ccdc"
            " --row=5099"
            " --col=3191"
            " --verbose\n");
    printf ("Note: The ccdc must run from the directory"
            " where the input data are located.\n\n");
}
