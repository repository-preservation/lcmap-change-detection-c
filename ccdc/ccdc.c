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
#define T_CONST 4.89 /* Threshold for cloud, shadow, and snow detection */
#define MIN_YEARS 1 /* minimum year for model intialization */
#define T_SN 0.75        /* no change detection for permanent snow pixels */ 
#define T_CLR 0.25       /* Fmask fails threshold */
#define T_CG 15.0863     /* chi-square inversed T_cg (0.99) for noise removal */
#define T_MAX_CG 35.8882 /* chi-square inversed T_max_cg (1e-6) for 
                            last step noise removal */
#define CFMASK_CLEAR 0
#define CFMASK_WATER 1
#define CFMASK_CLOUD 2
#define CFMASK_SNOW 3
#define CFMASK_SHADOW 4 

const char scene_list_name[] = {"scene_list.txt"};
int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5}; /* This is band index */
int cmpfunc (const void * a, const void * b)
{
   return ( *(float*)a - *(float*)b );
}

char *sub_string
(
    const char *source,
    size_t start,
    size_t length
) 
{
    int i;
    char *target;

    target = malloc(length*sizeof(char));

    for(i = 0; i != length; ++i) 
    {
        target[i] = source[start + i];
    }
    target[i] = 0;
    return target;
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
20151203    Brian Davis      Added arguments for input and output file
                             locations and optional scene list file.
20150104    Song Guo         Numerous bug fixes.

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

Note: The commented out parts of code is the inputs using ESPA putputs,
      the current input part is only for reading inputs from Zhe's code
******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";       /* for printing error messages          */
    char msg_str[MAX_STR_LEN];       /* input data scene name                */
    char filename[MAX_STR_LEN];      /* input binary filenames               */
    int status;                      /* return value from function call      */
    Output_t *rec_cg = NULL;         /* output structure and metadata        */
    bool verbose;                    /* verbose flag for printing messages   */
    int i, k, m, b, k_new;           /* loop counters                        */
    char **scene_list = NULL;        /* 2-D array for list of scene IDs      */
    FILE *fd;                        /* file descriptor for file             */
                                     /* containing scene names               */
    int num_scenes = MAX_SCENE_LIST; /* number of input scenes defined       */
    int num_c = 8;                   /* max number of coefficients for model */
    int num_fc = 0;                  /* intialize NUM of Functional Curves   */
    int rec_fc;
    float v_start[NUM_LASSO_BANDS];
    float v_end[NUM_LASSO_BANDS];
    float v_slope[NUM_LASSO_BANDS];
    float v_dif[NUM_LASSO_BANDS];
    float **v_diff;
    int *sdate;
    Input_meta_t *meta;
    int row, col;
    int landsat_number;
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
    float **clry;
    int *cpx;
    float **cpy;
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
    float v_dif_norm = 0.0;
    int i_count;
    float **v_dif_mag;
    float **v_diff_mag;
    int i_conse, i_b;
    float *vec_mag;
    float *vec_magg;
    float v_dif_mean;
    float vec_magg_min;
    float **rec_v_dif;
    float **temp_v_dif;
    float adj_rmse[TOTAL_IMAGE_BANDS];
    float mini_rmse;
    int bl_tmask;
    int n_rmse;
    float tmpcg_rmse[NUM_LASSO_BANDS];
    int d_rt;
    float *d_yr;
    int id_last;
    float ts_pred_temp;
    FILE *fp_bin_out;
    int ids_old_len;
    int i_break;
    int i_ini;
    int ini_conse;
    float break_mag;
    int ids_len;
    unsigned char *fmask_buf;
    int **buf;
    FILE ***fp_bin;
    int clr_land_water_counter;
    int fmask_total_counter;

    char inDir[MAX_STR_LEN];    /* directory location of input data/files     */
    char outDir[MAX_STR_LEN];   /* directory location for output files        */
    char sceneListFileName[MAX_STR_LEN];   /* file name containing list of input sceneIDs */
    char sceneList[MAX_STR_LEN]; /* optional input argument for file of list of scenes */
    char tmpstr[MAX_STR_LEN];   /* char string for text manipulation */
    int len;                    /* return of strlen for manipulating strings */
    char outputBinary[MAX_STR_LEN]; /* directory and file name for output.bin */
    char hdr_name[MAX_STR_LEN]; /* string for header file to read, different for LC8 */


    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Initialize the input and output directory specification, */
    /* because they are optional.                               */
    strcpy(inDir, "");
    strcpy(outDir, "");

    /* Read the command-line arguments */
    status = get_args (argc, argv, &row, &col, inDir, outDir, sceneList, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    printf("row,col,verbose=%d,%d,%d\n",row,col,verbose);
    printf("Input Directory:  %s\n", inDir);
    printf("Output Directory: %s\n", outDir);
    printf("sceneList: %s\n", sceneList);

    /* allocate memory for scene_list */
    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                         sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR ("Allocating scene_list memory", FUNC_NAME, FAILURE);
    }
//#if 0
    /* check if scene_list.txt file exists, if not, create the scene_list
       from existing files in the current data working directory */
    //if (access(scene_list_name, F_OK) != -1) /* File exists */
    if (access(sceneList, F_OK) != 0) /* File does not exist */
    {
        strcpy (sceneListFileName, inDir);
        strcat(sceneListFileName, "/scene_list.txt");
        if (access(sceneListFileName, F_OK) != -1) /* File exists */
        {
            num_scenes = MAX_SCENE_LIST;
        }
        else /* Default File exists */
        {
            status = create_scene_list(hdr_name, num_scenes, sceneListFileName, scene_list);
            if(status != SUCCESS)
            RETURN_ERROR("Running create_scene_list file", FUNC_NAME, FAILURE);
        }
    }
    else /* File exists */
    {
        strcpy(sceneListFileName, sceneList);
    }

    //fd = fopen("scene_list.txt", "r");
    fd = fopen(sceneListFileName, "r");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_scenes; i++)
    {
        //if (fscanf(fd, "%s", scene_list[i]) == EOF)
        if (fscanf(fd, "%s", tmpstr) == EOF)
            break;
        strcpy(scene_list[i], inDir);
        strcat(scene_list[i], "/");
        strcat(scene_list[i], tmpstr);
    }
    num_scenes = i;

    /* allocate memory for fp_bin, need testing */
    fp_bin = (FILE ***) allocate_2d_array (TOTAL_BANDS, num_scenes,
                                         sizeof (FILE*));
    if (fp_bin == NULL)
    {
        RETURN_ERROR ("Allocating fp_bin memory", FUNC_NAME, FAILURE);
    }

//#endif
#if 0
    num_scenes = 455;
#endif

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

    fmask_buf = malloc(num_scenes * sizeof(unsigned char));
    if (fmask_buf == NULL)
    {
        RETURN_ERROR("ERROR allocating fmask_buf memory", FUNC_NAME, FAILURE);
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

    buf = (int **) allocate_2d_array (TOTAL_BANDS, num_scenes,
                                         sizeof (int));
    if (buf == NULL)
    {
        RETURN_ERROR ("Allocating buf memory", FUNC_NAME, FAILURE);
    }

    clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, num_scenes,
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

    /* allocate memory for v_dif_mag */ 
    v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, CONSE,
                sizeof (float));
    if (v_dif_mag == NULL)
    {
        RETURN_ERROR ("Allocating v_dif_mag memory", 
                                 FUNC_NAME, FAILURE);
    }

    /* Allocate memory for rec_v_dif */
    rec_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, num_scenes,
                                     sizeof (float));
    if (rec_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating rec_v_dif memory",FUNC_NAME, FAILURE);
    }

    /* Allocate memory for temp_v_dif */
    temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, num_scenes,
                                     sizeof (float));
    if (temp_v_dif == NULL)
    {
        RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
    }
//#if 0
    /* sort scene_list based on year & julian_day */
    printf("num_scenes %d\n", num_scenes);
    printf("scene_list[0]=%s\n", scene_list[0]);
    status = sort_scene_based_on_year_doy(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, FAILURE);
    }
//#endif
    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
    {
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
    }

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    status = read_envi_header(scene_list[0], meta);
    //status = read_envi_header("LC80460272013120LGN01_MTLstack", meta);
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
//#if 0
    /* Open input files 
       Note: may need switch num_scenes and TOTAL_IMAGE_BANDS in buf/fp_bin */
    //    FILE *fp_bin[TOTAL_BANDS][num_scenes];
    //    short int buf[TOTAL_IMAGE_BANDS][num_scenes];
    //    unsigned char fmask_buf[num_scenes];
    /* Open input files */
    for (k = 0; k < TOTAL_BANDS; k++)
    {
        for (i = 0; i < num_scenes; i++)
        {
            //landsat_number = atoi(sub_string(scene_list[i],2,1));
            len = strlen(scene_list[i]);
            landsat_number = atoi(sub_string(scene_list[i],(len-19),1));
            if (landsat_number != 8)
            {
                if (k == 5)
		  sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+2);
                else if (k == 6)
                    sprintf(filename, "%s_toa_band6.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+1);
            }
            else
            {
                if (k == 6)
                    sprintf(filename, "%s_toa_band10.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else 
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+2);
            }

            fp_bin[k][i] = open_raw_binary(filename,"rb");
            if (fp_bin[k][i] == NULL)
                printf("error open %d scene, %d bands files\n",i, k+1);
            if (k != 7)
            {
                fseek(fp_bin[k][i], (row * meta->samples + col)*sizeof(short int), SEEK_SET);
                if (read_raw_binary(fp_bin[k][i], 1, 1,
                    sizeof(short int), &buf[k][i]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
            }
            else 
            {
                fseek(fp_bin[k][i], (row * meta->samples + col)*sizeof(unsigned char), 
                    SEEK_SET);
                if (read_raw_binary(fp_bin[k][i], 1, 1,
                    sizeof(unsigned char), &fmask_buf[i]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
            }
            close_raw_binary(fp_bin[k][i]);
        }

	/* Eliminate scenes that have less than 20% clear-sky pixels */
	clr_land_water_counter = 0;
	fmask_total_counter = 0;
	if (fmask_buf[i] == 0 || fmask_buf[i] == 1)
	{
            clr_land_water_counter++;
	}
	if (fmask_buf[i] > 0 && fmask_buf[i] < 255)
	{
	    fmask_total_counter++;
	}

	if (((float) clr_land_water_counter / (float)fmask_total_counter) < 0.2)
	{
            for (k = i; k < num_scenes; k++)
	    {
                strcpy(scene_list[k], scene_list[k+1]);
	    }
	    num_scenes--;
	}
    }
    snprintf (msg_str, sizeof(msg_str), "CCDC read_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);

//#endif
#if 0
    /* Temporary code for testing purpose */
    //    unsigned char fmask_buf[num_scenes];
    //    int buf[TOTAL_BANDS + 1][num_scenes];
    fd = fopen("pixel_inputs.txt","r");
    if (fd == NULL)
        RETURN_ERROR ("Open input file", FUNC_NAME, FAILURE);

    for (i = 0; i < num_scenes; i++)
    {
        for (i_b = 0; i_b < TOTAL_BANDS + 1; i_b++)
        {
            if (i_b == 0)
            {
                fscanf(fd, "%d", &sdate[i]);
                clrx[i] = sdate[i];
            }
            else if (i_b == 8)
            {
                fscanf(fd, "%d\n", &buf[i_b-1][i]);
                fmask_buf[i] = (unsigned char)buf[i_b-1][i];
            }
            else
            {
                fscanf(fd,"%d", &buf[i_b-1][i]);
                clry[i_b-1][i] = (float)buf[i_b-1][i];
            }
        }
    }
#endif
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

    snprintf (msg_str, sizeof(msg_str), "CCDC check_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);

    /* CHANGE: need change back to 0-6 from 1-7 if original 
       inputs are used */

    /* pixel value ranges should follow physical rules */
    for (i = 0; i < num_scenes; i++)
    { 
        /* convert Kelvin to Celsius (for new espa data) */
        if (buf[6][i] != -9999)
            buf[6][i] = (int16)(buf[6][i] * 10 - 27315);

        if ((buf[0][i] > 0) && (buf[0][i] < 10000) &&
            (buf[1][i] > 0) && (buf[1][i] < 10000) &&
            (buf[2][i] > 0) && (buf[2][i] < 10000) &&
            (buf[3][i] > 0) && (buf[3][i] < 10000) &&
            (buf[4][i] > 0) && (buf[4][i] < 10000) &&
            (buf[5][i] > 0) && (buf[5][i] < 10000) &&
            (buf[6][i] > -9320) && (buf[6][i] < 7070))
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
        if (fmask_buf[i] < 2)
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
    if ((clr_sum + sn_sum) != 0)
        sn_pct =  (float) sn_sum / (float)(clr_sum + sn_sum); 
    else
        sn_pct = (float)sn_sum;

    printf("clr_sum,all_sum,sn_sum,clr_pct,sn_pct=%d,%d,%d,%f,%f\n",
            clr_sum,all_sum,sn_sum,clr_pct,sn_pct);

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
	        if (((fmask_buf[i] == CFMASK_SNOW) || (fmask_buf[i] < 2)) 
                     && id_range[i] == 1)
                {
                    clrx[n_sn] = sdate[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		    {
         	        clry[k][n_sn] = (float)buf[k][i];
		    }
                    n_sn++;
                }
            }  
            end = n_sn;

            free(fmask_buf);
            status = free_2d_array ((void **) buf);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: buf\n", 
                              FUNC_NAME, FAILURE);
            }

            if (n_sn < N_TIMES * MIN_NUM_C) /* not enough snow pixels */
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
            num_f++; /* not needed, as array index starts with zero */
#endif
            /* treat saturated and unsaturated pixels differently */
            for (i = 0; i < end; i++)
            {
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                {
         	    i_span = 0;
                    if (k != TOTAL_IMAGE_BANDS - 1) /* for optical bands */
                    {
                        if (clry[i][k] > 0.0 && clry[i][k] < 10000.0)
                        {
                            clrx[i_span] = clrx[i];
                            clry[i_span][k] = clry[i][k];
                            i_span++;
                        }

                        if (i_span < MIN_NUM_C * N_TIMES)
                            fit_cft[i][k] = 10000; /* fixed value for saturated pixels */
                        else
                        {
                            status = auto_ts_fit(clrx, clry, k, 0, i_span-1, MIN_NUM_C, 
                                     fit_cft, &rmse[k], temp_v_dif); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit1\n", 
                                       FUNC_NAME, EXIT_FAILURE);

                        } 
                    }
                    else /* for thermal band */
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
#if 0
            if (n_sn < N_TIMES * MIN_NUM_C)
            {
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
		    matlab_2d_float_median(clry, i_b, end, &fit_cft[i_b][0]);
		    rmse_from_square_root_mean(clry, fit_cft[i_b][0], i_b, end, &rmse[i_b]); 
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
#endif
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
            rec_cg[num_fc].category = 50 + MIN_NUM_C; /* snow pixel */
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
                if (id_range[i] == 1)
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		      clry[k][n_clr] = (float)buf[k][i];
                    n_clr++;
                }   
            }
            end = n_clr;

            free(fmask_buf);
            status = free_2d_array ((void **) buf);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: buf\n", 
                              FUNC_NAME, FAILURE);
            }

            /* start model fit for clear persistent pixels */
            printf ("Fmask failed, clear pixel = %f\n", 
                   100.0 * clr_pct); 

	    n_clr = 0;
	    float band2_median;
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

	    /* the first observation for TSFit */
            i_start = 1; /* the first observation for TSFit */
#if 0
                /* identified and move on for the next curve */
            num_fc++; /* not needed, as array index starts with zero */
#endif
            if (n_clr < N_TIMES * MIN_NUM_C)
            {
                RETURN_ERROR("Not enough good clear observations\n", 
                            FUNC_NAME, FAILURE);
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
            rec_cg[num_fc].category = 40 + MIN_NUM_C; 
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
            if ((fmask_buf[i] < 2) && (id_range[i] == 1))
            {
                clrx[n_clr] = sdate[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		{
		  clry[k][n_clr] = (float)buf[k][i];
		}
                n_clr++;
            }   
        }
        end = n_clr;
        printf("end_clr=%d\n",end);

        free(fmask_buf);
        status = free_2d_array ((void **) buf);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: buf\n", 
                          FUNC_NAME, FAILURE);
        }

        /* calculate median variogram */
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

        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
         printf("k,adj_rmse[k]=%d,%f\n",k,adj_rmse[k]);

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
        snprintf (msg_str, sizeof(msg_str), "CCDC init_time=%s\n", ctime (&now));
            LOG_MESSAGE (msg_str, FUNC_NAME);

        while (i <= end - CONSE)
        {
            /* span of "i" */
            i_span = i - i_start + 1;

            // bdavis debug
            //printf("end,CONSE=%d,%d\n",end,CONSE);
            //printf("i_start,i,i_span=%d,%d,%d\n",i_start,i,i_span);

            /* span of time (num of years) */
            time_span = (float)(clrx[i-1] - clrx[i_start-1]) / NUM_YEARS;

            /* basic requrirements: 1) enough observations; 2) enough time */
            if ((i_span >= N_TIMES * MIN_NUM_C) && (time_span >= (float)MIN_YEARS))
            {
                /* initializing model */
                if (bl_train == 0)
                {
#if 0
		  for (k = i_start-1; k <  i+CONSE-1; k++)
		    printf("i,clrx[i],clry[i]=%d,%d,%f,%f,%f,%f,%f,%f\n",i,clrx[k],clry[0][k],clry[1][k],
                           clry[2][k],clry[3][k],clry[4][k],clry[5][k],clry[6][k]);
#endif
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
                            // bdavis debug
                            // printf("m,rm_ids[m]=%d,%d\n",m,rm_ids[m]);
                            m++;
                        }
                        else
                            i_span++;  /* update i_span after noise removal */
                    }

                    rm_ids_len = m;
                    /* update i_span after noise removal */
                    // i_span = i - i_start +1 - rm_ids_len;

                    // bdavis debug
                    // printf("rm_ids_len, update i_span = %d,%d\n", rm_ids_len, i_span);
                    /* check if there is enough observation */
                    if (i_span < (N_TIMES * MIN_NUM_C))
                    {
                        /* move forward to the i+1th clear observation */
                        i++;
                        // bdavis debug
                        // printf("i++=%d\n",i);
                        /* not enough clear observations */
                        continue;
                    }

		    if (end == 0)
                        RETURN_ERROR("No available data point", FUNC_NAME, FAILURE);

                    /* allocate memory for cpx, cpy */
                    cpx = malloc(end * sizeof(int));
                    if (cpx == NULL)
                        RETURN_ERROR("ERROR allocating cpx memory", FUNC_NAME, FAILURE);

                    cpy = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, end,
                                     sizeof (float));
                    if (cpy == NULL)
                    {
                        RETURN_ERROR ("Allocating cpy memory", FUNC_NAME, FAILURE);
                    }

                    // bdavis debug
                    // printf("end1,rm_ids_len=%d,%d\n",end,rm_ids_len);
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
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			{
                            cpy[b][k_new] = clry[b][k];
			}
                        k_new++;
                    }
                    end = k_new;
                    // bdavis debug
                    printf("end2=%d\n",end);

                    /* record i before noise removal 
                       This is very important as if model is not initialized 
                       the multitemporal masking shall be done again instead 
                       of removing outliers in every masking */
                    i_rec=i;

                    /* update i afer noise removal (i_start stays the same) */
                    i=i_start + i_span - 1;
                    // bdavis debug
                    printf("update i = %d\n",i);

                    /* update span of time (num of years) */
                    time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                    /* check if there is enough time */
                    if (time_span < MIN_YEARS)
                    {
                        i = i_rec;   /* keep the original i */
                        /* move forward to the i+1th clear observation */
                        i++;        
                        // bdavid debug
                        // printf("i++2=%d\n",i);
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

                    /* Step 2: model fitting: initialize model testing variables
                               defining computed variables */
                    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    {
#if 0
			 for (k = i_start-1; k < i-1; k++)
			   printf("i_b,k,clrx[k],clry[b][k]9=%d,%d,%d,%f\n",i_b,k,clrx[k],clry[b][k]);
#endif
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
                        //   bdavis question is rmse[b]9 correct?
                        printf("b,k,fit_cft[b][k],rmse[b]9=%d,%d,%f,%f\n",b,k,fit_cft[b][k],rmse[b]);
#endif
                    }

                    v_dif_norm = 0.0;
                    for(b = 0; b < NUM_LASSO_BANDS; b++)
                    {
                        /* calculate mini rmse */
                        mini_rmse = max(adj_rmse[lasso_blist[b]], rmse[lasso_blist[b]]);

                        /* compare the first observation */
                        v_start[b] = rec_v_dif[lasso_blist[b]][0] 
                                / mini_rmse;

                        /* compare the last clear observation */
                        v_end[b] = rec_v_dif[lasso_blist[b]][i-i_start]
                                                / mini_rmse;

                        /* anormalized slope values */
                        v_slope[b] = fit_cft[lasso_blist[b]][1] *
                                        (clrx[i-1]-clrx[i_start-1])/mini_rmse;
                            
                        /* difference in model intialization */
                        v_dif[b] = fabs(v_slope[b]) + fabs(v_start[b]) + fabs(v_end[b]); 
#if 0
			printf("lasso_blist[b],v_start[b],v_end[b],v_slope[b]=%d,%f,%f,%f\n",
			       lasso_blist[b],v_start[b],v_end[b],v_slope[b]); 
#endif
                        v_dif_norm += v_dif[b] * v_dif[b];               
                    }
                    // bdavis debug
                    // printf("v_dif_norm=%f\n",v_dif_norm);

                    /* find stable start for each curve */
                    if (v_dif_norm > T_CG)
                    {
                        /* start from next clear obs */
                        i_start++;
                        // bdavis debug
                        // printf("i_start1=%d\n",i_start);
                        
                        /* move forward to the i+1th clear observation */
                        i++;
                        // bdavis debug
                        // printf("i_start++,i++3=%d,%d\n",i_start,i);

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
		        /* after the first curve, compare rmse to determine which 
                           curve to determine t_break */
                        for (k = 0; k < end; k++) 
                        {
                            if (clrx[k] >= rec_cg[num_fc-1].t_break)
                            {
                                i_break = k + 1;
                                break;
                            }
                        }
                    }

                    // bdavis debug
		    // printf("rec_cg[num_fc-1].t_break=%d\n",rec_cg[num_fc-1].t_break);
		    // printf("i_start,i_break = %d,%d\n",i_start,i_break);
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

                            // bdavis debug
			    printf("ini_conse, i_start, i_break, %d %d %d\n",ini_conse,i_start,i_break);
			    if (ini_conse == 0)
                            {
                                RETURN_ERROR ("No data point for model fit at "
                                      "the begining", FUNC_NAME, FAILURE);
                            }

                            /* allocate memory for model_v_dif */ 
                            v_diff_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, 
                                         ini_conse, sizeof (float));
                            if (v_diff_mag == NULL)
                            {
                                RETURN_ERROR ("Allocating v_diff_mag memory", 
                                      FUNC_NAME, FAILURE);
                            }
   
                            /* allocate memory for v_diff */ 
                            v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                        ini_conse, sizeof (float));
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
                            // bdavis debug
                            // printf("ini_conse,i_ini=%d,%d\n",ini_conse,i_ini);
                            vec_magg_min = 9999.0;
                            for (i_conse = 0; i_conse < ini_conse; i_conse++)
                            {
                                v_dif_norm = 0.0;
                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    /* absolute differences */
        			    auto_ts_predict(clrx, fit_cft, MIN_NUM_C, i_b, i_ini-i_conse,
                                                    i_ini-i_conse, &ts_pred_temp);
                                    v_dif_mag[i_b][i_conse] = (float)clry[i_b][i_ini-i_conse] - 
                                                       ts_pred_temp;
#if 0
				    printf("i_b,i_conse,v_dif_mag[i_b][i_conse],ts_pred_temp=%d,%d,%f,%f\n",
					   i_b,i_conse,v_dif_mag[i_b][i_conse],ts_pred_temp);
#endif
                                    /* normalize to z-score */
                                    for (b = 0; b < NUM_LASSO_BANDS; b++)
                                    {
                                        if (i_b == lasso_blist[b])
                                        {
                                            /* minimum rmse */ 
                                            mini_rmse = max(adj_rmse[i_b], rmse[i_b]);

                                            /* z-scores */
                                            v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] 
                                                                          / mini_rmse;
                                            v_dif_norm += v_diff[b][i_conse] * v_diff[b][i_conse];
#if 0
                                            printf("i_conse,i_b,i_ini,clry[i_b][i_ini-i_conse],ts_pred_temp,v_dif_mag[i_b][i_conse],mini_rmse,v_diff[b][i_conse]=%d,%d,%d,%f,%f,%f,%f,%f\n",i_conse,i_b,i_ini,clry[i_b][i_ini-i_conse],ts_pred_temp,v_dif_mag[i_b][i_conse],mini_rmse,v_diff[b][i_conse]);
#endif
                                        }
                                    }
                                }
                                vec_magg[i_conse] = v_dif_norm; 

                                // bdavis debug
                                // printf("i_conse,vec_magg[i_conse]1 = %d,%f\n",i_conse,vec_magg[i_conse]);

                                if (vec_magg_min > vec_magg[i_conse])
				{
                                    vec_magg_min =  vec_magg[i_conse];
				}
                            }

                            // bdavis debug
			    // printf("vec_magg_min,vec_magg[0]=%f,%f\n",vec_magg_min,vec_magg[0]);
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
                                        clry[b][k] = clry[b][k+1];
				    }
                                }
                                i--;
                                // bdavis debug
				// printf("i--,end=%d,%d\n",i,end);

                                /* update i_start if i_ini is not a confirmed break */
                                i_start = i_ini;
                                // bdavis debug
                                // printf("up-update i_start=%d\n",i_start);
			    }

                            /* free the memory */
                            free(vec_magg);
                            status = free_2d_array ((void **) v_diff);
                            if (status != SUCCESS)
			    {
                                RETURN_ERROR ("Freeing memory: v_diff\n", 
                                              FUNC_NAME, FAILURE);
			    }
#if 0
                                /* remove noise pixels */
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
                                        cpy[i_b][k_new] = clry[i_b][k];
                                    k_new++;
                                }
                                end = k_new;
                                i--;

                            /* update i_start if i_ini is not a confirmed break */
                            i_start = i_ini;
                            printf("up-update i_start=%d\n",i_start);
#endif
                        }
                    }

// bdavis debug
// printf("i_start,i_break3=%d,%d,%d\n",i_start,i_break);
                    /* enough to fit simple model and confirm a break */
                    if ((num_fc == rec_fc) && ((i_start - i_break) >= CONSE))
		    {
                        /* defining computed variables */
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
#if 0
			 for (k = i_start-1; k < i-1; k++)
			   printf("i_b,k,clrx[k],clry[i_b][k]8=%d,%d,%f,%f\n",i_b,k,clrx[k],clry[i_b][k]);
#endif
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
			/* recored fit category */
			rec_cg[num_fc].category = 10 + MIN_NUM_C;
			/* record change probability */
			rec_cg[num_fc].change_prob = 1.0;
			/* record time of curve start */
			rec_cg[num_fc].t_start = clrx[0];
			/* record number of observations */
			rec_cg[num_fc].num_obs = i_start - i_break;
                        /* record change magnitude */
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            matlab_2d_float_median(v_dif_mag, i_b, ini_conse, 
                                                  &v_dif_mean);
                            // bdavis debug
			    // printf("i_b, ini_conse, v_dif_mean =%d,%d,%f\n", i_b, ini_conse, v_dif_mean);
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
#if 0
		    else if (clrx[i_start-1] < rec_cg[num_fc-1].t_break)
		    {
         		int min_num;
	        	float min_time;
         		/* update or remove num_fc-1 curve */
		        num_fc--;
                        /* two situations */
                        if (num_fc == rec_fc)                           
                        {
         		    min_num = CONSE;
			    min_time = 0.0;
			    qa = 10;
                        }
                        else
                        {
         		    min_num = N_TIMES * MIN_NUM_C;
			    min_time = MIN_YEARS * NUM_YEARS;
			    qa = 10;
                        } 

			/* update if these requirements are met */
			if ((i_start-i_break) >= min_num &&
                            (float)(clrx[i_start-2] - clrx[i_break-1]) >= min_time) 
			{
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

                            /* record the curve start */
                            rec_cg[num_fc].t_start = clrx[i_break-1]; 
                            /* record the curve end */
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
			    /* recored fit category */
			    rec_cg[num_fc].category = qa + MIN_NUM_C;
			    /* record change probability */
			    rec_cg[num_fc].change_prob = 1.0;
                            /* record number of observations */
                            rec_cg[num_fc].num_obs = i_start - i_break; 
                            /* record change magnitude */
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                matlab_2d_float_median(v_diff_mag, i_b, ini_conse, 
                                                       &v_dif_mean);  
                                rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                            }
                            /* identified and move on for the next functional curve */
                            num_fc++;  
			    printf("num_fc1=%d\n",num_fc);                                    
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
                        status = free_2d_array ((void **) v_diff_mag);
                        if (status != SUCCESS)
		        {
                            RETURN_ERROR ("Freeing memory: v_diff_mag\n", 
                                              FUNC_NAME, FAILURE);
		        }
		    }
                    snprintf (msg_str, sizeof(msg_str), "CCDC end_init=%s\n", ctime (&now));
                       LOG_MESSAGE (msg_str, FUNC_NAME);

#endif
		} /* end of initializing model */
 
                /* allocate memory for v_diff */ 
                v_diff = (float **) allocate_2d_array(NUM_LASSO_BANDS,
                                  CONSE, sizeof (float));
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

                    // bdavis debug
                    // printf("i_start, i, up-update i_span=%d,%d,%d\n",
                    //   i_start, i, i_span);

                    /* determine the time series model */
                    update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C, 
                              num_c, &update_num_c);

                    // bdavis debug
                    // printf("update_num_c=%d\n",update_num_c);
                    /* initial model fit when there are not many obs */
                    //                if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C))
                    if (i_count == 0 || i_span <= (N_TIMES * MAX_NUM_C))
                    {
                       /* update i_count at each iteration */
                       i_count = clrx[i-1] - clrx[i_start-1];

                        // bdavis debug
		       // printf("update i_count=%d\n",i_count);
                       for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                       {
#if 0
			 for (k = i_start-1; k < i-1; k++)
			   printf("i_b,k,clrx[k],clry[i_b][k]5=%d,%d,%f,%f\n",i_b,k,clrx[k],clry[i_b][k]);
#endif
                           status = auto_ts_fit(clrx, clry, i_b, i_start-1, i-1, update_num_c, 
                                             fit_cft, &rmse[i_b], rec_v_dif); 
                           if (status != SUCCESS) 
			   { 
                               RETURN_ERROR ("Calling auto_ts_fit during continuous monitoring\n", 
                                     FUNC_NAME, FAILURE);
			   }
#if 0
                            for(k = 0; k < i-i_start+1; k++)
                             printf("i_b,k,fit_cft[i_b][k],rmse[i_b],rec_v_dif[i_b][k]5=%d,%d,%f,%f,%f\n",
                                    i_b,k,rmse[i_b],rec_v_dif[i_b][k],rec_v_dif[i_b][k]);
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
        		       auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i+i_conse, i+i_conse, 
                                    &ts_pred_temp);
                               v_dif_mag[i_b][i_conse] = (float)clry[i_b][i+i_conse] - ts_pred_temp; 
       
                               /* normalize to z-score */
                               for (b = 0; b < NUM_LASSO_BANDS; b++)
                               {
                                   if (i_b == lasso_blist[b])
                                   {
                                       /* minimum rmse */ 
                                       mini_rmse = max(adj_rmse[i_b], rmse[i_b]);
 
                                       /* z-scores */
                                       v_diff[b][i_conse] = v_dif_mag[i_b][i_conse] / mini_rmse;
                                       v_dif_norm += v_diff[b][i_conse] * v_diff[b][i_conse];
#if 0
printf("i_conse,i_b,clry[i_b][i+i_conse],ts_pred_temp,"
       "v_dif_mag[i_b][i_conse],mini_rmse,v_diff[b][i_conse]=%d,%d,%d,%f,%f,%f,%f,%f\n",
       i_conse,i_b,clry[i_b][i+i_conse],ts_pred_temp,
       v_dif_mag[i_b][i_conse],mini_rmse,v_diff[b][i_conse]);
#endif
                                   }
                               }
                            }
                            vec_mag[i_conse] = v_dif_norm;

                            // bdavis debug
                            // printf("i_conse,vec_mag[i_conse]2=%d,%f\n",i_conse,vec_mag[i_conse]);

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

                        // bdavis debug
                        // printf("length(ids)1,length(ids_old)1,i-i_start+1=%d,%d,%d\n",
                        //       ids_len,ids_old_len,i-i_start+1);
                    }
                    else
                    {
        		if ((float)(clrx[i-1] - clrx[i_start-1]) >= (1.33*(float)i_count))
                        {
                            /* update i_count at each iteration year */
                            i_count = clrx[i-1] - clrx[i_start-1];

                            // bdavis debug
			    // printf("update i_count2=%d\n",i_count);
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
                                 printf("i_b,rec_v_dif[i_b][k]2=%d,%f\n",i_b,rec_v_dif[i_b][k]);
#endif
                            }

                            /* record fitted coefficients */
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

                            // bdavis debug
                            // printf("length(ids)2,length(ids_old)2,i-i_start+1=%d,%d,%d\n",
                            //   ids_len,ids_old_len,i-i_start+1);
                        }

                        /* record time of curve end */
                        rec_cg[num_fc].t_end = clrx[i-1]; /* record time of curve end */

                        // bdavis debug
                        // printf("length(ids_len)3, length(ids_old_len)3, i- i_start +1=%d,%d,%d\n",
                        //       ids_len, ids_old_len, i - i_start +1);

                        /* use fixed number for RMSE computing */
                        n_rmse = N_TIMES * rec_cg[num_fc].category;

                        /* better days counting for RMSE calculating */
                        /* relative days distance */
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
#if 0
                            printf("m,d_rt,d_yr[m]=%d,%d,%f\n",m,d_rt,d_yr[m]);
#endif
                        }

                        // bdavis debug
                        // printf("length(IDs),length(IDsOld)=%d,%d\n",  ids_len, ids_old_len); 
                        /* sort the d_yr */
                        qsort(d_yr, ids_old_len, sizeof(float), cmpfunc);
#if 0
                        for(m = 0; m < ids_old_len; m++)
                            printf("m,ids_old[m],d_yr[m]=%d,%d,%f\n",m,ids_old[m],d_yr[m]);
#endif

                        for(b = 0; b < NUM_LASSO_BANDS; b++)
                            tmpcg_rmse[b] = 0.0;

                        /* temporarily changing RMSE */
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            matlab_2d_array_norm(rec_v_dif, lasso_blist[b], n_rmse,
                                             &tmpcg_rmse[lasso_blist[b]]);
                            tmpcg_rmse[lasso_blist[b]] /= sqrt(n_rmse - update_num_c);

                         // bdavis debug
                         //   for (m = 0; m < n_rmse; m++)
                         // printf("m,b,rec_v_dif[lasso_b][m]===%d,%d,%f\n",m,b,rec_v_dif[lasso_blist[b]][m]);
                         // printf("b,tmpcg_rmse[lasso_blist[b]]=%d,%f\n",b,tmpcg_rmse[lasso_blist[b]]); 

                        }

                        /* free allocated memories */
                        free(d_yr);

                        /* move the ith col to i-1th col */
                        for (b = 0; b < NUM_LASSO_BANDS; b++)
                        {
                            for (m = 0; m < CONSE-1; m++)
                            {
                                v_diff[b][m] = v_diff[b][m+1];
                                v_dif_mag[b][m] = v_dif_mag[b][m+1];
                                vec_mag[m] = vec_mag[m+1];
                            }

                            v_diff[b][CONSE-1] = 0.0;
                            v_dif_mag[b][CONSE-1] = 0.0;
                        }
                        vec_mag[CONSE-1] = 0.0;

                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        {
                            /* absolute difference for all bands */
        		    auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i+CONSE-1, 
                                 i+CONSE-1, &ts_pred_temp);
                            v_dif_mag[i_b][CONSE-1] = clry[i_b][i+CONSE-1] - ts_pred_temp;
#if 0
                            printf("===%d,%d,%f,%f\n",i_b,clry[i_b][i+CONSE-1],ts_pred_temp,
                              v_dif_mag[i_b][CONSE-1]);
#endif
                            /* normalized to z-scores */
                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                if (i_b == lasso_blist[b])
                                {
                                    /* mini rmse */
                                    mini_rmse = max(adj_rmse[i_b], tmpcg_rmse[i_b]);

                                    /* z-score */
                                    v_diff[b][CONSE-1] = v_dif_mag[i_b][CONSE-1] / mini_rmse;
                                    vec_mag[CONSE-1] += v_diff[b][CONSE-1] * v_diff[b][CONSE-1]; 
                                    // bdavis debug
				    // printf("b,v_diff[b][CONSE-1] = %d,%f\n",b,v_diff[b][CONSE-1]);
                                }         
                            }
                        }
// bdavis debug
// printf("m,vec_mag[CONSE-1]3=%d,%f\n",CONSE-1,vec_mag[CONSE-1]);
		    }

                    break_mag = 9999.0;
                    for (m = 0; m < CONSE; m++)
                    {
                        // bdavis debug
                        // printf("m,vec_mag[m]4=%d,%f\n",m,vec_mag[m]);
                        if (break_mag > vec_mag[m])
			{
                            break_mag = vec_mag[m];
			}
                    }

                    // bdavis debug
                    // printf("break_mag,vec_mag[0],T_MAX_CG=%f,%f,%f\n",break_mag,vec_mag[0],T_MAX_CG);
                    if (break_mag > T_CG)
                    {    
                        // bdavis debug
        		// printf("Change Magnitude = %.2f\n", break_mag); 

                        /* record break time */
                        rec_cg[num_fc].t_break = clrx[i];
                        rec_cg[num_fc].change_prob = 1.0;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
			{
                            matlab_2d_array_mean(v_dif_mag, 0, CONSE, 
                                &rec_cg[num_fc].magnitude[i_b]);
                            // bdavis debug
			    // printf("i_b, ini_conse, rec_cg[num_fc].magnitude[i_b] =%d,%d,%f\n", 
                            //      i_b, ini_conse, rec_cg[num_fc].magnitude[i_b]);
			}
                        /* identified and move on for the next functional curve */
                        num_fc++;
                        // bdavis debug
			// printf("num_fc3=%d\n",num_fc);                                    
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
                        // bdavis debug
                        // printf("i_start2=%d\n",i_start);
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
                            if (m < rm_ids_len && k == rm_ids[m])
                            {
                                m++;
			        continue;
                            }
                            cpx[k_new] = clrx[k];
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                cpy[i_b][k_new] = clry[i_b][k];
                            k_new++;
                        }
                        end = k_new;
#endif
                        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                        {
                            for (m = i; m < num_scenes - 1; m++)
			    {
                                clrx[m] = clrx[m+1];
                                clry[b][m] = clry[b][m+1];
			    }
                        }
                        end--; /* check if this is needed */

                        /* stay & check again after noise removal */
                        i--;
                        // bdavis debug
                        // printf("i--2=%d\n",i);
                    }
                    status = free_2d_array ((void **) v_diff);
                    if (status != SUCCESS)
		    {
                        RETURN_ERROR ("Freeing memory: v_diff\n", 
                                      FUNC_NAME, FAILURE);
		    }
		} /* end of continuous monitoring */ 
	    }  /* end of checking basic requrirements */ 
            /* move forward to the i+1th clear observation */
            i++;
            // bdavis debug
            // printf("i++4=%d\n",i);
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

            // bdavis debug
            // printf("id_last = %d\n",id_last);
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
                            // bdavis debug
			    // printf("i_b,rec_cg[num_fc].magnitude[i_b] =%d,%f\n", i_b, rec_cg[num_fc].magnitude[i_b]);
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
            }
            else
            {
                for (k = 0; k < num_scenes; k++) 
                {
                     if (clrx[k] > rec_cg[num_fc-1].t_break)
		     {
                         i_start  = k + 1;
			 break;
		     }
                }
                // bdavis debug
                // printf("num_fc,end,rec_cg[num_fc-1].t_break,rec_cg[num_fc-1].t_end=%d,%d,%d,%d\n",
		//       num_fc,end,rec_cg[num_fc-1].t_break,rec_cg[num_fc-1].t_end);
            }
            // bdavis debug
            // printf("update i_start2=%d\n",i_start);

            for (m = 0; m < num_scenes; m++)
	    {
                bl_ids[m] = 0;
	    }

	    //            get_ids_length(clrx, 0, num_scenes-1, &end);
            // bdavis debug
            // printf("i_start,end=%d,%d\n",i_start,end);
	    if ((end - i_start + 1) > CONSE)
	    {
                /* multitemporal cloud mask */
                status = auto_mask(clrx, clry, i_start-1, end-1,
                               (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS, 
                               adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
                if (status != SUCCESS)
                    RETURN_ERROR("ERROR calling auto_mask at the end of time series", 
                                  FUNC_NAME, FAILURE);

                /* Clears the IDs buffers */
                for (m = 0; m < num_scenes-1; m++)
        	{
                    ids[m] = 0;
	        }

                /* IDs to be removed */
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
                        // bdavis debug
                        // printf("m,rm_ids[m]=%d,%d\n",m,rm_ids[m]);
                        m++;
                    }
                    else
                        i_span++;  /* update i_span after noise removal */
                }
                rm_ids_len = m;

                // bdavis debug
                // printf("rm_ids_len, update i_span2 =%d,%d\n", rm_ids_len, i_span);
                /* remove noise pixels between i_start & i */
                m = 0;
                for (k = 0, k_new=0; k < end-i_start+1; k++)
                {
                    if (m < rm_ids_len && k == rm_ids[m])
                    {
                        m++;
			continue;
                    }
                    clrx[k_new] = clrx[k];
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                        clry[i_b][k_new] = clry[i_b][k];
                    k_new++;
                }
                end = k_new;
	    }

	    if ((end - i_start + 1) >= CONSE)
	    {
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    // bdavis
                    // not sure if this is the intention, but update_num_c is MAXINT
                    // without this assignment. which causes a crash below in the
                    // loop below assigning gec_cg[i_b].coefs[k]
                    update_num_c = MIN_NUM_C;
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
#if 0
		      printf("i_b,k,fit_cft[i_b][k]=%d,%d,%f\n",i_b,k,fit_cft[i_b][k]);
#endif
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
                    rec_cg[num_fc].magnitude[i_b] = 0.0; 
	        }
                num_fc++;
                // bdavis debug
	        // printf("num_fc5=%d\n",num_fc);                                    
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

    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      FAILURE);
    }

    status = free_2d_array ((void **) fp_bin);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fp_bin\n", FUNC_NAME,
                      FAILURE);
    }

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
    strcpy(outputBinary, outDir);
    strcat(outputBinary, "/output.bin");
    //fp_bin_out = fopen("output.bin", "wb");
    fp_bin_out = fopen(outputBinary, "wb");
    if (fp_bin_out == NULL)
    {
        RETURN_ERROR ("Opening output.bin file\n", FUNC_NAME,
                      FAILURE);
    }
    if (num_fc == 0)
        fwrite(rec_cg, sizeof(Output_t), 1, fp_bin_out);
    else
        fwrite(rec_cg, sizeof(Output_t), num_fc-1, fp_bin_out);
    fclose(fp_bin_out);   

    if (verbose)
    {
        if (num_fc == 0)
	{
            printf("rec_cg[0].t_start=%d\n",rec_cg[0].t_start);
            printf("rec_cg[0].t_end=%d\n",rec_cg[0].t_end);
            printf("rec_cg[0].t_break=%d\n",rec_cg[0].t_break);
            printf("rec_cg[0].pos.row=%d\n",rec_cg[0].pos.row);
            printf("rec_cg[0].pos.col=%d\n",rec_cg[0].pos.col);
            printf("rec_cg[0].num_obs=%d\n",rec_cg[0].num_obs);
            printf("rec_cg[0].category=%d\n",rec_cg[0].category);
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
		{
                    // bdavis debug
                    printf("i_b,k,rec_cg[0].coefs[i_b][k] = %d,%d,%f\n", 
                         i_b,k,rec_cg[0].coefs[i_b][k]); 
		}
                // bdavis debug
                printf("rec_cg[0].rmse[i_b] = %f\n",rec_cg[0].rmse[i_b]);
                printf("rec_cg[0].magnitude[i_b]=%f\n",rec_cg[0].magnitude[i_b]); 
            }
	}
	else
	{
            for (i = 0; i < num_fc; i++)
            {
                printf("i=%d\n",i);
                printf("rec_cg[%d].t_start=%d\n",i,rec_cg[i].t_start);
                printf("rec_cg[%d].t_end=%d\n",i,rec_cg[i].t_end);
                printf("rec_cg[%d].t_break=%d\n",i,rec_cg[i].t_break);
                printf("rec_cg[%d].pos.row=%d\n",i,rec_cg[i].pos.row);
                printf("rec_cg[%d].pos.col=%d\n",i,rec_cg[i].pos.col);
                printf("rec_cg[%d].num_obs=%d\n",i,rec_cg[i].num_obs);
                printf("rec_cg[%d].category=%d\n",i,rec_cg[i].category);
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < update_num_c; k++)
		    {
                        // bdavis
                        // I belive the indecies being printed were incorrect.
                        // changed i_b,k,i
                        // to      i,i_b,k
                        printf("i_b,k,rec_cg[%d].coefs[i_b][k] = %d,%d,%f\n", 
                             i,i_b,k,rec_cg[i].coefs[i_b][k]); 
                        //printf("i_b,k,rec_cg[%d].coefs[i_b][k] = %d,%d,%f\n", 
                        //     i_b,k,i,rec_cg[i].coefs[i_b][k]); 
		    }
                    printf("rec_cg[%d].rmse[i_b] = %f\n",i,rec_cg[i].rmse[i_b]);
                    printf("rec_cg[%d].magnitude[i_b]=%f\n",i,rec_cg[i].magnitude[i_b]); 
                }
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
            " --inDir=<input directory>"
            " --outDir=<output directory>"
            " [--sceneList=<file with list of sceneIDs>]"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --row=: input row number\n");
    printf ("    --col=: input col number\n");
    printf ("    --inDir=: input data directory location\n");
    printf ("    --outDir=: directory location for output files\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    --sceneList=: file name containing list of sceneIDs"
            " (default is all files in inDir)\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("ccdc"
            " --row=3845"
            " --col=2918"
            " --inDir=/data/user/in"
            " --outDir=/home/user/out"
            " --sceneList=/home/user/scene_list.txt"
            " --verbose\n");
    printf ("Note: Previously, the ccdc had to be run from the directory"
            " where the input data are located.\n");
    printf ("      Now, input and output directory locations specifications are required.\n\n");
    //printf ("Note: The ccdc must run from the directory"
    //        " where the input data are located.\n\n");
}
