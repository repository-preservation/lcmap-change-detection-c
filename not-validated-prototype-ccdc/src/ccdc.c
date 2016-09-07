#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "ccdc.h"
#include "matio.h"

#define NUM_LASSO_BANDS 5
#define TOTAL_IMAGE_BANDS 7
#define TOTAL_BANDS 8
#define MIN_NUM_C 4
#define MID_NUM_C 6
#define MAX_NUM_C 8
#define CONSE 6
#define N_TIMES 3 /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25 /* average number of days per year */
#define MAX_NUM_FC 30000 /* Values change with number of pixels run */
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
#define CFMASK_FILL  255 
#define IMAGE_FILL -9999

const char scene_list_name[] = {"scene_list.txt"};
int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5}; /* This is band index */

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
--------    ---------------  ----------------------------------------------
1/15/2013   Song Guo         Original Development
20151203    Brian Davis      Added arguments for input and output file
                             locations and optional scene list file.
20160104    Song Guo         Numerous bug fixes.
20160218    Song Guo         Changed to use Zhe's stacked ARD data as inputs
 
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
    int num_fc = -1;
    int rec_fc;
    float v_start[NUM_LASSO_BANDS];
    float v_end[NUM_LASSO_BANDS];
    float v_slope[NUM_LASSO_BANDS];
    float v_dif[NUM_LASSO_BANDS];
    float **v_diff;
    int *sdate;
    int *sdate_pix;
    Input_meta_t *meta;
    int row, col;
    int ncols;
    int landsat_number;
    int clear_sum;
    int water_sum;
    int shadow_sum;
    int cloud_sum;
    int fill_sum;
    int clr_sum;
    int sn_sum;
    int all_sum;
    float sn_pct;
    float clr_pct;
    int n_sn;
    int n_clr;
    int *id_clr;
    int *id_all;
    int *id_sn;
    int *clrx;
    float **clry;
    int *cpx;
    float **cpy;
    int i_start;
    int end;
    int end2;
    float **fit_cft;
    float *rmse;
    int i_span;
    //int update_num_c;
    int update_num_c = MIN_NUM_C;
    int bl_train;
    float time_span;
    int *bl_ids;
    int *id_range;
    int *ids;
    int *ids_old;
    int *rm_ids;
    int rm_ids_len;
    int i_rec;
    float v_dif_norm;
    int i_count;
    float **v_dif_mag;
    int i_conse, i_b;
    float *vec_mag;
    float *vec_magg;
    float v_dif_mean;
    float vec_magg_min;
    float **rec_v_dif;
    float **rec_v_dif_copy;
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
    unsigned char *fmask_pix_buf;
    int **buf;
    int **pix_buf;
    FILE *fp_bin;
    int clr_land_water_counter;
    int fmask_total_counter;
    char in_path[MAX_STR_LEN];    /* directory location of input data/files     */
    char out_path[MAX_STR_LEN];   /* directory location for output files        */
    char scene_list_filename[MAX_STR_LEN];   /* file name containing list of input sceneIDs */
    char scene_list_file[MAX_STR_LEN]; /* optional input argument for file of list of scenes */
    char tmpstr[MAX_STR_LEN];   /* char string for text manipulation */
    char *short_scene;          /* char string for text manipulation */
    int len;                    /* return of strlen for manipulating strings */
    char output_mat[MAX_STR_LEN]; /* directory and file name for output.bin */
    char hdr_name[MAX_STR_LEN]; /* string for header file to read, different for LC8 */
    char directory[MAX_STR_LEN];
    char scene_name[MAX_STR_LEN];
    char errmsg[MAX_STR_LEN];
    int valid_num_scenes;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Initialize the input and output directory specification, */
    /* because they are optional.                               */
    strcpy(in_path, "");
    strcpy(out_path, "");

    /* Read the command-line arguments */
    status = get_args (argc, argv, &row, &ncols, in_path, out_path, scene_list_file, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    printf("row,ncols,verbose=%d,%d,%d\n",row,ncols,verbose);
    printf("Input Directory:  %s\n", in_path);
    printf("Output Directory: %s\n", out_path);
    printf("scene_list_file: %s\n", scene_list_file);

    /* allocate memory for scene_list */
    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                         sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR ("Allocating scene_list memory", FUNC_NAME, FAILURE);
    }

    /* check if scene_list.txt file exists, if not, create the scene_list
       from existing files in the current data working directory */
    if (access(scene_list_file, F_OK) != 0) /* File does not exist */
    {
        strcpy (scene_list_filename, in_path);
        strcat(scene_list_filename, "/scene_list.txt");
        if (access(scene_list_filename, F_OK) != -1) /* File exists */
        {
            num_scenes = MAX_SCENE_LIST;
        }
        else /* Default File exists */
        {
	    status = create_scene_list("L*", &num_scenes, scene_list_filename);
            if(status != SUCCESS)
            RETURN_ERROR("Running create_scene_list file", FUNC_NAME, FAILURE);
        }
    }
    else /* File exists */
    {
        strcpy(scene_list_filename, scene_list_file);
    }

    printf("num_scenes=%d\n",num_scenes);
    printf("scene_list_filename=%s\n",scene_list_filename);
    fd = fopen(scene_list_filename, "r");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_scenes; i++)
    {
        if (fscanf(fd, "%s", tmpstr) == EOF)
            break;
        strcpy(scene_list[i], in_path);
        strcat(scene_list[i], "/");
        strcat(scene_list[i], tmpstr);
    }
    num_scenes = i;

    /* Allocate memory */
    sdate = malloc(num_scenes * sizeof(int));
    if (sdate == NULL)
    {
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
    }

    fmask_buf = malloc(num_scenes * ncols * sizeof(unsigned char));
    if (fmask_buf == NULL)
    {
        RETURN_ERROR("ERROR allocating fmask_buf memory", FUNC_NAME, FAILURE);
    }

    buf = (int **) allocate_2d_array (TOTAL_BANDS * ncols, num_scenes,
                                         sizeof (int));
    if (buf == NULL)
    {
        RETURN_ERROR ("Allocating buf memory", FUNC_NAME, FAILURE);
    }

    if (verbose)
    {
        printf("num_scenes = %d\n", num_scenes);
        printf("scene_list[0]=%s\n", scene_list[0]);
    }

    /* sort scene_list based on year & julian_day */
    status = sort_scene_based_on_year_doy_row(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, FAILURE);
    }

    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
    {
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
    }

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    len = strlen(scene_list[0]);
    short_scene= strndup(scene_list[0], len-5);
    split_directory_scenename(scene_list[0], directory, scene_name);
    if (strncmp(short_scene, ".", 1) == 0)
    {
        strncpy(tmpstr, short_scene + 2, len - 2);
        sprintf(filename, "%s/%s_MTLstack.hdr", tmpstr, scene_name);
    }
    else
        sprintf(filename, "%s/%s_MTLstack.hdr", short_scene, scene_name);
    printf("filename=%s\n",filename);
    status = read_envi_header(filename, meta);
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

    /* Read in matlab version ARD data */
    for (i = 0; i < num_scenes; i++)
    {
        len = strlen(scene_list[i]);
        strncpy(short_scene, scene_list[i], len-5);
        split_directory_scenename(scene_list[i], directory, scene_name);
        if (strncmp(short_scene, ".", 1) == 0)
        {
            strncpy(tmpstr, short_scene + 2, len - 2);
            sprintf(filename, "%s/%s_MTLstack", tmpstr, scene_name);
        }
        else
            sprintf(filename, "%s/%s_MTLstack", short_scene, scene_name);
        fp_bin = open_raw_binary(filename,"rb");
        if (fp_bin == NULL)
	{
	    sprintf(errmsg, "Opening %d scene name: %s\n", i, filename);
            RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
	}
        fseek(fp_bin, ((row - 1)* meta->samples) * 
              TOTAL_BANDS * sizeof(short int), SEEK_SET);
	for (col = 0; col < ncols; col++)
	{
	    for (b = 0; b < TOTAL_BANDS; b++)
	    {
                if (read_raw_binary(fp_bin, 1, 1,
                    sizeof(short int), &buf[col * TOTAL_BANDS + b][i]) != 0)
	        {
                    sprintf(errmsg, "error reading %d scene\n", i);
                    RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
	        }
	    }
	    fmask_buf[col * num_scenes + i] = 
                (unsigned char)buf[col * TOTAL_BANDS + 7][i]; 
	}
	close_raw_binary(fp_bin);
    }
    free(short_scene);
	
    /* Allocate memory for rec_cg */ 
    rec_cg = malloc(MAX_NUM_FC * sizeof(Output_t));
    if (rec_cg == NULL)
    {
        RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);
    }

    sdate_pix = malloc(num_scenes * sizeof(int));
    if (sdate_pix == NULL)
    {
        RETURN_ERROR("ERROR allocating sdate_pix memory", FUNC_NAME, FAILURE);
    }

     fmask_pix_buf = malloc(num_scenes * sizeof(unsigned char));
    if (fmask_pix_buf == NULL)
    {
        RETURN_ERROR("ERROR allocating fmask_pix_buf memory", FUNC_NAME, FAILURE);
    }

    pix_buf = (int **) allocate_2d_array (TOTAL_IMAGE_BANDS, MAX_SCENE_LIST, sizeof (int));
    if (pix_buf == NULL)
    {
        RETURN_ERROR ("Allocating pix_buf memory", FUNC_NAME, FAILURE);
    }

    /* Loop through pixels in a line */
    for (col = 0; col < ncols; col++)
    {
        printf("col number = %d\n", col + 1);

        //sdate_pix = malloc(num_scenes * sizeof(int));
        //if (sdate_pix == NULL)
        //{
        //    RETURN_ERROR("ERROR allocating sdate_pix memory", FUNC_NAME, FAILURE);
        //}

        //fmask_pix_buf = malloc(num_scenes * sizeof(unsigned char));
        //if (fmask_pix_buf == NULL)
        //{
        //    RETURN_ERROR("ERROR allocating fmask_pix_buf memory", FUNC_NAME, FAILURE);
        //}

        //pix_buf = (int **) allocate_2d_array (TOTAL_BANDS, num_scenes,
        //                                 sizeof (int));
        //if (pix_buf == NULL)
        //{
        //    RETURN_ERROR ("Allocating pix_buf memory", FUNC_NAME, FAILURE);
        //}

        /* Removed fill pixels */
        clear_sum = 0;
        water_sum = 0;
        shadow_sum = 0;
        cloud_sum = 0;
        fill_sum = 0;
	clr_sum = 0;
	sn_sum = 0;
	all_sum = 0;
        valid_num_scenes = 0;
        for (i = 0; i < num_scenes; i++)
        {
            /* a bitmask should probably be set up to do these....... */
            switch (fmask_buf[col * num_scenes + i])
            {
                case CFMASK_CLEAR:
                    clear_sum++;
                    clr_sum++;
                    break;
                case CFMASK_WATER:
                    water_sum++;
                    clr_sum++;
                    break;
                case CFMASK_SHADOW:
                    shadow_sum++;
                    break;
                case CFMASK_SNOW:
                    sn_sum++;
                    break;
                case CFMASK_CLOUD:
                    cloud_sum++;
                    break;
                case CFMASK_FILL:
                    fill_sum++;
                    break;
                default:
                    printf ("Unknown fmask value %d\n", fmask_buf[col * num_scenes + i]);
                    break;
            }

            if (fmask_buf[col * num_scenes + i] < CFMASK_FILL)
            {
                all_sum++;
                fmask_pix_buf[valid_num_scenes] = fmask_buf[col * num_scenes + i];
                sdate_pix[valid_num_scenes] = sdate[i];
	        for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                    pix_buf[b][valid_num_scenes] = buf[col * TOTAL_BANDS + b][i];
                valid_num_scenes++;
	    }
        }

        /* Allocate memory */
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

        id_clr = (int *)calloc(valid_num_scenes, sizeof(int));
        if (id_clr == NULL)
        {
            RETURN_ERROR("ERROR allocating id_clr memory", FUNC_NAME, FAILURE);
        }

        id_all = (int *)calloc(valid_num_scenes, sizeof(int));
        if (id_all == NULL)
        {
            RETURN_ERROR("ERROR allocating id_all memory", FUNC_NAME, FAILURE);
        }

        id_sn = (int *)calloc(valid_num_scenes, sizeof(int));
        if (id_sn == NULL)
        {
            RETURN_ERROR("ERROR allocating id_sn memory", FUNC_NAME, FAILURE);
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

        //clry = (float **) allocate_2d_array ((TOTAL_IMAGE_BANDS + 1), valid_num_scenes, sizeof (float));
        clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes, sizeof (float));
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

        //vec_mag = (float *)calloc(NUM_LASSO_BANDS, sizeof(float));
        vec_mag = (float *)calloc(CONSE, sizeof(float));
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

        /* Allocate memory for temp_v_dif */
        temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                     sizeof (float));
        if (temp_v_dif == NULL)
        {
            RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
        }

        /* percent of clear pixels, clear (0) or water (1) */    
        clr_pct = (float) clr_sum / (float) all_sum;

        /* percent of snow observations, (3) */
        if ((clr_sum + sn_sum) != 0)
            sn_pct =  (float) sn_sum / (float)(clr_sum + sn_sum); 
        else
            sn_pct = (float)sn_sum;

    if (verbose)
    {
        printf(" clear_sum =   %d\n", clear_sum);
        printf(" water_sum =   %d\n", water_sum);
        printf(" shadow_sum =  %d\n", shadow_sum);
        printf(" snow_sum =    %d\n", sn_sum);
        printf(" cloud_sum =   %d\n", cloud_sum);
        printf(" fill_sum =    %d\n", fill_sum);
        printf(" clear+water = %d\n", clr_sum);
        printf(" clear_pct   = %f\n", clr_pct);
        printf(" snow_pct    = %f\n", sn_pct);
    }

        /* pixel value ranges should follow physical rules */
        for (i = 0; i < valid_num_scenes; i++)
        {  
            /* convert Kelvin to Celsius (for new espa data) */
	    if (pix_buf[6][i] != -9999)
                pix_buf[6][i] = (int16)(pix_buf[6][i] * 10 - 27315);

            if ((pix_buf[0][i] > 0) && (pix_buf[0][i] < 10000) &&
                (pix_buf[1][i] > 0) && (pix_buf[1][i] < 10000) &&
                (pix_buf[2][i] > 0) && (pix_buf[2][i] < 10000) &&
                (pix_buf[3][i] > 0) && (pix_buf[3][i] < 10000) &&
                (pix_buf[4][i] > 0) && (pix_buf[4][i] < 10000) &&
                (pix_buf[5][i] > 0) && (pix_buf[5][i] < 10000) &&
                (pix_buf[6][i] > -9320) && (pix_buf[6][i] < 7070))
	    {
                id_range[i] = 1;
	    }
            else
	    {
                id_range[i] = 0;
	    }
        }

	/* Initilialize rec_cg */
	rec_cg[0].t_break = 0;

        /* fit permanent snow observations */
        if (clr_pct < T_CLR)
        {
            if (sn_pct > T_SN)
            {
                n_sn = 0;
                /* snow observations are "good" now */
                for (i = 0; i < valid_num_scenes; i++)
                { 
	            //if (((fmask_pix_buf[i] == CFMASK_SNOW) || (fmask_pix_buf[i] < 2)) 
                    //     && id_range[i] == 1)
	            if ((fmask_pix_buf[i] == CFMASK_SNOW) || (fmask_pix_buf[i] < 2))
                       
                    {
                        clrx[n_sn] = sdate_pix[i];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		        {
         	            clry[k][n_sn] = (float)pix_buf[k][i];
		        }
                        n_sn++;
                    }
                }  
                end = n_sn;

                if (n_sn < N_TIMES * MIN_NUM_C) /* not enough snow pixels */
                {
        	    continue;
                }

        	/* Remove repeated ids */
                matlab_unique(clrx, clry, n_sn, &end);

                /* start model fit for snow persistent pixels */
                if (verbose)
                    printf ("Fit permanent snow observations, now pixel = %f\n", 100.0 * sn_pct); 

                /* the first observation for TSFit */
                i_start = 1; /* the first observation for TSFit */

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

                /* treat saturated and unsaturated pixels differently */
                //for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                {
         	    i_span = 0;
                    //if (k != TOTAL_IMAGE_BANDS - 1) /* for optical bands */
                    if (b != TOTAL_IMAGE_BANDS - 1) /* for optical bands */
                    {
                        //for (i = 0; i < end; i++)
                        for (k = 0; k < end; k++)
                        {
                            //if (clry[k][i] > 0.0 && clry[k][i] < 10000.0)
                            if (clry[b][k] > 0.0 && clry[b][k] < 10000.0)
                            {
                                //clrx[i_span] = clrx[i];
                                //clry[k][i_span] = clry[k][i];
                                //cpx[i_span] = clrx[i];
                                //cpy[k][i_span] = clry[k][i];
                                cpx[i_span] = clrx[k];
                                cpy[b][i_span] = clry[b][k];
                                i_span++;
                            }
                        }

                        if (i_span < MIN_NUM_C * N_TIMES)
                        {
                            for (k = 0; k < MIN_NUM_C; k++)
                                //fit_cft[k][i] = 10000.0; /* fixed value for saturated pixels */
                                fit_cft[b][k] = 10000.0; /* fixed value for saturated pixels */
                        }
                        else
                        {
                            //status = auto_ts_fit(cpx, cpy, k, 0, i_span-1, MIN_NUM_C, 
                            //         fit_cft, &rmse[k], temp_v_dif); 
                            status = auto_ts_fit(cpx, cpy, b, 0, i_span-1, MIN_NUM_C, 
                                     fit_cft, &rmse[b], temp_v_dif); 
                            if (status != SUCCESS)  
                                RETURN_ERROR ("Calling auto_ts_fit1\n", FUNC_NAME, EXIT_FAILURE);
                        }
                    }
                    else /* for thermal band */
                    {
                        //for (i = 0; i < end; i++)
                        for (k = 0; k < end; k++)
                        {
                            //if (clry[k][i] > -9300.0 && clry[k][i] < 7070.0)
                            if (clry[b][k] > -9300.0 && clry[b][k] < 7070.0)
                            {
                                //clrx[i_span] = clrx[i];
                                //clry[k][i_span] = clry[k][i];
                                //cpx[i_span] = clrx[i];
                                //cpy[k][i_span] = clry[k][i];
                                cpx[i_span] = clrx[k];
                                cpy[b][i_span] = clry[b][k];
                                i_span++;
                            }
                        }
                            
                        //status = auto_ts_fit(clrx, clry, k, 0, i_span-1, MIN_NUM_C, fit_cft, 
                        //         &rmse[k], temp_v_dif); 
                        //status = auto_ts_fit(cpx, cpy, k, 0, i_span-1, MIN_NUM_C, fit_cft, 
                        //         &rmse[k], temp_v_dif); 
                        status = auto_ts_fit(cpx, cpy, b, 0, i_span-1, MIN_NUM_C, fit_cft, 
                                 &rmse[b], temp_v_dif); 
                        if (status != SUCCESS)  
                            RETURN_ERROR ("Calling auto_ts_fit2\n", FUNC_NAME, EXIT_FAILURE);
                    }
                }

                /* free memories */
                free(cpx);
                status = free_2d_array ((void **) cpy);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: cpy\n", FUNC_NAME, FAILURE);
                }

                /* NUM of Fitted Curves (num_fc) */
                num_fc++;
		printf("    num_fc1 = %d\n",num_fc);   
           	if (num_fc >= MAX_NUM_FC)
	        {
                    /* Reallocate memory for rec_cg */ 
                    rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
	            {
                        RETURN_ERROR("ERROR allocating rec_cg memory", 
                                     FUNC_NAME, FAILURE);
		    }
                }
                /* update information at each iteration */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[i_start-1]; 
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end-1]; 
                /* no break at the moment */
                rec_cg[num_fc].t_break = 0; 
		/* record position of pixel */
		rec_cg[num_fc].pos = (row-1) * ncols + col + 1;
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
		    {
                        /* record fitted coefficients */
                        rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k];
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
            }
            else
            {
        	n_clr = 0;
                /* no change detection for clear observations, backup algorithm */
	        for (i = 0; i < valid_num_scenes; i++)
                { 
                    if (id_range[i] == 1)
                    {
		        clrx[n_clr] = sdate_pix[i];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		            clry[k][n_clr] = (float)pix_buf[k][i];
                        n_clr++;
                    }   
                }
                end = n_clr;

        	/* Remove repeated ids */
                matlab_unique(clrx, clry, n_clr, &end);

                /* start model fit for clear persistent pixels */
                if (verbose)
                    printf ("Fmask failed, clear pixel = %f\n", 100.0 * clr_pct); 

                // v6
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

                for (i = 0; i < end; i++)
                {
                    cpx[i] = clrx[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                    {
                        cpy[k][i] = clry[k][i];
                    }
                }
                // v6

        	n_clr = 0;
	        float band2_median;
                quick_sort_float(clry[1], 0, end - 1);
                matlab_2d_float_median(clry, 1, end, &band2_median);

                for (i = 0; i < end; i++)
                {
        	    //if (clry[1][i] < (band2_median + 400.0))
        	    // v6 
        	    if (cpy[1][i] < (band2_median + 400.0))
	            {
        	        //clrx[n_clr] = clrx[i];
        	        // v6 
        	        clrx[n_clr] = cpx[i];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                        { 
         	            //clry[k][n_clr] = clry[k][i]; 
         	            // v6
         	            clry[k][n_clr] = cpy[k][i]; 
		        }
    		        n_clr++;
                    }
	        }
                end = n_clr;

        	/* the first observation for TSFit */
                i_start = 1; /* the first observation for TSFit */

                if (n_clr < N_TIMES * MIN_NUM_C)
                {
		    continue;
	        }
                else
                {
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        //status = auto_ts_fit(clrx, clry, i_b, 0, end-1, MIN_NUM_C, 
                        //                     fit_cft, &rmse[k], temp_v_dif); 
                        // v6
                        status = auto_ts_fit(clrx, clry, i_b, 0, end-1, MIN_NUM_C, 
                                             fit_cft, &rmse[i_b], temp_v_dif); 
                        if (status != SUCCESS)
		        {  
                            RETURN_ERROR ("Calling auto_ts_fit for clear persistent pixels\n", 
                                          FUNC_NAME, FAILURE);
		        }
		    }
                }

                // v6
                /* Free the memories */
                free(cpx);
                status = free_2d_array ((void **) cpy);
                if (status != SUCCESS)
                {
                    RETURN_ERROR ("Freeing memory: cpy\n", FUNC_NAME, FAILURE);
                }
                // v6

                /* NUM of Fitted Curves (num_fc) */
                num_fc++;   
		printf("    num_fc2 = %d\n",num_fc);   
           	if (num_fc >= MAX_NUM_FC)
	        {
                    /* Reallocate memory for rec_cg */ 
                    rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                    if (rec_cg == NULL)
	            {
                        RETURN_ERROR("ERROR allocating rec_cg memory", 
                                     FUNC_NAME, FAILURE);
		    }
                }
                /* update information at each iteration */
                /* record time of curve start */
                rec_cg[num_fc].t_start = clrx[i_start-1]; 
                /* record time of curve end */
                rec_cg[num_fc].t_end = clrx[end-1]; 
                /* no break at the moment */
                rec_cg[num_fc].t_break = 0; 
		rec_cg[num_fc].pos = (row-1) * ncols + col + 1;
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
		    {
                        /* record fitted coefficients */
                        rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k];
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
	    for (i = 0; i < valid_num_scenes; i++)
            { 
	        if ((fmask_pix_buf[i] < 2) && (id_range[i] == 1))
                {
	            clrx[n_clr] = sdate_pix[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		    {
		        clry[k][n_clr] = (float)pix_buf[k][i];
		    }
                    n_clr++;
                }   
            }
            end = n_clr;

            /* Remove repeated ids */
            matlab_unique(clrx, clry, n_clr, &end);

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

            /* NUM of Fitted Curves (num_fc) */
            num_fc++;

	    printf("    num_fc3 = %d\n",num_fc);   

	    if (num_fc >= MAX_NUM_FC)
	    {
                /* Reallocate memory for rec_cg */ 
		rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
                if (rec_cg == NULL)
		{
                    RETURN_ERROR("ERROR allocating rec_cg memory", 
                                 FUNC_NAME, FAILURE);
         	}
            }   

            /* record the num_fc at the beginning of each pixel */
            rec_fc = num_fc;

            /* record the start of Tmask (0=>initial;1=>done) */
            bl_tmask = 0;

            /* while loop - process till the last clear observation - CONSE */
            snprintf (msg_str, sizeof(msg_str), "CCDC init_time=%s\n", ctime (&now));
                LOG_MESSAGE (msg_str, FUNC_NAME);

	    while (i <= (end - CONSE))
            {
                /* span of "i" */
                i_span = i - i_start + 1;

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
                        /* Clears the IDs buffers */
		        for (k = 0; k < valid_num_scenes; k++)
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
                                m++;
                            }
                            else
                                i_span++;  /* update i_span after noise removal */
                        }

                        rm_ids_len = m;
                        /* update i_span after noise removal */
                        // i_span = i - i_start +1 - rm_ids_len;

                        /* check if there is enough observation */
                        if (i_span < (N_TIMES * MIN_NUM_C))
                        {
                            /* move forward to the i+1th clear observation */
                            i++;
                            /* not enough clear observations */
                            continue;
                        }
			else
		        {
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
                        end2 = k_new;
			end = end2;

                        /* record i before noise removal 
                           This is very important as if model is not initialized 
                           the multitemporal masking shall be done again instead 
                           of removing outliers in every masking */
                        i_rec=i;

                        /* update i afer noise removal (i_start stays the same) */
                        i=i_start + i_span - 1;

                        /* update span of time (num of years) */
                        time_span=(cpx[i-1] - cpx[i_start-1]) / NUM_YEARS;

                        /* check if there is enough time */
                        if (time_span < MIN_YEARS)
                        {
                            i = i_rec;   /* keep the original i */
                            /* move forward to the i+1th clear observation */
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
                        /* remove noise in original arrays */
                        for (k = 0; k < end2; k++)
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
                            /* Initial model fit */
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
                            v_dif_norm += v_dif[b] * v_dif[b];               
                        }

                        /* find stable start for each curve */
                        if (v_dif_norm > T_CG)
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
                                            }
                                        }
                                    }
                                    vec_magg[i_conse] = v_dif_norm; 

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
                                            clry[b][k] = clry[b][k+1];
				        }
                                    }
                                    i--;
				    end--;
				}

                                /* update i_start if i_ini is not a confirmed break */
                                i_start = i_ini + 1;

                                /* free the memory */
                                free(vec_magg);
                                status = free_2d_array ((void **) v_diff);
                                if (status != SUCCESS)
			        {
                                    RETURN_ERROR ("Freeing memory: v_diff\n", 
                                                  FUNC_NAME, FAILURE);
			        }
                            }
                        }

                        /* enough to fit simple model and confirm a break */
                        if ((num_fc == rec_fc) && ((i_start - i_break) >= CONSE))
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
                	    /* record position of pixel */
		            rec_cg[num_fc].pos = (row-1) * ncols + col + 1;
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                for (k = 0; k < MIN_NUM_C; k++)
			        {
                                    /* record fitted coefficients */
                                    rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k]; 
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
         		        quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                                matlab_2d_float_median(v_dif_mag, i_b, ini_conse, 
                                                      &v_dif_mean);
                                rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                            }
                            /* identified and move on for the next functional curve */
                            num_fc++; 

		            printf("    num_fc4 = %d\n",num_fc);   

                            if (num_fc >= MAX_NUM_FC)
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
			}
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
                        for (k = 0; k < valid_num_scenes; k++)
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

                        /* determine the time series model */
                        update_cft(i_span, N_TIMES, MIN_NUM_C, MID_NUM_C, MAX_NUM_C, 
                                  num_c, &update_num_c);

                        /* initial model fit when there are not many obs */
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
                           }

                           /* updating information for the first iteration */
                           /* record time of curve start */
                           rec_cg[num_fc].t_start = clrx[i_start-1]; 
                           /* record time of curve end */
                           rec_cg[num_fc].t_end = clrx[i-1]; 
                           /* no break at the moment */
                           rec_cg[num_fc].t_break = 0; 
	 	           /* record position of pixel */
		           rec_cg[num_fc].pos = (row-1) * ncols + col + 1;

                           for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                           {
                               for (k = 0; k < MAX_NUM_C; k++)
			       {
                                   /* record fitted coefficients */
                                   rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k]; 
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
                                       }
                                   }
			       }
                               vec_mag[i_conse] = v_dif_norm;
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
                        }
                        else
                        {
            		    if ((float)(clrx[i-1] - clrx[i_start-1]) >= (1.33*(float)i_count))
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

                                }

                                /* record fitted coefficients */
                                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                {
                                    for (k = 0; k < MAX_NUM_C; k++)
			            {
                                        /* record fitted coefficients */
                                        rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k];
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
                            }

                            /* record time of curve end */
                            rec_cg[num_fc].t_end = clrx[i-1]; /* record time of curve end */

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
                            }

			    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
			    {
         		        for (m = 0; m < ids_old_len; m++)
         			    rec_v_dif_copy[b][m] = rec_v_dif[b][m];
			    }

                            /* sort the rec_v_dif based on d_yr */
        		    quick_sort_2d_float(d_yr, rec_v_dif_copy, 0, ids_old_len-1); 

                            for(b = 0; b < NUM_LASSO_BANDS; b++)
                                tmpcg_rmse[b] = 0.0;

                            /* temporarily changing RMSE */
                            for (b = 0; b < NUM_LASSO_BANDS; b++)
                            {
                                matlab_2d_array_norm(rec_v_dif_copy, lasso_blist[b], n_rmse,
                                             &tmpcg_rmse[b]);
                                tmpcg_rmse[b] /= sqrt(n_rmse - rec_cg[num_fc].category);
                            }

                            /* free allocated memories */
                            free(d_yr);

                            /* move the ith col to i-1th col */
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
                                /* absolute difference for all bands */
        		        auto_ts_predict(clrx, fit_cft, update_num_c, i_b, i+CONSE-1, 
                                     i+CONSE-1, &ts_pred_temp);
                                v_dif_mag[i_b][CONSE-1] = clry[i_b][i+CONSE-1] - ts_pred_temp;

                                /* normalized to z-scores */
                                for (b = 0; b < NUM_LASSO_BANDS; b++)
                                {
                                    if (i_b == lasso_blist[b])
                                    {
                                        /* mini rmse */
                                        mini_rmse = max(adj_rmse[i_b], tmpcg_rmse[b]);

                                        /* z-score */
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
                            /* record break time */
                            rec_cg[num_fc].t_break = clrx[i];
                            rec_cg[num_fc].change_prob = 1.0;
                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
	        	    {
         		        quick_sort_float(v_dif_mag[i_b], 0, CONSE-1);

                                matlab_2d_float_median(v_dif_mag, i_b, CONSE, 
                                    &rec_cg[num_fc].magnitude[i_b]);
			    }
                            /* identified and move on for the next functional curve */
                            num_fc++;

		            printf("    num_fc5 = %d\n",num_fc);   

                            if (num_fc >= MAX_NUM_FC)
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
                            /* start training again */
                            bl_train = 0;
                        }
                        else if (vec_mag[0] > T_MAX_CG)
                        {
                            /* remove noise */
                            for (m = i; m < end -1; m++)
                            {
                                clrx[m] = clrx[m+1];
                                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                                    clry[b][m] = clry[b][m+1];
                            }
                            end--; /* check if this is needed */

                            /* stay & check again after noise removal */
                            i--;
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
            } /* end of "while (i <= end - CONSE) */

            /* Two ways for processing the end of the time series */ 
            if (bl_train == 1)
            {
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

                /* update change probability */
                //rec_cg[num_fc].change_prob = (CONSE - id_last) / CONSE; 
                // v6
                rec_cg[num_fc].change_prob = (float)(CONSE - id_last) / (float)CONSE; 
                /* update end time of the curve */
                //rec_cg[num_fc].t_end = clrx[end - CONSE + id_last];
                // v6 
                rec_cg[num_fc].t_end = clrx[end - CONSE + id_last - 1];

                /* mean value fit for the rest of the pixels < CONSE & > 1 */
                if (CONSE > id_last)
                {
                    /* update time of the probable change */
                    rec_cg[num_fc].t_break = clrx[end-CONSE+id_last+1];
                    /* update magnitude of change */
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
		    {
         	        //quick_sort_float(v_dif_mag[i_b], id_last, CONSE-2);
                        //matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, CONSE-1, 
                        //                 &rec_cg[num_fc].magnitude[i_b]);
                        // v6 
                        if (id_last == CONSE - 1)
                        {
                            rec_cg[num_fc].magnitude[i_b] = v_dif_mag[i_b][id_last];
                        }
                        else
                        {
                            quick_sort_float(v_dif_mag[i_b], id_last, CONSE-1);
                            matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, CONSE-1,
                                         &rec_cg[num_fc].magnitude[i_b]);
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
                }
                else
                {
    	            for (k = 0; k < valid_num_scenes; k++) 
                    {
                        if (clrx[k] >= rec_cg[num_fc-1].t_break)
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
                    /* multitemporal cloud mask */
                    status = auto_mask(clrx, clry, i_start-1, end-1,
                               (float)(clrx[end-1]-clrx[i_start-1]) / NUM_YEARS, 
                               adj_rmse[1], adj_rmse[4], T_CONST, bl_ids);
                    if (status != SUCCESS)
                        RETURN_ERROR("ERROR calling auto_mask at the end of time series", 
                                      FUNC_NAME, FAILURE);

                    /* Clears the IDs buffers */
        	    for (m = 0; m < valid_num_scenes-1; m++)
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
                            m++;
                        }
                        else
                            i_span++;  /* update i_span after noise removal */
                    }
                    rm_ids_len = m;

                    /* remove noise pixels between i_start & i */
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

        	    /* record time of curve start */
                    rec_cg[num_fc].t_start = clrx[i_start-1];
                    /* record time of curve end */
                    rec_cg[num_fc].t_end = clrx[end-1];
                    /* record break time */
                    rec_cg[num_fc].t_break = 0;
		    /* record position of pixel */
		    rec_cg[num_fc].pos = (row-1) * ncols + col + 1;

                    /* record fitted coefficients */
                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        for (k = 0; k < MAX_NUM_C; k++)
		        {
                            rec_cg[num_fc].coefs[k][i_b] = fit_cft[i_b][k];
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
	        }
            }
        }
    
        /* Free memory allocation */
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
        free(sdate_pix);
        free(fmask_pix_buf);
        status = free_2d_array ((void **) pix_buf);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: pix_buf\n", FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) rec_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif\n", FUNC_NAME, FAILURE);
        }
        status = free_2d_array ((void **) rec_v_dif_copy);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: rec_v_dif_copy\n", FUNC_NAME, FAILURE);
        }
        status = free_2d_array ((void **) temp_v_dif);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: temp_v_dif\n", FUNC_NAME, FAILURE);
        }
        status = free_2d_array ((void **) v_dif_mag);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: v_dif_mag\n", FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) clry);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME, FAILURE);
        }

        status = free_2d_array ((void **) fit_cft);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME, FAILURE);
        }
    }

    mat_t *matfp;
    matvar_t* mat_struct;
    char *structname = "rec_cg";
    struct stat sb;
    //size_t structdim[10] = {1, num_fc}; // create 1x10 struct
    // v7 
    size_t structdim[10] = {1, num_fc + 1};
    // v7 
    const char *fieldnames[10] = { "t_start","t_end","t_break","coefs","rmse","pos",
				  "change_prob","num_obs","category","magnitude" };
    size_t dim1d[2]= {1, 1};
    size_t dim4d[2]= {MAX_NUM_C, TOTAL_IMAGE_BANDS};
    size_t dim5d[1]= {TOTAL_IMAGE_BANDS};
    matvar_t *variable1, *variable2, *variable3, *variable4, *variable5,
             *variable6, *variable7, *variable8, *variable9, *variable10;

#if 0
    /* Open the output file */
    strcpy(output_binary, out_path);
    strcat(output_binary, "/output.bin");
    printf("output_binary=%s\n",output_binary);
    if (access(output_binary, F_OK) != 0) /* File does not exist */
        fp_bin_out = fopen(output_binary, "wb");
    else
        fp_bin_out = fopen(output_binary, "ab");
    if (fp_bin_out == NULL)
    {
        RETURN_ERROR ("Opening output.bin file\n", FUNC_NAME,
                      FAILURE);
    }
#endif 
    strcat(out_path, "/TSFitMap"); 
    sprintf(output_mat, "%s/record_change%d.mat", out_path, row);

    /* Check the existence of the TSFitMap directory, create one if not existed */
    if (!(stat(out_path, &sb) == 0 && S_ISDIR(sb.st_mode)))
    {
        status = mkdir(out_path, 0777);
        if (status != SUCCESS)
            RETURN_ERROR ("Creating TSFitMap directory", FUNC_NAME, FAILURE);
    }
    
    /* Check the existence of the output file */
#if 0
    if (access(output_mat, F_OK) != 0) /* File does not exist */
        matfp = Mat_CreateVer(output_mat, NULL, MAT_FT_MAT5);
        //        matfp = Mat_CreateVer(output_mat, NULL, MAT_FT_MAT73);
    else
        matfp = Mat_Open(output_mat,MAT_ACC_RDWR);
    if (matfp == NULL) 
        RETURN_ERROR ("Opening MAT file", FUNC_NAME, FAILURE);
#endif
    matfp = Mat_CreateVer(output_mat, NULL, MAT_FT_MAT5);
    if (matfp == NULL) 
        RETURN_ERROR ("Opening MAT file", FUNC_NAME, FAILURE);

    mat_struct = Mat_VarCreateStruct(structname, 2, structdim, fieldnames, 10);
    if (mat_struct == NULL)
    { 
        Mat_Close(matfp);
        RETURN_ERROR ("Creating MAT structure\n", FUNC_NAME, FAILURE);
    }

    /* Output rec_cg structure to the output file 
       note: can use fread to read out the structure from the output file */
    for (i = 0; i <= num_fc; i++)
    {
        variable1 = Mat_VarCreate(fieldnames[0], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].t_start, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[0], i, variable1); 

        variable2 = Mat_VarCreate(fieldnames[1], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].t_end, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[1], i, variable2); 

        variable3 = Mat_VarCreate(fieldnames[2], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].t_break, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[2], i, variable3); 

        variable4 = Mat_VarCreate(fieldnames[3], MAT_C_SINGLE, MAT_T_SINGLE, 2, dim4d, 
                    &rec_cg[i].coefs, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[3], i, variable4); 

        variable5 = Mat_VarCreate(fieldnames[4], MAT_C_SINGLE, MAT_T_SINGLE, 1, dim5d, 
                    &rec_cg[i].rmse, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[4], i, variable5); 

        variable6 = Mat_VarCreate(fieldnames[5], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].pos, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[5], i, variable6); 

        variable7 = Mat_VarCreate(fieldnames[6], MAT_C_SINGLE, MAT_T_SINGLE, 2, dim1d, 
                    &rec_cg[i].change_prob, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[6], i, variable7); 

        variable8 = Mat_VarCreate(fieldnames[7], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].num_obs, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[7], i, variable8); 

        variable9 = Mat_VarCreate(fieldnames[8], MAT_C_INT32, MAT_T_INT32, 2, dim1d, 
                    &rec_cg[i].category, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[8], i, variable9); 

        variable10 = Mat_VarCreate(fieldnames[9], MAT_C_SINGLE, MAT_T_SINGLE, 1, dim5d, 
                    &rec_cg[i].magnitude, 0);
        Mat_VarSetStructFieldByName(mat_struct, fieldnames[9], i, variable10); 
    }
    Mat_VarWrite(matfp, mat_struct, 0);
    Mat_VarFree(mat_struct);

    if (verbose)
    {
        for (i = 0; i <= num_fc; i++)
        {
            printf("num_fc=%d\n",i);
            printf("rec_cg[%d].t_start=%d\n",i,rec_cg[i].t_start);
            printf("rec_cg[%d].t_end=%d\n",i,rec_cg[i].t_end);
            printf("rec_cg[%d].t_break=%d\n",i,rec_cg[i].t_break);
            printf("rec_cg[%d].pos=%d\n",i,rec_cg[i].pos);
            printf("rec_cg[%d].num_obs=%d\n",i,rec_cg[i].num_obs);
            printf("rec_cg[%d].category=%d\n",i,rec_cg[i].category);
            printf("rec_cg[%d].change_prob=%f\n",i,rec_cg[i].change_prob);
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                printf("rec_cg[%d].coefs[%d][0..7] = %f,%f,%f,%f,%f,%f,%f,%f\n", 
		       i,i_b,k,rec_cg[i].coefs[i_b][0],i,i_b,k,rec_cg[i].coefs[1][i_b],
		       rec_cg[i].coefs[2][i_b],i,i_b,k,rec_cg[i].coefs[3][i_b],
		       rec_cg[i].coefs[4][i_b],i,i_b,k,rec_cg[i].coefs[5][i_b],
		       rec_cg[i].coefs[6][i_b],i,i_b,k,rec_cg[i].coefs[7][i_b]); 
                printf("rec_cg[%d].rmse[%d] = %f\n",i,i_b,rec_cg[i].rmse[i_b]);
                printf("rec_cg[%d].magnitude[%d]=%f\n",i,i_b,rec_cg[i].magnitude[i_b]); 
            }
        }
    }

    /* Free rec_cg memory*/
    free(rec_cg);

    /* close the output file */
    //    fclose(fp_bin_out);  
    free(variable1);   
    free(variable2);   
    free(variable3);   
    free(variable4);   
    free(variable5);   
    free(variable6);   
    free(variable7);   
    free(variable8);   
    free(variable9);   
    free(variable10);   
    Mat_Close(matfp);
 
    /* free memory outside column loop */
    free(meta);
    free(sdate);
    free(fmask_buf);

    status = free_2d_array ((void **) buf);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: buf\n", 
                      FUNC_NAME, FAILURE);
    }

    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      FAILURE);
    }

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
            " --ncols=<input number of columns>"
            " --in_path=<input directory>"
            " --out_path=<output directory>"
            " [--scene_list_file=<file with list of sceneIDs>]"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --row=: input row number\n");
    printf ("    --ncols=: input number of columns\n");
    printf ("    --in_path=: input data directory location\n");
    printf ("    --out_path=: directory location for output files\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("ccdc"
            " --row=23"
            " --ncols=5000"
            " --in_path=/data/user/in"
            " --out_path=/home/user/out"
            " --scene_list_file=/home/user/scene_list.txt"
            " --verbose\n");
    printf ("Note: Default is ccdc must run from the directory"
            " where the input data are located.\n\n");
}
