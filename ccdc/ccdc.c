#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

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
#define N_TIMES 3         /* number of clear observations/coefficients*/
#define NUM_YEARS 365.25  /* average number of days per year          */
#define NUM_FC 10         /* Values change with number of pixels run  */
#define T_CONST 4.89      /* Threshold for cloud, shadow, and snow detection */
#define MIN_YEARS 1       /* minimum year for model intialization     */
#define T_SN 0.75         /* no change detection for permanent snow pixels */ 
#define T_CLR 0.25        /* Fmask fails threshold                    */
#define T_CG 15.0863      /* chi-square inversed T_cg (0.99) for noise removal */
#define T_MAX_CG 35.8882  /* chi-square inversed T_max_cg (1e-6) for 
                             last step noise removal                  */
#define CFMASK_CLEAR   0
#define CFMASK_WATER   1
#define CFMASK_SHADOW  2
#define CFMASK_SNOW    3
#define CFMASK_CLOUD   4
#define CFMASK_FILL  255 
#define IMAGE_FILL -9999
#define CFMASK_BAND    7

const char scene_list_name[] = {"scene_list.txt"};
int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5}; /* This is band index */

char *sub_string
(
    const char *source,
    size_t start,
    size_t length
) 
{
    size_t i;
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

RETURN VALUE: Type = int

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
20160104    Song Guo         Numerous bug fixes.
20160126    Brian Davis      Fixed all problems with calculating and
                             checking of percentages clear, etc.
                             Also printing out that information.
                             Removed most debug prints and ifdef 0s.
20160210    Brian Davis      Added stdin input option.
                             Added stdout output option.
                             Added swath overlap filtering after scene
                             list is sorted. This needs to change to a
                             per-pixel basis instead of per-scene.
                             Resolved alloc and free associations, and
                             moved some around.
                             Added debug flag for printfs.
                             Added some comments and general cleanup.
20160216    Brian Davis      Moved the swath overlap filtering down
                             to the reading of surface reflaceance
                             files into buf, so that it becomes pixel
                             based overlap checking.  The amount of
                             overlap in two scenes in a swath is very
                             small, we do not want to throw away the
                             entire "scene".  However, this requires
                             checking the fmask values for fill.  Only
                             if both upper scene and lower scene cfmask
                             values are non-fill is the x/y location in
                             an area of "pixel overlap", otherwise it is
                             "scene overlap" if one of them is fill.

When reading fmask values to detmermine fill vs. usability,
those "scenes" have been sorted, so one could test for:
If  current fmask and prev fmask are not FILL AND
    current path   =  prev path               AND
    current row    =  prev row -1             AND
    current year   =  prev year               AND
    current jdate  =  prev jdate             
then assume "pixel overlap".

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

*******************************************************************************/
int main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";       /* for printing error messages           */
    char msg_str[MAX_STR_LEN];       /* input data scene name                 */
    char filename[MAX_STR_LEN];      /* input binary filenames                */
    int status;                      /* return value from function call       */
    Output_t *rec_cg = NULL;         /* output structure and metadata         */
    bool verbose;                    /* verbose flag for printing messages    */
    int i, k, m, b, k_new;           /* loop counters                         */
    char **scene_list = NULL;        /* 2-D array for list of scene IDs       */
    char **valid_scene_list = NULL;  /* 2-D array for list of filtered        */
                                     /* scene IDs                             */
    FILE *fd;                        /* file descriptor for file              */
                                     /* containing scene names                */
    int num_scenes = MAX_SCENE_LIST; /* number of input scenes defined        */
    int num_c = 8;                   /* max number of coefficients for model  */
    int num_fc = 0;                  /* intialize NUM of Functional Curves    */
    int rec_fc;
    float v_start[NUM_LASSO_BANDS];
    float v_end[NUM_LASSO_BANDS];
    float v_slope[NUM_LASSO_BANDS];
    float v_dif[NUM_LASSO_BANDS];
    float **v_diff;
    int *sdate;
    int *updated_sdate_array;
    Input_meta_t *meta;
    int row, col;
    int landsat_number;
    int clr_sum = 0;                 /* total number of clear cfmask pixels   */
    int sn_sum = 0;                  /* total number of snow  cfmask pixels   */
    int all_sum = 0;                 /* total of all cfmask pixels            */
    float sn_pct;                    /* percent snow cfmask pixels            */
    float clr_pct;                   /* percent clear cfmask pixels           */
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
    int update_num_c = 8;
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
    unsigned char *fmask_buf;       /* cfmask pixel value array.              */
    unsigned char *updated_fmask_buf;/*sub-set of fmask buf, valid pixels only*/
    int **buf;
    FILE ***fp_bin;

    char in_path[MAX_STR_LEN];      /* directory location of input data/files */
    char out_path[MAX_STR_LEN];     /* directory location for output files    */
    char data_type[MAX_STR_LEN];    /* tifs, bip, stdin. Future: bsq, "rods". */
    char scene_list_filename[MAX_STR_LEN]; /* file name containing list of input sceneIDs */
    char scene_list_file[MAX_STR_LEN]; /* optional input argument for file of list of scenes */
    char tmpstr[MAX_STR_LEN];       /* char string for text manipulation      */
    char short_scene[MAX_STR_LEN];  /* char string for text manipulation      */
    int len;                        /* return of strlen for manipulating strings */
    char output_binary[MAX_STR_LEN];/* directory and file name for output.bin */
    char directory[MAX_STR_LEN];
    char scene_name[MAX_STR_LEN];

    int inputs_specified;        /* the number input scenes defined.          */
    int valid_num_scenes;        /* number of scenes after cfmask counts and  */
                                 /* swath overlap eliminated                  */
    int water_sum = 0;           /* counter for cfmask water pixels.          */
    int shadow_sum = 0;          /* counter for cfmask shadow pixels.         */
    int cloud_sum = 0;           /* counter for cfmask cloud pixels.          */
    int fill_sum = 0;            /* counter for cfmask fill pixels.           */
    bool debug = 1;              /* This replaces the "ifdef 0" convention.   */
    bool std_in = 0;             /* For doing lots of ifs.  "stdin"           */
    bool std_out = 0;            /* and "stdout" are reserved words.          */
    bool end_of_file = 0;        /* To identify end of stdin.                 */
    int j;                       /* Loop counter.                             */
    int wrs_path;                /* Worldwide Reference System row            */
    int wrs_row = 0;             /* for the current swath, this               */
    int year;                    /* group of variables is for                 */
    int jday;                    /* filtering out swath overlap, and          */
    int prev_wrs_path = 0;       /* using the first of two scenes in a        */
    int prev_wrs_row = 0;        /* swath, because it is recommended          */
    int prev_year = 0;           /* to use the meta  data from the            */
    int prev_jday = 0;           /* first for things like sun anle,           */
    unsigned char prev_fmask_buf;/* etc. However, always removing a specific  */
    int valid_scene_count = 0;   /* x/y location specified is not valid be-   */
    int swath_overlap_count = 0; /* it may or may not be in an overlap area.  */
    time_t now;
    time (&now);

    if (verbose)
    {
        snprintf (msg_str, sizeof(msg_str), "CCDC version 05.00 start_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);
    }

    /******************************************************************/
    /*                                                                */
    /* Initialize the input and output directory specification.       */
    /* Because they are optional, this prevents fails of strcmp.      */
    /*                                                                */
    /******************************************************************/

    strcpy(in_path, "");
    strcpy(out_path, "");

    /* Read the command-line arguments */
    status = get_args (argc, argv, &row, &col, in_path, out_path, data_type,
                       scene_list_file, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /******************************************************************/
    /*                                                                */
    /* Check for stdin and stdout, and then allocate memory here, for */
    /* pointers used in both stdin and file-system I/O branches,  and */
    /* for pointers used in the ccdc algorithm portion for both cases.*/
    /*                                                                */
    /******************************************************************/

    if (strcmp(data_type, "stdin") == 0)

    {
        std_in = true;
        valid_num_scenes = num_scenes; // worst case
        // bdavis 
        // therefore, need to move down all mallocs possible....
    }
    else
    {
        valid_num_scenes = num_scenes; // worst case
    }

    if (strcmp(out_path, "stdout") == 0)
    {
        std_out = true;
        verbose = false;
        debug = false;
    }

    if (verbose)
    {
        printf("row,col,verbose=%d,%d,%d\n",row,col,verbose);
        printf("Input Directory:  %s\n", in_path);
        printf("Output Directory: %s\n", out_path);
        printf("scene_list_file: %s\n", scene_list_file);
        printf("output_binary: %s\n", output_binary);
    }

    updated_fmask_buf = malloc(valid_num_scenes * sizeof(unsigned char));
    if (updated_fmask_buf == NULL)
    {
        RETURN_ERROR("ERROR allocating updated_fmask_buf memory", FUNC_NAME, FAILURE);
    }

    updated_sdate_array = malloc(valid_num_scenes * sizeof(int));
    if (updated_sdate_array == NULL)
    {
        RETURN_ERROR("ERROR allocating updated_sdate memory", FUNC_NAME, FAILURE);
    }

    /******************************************************************/
    /*                                                                */
    /* allocate memory here, for pointers only used in the stdin I/O  */
    /* branch.                                                        */
    /*                                                                */
    /******************************************************************/

    if (std_in)

    {

        buf = (int **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes, sizeof (int));
        if (buf == NULL)
        {
            RETURN_ERROR ("Allocating buf memory", FUNC_NAME, FAILURE);
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

        clrx = malloc(valid_num_scenes * sizeof(int));
        if (clrx == NULL)
        {
            RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);
        }

        clry = (float **) allocate_2d_array (TOTAL_IMAGE_BANDS, valid_num_scenes,
                                             sizeof (float));
        if (clry == NULL)
        {
            RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
        }

        id_range = (int *)calloc(valid_num_scenes, sizeof(int));
        if (id_range == NULL)
        {
            RETURN_ERROR("ERROR allocating id_range memory", FUNC_NAME, FAILURE);
        }

        /* Allocate memory for temp_v_dif */
        temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));
        if (temp_v_dif == NULL)
        {
            RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
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
    
        /**************************************************************/
        /*                                                            */
        /* For stdin:                                                 */
        /* this assumes order of: julian date value, then 6 SR and 1  */
        /* thermal band values, then cfmask band values, each         */
        /* set/group together, culiminated with a newline per scene,  */
        /* for number of scenes. For example:                         */
        /* 2456445 94 156 164 758 807 492 2809 0                      */
        /*                                                            */
        /**************************************************************/

        i = 0;
        while (!end_of_file)
        {
            if (scanf("%d", &updated_sdate_array[i]) == EOF)
                {
                end_of_file = true;
                break;
                }
            else
            {
                if (debug)
                    printf( "You entered: %d\n", updated_sdate_array[i]);
            }

            for (j = 0; j < TOTAL_IMAGE_BANDS; j++)
            {
                scanf("%d", &buf[j][i]);
                if (debug)
                    printf( "You entered: %d\n", buf[j][i]);
            }

            scanf("%hhu", &updated_fmask_buf[i]);
            if (debug)
                printf( "You entered: %u\n", updated_fmask_buf[i]);

            // bdavis 
            /* this needs to be put in a function call cause it exists in 2 places */
            switch (updated_fmask_buf[i])
            {
                case CFMASK_CLEAR:
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
                    printf ("Unknown fmask value %d", updated_fmask_buf[i]);
                    break;
            }

            i++;
            all_sum++;
        }
   
        valid_num_scenes = i;
        printf ("\n");

    }

    else

    {

        /* allocate memory for scene_list */
        scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                             sizeof (char));
        if (scene_list == NULL)
        {
            RETURN_ERROR ("Allocating scene_list memory", FUNC_NAME, FAILURE);
        }
        valid_scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                             sizeof (char));
        if (valid_scene_list == NULL)
        {
            RETURN_ERROR ("Allocating valid_scene_list memory", FUNC_NAME, FAILURE);
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
    
        fd = fopen(scene_list_filename, "r");
        if (fd == NULL)
        {
            RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
        }

        /******************************************************************/
        /*                                                                */
        /* Fill the scene list array with full path names.                */
        /*                                                                */
        /******************************************************************/

        for (i = 0; i < num_scenes; i++)
        {
            if (fscanf(fd, "%s", tmpstr) == EOF)
                break;
            strcpy(scene_list[i], in_path);
            strcat(scene_list[i], "/");
            strcat(scene_list[i], tmpstr);
        }
        num_scenes = i;
        inputs_specified = num_scenes;
        
        /* Now that we konw the actual number of scenes, allocate     */
        /* memory for date array.                                     */

        sdate = malloc(num_scenes * sizeof(int));
        if (sdate == NULL)
        {
            RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
        }
    
        /* sort scene_list based on year & julian_day */ // then do the swath filter, but read it above first
        if (verbose)
        {
            printf("num_scenes %d\n", num_scenes);
            printf("scene_list[0]=%s\n", scene_list[0]);
        }
        status = sort_scene_based_on_year_doy_row(scene_list, num_scenes, sdate);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                          FUNC_NAME, FAILURE);
        }
    
        /******************************************************************/
        /*                                                                */
        /* Do all of the me memory allocations for buffers/arrays which   */
        /* required to do the reading of files for pixel value assessment */
        /* and storage. Fill value pixels and "pixel overlap" pixels      */
        /* will not be saved and used to update counters.                 */
        /* Also, do the remaining allocations for all buffers/arrays      */
        /* necessary for the ccdc algorithm, because after filling the    */
        /* band data arrays here, that will be the next step.             */
        /*                                                                */
        /******************************************************************/

        /* allocate memory for fp_bin */
        fp_bin = (FILE ***) allocate_2d_array (TOTAL_BANDS, num_scenes,
                                             sizeof (FILE*));
        if (fp_bin == NULL)
        {
            RETURN_ERROR ("Allocating fp_bin memory", FUNC_NAME, FAILURE);
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
        rec_v_dif_copy = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, valid_num_scenes,
                                         sizeof (float));
        if (rec_v_dif_copy == NULL)
        {
            RETURN_ERROR ("Allocating rec_v_dif_copy memory",FUNC_NAME, FAILURE);
        }
    
        /* Allocate memory for temp_v_dif */
        temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, num_scenes,
                                         sizeof (float));
        if (temp_v_dif == NULL)
        {
            RETURN_ERROR ("Allocating temp_v_dif memory",FUNC_NAME, FAILURE);
        }
    
        buf = (int **) allocate_2d_array (TOTAL_BANDS, num_scenes, sizeof (int));
        if (buf == NULL)
        {
            RETURN_ERROR ("Allocating buf memory", FUNC_NAME, FAILURE);
        }
    
        /* Create the Input metadata structure */
        meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
        if (meta == NULL) 
        {
            RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
        }
    
        /* Get the metadata, all scene metadata are the same for stacked scenes */
        if (strcmp(data_type, "tifs") == 0)
        {
            //strcpy(filename, scene_list[0]);
            sprintf(filename, "%s_sr_band1.hdr", scene_list[0]);
        }
        else
        {
            len = strlen(scene_list[0]);
            strncpy(short_scene, scene_list[0], len-5);
            split_directory_scenename(scene_list[0], directory, scene_name);
            if (strncmp(short_scene, ".", 1) == 0)
            {
                strncpy(tmpstr, short_scene + 2, len - 2);
                sprintf(filename, "%s/%s_MTLstack.hdr", tmpstr, scene_name);
            }
            else
                sprintf(filename, "%s/%s_MTLstack.hdr", short_scene, scene_name);
        }

        status = read_envi_header(filename, meta);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling read_envi_header", 
                          FUNC_NAME, FAILURE);
        }
    
        /******************************************************************/
        /*                                                                */
        /* Read the cfmask file first, determine which pixels are valid.  */
        /*   0: clear                                                     */
        /*   1: water                                                     */
        /*   2: cloudShadow                                               */
        /*   3: snow                                                      */
        /*   4: cloud                                                     */
        /* 255: CFMASK_FILL                                               */
        /*                                                                */
        /* While reading, update the counters for clear, water, cloud,    */
        /* and snow, essentially skipping over fill values,               */
        /* and update the number of "scenes" so that only the data bands  */
        /* buffer and date information array are filled with values from  */
        /* from valid "scenes".                                           */
        /*                                                                */
        /* Previously, all data was read, then cfmask was checked after   */
        /* the last loop, and was accidentally checked because the cfmask */
        /* file just happend to be the last file read.  However, when     */
        /* checking for valid values of 50 percent or greater the         */
        /* calculation was done incorecctly, was always zero, so          */
        /* all cases were passed along to the ccdc algorithm, even if     */
        /* there were less than 50 percent valid pixels, and sometimes,   */
        /* even if there was no valid data, or not enough valid data,     */
        /* which caused data-dependent crashes.                           */
        /*                                                                */
        /* Furthermore, FILL_VALUE for cfmask is defined in ccdc.h,       */
        /* but nowhere for the image bands.  Both should really be read   */
        /* from the .xml metadata file provided by ESPA, but it is not    */
        /* populated along with other ARD files, currently. Optionally,   */
        /* user environment variables could be defined and parsed.        */
        /*                                                                */
        /******************************************************************/
    
        valid_num_scenes = 0;
        prev_fmask_buf = 254;

        for (i = 0; i < num_scenes; i++)
        {
    
            /**************************************************************/
            /*                                                            */
            /* Determine the cfmask file name to read.                    */
            /*                                                            */
            /**************************************************************/
    
            len = strlen(scene_list[i]);
            landsat_number = atoi(sub_string(scene_list[i],(len-19),1));
            wrs_path = atoi(sub_string(scene_list[i],(len-18),3));
            wrs_row =  atoi(sub_string(scene_list[i],(len-15),3));
            year = atoi(sub_string(scene_list[i],(len-12),4));
            jday = atoi(sub_string(scene_list[i],(len- 8),3));
            sprintf(filename, "%s_cfmask.img", scene_list[i]);
    
            /**************************************************************/
            /*                                                            */
            /* Open the cfmask file, fseek and read.                      */
            /* if the path, year, and jdate of adjacent sorted scenes are */
            /* the same, and the rows are different by 1, and both fmask  */
            /* values are NOT fill, then this is a case of swath (pixel)  */
            /* overlap, so use the lesser row number of the two, and      */
            /* throw away the grerater row of the two, and update the     */
            /* valid scene list accordingly.   If only one of the two     */
            /* fmask values are FILL, then we are not in an area of swath */
            /* overlap, and the later check for fill will elinimate the   */
            /* uneccesary scene pixels.                                   */
            /*                                                            */
            /**************************************************************/
    
            fp_bin[CFMASK_BAND][i] = open_raw_binary(filename,"rb");
            if (fp_bin[CFMASK_BAND][i] == NULL)
                printf("error open %d scene, %d bands files\n",i, CFMASK_BAND+1);
    
            fseek(fp_bin[CFMASK_BAND][i], (row * meta->samples + col)*sizeof(unsigned char), 
                SEEK_SET);
    
            if (read_raw_binary(fp_bin[CFMASK_BAND][i], 1, 1,
                sizeof(unsigned char), &fmask_buf[i]) != 0)
                printf("error reading %d scene, %d bands\n",i, CFMASK_BAND+1);
    
            if ((wrs_path == prev_wrs_path) && (wrs_row == (prev_wrs_row - 1)) && (year == prev_year) && (jday == prev_jday) && (fmask_buf[i] != CFMASK_FILL) && (prev_fmask_buf != CFMASK_FILL))

            {
                swath_overlap_count++;
                strcpy(valid_scene_list[valid_scene_count - 1], scene_list[i]);
                if (debug)
                {
                    printf("i = %d swath overlap %s\n", i, scene_list[i -1]);
                }
            }
            else
            {
                strcpy(valid_scene_list[valid_scene_count], scene_list[i]);
                valid_scene_count++;

                /**************************************************************/
                /*                                                            */
                /* If clear (0) water (1) or snow (3) update the counters and */
                /* then read the image bands for this scene. Otherwise, save  */
                /* time space and energy and skip this scene.                 */
                /*                                                            */
                /**************************************************************/

                // a bitmask could probably be set up to do these.......
                //
                switch (fmask_buf[i])
                {
                    case CFMASK_CLEAR:
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
                        printf ("Unknown fmask value %d", fmask_buf[i]);
                        break;
                }

                if (fmask_buf[i] < CFMASK_FILL)
                {
                    /**********************************************************/
                    /*                                                        */
                    /* valid pixel according to cfmask, so read the image     */
                    /* bands and update the number of valid scenes, the date  */
                    /* array, and the updated cfmask array.                   */
                    /*                                                        */
                    /**********************************************************/
                    all_sum++;

                    if (debug)
                    {
                        printf("%d ", (int)fmask_buf[i]);
                    }

                    for (k = 0; k < TOTAL_BANDS - 1; k++)
                    {
                        len = strlen(valid_scene_list[valid_num_scenes]);
                        landsat_number = atoi(sub_string(valid_scene_list[valid_num_scenes],(len-19),1));
                        if (landsat_number != 8)
                        {
                            if (k == 5)
                	        sprintf(filename, "%s_sr_band%d.img", valid_scene_list[valid_num_scenes], k+2);
                            else if (k == 6)
                                sprintf(filename, "%s_toa_band6.img", valid_scene_list[valid_num_scenes]);
                            else
                                sprintf(filename, "%s_sr_band%d.img", valid_scene_list[valid_num_scenes], k+1);
                        }
                        else
                        {
                            if (k == 6)
                                sprintf(filename, "%s_toa_band10.img", valid_scene_list[valid_num_scenes]);
                            else 
                                sprintf(filename, "%s_sr_band%d.img", valid_scene_list[valid_num_scenes], k+2);
                        }

                        fp_bin[k][valid_num_scenes] = open_raw_binary(filename,"rb");
                        if (fp_bin[k][valid_num_scenes] == NULL)
                            printf("error open %d scene, %d bands files\n",valid_num_scenes, k+1);
                        if (k != 7)

                        fseek(fp_bin[k][valid_num_scenes], (row * meta->samples + col)*sizeof(short int), SEEK_SET);
                        if (read_raw_binary(fp_bin[k][valid_num_scenes], 1, 1,
                            sizeof(short int), &buf[k][valid_num_scenes]) != 0)
                            printf("error reading %d scene, %d bands\n", valid_num_scenes, (k + 1));
                        
                        close_raw_binary(fp_bin[k][valid_num_scenes]);

                        if (debug)
                        {
                            printf("%d ", (short int)buf[k][valid_num_scenes]);
                        }

                    }

                    // update stuff
                    updated_fmask_buf[valid_num_scenes] = fmask_buf[i];
                    if (debug)
                    {
                        printf("%d ", sdate[i]);
                    }
                    updated_sdate_array[valid_num_scenes] = sdate[i];
                    valid_num_scenes++;
                }

                close_raw_binary(fp_bin[CFMASK_BAND][i]);
            }
        prev_wrs_path = wrs_path;
        prev_wrs_row  = wrs_row;
        prev_year = year;
        prev_jday = jday;
        prev_fmask_buf = fmask_buf[i];
        }

    } // end of elseif stdin bracket

    if (verbose)
    {
        snprintf (msg_str, sizeof(msg_str), "CCDC read_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);

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
        printf ("DEBUG: Number of pixel overlap scenes removed: %d\n", swath_overlap_count);
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
        printf("  Number inputs specified      = %d\n", inputs_specified);
        printf("  Number of non-overlap pixels = %d\n", (inputs_specified - swath_overlap_count));
        printf("  Number of fill (255)  pixels = %d\n", fill_sum);
        printf("  Number of non-fill    pixels = %d\n", all_sum);
        printf("  Number of clear  (0)  pixels = %d\n", clr_sum);
        printf("  Number of water  (1)  pixels = %d\n", water_sum);
        printf("  Number of shadow (2)  pixels = %d\n", shadow_sum);
        printf("  Number of snow   (3)  pixels = %d\n", sn_sum);
        printf("  Number of cloud  (4)  pixels = %d\n", cloud_sum);
        printf("  Number of clear+water pixels = %d\n", clr_sum);
        printf("  Percent of clear pixels      = %f (of non-fill pixels)\n", clr_pct);
        printf("  Percent of clear pixels      = %f (of non-fill, non-cloud, non-shadow pixels)\n",
               (float) clr_sum / (float) (clr_sum + sn_sum));
        printf("  Percent of snow  pixels      = %f (of non-fill, non-cloud, non-shadow pixels)\n", sn_pct);
    }

    // if clr pct less than 50, return error, this syntax is invalid
    //if (clr_sum < (int) (0.5 * valid_num_scenes))
    // bdavis 
    // 20160211, per zhezhu, not that we are processing from gridded
    // inputs, everything is pixel-based, not scene-based, so the
    // algorithm itslef can determine whether there are enough pixels
    // to make a judgement.
    /*
    if (clr_sum < (int) (0.5 * all_sum))
        {
            RETURN_ERROR ("Not enough clear-sky pixels", FUNC_NAME, FAILURE);
        }
     */

    // bdavis 
    // not sure what this means.....
    /* CHANGE: need change back to 0-6 from 1-7 if original 
       inputs are used ???????? */
    /* pixel value ranges should follow physical rules */

    for (i = 0; i < valid_num_scenes; i++)
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
            n_sn = 0;
            /* snow observations are "good" now */
            for (i = 0; i < valid_num_scenes; i++)
            { 
	        if (((updated_fmask_buf[i] == CFMASK_SNOW) || (updated_fmask_buf[i] < 2)) 
                     && id_range[i] == 1)
                {
                    clrx[n_sn] = updated_sdate_array[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		    {
         	        clry[k][n_sn] = (float)buf[k][i];
		    }
                    n_sn++;
                }
            }  
            end = n_sn;

            /* Remove repeated ids */
            matlab_unique(clrx, clry, n_sn, &end);

            free(updated_fmask_buf);
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
                for (k = 0; k < MAX_NUM_C; k++)
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
            n_clr = 0;
            for (i = 0; i < valid_num_scenes; i++)
            { 
                if (id_range[i] == 1)
                {
                    clrx[n_clr] = updated_sdate_array[i];
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		      clry[k][n_clr] = (float)buf[k][i];
                    n_clr++;
                }   
            }
            end = n_clr;

            /* Remove repeated ids */
            matlab_unique(clrx, clry, n_clr, &end);

            free(updated_fmask_buf);
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
            quick_sort_float(clry[1], 0, end - 1);
            matlab_2d_float_median(clry, 1, end, &band2_median);

            n_clr = 0;
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
                for (k = 0; k < MAX_NUM_C; k++)
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
        if (verbose)
        {
            printf("Seasonal Snow (Snow < %f)\n", 100.0 * sn_pct);
            printf("Fmask works, clear pixels (land/water) = %f\n", 100.0 * clr_pct);
        }

        n_clr = 0;
        for (i = 0; i < valid_num_scenes; i++)
        { 
            if ((updated_fmask_buf[i] < 2) && (id_range[i] == 1))
            {
                clrx[n_clr] = updated_sdate_array[i];
                for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
		{
		    clry[k][n_clr] = (float)buf[k][i];
		}
                n_clr++;
            }   
        }
        end = n_clr;
        if (debug)
        {
            printf("end_clr=%d\n",end);
        }

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

        if (!std_in)
        {
            for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                printf("k,adj_rmse[k]=%d,%f\n",k,adj_rmse[k]);
        }

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
        if (verbose)
        {
            snprintf (msg_str, sizeof(msg_str), "CCDC init_time=%s\n", ctime (&now));
            LOG_MESSAGE (msg_str, FUNC_NAME);
        }

        while (i <= end - CONSE)
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
                    end = k_new;

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

                    /* remove noise in original arrays */
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

                                /* update i_start if i_ini is not a confirmed break */
                                i_start = i_ini;
			    }

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
                            quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                            matlab_2d_float_median(v_dif_mag, i_b, ini_conse, 
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
                    // if (i_count == 0 || ids_old_len < (N_TIMES * MAX_NUM_C))
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
                       /* record postion of the pixel */
                       rec_cg[num_fc].pos.row = row; 
                       /* record postion of the pixel */
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

                        if (verbose)
                        {
                            printf("Change Magnitude = %.2f\n", break_mag - T_CG);
                        }

                        /* record break time */
                        rec_cg[num_fc].t_break = clrx[i];
                        rec_cg[num_fc].change_prob = 1.0;
                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
			{
                            quick_sort_float(v_dif_mag[i_b], 0, CONSE-1);
                            matlab_2d_array_mean(v_dif_mag, 0, CONSE, 
                                &rec_cg[num_fc].magnitude[i_b]);
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
            end = valid_num_scenes - 1;
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
                    quick_sort_float(v_dif_mag[i_b], id_last, CONSE-2);
                    matlab_2d_partial_mean(v_dif_mag, i_b, id_last, CONSE-1, 
                                         &rec_cg[num_fc].magnitude[i_b]);
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
                    for (k = 0; k < MAX_NUM_C; k++)
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
                    rec_cg[num_fc].magnitude[i_b] = 0.0; 
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
    }

    free(updated_fmask_buf);
    status = free_2d_array ((void **) buf);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: buf\n", 
                      FUNC_NAME, FAILURE);
    }

    /* Free memory allocation */
    free(updated_sdate_array);
    free(clrx);
    free(rmse);
    free(vec_mag);
    free(id_range);
    status = free_2d_array ((void **) rec_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif\n", 
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) rec_v_dif_copy);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: rec_v_dif_copy\n",
                      FUNC_NAME, FAILURE);
    }
    status = free_2d_array ((void **) v_dif_mag);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
                   FUNC_NAME, FAILURE);
    }

    if (!std_in)
    {
        free(fmask_buf);
        free(meta);
        free(sdate);
        free(id_clr);
        free(id_all);
        free(id_sn);
        status = free_2d_array ((void **) fp_bin);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: fp_bin\n", FUNC_NAME,
                          FAILURE);
        }
        status = free_2d_array ((void **) scene_list);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      FAILURE);
        }
        status = free_2d_array ((void **) valid_scene_list);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Freeing memory: valid_scene_list\n", FUNC_NAME,
                      FAILURE);
        }
    }

    free(ids);
    free(ids_old);
    free(rm_ids);
    free(bl_ids);
    status = free_2d_array ((void **) temp_v_dif);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_v_dif\n", 
                      FUNC_NAME, FAILURE);
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
    if (!std_out)
    {
        strcpy(output_binary, out_path);
        strcat(output_binary, "/output.bin");
        if (access(output_binary, F_OK) != 0) /* File does not exist */
            fp_bin_out = fopen(output_binary, "wb");
        else
            fp_bin_out = fopen(output_binary, "ab");
        if (fp_bin_out == NULL)
        {
            RETURN_ERROR ("Opening output.bin file\n", FUNC_NAME,
                          FAILURE);
        }
        if (num_fc == 0)
        {
            status = fwrite(rec_cg, sizeof(Output_t), 1, fp_bin_out);
            if (status != 1)
            {
                RETURN_ERROR ("Writing output.bin file\n", FUNC_NAME, FAILURE);
            }
        }
        else
        {
            status = fwrite(rec_cg, sizeof(Output_t), num_fc-1, fp_bin_out);
            if (status != num_fc)
            {
                RETURN_ERROR ("Writing output.bin file\n", FUNC_NAME, FAILURE);
            }
        }
        fclose(fp_bin_out);   
    }

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
                    if (debug)
                    {
                        printf("i_b,k,rec_cg[0].coefs[i_b][k] = %d,%d,%f\n", 
                                i_b,k,rec_cg[0].coefs[i_b][k]); 
                    }
		}
                if (debug)
                {
                    printf("rec_cg[0].rmse[%d] = %f\n",i_b,rec_cg[0].rmse[i_b]);
                    printf("rec_cg[0].magnitude[%d]=%f\n",i_b,rec_cg[0].magnitude[i_b]); 
                }
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
                        printf("i_b,k,rec_cg[%d].coefs[i_b][k] = %d,%d,%f\n", 
                             i,i_b,k,rec_cg[i].coefs[i_b][k]); 
		    }
                    printf("rec_cg[%d].rmse[i_b] = %f\n",i,rec_cg[i].rmse[i_b]);
                    printf("rec_cg[%d].magnitude[i_b]=%f\n",i,rec_cg[i].magnitude[i_b]); 
                }
	    }
        }
    }

    if (std_out)
    {
        if (num_fc == 0)
	{
            printf("%d\n",rec_cg[0].t_start);
            printf("%d\n",rec_cg[0].t_end);
            printf("%d\n",rec_cg[0].t_break);
            printf("%d\n",rec_cg[0].pos.row);
            printf("%d\n",rec_cg[0].pos.col);
            printf("%d\n",rec_cg[0].num_obs);
            printf("%d\n",rec_cg[0].category);
            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
            {
                for (k = 0; k < update_num_c; k++)
		{
                    if ((debug) || (std_out))
                    {
                        printf("%f\n", 
                                rec_cg[0].coefs[i_b][k]); 
                    }
		}
                if ((debug) || (std_out))
                {
                    printf("%f\n",rec_cg[0].rmse[i_b]);
                    printf("%f\n",rec_cg[0].magnitude[i_b]); 
                }
            }
	}
	else
	{
            for (i = 0; i < num_fc; i++)
            {
                //printf("i=%d\n",i);
                printf("%d\n",rec_cg[i].t_start);
                printf("%d\n",rec_cg[i].t_end);
                printf("%d\n",rec_cg[i].t_break);
                printf("%d\n",rec_cg[i].pos.row);
                printf("%d\n",rec_cg[i].pos.col);
                printf("%d\n",rec_cg[i].num_obs);
                printf("%d\n",rec_cg[i].category);
                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < update_num_c; k++)
		    {
                        // bdavis
                        // I belive the indecies being printed were incorrect.
                        // changed i_b,k,i
                        // to      i,i_b,k
                        printf("%f\n", 
                             rec_cg[i].coefs[i_b][k]); 
                        //printf("i_b,k,rec_cg[%d].coefs[i_b][k] = %d,%d,%f\n", 
                        //     i_b,k,i,rec_cg[i].coefs[i_b][k]); 
		    }
                    printf("%f\n",rec_cg[i].rmse[i_b]);
                    printf("%f\n",rec_cg[i].magnitude[i_b]); 
                }
	    }
        }
    }

    /* Free rec_cg memory*/
    free(rec_cg);

    time (&now);

    if (verbose)
    {
        snprintf (msg_str, sizeof(msg_str), "CCDC end_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);
    }

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
    printf ("Version 05.00\n");
    printf ("\n");
    printf ("usage:\n");
    printf ("ccdc"
            " --row=<input row number>"
            " --col=<input col number>"
            " --in-path=<input directory>"
            " --out-path=<output directory>"
            " --data-type=<tifs|bip|stdin>"
            " [--scene-list-file=<file with list of sceneIDs>]"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --row=: input row number\n");
    printf ("    --col=: input col number\n");
    printf ("    --in-path=: input data directory location\n");
    printf ("    --out-path=: directory location for output files\n");
    printf ("    --data-type=: type of input data files to ingest\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    --scene-list-file=: file name containing list of sceneIDs"
            " (default is all files in in-path)\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("ccdc --help will print the usage statement\n");
    printf ("\n\n");
    printf ("Example:\n");
    printf ("ccdc"
            " --row=3845"
            " --col=2918"
            " --in-path=/data/user/in"
            " --out-path=/home/user/out"
            " --data-type=bip"
            " --scene-list-file=/home/user/scene_list.txt"
            " --verbose\n\n");
    printf ("An example of how to pipe input from stdin and output to stdout:\n");
    printf ("ccdc"
            " --row=3845"
            " --col=2918"
            " --in-path=/data/user/in"
            " --out-path=stdout"
            " --data-type=stdin"
            " --verbose < pixel_value_text_file.txt > coeffs_results_text_file.txt\n\n");
    printf ("The stdout option eliminates the creation of the output binary file, \n");
    printf ("coeffs are just printed to stdout.  It could be "
            "re-directed to a text file, or piped to another program.\n");
    printf ("\nNote: Previously, the ccdc had to be run from the directory"
            " where the input data are located.\n");
    printf ("      Now, input and output directory locations specifications are required.\n\n");
}
