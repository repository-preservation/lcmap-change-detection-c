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
//#include "matio.h"
#include "defines.h"

const char scene_list_name[] = {"scene_list.txt"};  /* default, if none specified */
int lasso_blist[NUM_LASSO_BANDS] = {1, 2, 3, 4, 5}; /* This is LASSO band index */


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
---------   -------------    -------------------------------------
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
20160218    Song Guo         Changed to use Zhe's stacked ARD data as inputs
20160421    Brian Davis      Moved the reading from stdin to a function
                             in input.c.  Added and formatted comments.
                             Moved the accumlators of fmask pixel
                             categories to function in input.c.
                             Moved all defines to a new .h, defines.h.
                             There were multiple conflicts in multiple
                             places.  Moved sub_string function to utilities.c
                             Removed various obsolete/unused variables.
                             Moved reading of individual .tif files to
                             input.c
                             Added capability to read ENVI BIP format files
                             (input.c).
                             Latest bug fixes and algorithm updates from
                             Song.
                             Initialized update_num_c to 8.
                             calling matlab_2d_float_median insead of
                                     matlab_2d_float_mean
                             calling matlab_float_2d_partial_median instead of
                                     matlab_2d_partial_mean
20160526    Brian Davis      Added reading of bip and bip_lines
                             More consolodation of things in input.c instead
                             of ccdc.c
                             Changed initialization of num_fc from 0 to 1 (sguo)
                             removed intialization of update_num_c to 8.
                             removed intialization of v_dif_nrm to 0.0.
20160901    Brian Davis      Initialized update_num_c to 8.

When reading fmask values to detmermine fill vs. usability,
those "scenes" have been sorted, so one could test for:
If  current fmask and prev fmask are not FILL AND
    current path   =  prev path               AND
    current row    =  prev row -1             AND
    current year   =  prev year               AND
    current jdate  =  prev jdate             
then assume "pixel swath overlap".

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
    char FUNC_NAME[] = "main";       /* For printing error messages           */
    char msg_str[MAX_STR_LEN];       /* Input data scene name                 */
    int status;                      /* Return value from function call       */
    Output_t *rec_cg = NULL;         /* Output structure and metadata         */
    bool verbose;                    /* Verbose flag for printing messages    */
    int i, k, m, b, k_new;           /* Loop counters                         */
    char **scene_list = NULL;        /* 2-D array for list of scene IDs       */
    char **valid_scene_list = NULL;  /* 2-D array for list of filtered        */
                                     /* scene IDs                             */
    FILE *fd;                        /* File descriptor for file              */
                                     /* containing scene names                */
    int num_scenes = MAX_SCENE_LIST; /* Number of input scenes defined        */
    int num_c = 8;                   /* Max number of coefficients for model  */
    int num_fc = -1;                  /* Intialize NUM of Functional Curves    */
    int rec_fc;                      /* Record num. of functional curves      */
    float v_start[NUM_LASSO_BANDS];  /* Vector for start of observation(s)    */
    float v_end[NUM_LASSO_BANDS];    /* Vector for end of observastion(s)     */
    float v_slope[NUM_LASSO_BANDS];  /* Vector for anormalized slope values   */
    float v_dif[NUM_LASSO_BANDS];    /* Vector for difference values          */
    float **v_diff;
    int *sdate;                      /* Pointer to list of acquisition dates  */
    int *updated_sdate_array;        /* Sdate array after cfmask filtering    */
    int *sdate_pix;                  /* updated list of sdate values          */
    Input_meta_t *meta;              /* Structure for ENVI metadata hdr info  */
    int row, col;                    /* The input indecies of the data frame. */
    int nrows;                       /* the number of rows in "tile" of data  */
    int ncols;                       /* the number of columns in "row" of data*/
    int clr_sum = 0;                 /* Total number of clear cfmask pixels   */
    int sn_sum = 0;                  /* Total number of snow  cfmask pixels   */
    int all_sum = 0;                 /* Total of all cfmask pixels            */
    float sn_pct;                    /* Percent snow cfmask pixels            */
    float clr_pct;                   /* Percent clear cfmask pixels           */
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
    int i_rec;                       /* start of model before noise removal   */
    float v_dif_norm;                /* vector for normalized difference vals */
    int i_count;                     /* Count difference of i each iteration  */
    float **v_dif_mag;               /* vector for magnitude of differences.  */
    int i_conse, i_b;
    float *vec_mag;                  /* permanent storage for magnitude vector*/
    float *vec_magg;                 /* temp storage for magnitude vecor      */
    float v_dif_mean;
    float vec_magg_min;
    float **rec_v_dif;
    float **rec_v_dif_copy;
    float **temp_v_dif;              /* for the thermal band.......           */
    float adj_rmse[TOTAL_IMAGE_BANDS];/* Adjusted RMSE for all bands          */
    float mini_rmse;                 /* Mimimum RMSE                          */
    int bl_tmask;
    int n_rmse;                      /* number of RMSE values                 */
    float tmpcg_rmse[NUM_LASSO_BANDS]; /* to temporarily change RMSE          */
    int d_rt;
    float *d_yr;
    int id_last;                     /* The last stable id.                   */
    float ts_pred_temp;
    FILE *fp_bin_out;                /* Binary output file name.              */
    int ids_old_len;
    int i_break;                     /* for recording break points, i is index*/
    int i_ini;                       /* for recording begin of time, i is index*/
    int ini_conse;                   /* Initial CONSE.                        */
    float break_mag;
    int ids_len;                     /* number of ids, incremented continuously*/
    unsigned char *fmask_buf;       /* cfmask pixel value array.              */
    unsigned char *updated_fmask_buf;/*sub-set of fmask buf, valid pixels only*/
    //int **buf;                      /* This is the image bands buffer.        */
    int *buf;                       /* This is the image bands buffer.        */
    FILE ***fp_tifs;                /* Array of file pointers of multiple     */
                                    /*     band files for specific dates.     */
    FILE **fp_bip;                  /* Array of file pointers of BIP files    */
    char in_path[MAX_STR_LEN];      /* directory location of input data/files */
    char out_path[MAX_STR_LEN];     /* directory location for output files    */
    char data_type[MAX_STR_LEN];    /* tifs, bip. Future: bsq, "rods".        */
    char scene_list_filename[MAX_STR_LEN]; /* file name containing list of input sceneIDs */
    char scene_list_file[MAX_STR_LEN]; /* optional input argument for file of list of scenes */
    char tmpstr[MAX_STR_LEN];       /* char string for text manipulation      */
    char output_binary[MAX_STR_LEN];/* directory and file name for output.bin */
    int inputs_specified;        /* the number input scenes defined.          */
    int valid_num_scenes;        /* number of scenes after cfmask counts and  */
                                 /* swath overlap eliminated                  */
    int clear_sum = 0;           /* counter for cfmask clear pixels.          */
    int water_sum = 0;           /* counter for cfmask water pixels.          */
    int shadow_sum = 0;          /* counter for cfmask shadow pixels.         */
    int cloud_sum = 0;           /* counter for cfmask cloud pixels.          */
    int fill_sum = 0;            /* counter for cfmask fill pixels.           */
    bool debug = 1;              /* This replaces the "ifdef 0" convention.   */
    bool std_in = 0;             /* For doing lots of ifs.  "stdin"           */
    bool std_out = 0;            /* and "stdout" are reserved words.          */
    int prev_wrs_path = 0;       /* using the first of two scenes in a        */
    int prev_wrs_row = 0;        /* swath, because it is recommended          */
    int prev_year = 0;           /* to use the meta  data from the            */
    int prev_jday = 0;           /* first for things like sun anle,           */
    unsigned char prev_fmask_buf;/* etc. However, always removing a specific  */
    int valid_scene_count = 0;   /* x/y location specified is not valid be-   */
    int swath_overlap_count = 0; /* it may or may not be in an overlap area.  */
    time_t now;                  /* For logging the start, stop, and some     */
    time (&now);                 /*     intermediate times.                   */
    int row_inx, col_inx;        /* for looping through 2D array of pixels    */
    int row_count, col_count;    /* to keep track of row/col relative to start*/
    int offset;                  /* for determining buf pointer location of pix*/
    int scene_inx;               /* for looping through scenes                */

    if (verbose)
    {
        snprintf (msg_str, sizeof(msg_str), "CCDC version 05.02 start_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);
    }

    /******************************************************************/
    /*                                                                */
    /* Initialize the input and output directory specification.       */
    /* Because they are optional, this prevents fails of strcmp later.*/
    /*                                                                */
    /******************************************************************/

    strcpy(in_path, "");
    strcpy(out_path, "");

    /******************************************************************/
    /*                                                                */
    /* Read the command-line arguments.                               */
    /*                                                                */
    /******************************************************************/

    status = get_args (argc, argv, &row, &col, &nrows, &ncols, in_path,
                       out_path, data_type, scene_list_file, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /******************************************************************/
    /*                                                                */
    /* Check for stdin and stdout, and then allocate memory here, for */
    /* pointers used in both stdin and file-system I/O branches,  and */
    /* for pointers used in the ccdc algorithm portion for both cases.*/
    /* Not all are used everywhere....                                */
    /*                                                                */
    /******************************************************************/

    if (strcmp(in_path, "stdin") == 0)

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

    updated_fmask_buf = malloc(valid_num_scenes * ncols * nrows * sizeof(unsigned char));
    if (updated_fmask_buf == NULL)
    {
        RETURN_ERROR("ERROR allocating updated_fmask_buf memory", FUNC_NAME, FAILURE);
    }

    updated_sdate_array = malloc(valid_num_scenes * ncols * nrows * sizeof(int));
    if (updated_sdate_array == NULL)
    {
        RETURN_ERROR("ERROR allocating updated_sdate memory", FUNC_NAME, FAILURE);
    }

    /******************************************************************/
    /*                                                                */
    /* Allocate memory here, for pointers only used in the stdin I/O  */
    /* branch.                                                        */
    /*                                                                */
    /******************************************************************/

    if (std_in)

    {

        //buf = (int **) allocate_2d_array (TOTAL_IMAGE_BANDS * ncols * nrows,
        //                                  valid_num_scenes, sizeof (int));
        buf = malloc (TOTAL_IMAGE_BANDS * ncols * nrows * valid_num_scenes * sizeof (int));
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
    
        v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, CONSE,
                    sizeof (float));
        if (v_dif_mag == NULL)
        {
            RETURN_ERROR ("Allocating v_dif_mag memory", 
                                     FUNC_NAME, FAILURE);
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
    
        /**************************************************************/
        /*                                                            */
        /* For stdin:                                                 */
        /* This assumes order of: julian date value, then 6 SR and 1  */
        /* thermal band values, then cfmask band values, each         */
        /* set/group together, culiminated with a newline per scene,  */
        /* for number of scenes. For example:                         */
        /* 2456445 94 156 164 758 807 492 2809 0                      */
        /*                                                            */
        /**************************************************************/

        status = read_stdin (updated_sdate_array, buf, updated_fmask_buf,
                             TOTAL_IMAGE_BANDS, &clr_sum, &water_sum,
                             &shadow_sum, &sn_sum, &cloud_sum, &fill_sum,
                             &all_sum, &valid_num_scenes, debug);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("reading stdin",FUNC_NAME, FAILURE);
        }
        inputs_specified = all_sum;

    }   // end of if std_in

    else

    {   // start of not std_in

        /**************************************************************/
        /*                                                            */
        /* Allocate memory for scene_list.                            */
        /*                                                            */
        /**************************************************************/

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

        /**************************************************************/
        /*                                                            */
        /* Check if scene_list.txt file exists, if not, create the    */
        /* scene_list from existing files in the current data working */
        /* directory.                                                 */
        /*                                                            */
        /**************************************************************/

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

        /**************************************************************/
        /*                                                            */
        /* Fill the scene list array with full path names.            */
        /*                                                            */
        /**************************************************************/

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
        
        /**************************************************************/
        /*                                                            */
        /* Now that we konw the actual number of scenes, allocate     */
        /* memory for date array.                                     */
        /*                                                            */
        /**************************************************************/

        sdate = malloc(num_scenes * nrows * ncols * sizeof(int));
        if (sdate == NULL)
        {
            RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);
        }
    
        /**************************************************************/
        /*                                                            */
        /* Sort scene_list based on year & julian_day, then do the    */
        /* swath filter, but read it above first.                     */
        /*                                                            */
        /**************************************************************/

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
    
        /**************************************************************/
        /*                                                            */
        /* Allocate memory for fp_tifs and all the pointers required   */
        /* for the branch which reads files from the filesystem.      */
        /*                                                            */
        /**************************************************************/

        if (strcmp(data_type, "tifs") == 0)
        {
            fp_tifs = (FILE ***) allocate_2d_array (TOTAL_BANDS, num_scenes,
                                                 sizeof (FILE*));
            if (fp_tifs == NULL)
            {
                RETURN_ERROR ("Allocating fp_tifs memory", FUNC_NAME, FAILURE);
            }
        }
        else if (strcmp(data_type, "bip") == 0)
        {
            fp_bip = (FILE **)malloc(num_scenes * sizeof (FILE*));
            if (fp_bip == NULL)
            {
                RETURN_ERROR ("Allocating fp_bip memory", FUNC_NAME, FAILURE);
            }
        }
    
        /**************************************************************/
        /*                                                            */
        /* Do all of the memory allocations for buffers/arrays which  */
        /* are required to do the reading of files for pixel value    */
        /* assessment and storage. Fill value and "swath overlap"     */
        /* pixels will not be saved nor used to update counters.      */
        /*                                                            */
        /******************************************************************/

        fmask_buf = malloc(num_scenes * ncols * nrows * sizeof(unsigned char));
        if (fmask_buf == NULL)
        {
            RETURN_ERROR("ERROR allocating fmask_buf memory", FUNC_NAME, FAILURE);
        }
    
        //buf = (int **) allocate_2d_array ((TOTAL_IMAGE_BANDS * nrows * ncols), num_scenes, sizeof (int));
        buf = malloc ((TOTAL_IMAGE_BANDS * nrows * ncols * num_scenes * sizeof (int)));
        if (buf == NULL)
        {
            RETURN_ERROR ("Allocating buf memory", FUNC_NAME, FAILURE);
        }
    
    
        /**************************************************************/
        /*                                                            */
        /* Create the Input metadata structure.                       */
        /*                                                            */
        /**************************************************************/

        meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
        if (meta == NULL) 
        {
            RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);
        }
    
/**************************************************/
/*                                                */
/* put all these in a wrapper call like get data  */
/*                                                */
/**************************************************/

        /**************************************************************/
        /*                                                            */
        /* Get the metadata, all scene metadata are the same for      */
        /* stacked scenes.                                            */
        /*                                                            */
        /**************************************************************/

        status = read_envi_header(data_type, scene_list[0], meta);
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
        /* This will determine whether or not to use this pixel, which    */
        /* then requires subsequent reading of image bands.               */
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
            status = read_cfmask(i, data_type, scene_list, row, col, nrows, ncols,
                                 meta->samples, fp_tifs, fp_bip, fmask_buf,
                                 &prev_wrs_path, &prev_wrs_row, &prev_year,
                                 &prev_jday, &prev_fmask_buf, &valid_scene_count,
                                 &swath_overlap_count, valid_scene_list,
                                 &clr_sum, &water_sum, &shadow_sum, &sn_sum,
                                 &cloud_sum, &fill_sum, &all_sum, updated_fmask_buf,
                                 updated_sdate_array, sdate, &valid_num_scenes); // offsets

            if (fmask_buf[i] < CFMASK_FILL)
            {
                /******************************************************/
                /*                                                    */
                /* Valid pixel according to cfmask, so read the image */
                /* bands and update image values buffer.              */
                /*                                                    */
                /******************************************************/

                if (debug)
                {
                    printf("%d %d %d ", i, valid_scene_count -1, updated_sdate_array[valid_scene_count - 1]);
                }
                if (strcmp(data_type, "tifs") == 0)
                {
                    status = read_tifs(valid_scene_list[valid_scene_count - 1],
                                       fp_tifs, (valid_scene_count - 1), row, col,
                                       nrows, ncols, meta->samples, debug, buf); // offsets
                    if (status != SUCCESS)
                    {
                        RETURN_ERROR ("Calling read_tifs", 
                                      FUNC_NAME, FAILURE);
                    }
                }
                else if ((strcmp(data_type, "bip")       == 0) ||
                         (strcmp(data_type, "bip_lines") == 0))
                {
                    printf (" reading bip ");
                    status = read_bip(valid_scene_list[valid_scene_count - 1],
                                      fp_bip, (valid_scene_count - 1), row, col,
                                      nrows, ncols, meta->samples, buf); // offsets
                }
//                else if (strcmp(data_type, "bip_lines") == 0)
//                {
//                    printf ("reading bip lines");
//                    status = read_bip_lines(valid_scene_list[valid_scene_count - 1],
//                                            fp_bip, (valid_scene_count - 1), row, col,
//                                            meta->samples, buf);
//                }

                if (debug)
                {
                    printf("%d\n", updated_fmask_buf[valid_scene_count - 1]);
                }

            }

        }

    } // end of elseif stdin bracket, meaning not stdin, read cfmask and image files

    if ((verbose) && (!std_in))
    {
        /**************************************************************/
        /*                                                            */
        /* Print some info to show how the input metadata works.      */
        /*                                                            */
        /**************************************************************/

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

        /**************************************************************/
        /*                                                            */
        /* Log the time to denote end of I/O, algorithm to follow.    */
        /*                                                            */
        /**************************************************************/

        snprintf (msg_str, sizeof(msg_str), "CCDC read_time=%s\n", ctime (&now));
        LOG_MESSAGE (msg_str, FUNC_NAME);

    }


    /******************************************************************/
    /*                                                                */
    /* Allocate memory for rec_cg.                                    */ 
    /*                                                                */
    /******************************************************************/

    rec_cg = malloc(MAX_NUM_FC * sizeof(Output_t));
    if (rec_cg == NULL)
    {
        RETURN_ERROR("ERROR allocating rec_cg memory", FUNC_NAME, FAILURE);
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
    
        v_dif_mag = (float **) allocate_2d_array(TOTAL_IMAGE_BANDS, CONSE,
                    sizeof (float));
        if (v_dif_mag == NULL)
        {
            RETURN_ERROR ("Allocating v_dif_mag memory", 
                                     FUNC_NAME, FAILURE);
        }
    
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
    
        temp_v_dif = (float **)allocate_2d_array(TOTAL_IMAGE_BANDS, num_scenes,
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

        /**************************************************************/
        /*                                                            */
        /* Do the remaining allocations for all buffers/arrays        */
        /* necessary for the ccdc algorithm, because after filling    */
        /* the band data arrays here, that will be the next step.     */
        /*                                                            */
        /**************************************************************/


        /**************************************************************/
        /*                                                            */
        /* Clear cfmask pixel value accumulators, then call the       */
        /* function for totalling cfmask values.                      */
        /*                                                            */
        /* All this was done in input.c.  Assume valid 2d arrays for  */
        /* cfmask and buf already exist. However, we need to          */
        /* stats to determine whether or not this pixel gets          */
        /* processed.                                                 */
        /*                                                            */
        /* Correct solution is to save these values in read_cfmask.   */
        /*                                                            */
        /**************************************************************/

//        clear_sum = 0;
//        water_sum = 0;
//        shadow_sum = 0;
//        cloud_sum = 0;
//        fill_sum = 0;
//        clr_sum = 0;
//        sn_sum = 0;
//        all_sum = 0;
//        valid_num_scenes = 0;
//
//        for (i = 0; i < valid_num_scenes; i++)
//        {
//            status = assign_cfmask_values (fmask_buf[col * valid_num_scenes + i], &clear_sum, // offsets
//                                           &water_sum, &shadow_sum, &sn_sum,
//                                           &cloud_sum, &fill_sum, &all_sum);
//            if (status != SUCCESS)
//            {
//                RETURN_ERROR ("Calling assign_cfmask_values", FUNC_NAME, FAILURE);
//            }

// assume fmask buf and buf are already filled with only valid pixels
//            if (fmask_buf[col * num_scenes + i] < CFMASK_FILL)
//            {
//                all_sum++;
//                fmask_pix_buf[valid_num_scenes] = fmask_buf[col * num_scenes + i];
//                sdate_pix[valid_num_scenes] = sdate[i];
// fix this row col index offsets
//                for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
//                    pix_buf[b][valid_num_scenes] = buf[col * TOTAL_BANDS + b][i];
//                valid_num_scenes++;
//            }
//        }

        /******************************************************************/
        /*                                                                */
        /* Percent of clear pixels: clear (0) or water (1).               */    
        /*                                                                */
        /******************************************************************/

        clr_pct = (float) clr_sum / (float) all_sum;

        /******************************************************************/
        /*                                                                */
        /* percent of snow observations (3).                              */
        /*                                                                */
        /******************************************************************/

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
                   (float) clr_sum / (float) (clr_sum + sn_sum) * 100);
            printf("  Percent of snow  pixels      = %f (of non-fill, non-cloud, non-shadow pixels)\n", (sn_pct * 100));
        }

        for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
        { 
            /**************************************************************/
            /*                                                            */
            /* pixel value ranges should follow physical rules.i          */
            /* Convert Kelvin to Celsius (for new espa data).             */
            /*                                                            */
            /**************************************************************/
// offsets
            //if (buf[6][i] != -9999)
            //offset = (((row_inx * meta->samples * valid_num_scenes * TOTAL_IMAGE_BANDS) +
            //            (col_inx * valid_num_scenes * TOTAL_IMAGE_BANDS)) +
            //             scene_inx * TOTAL_IMAGE_BANDS);
            offset = (((row_count * meta->samples * valid_num_scenes * TOTAL_IMAGE_BANDS) +
                        (col_count * valid_num_scenes * TOTAL_IMAGE_BANDS)) +
                         scene_inx * TOTAL_IMAGE_BANDS);
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
                id_range[i] = 1;
            }
            else
            {
                id_range[i] = 0;
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
                    offset = (((row_inx * meta->samples) + 
                               (col_inx * TOTAL_BANDS))  +
                                scene_inx);
    	            if (((updated_fmask_buf[scene_inx] == CFMASK_SNOW) || 
                         (updated_fmask_buf[scene_inx] < 2))           &&
                          id_range[scene_inx] == 1)
                    {
                        clrx[n_sn] = updated_sdate_array[scene_inx];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                        {
                 	    clry[k][n_sn] = (float)buf[offset + k]; // offsets
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

                num_fc++;
                if (num_fc >= MAX_NUM_FC)
                {
                    /**************************************************/
                    /*                                                */
                    /* Reallocate memory for rec_cg                   */
                    /*                                                */
                    /**************************************************/
                    rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
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

                rec_cg[num_fc].t_start = clrx[i_start-1]; 
                rec_cg[num_fc].t_end = clrx[end-1]; 

                /**********************************************************/
                /*                                                        */
                /* No break at the moment.                                */
                /*                                                        */
                /**********************************************************/

                rec_cg[num_fc].t_break = 0; 

                /**********************************************************/
                /*                                                        */
                /* Record postion of the pixel.                           */
                /*                                                        */
                /**********************************************************/

                rec_cg[num_fc].pos = (row) * ncols + col_inx + 1; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Record fitted coefficients.                    */
                        /*                                                */
                        /**************************************************/

                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record rmse of the pixel.                          */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                }

                /**********************************************************/
                /*                                                        */
                /* Record change probability, number of observations.     */
                /*                                                        */
                /**********************************************************/

                rec_cg[num_fc].change_prob = 0.0; 
                rec_cg[num_fc].num_obs = n_sn; 
                rec_cg[num_fc].category = 50 + MIN_NUM_C; /* snow pixel */

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */ 
                    /*                                                    */
                    /******************************************************/

                    rec_cg[num_fc].magnitude[i_b] = 0.0; 
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
                    offset = (((row_inx * meta->samples) + 
                               (col_inx * TOTAL_BANDS))  +
                                scene_inx);
                    if (id_range[scene_inx] == 1)
                    {
                        clrx[n_clr] = updated_sdate_array[scene_inx];
                        for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                            clry[k][n_clr] = (float)buf[offset + k]; // offsets
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

                num_fc++;
                if (num_fc >= MAX_NUM_FC)
                {
                    /**************************************************/
                    /*                                                */
                    /* Reallocate memory for rec_cg                   */
                    /*                                                */
                    /**************************************************/
                    rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
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

                rec_cg[num_fc].t_start = clrx[i_start-1]; 
                rec_cg[num_fc].t_end = clrx[end-1]; 

                /**********************************************************/
                /*                                                        */
                /* No break at the moment.                                */
                /*                                                        */
                /**********************************************************/

                rec_cg[num_fc].t_break = 0; 

                /**********************************************************/
                /*                                                        */
                /* Record postion of the pixel.                           */
                /*                                                        */
                /**********************************************************/

                rec_cg[num_fc].pos = (row-1) * ncols + col_inx + 1; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    for (k = 0; k < MAX_NUM_C; k++)
                    {
                        /**************************************************/
                        /*                                                */
                        /* Record fitted coefficients.                    */
                        /*                                                */
                        /**************************************************/

                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record rmse of the pixel.                          */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                }

                /**********************************************************/
                /*                                                        */
                /* Record change probability, number of observations,     */
                /* fit category.                                          */
                /*                                                        */
                /**********************************************************/
                rec_cg[num_fc].change_prob = 0.0; 
                rec_cg[num_fc].num_obs = n_clr; 
                rec_cg[num_fc].category = 40 + MIN_NUM_C; 

                for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                {
                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */ 
                    /*                                                    */
                    /******************************************************/
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
            for (scene_inx = 0; scene_inx < valid_num_scenes; scene_inx++)
            { 
                offset = (((row_inx * meta->samples) + 
                           (col_inx * TOTAL_BANDS))  +
                            scene_inx);
                if ((updated_fmask_buf[scene_inx] < 2) && (id_range[scene_inx] == 1)) // offsets
                {
                    clrx[n_clr] = updated_sdate_array[scene_inx]; // offsets
                    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
                    {
                        clry[k][n_clr] = (float)buf[offset + k]; // offsets
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

            if ((!std_out) && (debug))
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

            num_fc++;
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

            /**************************************************************/
            /*                                                            */
            /* Record the num_fc at the beginning of each pixel.          */
            /*                                                            */
            /**************************************************************/
            rec_fc = num_fc;

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

                                    if (num_fc == rec_fc)
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
                                            if (clrx[k] >= rec_cg[num_fc-1].t_break)
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

                                    if ((num_fc == rec_fc) && ((i_start - i_break) >= CONSE))
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

                                        rec_cg[num_fc].t_end = clrx[i_start-2]; 
                                        rec_cg[num_fc].pos = (row-1) * ncols + col_inx + 1; 

                                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                        {
                                            for (k = 0; k < MIN_NUM_C; k++)
                                            {
                                                /******************************/
                                                /*                            */
                                                /* Record fitted coefficients.*/
                                                /*                            */
                                                /******************************/

                                                rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                             }

                                            /**********************************/
                                            /*                                */
                                            /* Record rmse of the pixel.      */                         
                                            /*                                */
                                            /**********************************/

                                            rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                        }

                                        /**************************************/
                                        /*                                    */
                                        /* Record break time, fit category,   */
                                        /* change probability, time of curve  */
                                        /* start, number of observations,     */
                                        /* change magnitude.                  */
                                        /*                                    */
                                        /**************************************/

                                        rec_cg[num_fc].t_break = clrx[i_start -1];
                                        rec_cg[num_fc].category = 10 + MIN_NUM_C;
                                        rec_cg[num_fc].change_prob = 1.0;
                                        rec_cg[num_fc].t_start = clrx[0];
                                        rec_cg[num_fc].num_obs = i_start - i_break;

                                        for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                                        {
                                            quick_sort_float(v_dif_mag[i_b], 0, ini_conse-1);
                                            matlab_2d_float_median(v_dif_mag, i_b, ini_conse, 
                                                                  &v_dif_mean);
                                            rec_cg[num_fc].magnitude[i_b] = -v_dif_mean; 
                                        }

                                        /**************************************/
                                        /*                                    */
                                        /* Identified and move on for the     */
                                        /* nex tfunctional curve.             */
                                        /*                                    */
                                        /**************************************/

                                        num_fc++;  
                                        if (num_fc >= MAX_NUM_FC)
                                        {
                                            /**********************************/
                                            /*                                */
                                            /* Reallocate memory for rec_cg.  */ 
                                            /*                                */
                                            /**********************************/

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

                            rec_cg[num_fc].t_start = clrx[i_start-1]; 
                            rec_cg[num_fc].t_end = clrx[i-1]; 

                            /**********************************************/
                            /*                                            */
                            /* No break at the moment.                    */
                            /*                                            */
                            /**********************************************/

                            rec_cg[num_fc].t_break = 0; 

                            /**********************************************/
                            /*                                            */
                            /* Record postion of the pixel.               */
                            /*                                            */
                            /**********************************************/

                            rec_cg[num_fc].pos = (row-1) * ncols + col_inx + 1; // offsets

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                for (k = 0; k < MAX_NUM_C; k++)
                                {
                                    /**************************************/
                                    /*                                    */
                                    /* Record fitted coefficients.        */
                                    /*                                    */
                                    /**************************************/

                                    rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k]; 
                                }
                                /******************************************/
                                /*                                        */
                                /* Record rmse of the pixel.              */
                                /*                                        */
                                /******************************************/

                                rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                            }
                            /**********************************************/
                            /*                                            */
                            /* Record change probability, number of       */
                            /* observations, fit category.                */
                            /*                                            */
                            /**********************************************/

                            rec_cg[num_fc].change_prob = 0.0; 
                            rec_cg[num_fc].num_obs = i-i_start+1; 
                            rec_cg[num_fc].category = 0 + update_num_c; 

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                /******************************************/
                                /*                                        */
                                /* Record change magnitude.               */
                                /*                                        */
                                /******************************************/
                                rec_cg[num_fc].magnitude[i_b] = 0.0; 
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

                                        rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                                    } 
                                    /**************************************/
                                    /*                                    */
                                    /* Record rmse of the pixel.          */
                                    /*                                    */
                                    /**************************************/

                                    rec_cg[num_fc].rmse[i_b] = rmse[i_b]; 
                                }
                                /******************************************/
                                /*                                        */
                                /* Record number of observations, fit     */
                                /* category.                              */
                                /*                                        */
                                /******************************************/

                                rec_cg[num_fc].num_obs = i-i_start+1; 
                                rec_cg[num_fc].category = 0 + update_num_c; 

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

                            rec_cg[num_fc].t_end = clrx[i-1];

                            /**********************************************/
                            /*                                            */
                            /* Use fixed number for RMSE computing.       */
                            /*                                            */
                            /**********************************************/

                            n_rmse = N_TIMES * rec_cg[num_fc].category;

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
                                tmpcg_rmse[b] /= sqrt(n_rmse - rec_cg[num_fc].category);
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

                            rec_cg[num_fc].t_break = clrx[i];
                            rec_cg[num_fc].change_prob = 1.0;

                            for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                            {
                                quick_sort_float(v_dif_mag[i_b], 0, CONSE-1);
                                matlab_2d_float_median(v_dif_mag, i_b, CONSE,
                                                       &rec_cg[num_fc].magnitude[i_b]);
                            }
                            /**********************************************/
                            /*                                            */
                            /* Identified and move on for the next        */
                            /* functional curve.                          */
                            /*                                            */
                            /**********************************************/

                            num_fc++;

                            if (num_fc >= MAX_NUM_FC)
                            {
                                /******************************************/
                                /*                                        */
                                /* Reallocate memory for rec_cg.          */ 
                                /*                                        */
                                /******************************************/

                                rec_cg = realloc(rec_cg, (num_fc + 1) * sizeof(Output_t));
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

                rec_cg[num_fc].change_prob = (CONSE - id_last) / CONSE; 
                rec_cg[num_fc].t_end = clrx[end - CONSE + id_last];

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

                    rec_cg[num_fc].t_break = clrx[end-CONSE+id_last+1];

                    /******************************************************/
                    /*                                                    */
                    /* Update magnitude of change.                        */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        quick_sort_float(v_dif_mag[i_b], id_last, CONSE-2);
                        matlab_float_2d_partial_median(v_dif_mag, i_b, id_last, CONSE-1,
                                                       &rec_cg[num_fc].magnitude[i_b]);
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
                if (num_fc == rec_fc)
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

                    rec_cg[num_fc].t_start = clrx[i_start-1];
                    rec_cg[num_fc].t_end = clrx[end-1];
                    rec_cg[num_fc].t_break = 0;
                    rec_cg[num_fc].pos = (row-1) * ncols + col_inx + 1; // offsets

                    /******************************************************/
                    /*                                                    */
                    /* Record fitted coefficients.                        */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        for (k = 0; k < MAX_NUM_C; k++)
                        {
                            rec_cg[num_fc].coefs[i_b][k] = fit_cft[i_b][k];
                        }
                        rec_cg[num_fc].rmse[i_b] = rmse[i_b];
                    }

                    /******************************************************/
                    /*                                                    */
                    /* Record change probability, number of observations, */
                    /* fit category.                                      */
                    /*                                                    */
                    /******************************************************/

                    rec_cg[num_fc].change_prob = 0.0;
                    rec_cg[num_fc].num_obs = i_span;
                    rec_cg[num_fc].category = 20 + MIN_NUM_C; /* simple model fit at the end */

                    /******************************************************/
                    /*                                                    */
                    /* Record change magnitude.                           */
                    /*                                                    */
                    /******************************************************/

                    for (i_b = 0; i_b < TOTAL_IMAGE_BANDS; i_b++)
                    {
                        rec_cg[num_fc].magnitude[i_b] = 0.0; 
                    }
                }
            }
        }
      } // for ncols
    } // for nrows

    free(updated_fmask_buf);
    //status = free_2d_array ((void **) buf);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: buf\n", 
    //                  FUNC_NAME, FAILURE);
    //}
    free(buf);

    /******************************************************************/
    /*                                                                */
    /* Free memory allocations for this section.                      */
    /*                                                                */
    /******************************************************************/

    free(updated_sdate_array);
    //free(clrx);
    free(rmse);
    free(vec_mag);
    //free(vec_magg);
    //free(id_range);
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
        if (strcmp(data_type, "tifs") == 0)
        {
            status = free_2d_array ((void **) fp_tifs);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: fp_tifs\n", FUNC_NAME,
                              FAILURE);
            }
        }
        else if (strcmp(data_type, "bip") == 0)
        {
            free(fp_bip);
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

    //free(ids);
    //free(ids_old);
    //free(rm_ids);
    //free(bl_ids);
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

    /******************************************************************/
    /*                                                                */
    /* Output rec_cg structure to the output file.                    */ 
    /* Note: can use fread to read out the structure from the output  */
    /* file.                                                          */
    /* If output was stdout, skip this step.                          */
    /*                                                                */
    /******************************************************************/

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
            if ( status != (num_fc -1) )
            {
                RETURN_ERROR ("Writing output.bin file\n", FUNC_NAME, FAILURE);
            }
        }
        fclose(fp_bin_out);   
    }

    /******************************************************************/
    /*                                                                */
    /* If one wants to capture the contents of the binary output      */
    /* file, but as test for debugging or information, it gets output */
    /* as text here.                                                  */
    /*                                                                */
    /******************************************************************/

    if (verbose)
    {
        if (num_fc == 0)
	{
            printf("rec_cg[0].t_start=%d\n",rec_cg[0].t_start);
            printf("rec_cg[0].t_end=%d\n",rec_cg[0].t_end);
            printf("rec_cg[0].t_break=%d\n",rec_cg[0].t_break);
            printf("rec_cg[0].pos=%d\n",rec_cg[0].pos);
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
                printf("rec_cg[%d].pos=%d\n",i,rec_cg[i].pos);
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

    /******************************************************************/
    /*                                                                */
    /* If output is to be stdout, then just output only values here,  */
    /* with none of the verbose labels above.                         */
    /*                                                                */
    /******************************************************************/

    if (std_out)
    {
        if (num_fc == 0)
	{
            printf("%d\n",rec_cg[0].t_start);
            printf("%d\n",rec_cg[0].t_end);
            printf("%d\n",rec_cg[0].t_break);
            printf("%d\n",rec_cg[0].pos);
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
                printf("%d\n",rec_cg[i].pos);
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

    /******************************************************************/
    /*                                                                */
    /* Free rec_cg memory, for the final time.                        */
    /*                                                                */
    /******************************************************************/

    free(rec_cg);

    /******************************************************************/
    /*                                                                */
    /* Obtain the current time and log the final completion time.     */
    /*                                                                */
    /******************************************************************/

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
20160304    Brian Davis      Updated to reflect stdin and stdout.

******************************************************************************/
void
usage ()
{
    printf ("\n");
    printf ("Continuous Change Detection and Classification\n");
    printf ("Version 05.02\n");
    printf ("\n");
    printf ("usage:\n");
    printf ("ccdc"
            " --row=<input row number>"
            " --col=<input col number>"
            " --num-rows=<number of rows>"
            " --num-cols=<number of columns>"
            " [--in-path=<input directory>]"
            " [--out-path=<output directory>]"
            " [--data-type=<tifs|bip]"
            " [--scene-list-file=<file with list of sceneIDs>]"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --row=: input row number\n");
    printf ("    --col=: input col number\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    --num-rows=: input number of rows\n");
    printf ("    --num-cols=: input number of columns\n");
    printf ("    --in-path=: input data directory location\n");
    printf ("    --out-path=: directory location for output files\n");
    printf ("    --data-type=: type of input data files to ingest\n");
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
            " --num-rows=1"
            " --num-cols=1"
            " --in-path=/data/user/in"
            " --out-path=/home/user/out"
            " --data-type=bip"
            " --scene-list-file=/home/user/scene_list.txt"
            " --verbose\n\n");
    printf ("An example of how to pipe input from stdin and output to stdout:\n");
    printf ("ccdc"
            " --row=3845"
            " --col=2918"
            " --in-path=stdin"
            " --out-path=stdout"
            " --verbose < pixel_value_text_file.txt > coeffs_results_text_file.txt\n\n");
    printf ("The stdout option eliminates the creation of the output binary file, \n");
    printf ("coeffs are just printed to stdout.  It could be "
            "re-directed to a text file, or piped to another program.\n");
    printf ("\nNote: Previously, the ccdc had to be run from the directory"
            " where the input data are located.\n");
    printf ("      Now, input and output directory locations specifications are used.\n");
    printf ("      If in-path or out-path are not specified, current working directory is assumed.\n");
    printf ("      If scene-file-name is not specified, all scenes in in-path are processed.\n\n");
}
