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
20160914    Brian Davis      Made all refereces to rec_cg[].coefs dimensions
                             [bands][coefs], they had been transposed.
20160914    Brian Davis      fixed some kind of a problem reading/assigning....
20160915    Brian Davis      Change calculation of position to use
                             meta->num_samples for row size instead of ncols,
                             which can vary depending on what was specified
                             for --num-cols argument.
20160923    Brian Davis      moved the main processing to a c-type-function
                             call in the_ccd.c.

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
    Output_t *rec_cg = NULL;         /* Output structure of vals and metadata */
    bool verbose;                    /* Verbose flag for printing messages    */
    int i, k;                        /* Loop counters                         */
    char **scene_list = NULL;        /* 2-D array for list of scene IDs       */
    char **valid_scene_list = NULL;  /* 2-D array for list of filtered        */
                                     /* scene IDs                             */
    FILE *fd;                        /* File descriptor for file              */
                                     /* containing scene names                */
    int num_scenes = MAX_SCENE_LIST; /* Number of input scenes defined        */
    int num_fc = -1;                 /* Intialize NUM of Functional Curves    */
                                     /*num_fc is cummulative,for all rows/cols*/
                                     /* for dumping/writing all outputs.      */
    int *sdate;                      /* Pointer to list of acquisition dates  */
    int *updated_sdate_array;        /* Sdate array after cfmask filtering    */
    Input_meta_t *meta;              /* Structure for ENVI metadata hdr info  */
    int row, col;                    /* The input indecies of the data frame. */
    int nrows;                       /* the number of rows in "tile" of data  */
    int ncols;                       /* the number of columns in "row" of data*/
    int clr_sum = 0;                 /* Total number of clear cfmask pixels   */
    int sn_sum = 0;                  /* Total number of snow  cfmask pixels   */
    int all_sum = 0;                 /* Total of all cfmask pixels            */
    int *clrx;                       /* clear pixel curve in X direction ?    */
    float **clry;                    /* clear pixel curve in Y direction ?    */
    float **fit_cft;                 /* Fitted coefficients 2-D array.        */
    float *rmse;                     /* Root Mean Squared Error array.        */
    int update_num_c=MIN_NUM_C;      /* Number of coefficients to update      */
    int *bl_ids;
    int *id_range;
    int *ids;
    int *ids_old;
    int *rm_ids;
    float **v_dif_mag;               /* vector for magnitude of differences.  */
    int i_b;
    float *vec_mag;                  /* permanent storage for magnitude vector*/
    float **rec_v_dif;
    float **rec_v_dif_copy;
    float **temp_v_dif;              /* for the thermal band.......           */
    FILE *fp_bin_out;                /* Binary output file name.              */
    unsigned char *fmask_buf;       /* cfmask pixel value array.              */
    unsigned char *updated_fmask_buf;/*sub-set of fmask buf, valid pixels only*/
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
                                 updated_sdate_array, sdate, &valid_num_scenes);

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
                                       nrows, ncols, meta->samples, debug, buf);
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
                                      nrows, ncols, meta->samples, buf);
                }

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

    /**************************************************************/
    /*                                                            */
    /* Do the remaining allocations for all buffers/arrays        */
    /* necessary for the ccdc algorithm, because after adjusting  */
    /* the band data arrays, that will be the next step.          */
    /*                                                            */
    /* This section can all be put in a function call to enable   */
    /* topic/c-function-type to be called from the lcmap api      */
    /* after it accesses the IW+DS.                               */
    /*                                                            */
    /*                                                            */
    /**************************************************************/

    status = the_ccd (valid_num_scenes, row, col, nrows, ncols, meta->samples,
                      updated_sdate_array, updated_fmask_buf, buf,
                      clr_sum, sn_sum, water_sum, all_sum, verbose,
                      rec_cg, &num_fc);
    if (status != SUCCESS)
        RETURN_ERROR("ERROR calling the_ccd", FUNC_NAME, FAILURE);


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
    //free(rmse);
    //free(vec_mag);
    //free(vec_magg);
    //free(id_range);
    //status = free_2d_array ((void **) rec_v_dif);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: rec_v_dif\n", 
    //                  FUNC_NAME, FAILURE);
    //}
    //status = free_2d_array ((void **) rec_v_dif_copy);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: rec_v_dif_copy\n",
    //                  FUNC_NAME, FAILURE);
    //}
    //status = free_2d_array ((void **) v_dif_mag);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: v_dif_mag\n", 
    //               FUNC_NAME, FAILURE);
    //}

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
    //status = free_2d_array ((void **) temp_v_dif);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: temp_v_dif\n", 
    //                  FUNC_NAME, FAILURE);
    //}

    //status = free_2d_array ((void **) clry);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
    //                  FAILURE);
    //}

    //status = free_2d_array ((void **) fit_cft);
    //if (status != SUCCESS)
    //{
    //    RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
    //                  FAILURE);
    //}

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
