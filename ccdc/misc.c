#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>
#include <dirent.h>
#include <fnmatch.h>
#include <string.h>
#include <limits.h>
#include <errno.h>

#include "2d_array.h"
#include "const.h"
#include "utilities.h"
#include "ccdc.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

#include "defines.h"

/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
FAILURE         Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/5/2015    Song Guo         Original Development
20151203    Brian Davis      Added arguments for input and output
                             directories, and scene list file.
20160401    Brian Davis      Added data type argument.
20160536    Brian Davis      Added arguments for number of rows and
                             number of columns.

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
  2. chi2inv(T_cg, num_bands) = chi2inv(0.99, 5) = 15.0863 
  3. chi2inv(T_max_cg, num_bands) = chi2inv(1-1e-6, 5) = 35.8882 
******************************************************************************/

int get_args
(
    int  argc,             /* I: number of cmd-line args                    */
    char *argv[],          /* I: string of cmd-line args                    */
    int  *row,             /* O: row number for the pixel                   */
    int  *col,             /* O: col number for the pixel                   */
    int  *nrows,           /* O: number of rows (Y) for line reading        */
    int  *ncols,           /* O: number of columns (X) for tile reading     */
    char *in_path,         /* O: directory locaiton for input data          */
    char *out_path,        /* O: directory location for output files        */
    char *data_type,       /* O: data type:tif,bip,stdin.Future: bsq,"rods".*/
    char *scene_list_file, /* O: opitonal file name of list of sceneIDs     */
    bool *verbose          /* O: verbose flag                               */
)
{
    int c;                         /* current argument index                */
    int option_index;              /* index for the command-line option     */
    static int verbose_flag = 0;   /* verbose flag                          */
    char errmsg[MAX_STR_LEN];      /* error message                         */
    char FUNC_NAME[] = "get_args"; /* function name                         */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"row", required_argument, 0, 'r'},
        {"col", required_argument, 0, 'c'},
        {"num-cols", required_argument, 0, 'n'},
        {"num-rows", required_argument, 0, 'n'},
        {"in-path", required_argument, 0, 'i'},
        {"out-path", required_argument, 0, 'o'},
        {"data-type", required_argument, 0, 'd'},
        {"scene-list-file", required_argument, 0, 's'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /******************************************************************/
    /*                                                                */
    /* Turn off getopt_long error msgs as we'll print our own         */
    /*                                                                */
    /******************************************************************/

    opterr = 0;

    /******************************************************************/
    /*                                                                */
    /* Initialize nrows and nocls to 1, the default.  If not          */
    /* specified, we will process 1 row and 1 column (1 pixel) only.  */
    /*                                                                */
    /******************************************************************/

    *nrows = 1;
    *ncols = 1;

    /******************************************************************/
    /*                                                                */
    /* Loop through all the cmd-line options                          */
    /*                                                                */
    /******************************************************************/

    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {

            /**********************************************************/
            /*                                                        */
            /* Out of cmd-line options                                */
            /*                                                        */
            /**********************************************************/

            break;
        }

        switch (c)
        {
            case 0:

                /******************************************************/
                /*                                                    */
                /* If this option set a flag, do nothing else now.    */
                /*                                                    */
                /******************************************************/

                if (long_options[option_index].flag != 0)
		{
                    break;
		}
		sprintf (errmsg, "option %s\n", long_options[option_index].name);
                if (optarg)
		{
		    sprintf (errmsg, "option %s with arg %s\n", 
                             long_options[option_index].name, optarg);
		}
                RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
                break;

            case 'h':              /* help */
                usage ();
                exit (SUCCESS);
                break;

            case 'i':
                strcpy (in_path, optarg);
                break;

            case 'o':
                strcpy (out_path, optarg);
                break;

            case 's':
                strcpy (scene_list_file, optarg);
                break;

            case 'd':
                strcpy (data_type, optarg);
                break;

            case 'r':             
                *row = atoi (optarg);
                break;

            case 'c':             
                *col = atoi (optarg);
                break;

            case 'n':             
                if (strcmp(long_options[option_index].name, "num-cols") == 0)
                    *ncols = atoi (optarg);
                else if (strcmp(long_options[option_index].name, "num-rows") == 0)
                    *nrows = atoi (optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
                break;
        }
    }

    /******************************************************************/
    /*                                                                */
    /* Check the input values.                                        */
    /* If reading single pixel locations from tifs or bip, ncols and  */
    /* nrows must be 1 processing just 1 x/y pixel locaiton.          */
    /* If reading bip_lines, col must be initialized to 0 (first)     */
    /* and then ncols will be used for seeking and reading indecies.  */
    /* An entire row is assumed for now (start at col=0),             */
    /* row/nrows will be added later for entire tiles.                */
    /*                                                                */
    /******************************************************************/

    if (*row < 0)
    {
        sprintf (errmsg, "row number must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }
    if (*col < 0)
    {
        sprintf (errmsg, "col number must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    if (strcmp(data_type, "tifs") == 0) 
    {
        *ncols = 1;
        *nrows = 1;
    }
    else if (strcmp(data_type, "bip_lines") == 0)
    {
        *col = 1;
        if (*ncols < 0)
        {
            sprintf (errmsg, "number of columns must be > 0");
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }
    else if (strcmp(data_type, "tiles") == 0)
    {
        printf("data-type is tiles\n");
        if (*nrows < 0)
        {
            sprintf (errmsg, "number of rows must be > 0");
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
        if (*ncols < 0)
        {
            sprintf (errmsg, "number of columns must be > 0");
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }

    /******************************************************************/
    /*                                                                */
    /* Current valid input types are separate tif files, single pixel */
    /* location from bip envi files, an entire line from a bip envi   */
    /* file, a tile of data of some sort, or values streamed/piped to */
    /* stdin, one group per pixel/scene.                              */
    /*                                                                */
    /******************************************************************/

    else if ((strcmp(in_path, "stdin") != 0) && (strcmp(data_type, "bip") != 0))
    {
        sprintf (errmsg, "in-path must be stdin, or data-type must be one of: tifs; bip; bip-lines; tiles\n");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    /******************************************************************/
    /*                                                                */
    /* If in_path and out_path were not specified, assign local       */
    /* directory, so that pre-pending a directory/path later on will  */
    /* not cause an error, but instead result in ./<filename> .       */
    /*                                                                */
    /******************************************************************/

    if (strlen(in_path) == 0)
    {
        strcpy (in_path, ".");
    }
    if (strlen(out_path) == 0)
    {
        strcpy (out_path, ".");
    }



    if (strcmp(in_path, "stdin") != 0)
    {
        if ((strcmp(data_type, "tifs"     ) != 0) &&
            (strcmp(data_type, "bip"      ) != 0)   &&
            (strcmp(data_type, "bip-lines") != 0))
        {
            sprintf (errmsg, "data-type must be one of: tifs, bip, bip-lines");
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }


    /******************************************************************/
    /*                                                                */
    /* Check the verbose flag                                         */
    /*                                                                */
    /******************************************************************/

    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    /******************************************************************/
    /*                                                                */
    /* We should to do this only here, not back in main. After        */
    /* determining whether stdin/stdout are being used or not,        */
    /* limit the printed output to only that required for stdout, if .*/
    /* that is the case.                                              */
    /*                                                                */
    /******************************************************************/
    
    if ((*verbose) && (strcmp(out_path, "stdout") != 0))
    {
        printf ("row = %d\n", *row);
        printf ("col = %d\n", *col);
        printf ("ncols = %d\n", *ncols);
        printf ("in-path = %s\n", in_path);
        printf ("out-path = %s\n", out_path);
        printf ("scene-list-file = %s\n", scene_list_file);
        printf ("data-type = %s\n", data_type);
        printf ("verbose = %d\n", *verbose);
    }

    return (SUCCESS);
}


/*****************************************************************************
MODULE:  get_scenename

PURPOSE:  get scene name based on full filename even with path

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/21/2015   Song Guo         Original Development

NOTES:
*****************************************************************************/

void get_scenename
(
    const char *filename, /* I: Name of file to split                       */
    char *directory,      /* O: Directory portion of file name              */
    char *scene_name,     /* O: Scene name portion of the file name.        */
    char *appendix        /* O: Appendix portion of the file name           */
)
{
    char file_name[PATH_MAX];   /* Local copy of filename                   */
    char *ptr;                  /* String pointer                           */

    /******************************************************************/
    /*                                                                */
    /* Make a local copy of filename so it is not destroyed           */
    /*                                                                */
    /******************************************************************/

    strcpy (file_name, filename);

    /******************************************************************/
    /*                                                                */
    /* Check for a directory path                                     */
    /* Find ending '/'                                                */
    /*                                                                */
    /******************************************************************/

    ptr = (char *) strrchr (file_name, '/');
    if (ptr != NULL)
    {
        strcpy (directory, file_name);
        ptr = (char *) strrchr (directory, '/');
        ptr++;
        strcpy (file_name, ptr);
        *ptr = '\0';
    }
    else
    {
        strcpy (directory, "");
    }

    /******************************************************************/
    /*                                                                */
    /* Check for the first "_"                                        */
    /*                                                                */
    /******************************************************************/

    ptr = (char *) strchr (file_name, '_');
    if (ptr != NULL)
    {
        *(ptr++) = '\0';
        strcpy (scene_name, file_name);
        strcpy (appendix, ptr);
    }
    else
    {
        strcpy (scene_name, file_name);
        strcpy (appendix, "");
    }
}


/*****************************************************************************
MODULE:  create_scene_list

PURPOSE:  Create scene list from existing files under working data 
          directory and pop them into scene_list string array

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/21/2015   Song Guo         Original Development
20151203    Brian Davis      Added argument for scene list file name.

NOTES:
*****************************************************************************/

int create_scene_list
(
    const char *item,         /* I: string of file items be found           */
    int *num_scenes,          /* I/O: number of scenes                      */
    char *scene_list_filename /* I: file name of list of scene IDs          */
)
{
    DIR *dirp;
    DIR *dirsubp;
    struct dirent *dp;        /* structure for directory entries            */
    struct dirent *sub_dp;    /* structure for sub-directory entries        */
    FILE *fd;                 /* file descriptor for scene list file        */
    char FUNC_NAME[] = "create_scene_list"; /* function name for messages   */
    char directory[MAX_STR_LEN]; /* directory portion of path               */    
    char scene_name[MAX_STR_LEN]; /* scene name portion of path             */  
    char scene_fullname[MAX_STR_LEN]; /* scene name including ground station*/
    char appendix[MAX_STR_LEN]; /* possible extension of file name          */
    int scene_counter = 0;     /* to record number of scenes                */

    fd = fopen(scene_list_filename, "w");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, ERROR);
    }

    dirp = opendir(".");
    if (dirp != NULL)
    {
        while ((dp = readdir(dirp)) != NULL)
        {
            if(fnmatch(item, dp->d_name, 0) == 0)
            {
                get_scenename(dp->d_name,directory,scene_name,appendix);
                dirsubp = opendir(scene_name);
                if (dirsubp != NULL)
                {
                    while ((sub_dp = readdir(dirsubp)) != NULL)
                    {
                        if(fnmatch("L*_MTLstack", sub_dp->d_name, 0) == 0)
                        {
                            get_scenename(sub_dp->d_name,directory,scene_fullname,appendix);
                            fprintf(fd, "%s\n", scene_fullname);
                            scene_counter++;
                        }
                    }
                }
                closedir(dirsubp);
            }
        }
    }
    (void) closedir(dirp);
    fclose(fd);
    *num_scenes = scene_counter;

    return (SUCCESS);
}


/******************************************************************************
MODULE:  partition

PURPOSE:  partition used for the quick_sort routine

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value   

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/

int partition (int arr[], char *brr[], int crr[], int drr[], int left, int right)
{
    int i = left, j = right;
    int tmp, tmp2, tmp3;
    char temp[MAX_STR_LEN];
    int pivot = arr[(left + right) / 2];

    while (i <= j)
    {
        while (arr[i] < pivot)
	{
            i++;
	}
        while (arr[j] > pivot)
	{
            j--;
	}
        if (i <= j)
        {
            tmp = arr[i];
            strcpy(&temp[0], brr[i]);
            tmp2 = crr[i];
            tmp3 = drr[i];
            arr[i] = arr[j];
            strcpy(brr[i], brr[j]);
            crr[i] = crr[j];
            drr[i] = drr[j];
            arr[j] = tmp;
            strcpy(brr[j],&temp[0]);
            crr[j] = tmp2;
            drr[j] = tmp3;
            i++;
            j--;
        }
    }

    return i;
}


/******************************************************************************
MODULE:  quick_sort

PURPOSE:  sort the scene_list & sdate based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development
20160523    Brian Davis      Added drr

NOTES:
******************************************************************************/

void quick_sort (int arr[], char *brr[], int crr[], int drr[], int left, int right)
{
    int index = partition (arr, brr, crr, drr, left, right);

    if (left < index - 1)
    {
        quick_sort (arr, brr, crr, drr, left, index - 1);
    }
    if (index < right)
    {
        quick_sort (arr, brr, crr, drr, index, right);
    }
}


/******************************************************************************
MODULE:  partition_2d_float

PURPOSE:  partition used for the quick_sort routine

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value   

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/21/2016   Song Guo         Original Development

NOTES:
******************************************************************************/

int partition_2d_float (float arr[], float *brr[], int left, int right)
{
    int i = left, j = right;
    float tmp;
    float tmp2[TOTAL_IMAGE_BANDS];
    float pivot = arr[(left + right) / 2];
    int b;

    while (i <= j)
    {
        while (arr[i] < pivot)
	{
            i++;
	}
        while (arr[j] > pivot)
	{
            j--;
	}
        if (i <= j)
        {
            tmp = arr[i];
	    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                tmp2[b] = brr[b][i];
            arr[i] = arr[j];
	    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                brr[b][i] = brr[b][j];
            arr[j] = tmp;
	    for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                brr[b][j] = tmp2[b];
            i++;
            j--;
        }
    }

    return i;
}


/******************************************************************************
MODULE:  quick_sort_2d_float

PURPOSE:  sort the scene_list & sdate based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/21/2016   Song Guo         Original Development

NOTES:
******************************************************************************/

void quick_sort_2d_float (float arr[], float *brr[], int left, int right)
{
    int index = partition_2d_float (arr, brr, left, right);

    if (left < index - 1)
    {
        quick_sort_2d_float (arr, brr, left, index - 1);
    }
    if (index < right)
    {
        quick_sort_2d_float (arr, brr, index, right);
    }
}


/************************************************************************
FUNCTION: is_leap_year

PURPOSE:
Test if year given is a leap year.

RETURN VALUE:
Type = int
Value    Description
-----    -----------
TRUE     the year is a leap year
FALSE    the year is NOT a leap year

**************************************************************************/
int is_leap_year
( 
    int year        /*I: Year to test         */
)
{
    if (((year % 4) != 0) || (((year % 100) == 0) && ((year % 400) != 0)))
    {
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/******************************************************************************
MODULE:  convert_year_doy_to_jday_from_0000

PURPOSE:  convert day of year in a year to julian day counted from year 0000

RETURN VALUE: int
ERROR           Error for year less than 1973 as input
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development
20160513    Brian Davis      At some point, a fix was made to add doy only
                             after the last year.  Previously, (v04.03) it
                             was erroneously being added in the loop for each
                             of the years, resulting in creating future dates
                             and time travel......

NOTES:
******************************************************************************/
int convert_year_doy_to_jday_from_0000
(
    int year,      /* I: year */
    int doy,       /* I: day of the year */
    int *jday      /* O: julian date since year 0000 */
)
{
    char FUNC_NAME[] = "convert_year_doy_to_jday_from_0000";
    int i;
    int status;


    if (year < 1973)
    {
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, ERROR);
    }

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    if (year != LANDSAT_START_YEAR)
    {
        *jday = JULIAN_DATE_LAST_DAY_1972;
        for (i = LANDSAT_START_YEAR; i < year; i++)
        {
            status = is_leap_year(i);
            if (status == true)
            {
                *jday += LEAP_YEAR_DAYS;
            }
            else
            {
                *jday += NON_LEAP_YEAR_DAYS;
            }
        }
    }
    *jday += doy;

    return (SUCCESS);
}


/******************************************************************************
MODULE:  sort_scene_based_on_year_doy_row

PURPOSE:  Sort scene list based on year and julian day of year and row number

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           error return
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development
2/1/2016    Song Guo         Added row number to quick_sort call. 

NOTES:
******************************************************************************/

int sort_scene_based_on_year_doy_row
(
    char **scene_list,      /* I/O: scene_list, sorted as output             */
    int num_scenes,         /* I: number of scenes in the scene list         */
    int *sdate              /* O: year plus date since 0000                  */
)
{
    int i;                  /* loop counter                                  */
    int status;             /* return of function calls for errror handling  */
    int year, doy;          /* to keep track of year, and day within year    */
    int *yeardoy;           /* combined year day of year as one string       */
    int *row;               /* row of path/row for ordering from same swath  */
    char temp_string[8];    /* for string manipulation                       */
    char temp_string2[5];   /* for string manipulation                       */
    char temp_string3[4];   /* for string manipulation                       */
    char temp_string4[4];   /* for string manipulation                       */
    char errmsg[MAX_STR_LEN]; /* for printing error messages                 */
    char FUNC_NAME[] = "sort_scene_based_on_year_doy_row"; /* function name  */
    int len; /* length of string returned from strlen for string manipulation*/
    char temp[MAX_STR_LEN]; /* for string manipuulation                      */

    /******************************************************************/
    /*                                                                */
    /* Allocate memory for yeardoy                                    */
    /*                                                                */
    /******************************************************************/

    yeardoy = malloc(num_scenes * sizeof(int));
    if (yeardoy == NULL)
    {
        RETURN_ERROR("Allocating yeardoy memory", FUNC_NAME, ERROR);
    }

    row = malloc(num_scenes * sizeof(int));
    if (row == NULL)
    {
        RETURN_ERROR("Allocating row memory", FUNC_NAME, ERROR);
    }
 
    /******************************************************************/
    /*                                                                */
    /* Get year plus doy from scene name                              */
    /*                                                                */
    /******************************************************************/

    for (i = 0; i < num_scenes; i++)
    {
        len = strlen(scene_list[i]);
        strncpy(temp_string, scene_list[i]+(len-12), 7);
        yeardoy[i] = atoi(temp_string);
        strncpy(temp_string2, scene_list[i]+(len-12), 4);
        year = atoi(temp_string2);
        strncpy(temp_string3, scene_list[i]+(len-8), 3);
        doy = atoi(temp_string3);
        strncpy(temp_string4, scene_list[i]+(len-15), 3);
        row[i] = atoi(temp_string4);
        status = convert_year_doy_to_jday_from_0000(year, doy, &sdate[i]);
        if (status != SUCCESS)
        {
            sprintf(errmsg, "Converting year %d doy %d", year, doy);
            RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
        }
    }

    /******************************************************************/
    /*                                                                */
    /* Sort the scene_list & sdate based on yeardoy                   */
    /*                                                                */
    /******************************************************************/

    quick_sort(yeardoy, scene_list, sdate, row, 0, num_scenes - 1);

    /******************************************************************/
    /*                                                                */
    /* If two identical date data exist (this is only the case        */
    /* for ARD data, then put the smaller row number data in          */
    /* the front.                                                     */
    /*                                                                */
    /******************************************************************/

    for (i = 0; i < num_scenes - 1; i++)
    {
        if (sdate[i+1] == sdate[i])
	{
            if (row[i] > row[i+1])
	    {
                strcpy(temp, scene_list[i]);
                strcpy(scene_list[i], scene_list[i+1]);
                strcpy(scene_list[i+1], temp);
	    }
	}
    }

    /******************************************************************/
    /*                                                                */
    /* Free memory                                                    */
    /*                                                                */
    /******************************************************************/

    free(yeardoy);
    free(row);

    return (SUCCESS);

}


/******************************************************************************
MODULE:  update_cft

PURPOSE:  determine the number of coefficient use in the time series model

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void update_cft
(
    int i_span,
    int n_times,
    int min_num_c,
    int mid_num_c,
    int max_num_c,
    int num_c,
    int *update_number_c
)
{
    /* start with 4 coefficients model */ 
    if (i_span < mid_num_c * n_times)
    {
        *update_number_c = min(min_num_c, num_c); 
    }
    /* start with 6 coefficients model */
    else if (i_span < max_num_c * n_times) 
    {
        *update_number_c = min(mid_num_c, num_c);  
    }
    /* start with 8 coefficients model */ 
    else
    {
        *update_number_c = min(max_num_c, num_c); 
    }

}

/******************************************************************************
MODULE:  partition_float

PURPOSE:  partition the sorted list

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value   

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int partition_float (float arr[], int left, int right)
{
    int i = left, j = right;
    float tmp;
    float pivot = arr[(left + right) / 2];

    while (i <= j)
    {
        while (arr[i] < pivot)
	{
            i++;
	}
        while (arr[j] > pivot)
	{
            j--;
	}
        if (i <= j)
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  quick_sort_float

PURPOSE:  sort the scene_list based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void quick_sort_float(float arr[], int left, int right)
{
    int index = partition_float (arr, left, right);

    if (left < index - 1)
    {
        quick_sort_float (arr, left, index - 1);
    }
    if (index < right)
    {
        quick_sort_float (arr, index, right);
    }
}

/******************************************************************************
MODULE:  median_variogram

PURPOSE:  simulate matlab medium function for 2d array case

RETURN VALUE: int
ERROR in allocating memory
SUCCESS non-error encounted

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
5/19/2015   Song Guo         Original Development
                             Changed:
                             output_array[i] = (var[m] + var[m+1]) / 2.0
                             to
                             output_array[i] = (var[m-1] + var[m]) / 2.0

******************************************************************************/
int median_variogram
(
    float **array,      /* I: input array                                    */
    int dim1_len,       /* I: dimension 1 length in input array              */
    int dim2_start,     /* I: dimension 2 start index                        */
    int dim2_end,       /* I: dimension 2 end index                          */
    float *output_array /* O: output array                                   */
)
{
    int i, j;           /* loop indecies                                     */
    float *var;         /* pointer for allocation variable memory            */
    int dim2_len = dim2_end - dim2_start + 1; /* perhaps should get defined  */
    int m = dim2_len / 2 - 1;                 /* later?                      */
    char FUNC_NAME[] = "median_variogram"; /* for error messages             */

    if (dim2_len == 1)
    {
        for (i = 0; i < dim1_len; i++)
        {
            output_array[i] = array[i][dim2_start];
            return (SUCCESS);
        }
    }

    var = malloc((dim2_len-1) * sizeof(float));
    if (var == NULL)
    {
        RETURN_ERROR ("Allocating var memory", FUNC_NAME, ERROR);
    }

    for (i = 0; i < dim1_len; i++)
    {
        for (j = dim2_start; j < dim2_end; j++)
        {
            var[j] = abs(array[i][j+1] - array[i][j]);
        }
        quick_sort_float(var, dim2_start, dim2_end-1);
        if ((dim2_len-1) % 2 == 0)
	{
            output_array[i] = (var[m] + var[m+1]) / 2.0;
	}
        else
            output_array[i] = var[m];
    }

    free(var);

    return (SUCCESS);
}


/******************************************************************************
NAME:           split_directory_scenename

PURPOSE:
Split the specified filename into a directory, a root file name, and an
extension.  The last character of the directory path will be a '/'.

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/16/2016   Song Guo         Original development
******************************************************************************/
#include <string.h>
#include <limits.h>           /* For PATH_MAX                                */
#include "error.h"
#include "input.h"
#include "const.h"

void split_directory_scenename
(
    const char *filename,     /* I: Name of scene with path to split         */
    char *directory,          /* O: Directory portion of file name           */
    char *scene_name          /* O: Scene name portion of the file name      */
)
{
    char file_name[PATH_MAX]; /* Local copy of filename                      */
    char *ptr = NULL;         /* String pointer                              */

    /******************************************************************/
    /*                                                                */
    /* Make a local copy of filename so it is not destroyed */
    /*                                                                */
    /******************************************************************/

    strcpy(file_name, filename);

    /******************************************************************/
    /*                                                                */
    /* Check for a directory path                                     */
    /* Find ending '/'                                                */
    /*                                                                */
    /******************************************************************/

    ptr = (char *) strrchr(file_name, '/');
    if (ptr != NULL)
    {
        *(ptr++) = '\0';
        strcpy (directory, file_name);
    }
    else
        strcpy (directory, "");

    /******************************************************************/
    /*                                                                */
    /* Obtain the scene name                                          */
    /*                                                                */
    /******************************************************************/

    strcpy(scene_name, ptr);

    return;
}


/******************************************************************************
MODULE:  rmse_from_square_root_mean

PURPOSE:  simulate matlab calculate rmse from square root mean

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
6/26/2015   Song Guo         Original Development

******************************************************************************/
void rmse_from_square_root_mean
(
    float **array,      /* I: input array */
    float fit_cft,      /* I: input fit_cft value */
    int dim1_index,     /* I: dimension 1 index in input array */
    int dim2_len,       /* I: dimension 2 length */
    float *rmse         /* O: output rmse */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < dim2_len; i++)
    {
         sum += (array[dim1_index][i] - fit_cft) * 
                (array[dim1_index][i] - fit_cft); 
    } 
    *rmse = sqrt(sum / dim2_len); 
}

/******************************************************************************
MODULE:  partition_int

PURPOSE:  partition used for the quick_sort_int routine

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value   

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int partition_int (int arr[], int left, int right)
{
    int i = left, j = right;
    int tmp;
    int pivot = arr[(left + right) / 2];

    while (i <= j)
    {
        while (arr[i] < pivot)
	{
            i++;
	}
        while (arr[j] > pivot)
	{
            j--;
	}
        if (i <= j)
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  quick_sort_int

PURPOSE:  sort the scene_list based on integer yeardoy 

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void quick_sort_int (int arr[], int left, int right)
{
    int index = partition_int (arr, left, right);

    if (left < index - 1)
    {
        quick_sort_int (arr, left, index - 1);
    }
    if (index < right)
    {
        quick_sort_int (arr, index, right);
    }
}

/******************************************************************************
MODULE:  partial_square_root_mean

PURPOSE:  simulate square root mean function of paritail of a 2d array

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void partial_square_root_mean
(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */   
    int start,           /* I: number of start elements in 2nd dim */
    int end,             /* I: number of end elements in 2nd dim */
    float **fit_ctf,     /* I: */
    float  *rmse         /* O: output rmse value */
)
{
    int i;
    float sum = 0;
    float rmse_square;

    for (i = start; i <= end; i++)
    {
        sum += ((array[dim1_index][i] - fit_ctf[dim1_index][0]) *
                (array[dim1_index][i] - fit_ctf[dim1_index][0]));
    }
    rmse_square = sum / (float)(end-start+1);
    *rmse = sqrt(rmse_square);
}

/******************************************************************************
MODULE:  matlab_2d_array_norm

PURPOSE:  simulate matlab norm function for 1 dimension in 2d array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/

void matlab_2d_array_norm
(
    float **array,       /* I: input array                                   */
    int dim1_index,      /* I: 1st dimension index                           */
    int dim2_len,        /* I: number of input elements in 2nd dim           */
    float  *output_norm  /* O: output norm value                             */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < dim2_len; i++)
    {
        sum += array[dim1_index][i] * array[dim1_index][i];
    }
    *output_norm = sqrt(sum);
}


/******************************************************************************
MODULE:  matlab_2d_float_median

PURPOSE:  simulate matlab median function for 1 dimesion in 2d array float point
          number case only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
6/26/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_2d_float_median
(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */   
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_median /* O: output norm value */
)
{
    int m = dim2_len / 2;

    if (dim2_len % 2 == 0)
    {
        *output_median = (array[dim1_index][m-1] + array[dim1_index][m]) / 2.0;
    }
    else
    {
        *output_median = array[dim1_index][m];
    }
}

/******************************************************************************
MODULE:  matlab_2d_array_mean

PURPOSE:  simulate matlab mean function for 1 dimension in 2d array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_2d_array_mean
(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */   
    int dim2_len,        /* I: number of input elements in 2nd dim */
    float *output_mean   /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < dim2_len; i++)
    {
        sum += array[dim1_index][i];
    }
    *output_mean = sum / (float)(dim2_len);
}


/******************************************************************************
MODULE:  matlab_float_2d_partial_median

PURPOSE: simulate matlab mean function for partial part of 1 dimesion in 2d
         array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES:
******************************************************************************/

void matlab_float_2d_partial_median
(
    float **array,         /* I: input array */
    int dim1_index,        /* I: 1st dimension index */
    int start,             /* I: number of start elements in 2nd dim */
    int end,               /* I: number of end elements in 2nd dim */
    float *output_median   /* O: output median value */
)
{
    int m = (end - start) / 2;

    if (m % 2 == 0)
    {
        *output_median = (array[dim1_index][m-1] + array[dim1_index][m]) / 2.0;
    }
    else
    {
        *output_median = array[dim1_index][m];
    }
}


/******************************************************************************
MODULE:  matlab_2d_partial_mean

PURPOSE:  simulate matlab mean function for partial part of 1 dimesion in 2d 
          array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_2d_partial_mean
(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = start; i <= end; i++)
    {
        sum += array[dim1_index][i];
    }
    *output_mean = sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  matlab_2d_partial_square_mean

PURPOSE:  simulate matlab square mean function for partial part of 1 dimension 
          in 2d array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/9/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_2d_partial_square_mean
(
    float **array,       /* I: input array */
    int dim1_index,      /* I: 1st dimension index */   
    int start,           /* I: number of start elements in 2nd dim */
    int end,             /* I: number of end elements in 2nd dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    float sum = 0;

    for (i = start; i <= end; i++)
    {
        sum += array[dim1_index][i] * array[dim1_index][i];
    }
    *output_mean = sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  get_ids_length

PURPOSE:  get total number of non-zero elements of id array

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/11/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void get_ids_length
(
    int *id_array,        /* I: input array */
    int start,            /* I: array start index */
    int end,              /* I: array end index */
    int *id_len           /* O: number of non-zero number in the array */
)
{
    int i;
    int length = 0;

    for (i = start; i <= end; i++)
    {
        if (id_array[i] != 0)
	{
            length++;
	}
    }

    *id_len = length;
}


/******************************************************************************
MODULE:  matlab_unique

PURPOSE:  use the first value only for identical date data

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
----------  ---------------  -------------------------------------
 1/11/2016  Song Guo         Original Development
05/23/2016  Song Guo         Changed:
                                 if (clrx[k] == clrx[k-1])
                             to
                                 if (clrx[k] != clrx[k-1])
                             eliminated continue branch.
                             moved for loop index outside 2nd for loop.

NOTES:
******************************************************************************/

void matlab_unique
(
    int *clrx,
    float **clry,
    int nums,
    int *new_nums
)
{
    int k, b;
    int k_new = 0;

    for (k = 1, k_new = 1; k < nums; k++)
    {
        if (clrx[k] != clrx[k-1])
        {
            clrx[k_new] = clrx[k];
            for (b = 0; b < TOTAL_IMAGE_BANDS; b++)
                clry[b][k_new] = clry[b][k];
            k_new++;
        }
    }
    *new_nums = k_new;
}


/******************************************************************************
MODULE:  dofit

PURPOSE: Declare data type and allocate memory and do multiple linear robust
         fit used for auto_robust_fit

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov)
{
  gsl_multifit_robust_workspace * work 
    = gsl_multifit_robust_alloc (T, X->size1, X->size2);

  gsl_multifit_robust (X, y, c, cov, work);
  gsl_multifit_robust_free (work);
}

/******************************************************************************
MODULE:  auto_robust_fit

PURPOSE:  Robust fit for one band

RETURN VALUE: None

HISTORY:
Date       Programmer       Reason
--------   ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development
20160525   Brian Davis      p delcaration from const size_t to int for type
                            mismatch.

NOTES:
******************************************************************************/
void auto_robust_fit
(
    float **clrx,
    float **clry,
    int nums,
    int start,
    int band_index,
    float *coefs
)
{
    int i, j;
    const int p = 5; /* linear fit */
    gsl_matrix *x, *cov;
    gsl_vector *y, *c;

    /******************************************************************/
    /*                                                                */
    /* Defines the inputs/outputs for robust fitting                  */
    /*                                                                */
    /******************************************************************/

    x = gsl_matrix_alloc (nums, p);
    y = gsl_vector_alloc (nums);

    c = gsl_vector_alloc (p);
    cov = gsl_matrix_alloc (p, p);

    /******************************************************************/
    /*                                                                */
    /* construct design matrix x for linear fit                       */
    /*                                                                */
    /******************************************************************/

    for (i = 0; i < nums; ++i)
    {
        for (j = 0; j < p; j++)
        {
            if (j == 0)
	    {
                gsl_matrix_set (x, i, j, 1.0);
	    }
            else
	    {
                gsl_matrix_set (x, i, j, clrx[i][j-1]);
	    }
        }
        gsl_vector_set(y,i,clry[band_index][i+start]);
    }

    /******************************************************************/
    /*                                                                */
    /* perform robust fit                                             */
    /*                                                                */
    /******************************************************************/

    dofit(gsl_multifit_robust_bisquare, x, y, c, cov);

    for (j = 0; j < c->size; j++)
    {
        coefs[j] = gsl_vector_get(c, j);
    }

    /******************************************************************/
    /*                                                                */
    /* Free the memories                                              */
    /*                                                                */
    /******************************************************************/

    gsl_matrix_free (x);
    gsl_vector_free (y);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
}


/******************************************************************************
MODULE:  auto_mask

PURPOSE:  Multitemporal cloud, cloud shadow, & snow masks (global version)

RETURN VALUE:
Type = int
ERROR error out due to memory allocation
SUCCESS no error encounted 

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015    Song Guo         Original Development
20160104    Song Guo         Numerous bug fixes.
20160513    Brian Davis      Added SUCCESS argument to int return.

NOTES:
******************************************************************************/
int auto_mask
(
    int *clrx,
    float **clry,
    int start,
    int end,
    float years,
    float t_b1,
    float t_b2,
    float n_t,
    int *bl_ids
)
{
    char FUNC_NAME[] = "auto_mask";
    int year;
    float w, w2;
    int i;
    float **x;
    float pred_b2, pred_b5;
    int nums;
    float coefs[ROBUST_COEFFS];
    float coefs2[ROBUST_COEFFS];

    nums = end - start + 1;
    /* Allocate memory */
    x = (float **)allocate_2d_array(nums, ROBUST_COEFFS - 1, sizeof(float));
    if (x == NULL)
    {
        RETURN_ERROR("ERROR allocating x memory", FUNC_NAME, ERROR);
    }

    year = ceil(years);
    w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    w2 = w / (float)year;

    for (i = 0; i < nums; i++)
    {
        x[i][0] = cos(w * (float)clrx[i+start]);
        x[i][1] = sin(w * (float)clrx[i+start]);
        x[i][2] = cos(w2 * (float)clrx[i+start]);
        x[i][3] = sin(w2 * (float)clrx[i+start]);
    }

    /******************************************************************/
    /*                                                                */
    /* Do robust fitting for band 2 */
    /*                                                                */
    /******************************************************************/

    auto_robust_fit(x, clry, nums, start, 1, coefs);

    /******************************************************************/
    /*                                                                */
    /* Do robust fitting for band 5 */
    /*                                                                */
    /******************************************************************/

    auto_robust_fit(x, clry, nums, start, 4, coefs2);

    /******************************************************************/
    /*                                                                */
    /* predict band 2 * band 5 refs, bl_ids value of 0 is clear and   */
    /* 1 otherwise                                                    */
    /*                                                                */
    /******************************************************************/

    for (i = 0; i < nums; i++)
    {
        pred_b2 = coefs[0] + coefs[1] * cos((float)clrx[i+start] * w ) + 
                  coefs[2] * sin((float)clrx[i+start] * w ) + coefs[3] * 
                  cos((float)clrx[i+start] * w2 ) 
                  + coefs[4] * sin((float)clrx[i+start] * w2);
        pred_b5 = coefs2[0] + coefs2[1] * cos((float)clrx[i+start] * w ) + 
                  coefs2[2] * sin((float)clrx[i+start] * w ) + coefs2[3] * 
                  cos((float)clrx[i+start] * w2 ) 
                  + coefs2[4] * sin((float)clrx[i+start] * w2);
        if (((clry[1][i+start]-pred_b2) > (n_t * t_b1)) || 
            ((clry[4][i+start]-pred_b5) < -(n_t * t_b2)))
	{
            bl_ids[i] = 1;
	}
        else
	{
            bl_ids[i] = 0;
	}
    }

    /* Free allocated memory */
    if (free_2d_array ((void **) x) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: x\n", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}

/******************************************************************************
MODULE:  auto_ts_predict

PURPOSE:  Using lasso regression fitting coefficients to predict new values

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015    Song Guo         Original Development
20160104    Song Guo         Numerous bug fixes.

NOTES:
******************************************************************************/
int auto_ts_predict
(
    int *clrx,
    float **coefs,
    int df,
    int band_index,
    int start,
    int end,
    float *pred_y
)
{
    char FUNC_NAME[] = "auto_ts_predict";
    int i;
    int nums = end - start + 1;
    float w, w2, w3;
    w = TWO_PI / AVE_DAYS_IN_A_YEAR;
    w2 = 2.0 * w;
    w3 = 3.0 * w;

    for (i = 0; i < nums; i++)
    { 
      if (df ==2)
      {
        pred_y[i]  = coefs[band_index][0] + coefs[band_index][1] * 
              (float)clrx[i+start];
      }
      else if (df == 4)
      {
        pred_y[i]  = coefs[band_index][0] + coefs[band_index][1] * 
              (float)clrx[i+start] + coefs[band_index][2] * 
              cos((float)clrx[i+start] * w ) + coefs[band_index][3] * 
              sin((float)clrx[i+start] * w );
      }
      else if (df == 6)
      {
        pred_y[i]  = coefs[band_index][0] + coefs[band_index][1] * 
              (float)clrx[i+start] + coefs[band_index][2] * 
              cos((float)clrx[i+start] * w ) + coefs[band_index][3] * 
              sin((float)clrx[i+start] * w ) + coefs[band_index][4] * 
              cos((float)clrx[i+start] * w2 ) + coefs[band_index][5] * 
              sin((float)clrx[i+start] * w2 );
      }
      else if (df == 8)
      {
        pred_y[i]  = coefs[band_index][0] + coefs[band_index][1] * 
              (float)clrx[i+start] + coefs[band_index][2] * 
              cos((float)clrx[i+start] * w ) + coefs[band_index][3] * 
              sin((float)clrx[i+start] * w ) + coefs[band_index][4] * 
              cos((float)clrx[i+start] * w2 ) + coefs[band_index][5] * 
              sin((float)clrx[i+start] * w2) + coefs[band_index][6] * 
              cos((float)clrx[i+start] * w3 ) +coefs[band_index][7] * 
              sin((float)clrx[i+start] * w3 );
      }
      else
      { 
        RETURN_ERROR("Unsupported df number", FUNC_NAME, ERROR);
      }
    }

    return (SUCCESS);
}


/******************************************************************************
MODULE:  c_glmnet

PURPOSE:  The R equivalant function in C for Lasso regression fitting 

RETURN VALUE:
Type = int
ERROR error in allocating memories
SUCCESS no error encounted

HISTORY:
Date        Programmer       Original development
--------    ---------------  -------------------------------------
11/23/2015   Song Guo         Geoscience Australia
12/16/2015   Song Guo         Fixed bug to return the 
                              correct order of coefficients

NOTES:
******************************************************************************/

int c_glmnet
(
    int no,		// number of observations (no)
    int ni,		// number of predictor variables (ni)
    double *x,		// input matrix, x[ni][no]
    double *y,		// response vaiable, of dimentions (no)
    int nlam,		// number of lambda values
    double *ulam,	// value of lambda values, of dimentions (nlam)
    double parm,	// the alpha variable

    int *lmu,		// lmu = actual number of lamda values (solutions)
    double cfs[nlam][ni+1]	// results = cfs[lmu][ni + 1]
)
{
    double w[no];	// weight(no), default to a sequence of 1's
    double vp[ni];	// penalty factor (ni), default to a sequence of 1's
    double cl[ni][2];	// lower and upper limits (ni, 2), default (-inf, +inf)
    
    int ka = (ni < 600)? 1:2;
    int jd[2] = {1, 0};
    int ne = ni + 1;			// dfmax = ni + 1
    int nx = min(ne * 2 + 20, ni);	// pmax = min(dfmax * 2 + 20, ni)

    double flmin = 1.0;	// if lambda is NULL, then flmin is derived from lambda.min.ration:
    			// (no < ni)? 0.01:0.0001; otherwise flmin = 1.0
    double thr = 1.0e-07;	// thresh in R, default 1.0e-07
    int isd = 1;		// derivide from standardize in R, default is True = 1
    int intr = 1;		// derived from intercept in R, default is True = 1
    int maxit = 10000;		// default is 10000

    double a0[nlam];	//   a0(lmu) = intercept values for each solution
    double ca[nlam][nx];// ca(nx,lmu) = compressed coefficient values for each solution
    int ia[nx];		//   ia(nx) = pointers to compressed coefficients
    int nin[nlam];	//   nin(lmu) = number of compressed coefficients for each solution
    double rsq[nlam];	//   rsq(lmu) = R**2 values for each solution
    double alm[nlam];	//   alm(lmu) = lamda values corresponding to each solution
    int nlp;		// actual number of passes over the data for all lamda values
    double b[nlam][ni]; // b(ni,lmu) = uncompressed coefficient values for each solution
    int jerr;		// error flag

    int i, j;

    for (i = 0; i < no; i++) 
    {
	w[i] = 1;
    }

    for (i = 0; i < ni; i++) 
    {
	vp[i] = 1;
	cl[i][0] = -INFINITY;
	cl[i][1] = INFINITY;
    }

    elnet_(&ka, &parm, &no, &ni, x, y, w, jd, vp, cl, &ne, &nx,
           &nlam, &flmin, ulam, &thr, &isd, &intr, &maxit,
           lmu, a0, &ca[0][0], ia, nin, rsq, alm, &nlp, &jerr);

    solns_(&ni, &nx, lmu, &ca[0][0], ia, nin, &b[0][0]);

    for (i = 0; i < *lmu; i++)
    {
        cfs[i][0] = a0[i];
        for (j = 0; j < ni; j++)
        {
            cfs[i][j + 1] = b[i][j];
        }
    }

    return jerr;
}


/******************************************************************************
MODULE:  auto_ts_fit

PURPOSE:  Lasso regression fitting with full outputs

RETURN VALUE:
Type = int
ERROR error in allocating memories
SUCCESS no error encounted

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
11/23/2015  Song Guo         Original Development
20151203    Brian Davis      Merge of versions.
20160104    Song Guo         Numerous bug fixes.
20160517    Brian Davis      Removed temporary hack of initializing lmu to 1.
                             There now seem to be no data-dependent crashes
                             with lmu initialized to some strange value(s).
                             Put df and nums correctly into error message
                             sprintf when malloc of y fails.
                             Incorporated fix to initialize coefs to 0.0.

NOTES:
******************************************************************************/
int auto_ts_fit
(
    int *clrx,
    float **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse,
    float **v_dif
)
{
    char FUNC_NAME[] = "auto_ts_fit";
    char errmsg[MAX_STR_LEN];
    float w;
    int i, j;
    double **x;
    double *y;
    int status;
    float *yhat;
    float v_dif_norm = 0.0;
    int nums = 0.0;
    int nlam = 1;		// number of lambda
    double ulam[1] = {20.0, };  // lambda = 20
    double alpha = 1.0;
    int lmu;
    double cfs[nlam][df];

    nums = end - start + 1;
    w = TWO_PI / 365.25;

    /* Allocate memory */
    if (df ==2 || df ==4 || df == 6 || df == 8)
    {
        x = (double **)allocate_2d_array(df - 1, nums, sizeof(double));
        if (x == NULL)
	{
            sprintf(errmsg, "Allocating x memory for %d - 1 times %d size of double", 
                    df, nums);
            RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
	}

	y = malloc(nums * sizeof(double));
        if (y == NULL)
	{
            sprintf(errmsg, "Allocating y memory %d %d", df, nums);
            RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
	}
    }
    else 
    {
        RETURN_ERROR("Unsupported df value", FUNC_NAME, ERROR);
    }

    yhat = (float *)malloc(nums * sizeof(float));
    if (yhat == NULL)
    {
        RETURN_ERROR("Allocating yhat memory", FUNC_NAME, ERROR);
    }

    switch (df)
    {

        case 2:            
            for (i = 0; i < nums; i++)
            {
                x[0][i] = (double)clrx[i+start];
		y[i] = (double)clry[band_index][i+start];
	    }

            status = c_glmnet(nums, df-1, &x[0][0], y, nlam, ulam, alpha, &lmu, cfs);
            if (status != SUCCESS) 
            {
                sprintf(errmsg, "Calling c_glmnet when df = %d", df);
                RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
            }

	    for (i = 0; i < LASSO_COEFFS; i++)
	        coefs[band_index][i] = 0.0;

            for (i = 0; i < lmu; i++) 
            {
	        for (j = 0; j < df; j++) 
                {
	            coefs[band_index][j]= cfs[i][j];
	        }
            }
            break;

       case 4:             
            for (i = 0; i < nums; i++)
            {
                x[0][i] = (double)clrx[i+start];
                x[1][i] = (double)cos(w * (float)clrx[i+start]);
                x[2][i] = (double)sin(w * (float)clrx[i+start]);
		y[i] = (double)clry[band_index][i+start];
	    }

            status = c_glmnet(nums, df-1, &x[0][0], y, nlam, ulam, alpha, &lmu, cfs);
            if (status != SUCCESS) 
            {
                sprintf(errmsg, "Calling c_glmnet when df = %d", df);
                RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
            }

	    for (i = 0; i < LASSO_COEFFS; i++)
	        coefs[band_index][i] = 0.0;

            for (i = 0; i < lmu; i++) 
            {
	        for (j = 0; j < df; j++) 
                {
	            coefs[band_index][j]= cfs[i][j];
	        }
            }
            break;

       case 6:             
            for (i = 0; i < nums; i++)
            {
                x[0][i] = (double)clrx[i+start];
                x[1][i] = (double)cos(w * (float)clrx[i+start]);
                x[2][i] = (double)sin(w * (float)clrx[i+start]);
                x[3][i] = (double)cos(2.0 * w * (float)clrx[i+start]);
                x[4][i] = (double)sin(2.0 * w * (float)clrx[i+start]);
		y[i] = (double)clry[band_index][i+start];
	    }

            status = c_glmnet(nums, df-1, &x[0][0], y, nlam, ulam, alpha, &lmu, cfs);
            if (status != SUCCESS) 
            {
                sprintf(errmsg, "Calling c_glmnet when df = %d", df);
                RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
            }

	    for (i = 0; i < LASSO_COEFFS; i++)
	        coefs[band_index][i] = 0.0;

            for (i = 0; i < lmu; i++) 
            {
	        for (j = 0; j < df; j++) 
                {
	            coefs[band_index][j]= cfs[i][j];
	        }
            }
            break;

       case 8:             
            for (i = 0; i < nums; i++)
            {
                x[0][i] = (double)clrx[i+start];
                x[1][i] = (double)cos(w * (float)clrx[i+start]);
                x[2][i] = (double)sin(w * (float)clrx[i+start]);
                x[3][i] = (double)cos(2.0 * w * (float)clrx[i+start]);
                x[4][i] = (double)sin(2.0 * w * (float)clrx[i+start]);
                x[5][i] = (double)cos(3.0 * w * (float)clrx[i+start]);
                x[6][i] = (double)sin(3.0 * w * (float)clrx[i+start]);
		y[i] = (double)clry[band_index][i+start];
	    }

            status = c_glmnet(nums, df-1, &x[0][0], y, nlam, ulam, alpha, &lmu, cfs);
            if (status != SUCCESS) 
            {
                sprintf(errmsg, "Calling c_glmnet when df = %d", df);
                RETURN_ERROR(errmsg, FUNC_NAME, ERROR);
            }

            for (i = 0; i < LASSO_COEFFS; i++)
                coefs[band_index][i] = 0.0;

	    for (i = 0; i < lmu; i++)
            {
                for (j = 0; j < df; j++)
                {
	            coefs[band_index][j]= cfs[i][j];
	        }
            }
            break;
    }

    /* predict lasso model results */
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        auto_ts_predict(clrx, coefs, df, band_index, start, end, yhat);
        for (i = 0; i < nums; i++)
        {
            v_dif[band_index][i] = clry[band_index][i+start] - yhat[i];
        }
        matlab_2d_array_norm(v_dif, band_index, nums, &v_dif_norm);
        *rmse = v_dif_norm / sqrt((float)(nums - df));
    }

    /* Free allocated memory */
    free(yhat);
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        if (free_2d_array ((void **) x) != SUCCESS)
	{
            RETURN_ERROR ("Freeing memory: x\n", FUNC_NAME, ERROR);
	}
	free(y);
    }

    return (SUCCESS);
}

