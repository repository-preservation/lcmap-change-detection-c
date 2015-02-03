
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>
#include <dirent.h>
#include <fnmatch.h>
#include <string.h>
#include <limits.h>

#include "const.h"
#include "utilities.h"
#include "input.h"

#ifndef max
#define max(a,b) (((a) (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

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
1/2/2013    Gail Schmidt     Original Development
3/15/2013   Song Guo         Changed to support Fmask
9/13/2013   Song Guo         Changed to use RETURN_ERROR
2/19/2014   Gail Schmidt     Modified to utilize the ESPA internal raw binary
                             file format

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
  2. chi2inv(T_cg, num_bands) = chi2inv(0.99, 5) = 15.0863 
  3. chi2inv(T_max_cg, num_bands) = chi2inv(1-1e-6, 5) = 35.8882 
******************************************************************************/
int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    float *min_rmse,          /* I: minimum rmse threshold value */
    int *row,                 /* I: row number for the pixel */
    int *col,                 /* I: col number for the pixel */
    float *t_cg,              /* I: chi-square inversed T_cg */
    float *t_max_cg,          /* I: chi-square inversed T_max_cg for 
                                    last step noise removal */
    int * conse,              /* I: number of points used for change detection */ 
    bool *verbose             /* O: verbose flag */
)
{
    int c;                          /* current argument index */
    int option_index;               /* index for the command-line option */
    static float min_rmse_default = 0.5;    /* default value of minimum RMSE */
    static float t_cg_default = 15.0863;    /* default value of chi2inv(0.99, 5) */
    static float t_max_cg_default = 35.8882;/* default value of t_max_cg */
    static int conse_default = 6;           /* default value of conse */ 
    static int verbose_flag = 0;    /* verbose flag */
    char errmsg[MAX_STR_LEN];       /* error message */
    char FUNC_NAME[] = "get_args";  /* function name */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"min_rmse", required_argument, 0, 'r'},
        {"row", required_argument, 0, 'w'},
        {"col", required_argument, 0, 'l'},
        {"t_cg", required_argument, 0, 't'},
        {"t_max_cg", required_argument, 0, 'm'},
        {"conse", required_argument, 0, 'c'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Assign the default values */
    *min_rmse = min_rmse_default;
    *t_cg = t_cg_default;
    *t_max_cg = t_max_cg_default;
    *conse = conse_default;

    /* Loop through all the cmd-line options */
    opterr = 0; /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {                       /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

            case 'h':              /* help */
                usage ();
                return FAILURE;
                break;

            case 'w':             
                *row = atoi (optarg);
                break;

            case 'l':             
                *col = atoi (optarg);
                break;

            case 'r':             
                *min_rmse = atof (optarg);
                break;

            case 't':             
                *t_cg = atof (optarg);
                break;

            case 'm':             
                *t_max_cg = atof (optarg);
                break;

            case 'c': 
                *conse = atoi(optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, FAILURE);
                break;
        }
    }

    /* Make sure this is some positive value */
    if (*conse <= 0)
    {
        sprintf (errmsg, "conse must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }


    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    if (*verbose)
    {
        printf ("row = %f\n", *row);
        printf ("col = %f\n", *col);
        printf ("min_rmse = %f\n", *min_rmse);
        printf ("t_cg = %f\n", *t_cg);
        printf ("t_max_cg = %f\n", *t_max_cg);
        printf ("conse = %d\n", *conse);
    }

    return SUCCESS;
}

/******************************************************************************
MODULE:  get_scenename

PURPOSE:  get scene name based on full filename even with path

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
1/21/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void get_scenename
(
    const char *filename, /* I: Name of file to split */
    char *directory,      /* O: Directory portion of file name */
    char *scene_name,     /* O: Scene name portion of the file name */
    char *appendix        /* O: Appendix portion of the file name */
)
{
    char file_name[PATH_MAX];   /* Local copy of filename */
    char *ptr;                  /* String pointer */

    /* Make a local copy of filename so it is not destroyed */
    strcpy (file_name, filename);

    /* Check for a directory path */
    /* Find ending '/' */
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
        strcpy (directory, "");

    /* Check for the first "_" */
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

/******************************************************************************
MODULE:  create_scene_list

PURPOSE:  Create scene list from existing files under working data 
          directory and pop them into scene_list string array

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
1/21/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int create_scene_list
(
    const char *item,         /* I: string of file items be found */
    int num_scenes,           /* O: number of scenes */
    char **scene_list,        /* O: scene_list used for ccdc processing */ 
)
{
    DIR *dirp;
    struct dirent *dp;
    FILE *fd;
    char FUNC_NAME[] = "create_scene_list"; /* function name */
    char directory[MAX_STR_LEN];     
    char scene_name[MAX_STR_LEN];   
    char appendix[MAX_STR_LEN];       

    fd = fopen("scene_list.txt", "w");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
    }

    if ((dirp = opendir(".")) == NULL) 
    {
        RETURN_ERROR("Opening current directory", FUNC_NAME, FAILURE);
    }

    do 
    {
        errno = 0;
        if ((dp = readdir(dirp)) != NULL) 
        {
            if(fnmatch(arg, dp->d_name, 0) == 0)
            {               
                get_scenename(dp->d_name,directory,scene_name,appendix);
                fprintf(fd, "%s\n", scene_name);
            }
        }
    } while (dp != NULL);

    if (errno != 0)
        RETURN_ERROR("ERROR reading directory", FUNC_NAME, FAILURE);
    (void) closedir(dirp);
    fclose(fd);

    return SUCCESS;
}

/******************************************************************************
MODULE:  sub_string

PURPOSE:  get sub part of a string

RETURN VALUE:
Type = char
Value           Description
-----           -----------
sub staring

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/22/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
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
MODULE:  partition

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
int partition (int arr[], char *brr[], int left, int right)
{
    int i = left, j = right;
    int tmp;
    char temp[MAX_STR_LEN];
    int pivot = arr[(left + right) / 2];

    while (i <= j)
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j)
        {
            tmp = arr[i];
            strcpy(&temp[0], brr[i]);
            arr[i] = arr[j];
            strcpy(brr[i], brr[j]);
            arr[j] = tmp;
            strcpy(brr[j],&temp[0]);
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  quick_sort

PURPOSE:  sort the scene_list based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void quick_sort (int arr[], char *brr[], int left, int right)
{
    int index = partition (arr, brr, left, right);

    if (left < index - 1)
     quick_sort (arr, brr, left, index - 1);
    if (index < right)
        quick_sort (arr, brr, index, right);
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
        return FALSE;
    else
        return TRUE;
}

/******************************************************************************
MODULE:  convert_year_doy_to_jday_from_0000

PURPOSE:  sort the scene_list based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

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
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, EXIT_FAILURE);
    }

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    *jday = 720624;
    for (i = 1973; i < year; i++)
    {

        if (i == 1973)
            *jday += doy;
        status = is_leap_year(i);
        if (status == true)
            *jday += 366;
        else
            *jday += 365;
        *jday += doy;
    }    

    return SUCCESS;
}

/******************************************************************************
MODULE:  convert_year_doy_to_jday_from_0000

PURPOSE:  sort the scene_list based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int convert_jday_from_0000_to_year_doy
(
    int jday,      /* I: julian date since year 0000 */
    int *year,     /* O: year */
    int *doy       /* O: day of the year */
)
{
    char FUNC_NAME[] = "convert_jday_from_0000_to_year_doy";
    int i;
    int status;

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    jday -= 720624;
    *year = 1973;

    if (jday < 0)
    {
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, EXIT_FAILURE);
    }

    while((is_leap_year(year) && jday > 366) || (!is_leap_year(year) && jday > 365))
    {
        status = is_leap_year(year);
        if (status == true)
        {
            jday -= 366;  
        }
        else
        {
            jday -= 365;  
        }
        *year++;
        *day = jday;
    }

    return SUCCESS;
}

/******************************************************************************
MODULE:  sort_scene_based_on_year_jday

PURPOSE:  Sort scene list based on year and julian day of year

RETURN VALUE:
Type = int
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int sort_scene_based_on_year_doy
(
    char **scene_list,      /* I/O: scene_list, sorted as output */
    int num_scenes,         /* I: number of scenes in the scene list */
    int *sdate              /* O: year plus date since 0000 */
)
{
    int i;
    int status;
    int year, doy;
    int *yeardoy;
    char FUNC_NAME[] = "sort_scene_based_on_year_doy"; /* function name */

    /* Allocate memory for yeardoy */
    yeardoy = malloc(num_scenes * sizeof(int));
    if (yeardoy == NULL)
        RETURN_ERROR("ERROR allocating memory", FUNC_NAME, FAILURE);

    /* Get year plus doy from scene name */
    for (i = 0; i < num_scenes; i++)
    {
        yeardoy[i] = atoi(sub_string(scene_list[i], 9, 7));
        year = atoi(sub_string(scene_list[i], 9, 4));
        doy = atoi(sub_string(scene_list[i], 13, 3));
        status = convert_year_doy_to_jday_from_0000(year, doy, &sdate[i]);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, EXIT_FAILURE);
        }
    }

    /* Sort the scene_list based on yeardoy */
    quick_sort(yeardoy, scene_list, 0, num_scenes - 1);

    /* Free memory */
    free(yeardoy);

    return SUCCESS;

}

/******************************************************************************
MODULE:  update_cft

PURPOSE:  determine the number of coefficientd use in the time series model

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
    if (i_span < mid_num_c * n_times)
        update_num_c = min(min_num_c, num_c);  /* start with 4 coefficients model */ 
    else if (i_span < max_num_c * n_times) 
        update_num_c = min(mid_num_c, num_c);  /* start with 6 coefficients model */ 
    else
        update_num_c = min(max_num_c, num_c);  /* start with 8 coefficients model */ 

}
