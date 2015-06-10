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
#include "input.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_randist.h>

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
    int *row,                 /* O: row number for the pixel */
    int *col,                 /* O: col number for the pixel */
    float *t_cg,              /* O: chi-square inversed T_cg */
    bool *verbose             /* O: verbose flag */
)
{
    int c;                          /* current argument index */
    int option_index;               /* index for the command-line option */
    static float t_cg_default = 15.0863;    /* default value of chi2inv(0.99, 5) */
    static int verbose_flag = 0;    /* verbose flag */
    char errmsg[MAX_STR_LEN];       /* error message */
    char FUNC_NAME[] = "get_args";  /* function name */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"row", required_argument, 0, 'w'},
        {"col", required_argument, 0, 'l'},
        {"t_cg", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Assign the default values */
    *t_cg = t_cg_default;

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

            case 't':             
                *t_cg = atof (optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, FAILURE);
                break;
        }
    }

    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    if (*verbose)
    {
        printf ("row = %d\n", *row);
        printf ("col = %d\n", *col);
        printf ("t_cg = %f\n", *t_cg);
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
    int num_scenes,           /* I/O: number of scenes */
    char **scene_list         /* O: scene_list used for ccdc processing */ 
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
            if(fnmatch(item, dp->d_name, 0) == 0)
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
int partition (int arr[], char *brr[], int crr[], int left, int right)
{
    int i = left, j = right;
    int tmp, tmp2;
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
            tmp2 = crr[i];
            arr[i] = arr[j];
            strcpy(brr[i], brr[j]);
            crr[i] = crr[j];
            arr[j] = tmp;
            strcpy(brr[j],&temp[0]);
            crr[j] = tmp2;
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

NOTES:
******************************************************************************/
void quick_sort (int arr[], char *brr[], int crr[], int left, int right)
{
 int index = partition (arr, brr, crr, left, right);

    if (left < index - 1)
        quick_sort (arr, brr, crr, left, index - 1);
    if (index < right)
        quick_sort (arr, brr, crr, index, right);
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
    int status;
    int local_year;

    /* 12-31-1972 is 720624 in julian day since year 0000 */
    jday -= 720624;
    local_year = 1973;

    if (jday < 0)
    {
        RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, EXIT_FAILURE);
    }

    while((is_leap_year(local_year) && jday > 366) || 
          (!is_leap_year(local_year) && jday > 365))
    {
        status = is_leap_year(local_year);
        if (status == true)
        {
            jday -= 366;  
        }
        else
        {
            jday -= 365;  
        }
        local_year++;
    }
    *year = local_year;
    *doy = jday;

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
            RETURN_ERROR ("Landsat data starts from 1973", FUNC_NAME, FAILURE);
        }
    }

    /* Sort the scene_list & sdate based on yeardoy */
    quick_sort(yeardoy, scene_list, sdate, 0, num_scenes - 1);

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
    /* start with 4 coefficients model */ 
    if (i_span < mid_num_c * n_times)
        *update_number_c = min(min_num_c, num_c); 
    /* start with 6 coefficients model */
    else if (i_span < max_num_c * n_times) 
        *update_number_c = min(mid_num_c, num_c);  
    /* start with 8 coefficients model */ 
    else
        *update_number_c = min(max_num_c, num_c); 

}

/******************************************************************************
MODULE:  partition_int16

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
int partition_int16 (int16 arr[], int left, int right)
{
    int i = left, j = right;
    int tmp;
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
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  quick_sort_int16

PURPOSE:  sort the scene_list based on yeardoy string

RETURN VALUE: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/23/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void quick_sort_int16 (int16 arr[], int left, int right)
{
    int index = partition_int16 (arr, left, right);

    if (left < index - 1)
        quick_sort_int16 (arr, left, index - 1);
    if (index < right)
        quick_sort_int16 (arr, index, right);
}

/******************************************************************************
MODULE:  median_variogram

PURPOSE:  simulate matlab medium function for 2d array case

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
5/19/2015   Song Guo         Original Development

******************************************************************************/
int median_variogram
(
    int16 **array,      /* I: input array */
    int dim1_start,     /* I: dimension 1 start index */
    int dim1_end,       /* I: dimension 1 end index */
    int dim2_len,       /* I: dimension 2 length in input array */
    float *output_array /* O: output array */
)
{
    int i, j;
    int16 *var; 
    int dim1_len = dim1_end - dim1_start + 1;
    int m = dim1_len / 2;
    char FUNC_NAME[] = "median_variogram";

    var = malloc((dim1_len-1) * sizeof(int16));
    if (var == NULL)
        RETURN_ERROR ("ERROR allocating memory", FUNC_NAME, FAILURE);

    for (j = 0; j < dim2_len; j++)
    {
        for (i = dim1_start; i < dim1_end; i++)
        {
            var[i] = abs(array[i+1][j] - array[i][j]);
        }
        quick_sort_int16(var, dim1_start, dim1_len-1);
        if (dim1_len % 2 == 0)
            output_array[j] = (float)(var[m-1] + var[m]) / 2;
        else
            output_array[j] = (float)var[m];
    }

    free(var);

    return SUCCESS;
}

/******************************************************************************
MODULE:  partition_int

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
int partition_int (int arr[], int left, int right)
{
    int i = left, j = right;
    int tmp;
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

PURPOSE:  sort the scene_list based on yeardoy string

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
        quick_sort_int (arr, left, index - 1);
    if (index < right)
        quick_sort_int (arr, index, right);
}

/******************************************************************************
MODULE:  array_intersection

PURPOSE:  simulate matlab intersect function for array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/5/2015   Song Guo         Original Development

NOTES: We only handle odd number input N case as it will be 2 * conse + 1 
       for CCDC run, if needed, we can add in even number case later
******************************************************************************/
void array_intersection
(
    int *array1,       /* I: input array 1 */
    int array_len1,    /* I: number of elements in input array1 */
    int *array2,       /* I: input array 2 */
    int array_len2,    /* I: number of elements in input array2 */
    int *output_array, /* O: output array */
    int *output_len    /* O: output array length */
)
{
    int i, j;
    int k = 0;

    if (array_len1 <= array_len2)
    {
        for (i = 0; i < array_len1; i++)
        {
            for (j = 0; j < array_len2; j++)
            {
                if (array1[i] == array2[j])
                {
                    output_array[k] = array1[i];
                    k++;
                    break;
                }
            }
        }
    }
    else
    {
        for (i = 0; i < array_len2; i++)
        {
            for (j = 0; j < array_len1; j++)
            {
                if (array2[i] == array1[j])
                {
                    output_array[k] = array2[i];
                    k++;
                    break;
                }
            }
        }
    }
    *output_len = k;           
    quick_sort_int(output_array, 0, k - 1); 
}

/******************************************************************************
MODULE:  matlab_norm

PURPOSE:  simulate matlab norm function for 1d array cases only

RETURN VALUE:
Type = void
Value           Description
-----           -----------


HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/7/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void matlab_norm
(
    float *array,        /* I: input array */
    int array_len,       /* I: number of elements in input array */
    float  *output_norm  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < array_len; i++)
        sum += array[i] * array[i];
    *output_norm = sqrt(sum);
}

/******************************************************************************
MODULE:  square_root_mean

PURPOSE:  simulate square root mean function

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
void square_root_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements */
    float **fit_ctf,     /* I: */
    float  *rmse         /* O: output rmse value */
)
{
    int i;
    int sum = 0;
    float rmse_square;

    for (i = 0; i < array_len1; i++)
    {
        sum += ((array[i][dim2_number] - fit_ctf[0][dim2_number]) *
                (array[i][dim2_number] - fit_ctf[0][dim2_number]));
    }
    rmse_square = (float)sum / (float)array_len1;
    *rmse = sqrt(rmse_square);
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
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float **fit_ctf,     /* I: */
    float  *rmse         /* O: output rmse value */
)
{
    int i;
    int sum = 0;
    float rmse_square;

    for (i = start; i <= end; i++)
    {
        sum += ((array[i][dim2_number] - fit_ctf[0][dim2_number]) *
                (array[i][dim2_number] - fit_ctf[0][dim2_number]));
    }
    rmse_square = (float)sum / (float)(end-start+1);
    *rmse = sqrt(rmse_square);
}

/******************************************************************************
MODULE:  matlab_2d_array_norm

PURPOSE:  simulate matlab mean function for 1 dimesion in 2d array cases only

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
    float **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float  *output_norm  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < array_len1; i++)
    {
        sum += array[i][dim2_number] * array[i][dim2_number];
    }
    *output_norm = sqrt(sum);
}

/******************************************************************************
MODULE:  matlab_2d_array_mean

PURPOSE:  simulate matlab mean function for 1 dimesion in 2d array cases only

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
    int dim2_number,     /* I: second dimension number used */   
    int array_len1,      /* I: number of input elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < array_len1; i++)
    {
        sum += array[i][dim2_number];
    }
    *output_mean = sum / (float)(array_len1);
}

/******************************************************************************
MODULE:  matlab_int_2d_partial_mean

PURPOSE:  simulate matlab mean function for parital part of 1 dimesion in 2d 
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
void matlab_int_2d_partial_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = start; i <= end; i++)
    {
        sum += (float)array[i][dim2_number];
    }
    *output_mean = sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  matlab_2d_partial_mean

PURPOSE:  simulate matlab mean function for parital part of 1 dimesion in 2d 
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
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    float sum = 0.0;

    for (i = start; i <= end; i++)
    {
        sum += array[i][dim2_number];
    }
    *output_mean = sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  matlab_2d_int_partial_mean

PURPOSE:  simulate matlab mean function for parital part of 1 dimesion in 2d 
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
void matlab_2d_int_partial_mean
(
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    int sum = 0;

    for (i = start; i <= end; i++)
    {
        sum += array[i][dim2_number];
    }
    *output_mean = (float)sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  matlab_2d_partial_square_mean

PURPOSE:  simulate matlab square mean function for parital part of 1 dimesion 
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
    int16 **array,       /* I: input array */
    int dim2_number,     /* I: second dimension number used */   
    int start,           /* I: number of start elements in 1st dim */
    int end,             /* I: number of end elements in 1st dim */
    float  *output_mean  /* O: output norm value */
)
{
    int i;
    int sum = 0;

    for (i = start; i <= end; i++)
    {
        sum += array[i][dim2_number] * array[i][dim2_number];
    }
    *output_mean = (float)sum / (float)(end - start + 1);
}

/******************************************************************************
MODULE:  get_ids_length

PURPOSE:  get total number of non-zero elements of id array

RETURN VALUE:
Type = int *
Value           Description
-----           -----------
total number of non-zero elements

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
            length++;
    }

    *id_len = length;
}

/******************************************************************************
MODULE:  get_array_length

PURPOSE:  get total number of non-zero elements of an array

RETURN VALUE:
Type = int *
Value           Description
-----           -----------
total number of non-zero elements

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/11/2015   Song Guo         Original Development

NOTES: 
******************************************************************************/
void get_array_length
(
    int *array,           /* I: input array */
    int array_len,        /* I: number of input elements in 1st dim */
    int *id_len           /* O: number of non-zero number in the array */
)
{
    int i;
    static int sum = 0;

    for (i = 0; i < array_len; i++)
    {
        if (array[i] != 0)
           sum ++;
    }

    *id_len = sum;
}

/******************************************************************************
MODULE:  partition_index

PURPOSE:  partition the sorted list with index

RETURN VALUE:
Type = int
Value           Description
-----           -----------
i               partitioned value   

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/14/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int partition_index (float arr[], int idx[], int left, int right)
{
    int i = left, j = right;
    float tmp;
    int index;
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
            index = idx[i];
            arr[i] = arr[j];
            idx[i] = idx[j];
            arr[j] = tmp;
            idx[j] = index;
            i++;
            j--;
        }
    }

    return i;
}

/******************************************************************************
MODULE:  quick_sort_index

PURPOSE:  sorted the array and return its index

RETURN VALUE:
Type = void

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
2/14/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void quick_sort_index (float arr[], int idx[], int left, int right)
{
    int index = partition_index (arr, idx, left, right);

    if (left < index - 1)
        quick_sort_index (arr, idx, left, index - 1);
    if (index < right)
        quick_sort_index (arr, idx, index, right);
}

/******************************************************************************
MODULE:  auto_robust_fit

PURPOSE:  Robust fit for one band

RETURN VALUE:
Type = int

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int dofit(const gsl_multifit_robust_type *T,
      const gsl_matrix *X, const gsl_vector *y,
      gsl_vector *c, gsl_matrix *cov)
{
  int s;
  gsl_multifit_robust_workspace * work 
    = gsl_multifit_robust_alloc (T, X->size1, X->size2);

  s = gsl_multifit_robust (X, y, c, cov, work);
  gsl_multifit_robust_free (work);

  return s;
}

int auto_robust_fit
(
    float **clrx,
    int16 **clry,
    int nums,
    int start,
    int band_index,
    float *coefs
)
{
    int i, j;
    const size_t p = 5; /* linear fit */
    gsl_matrix *x, *cov;
    gsl_vector *y, *c;

    /* Defines the inputs/outputs for robust fitting */
    x = gsl_matrix_alloc (nums, p);
    y = gsl_vector_alloc (nums);

    c = gsl_vector_alloc (p);
    cov = gsl_matrix_alloc (p, p);

    /* construct design matrix x for linear fit */
    for (i = 0; i < nums; ++i)
    {
        for (j = 0; j < p; j++)
        {
            if (j == 0)
                gsl_matrix_set (x, i, j, 1.0);
            else
                gsl_matrix_set (x, i, j, clrx[j-1][i]);
        }
        gsl_vector_set(y,i,clry[i+start][band_index]);
    }

    /* perform robust fit */
    dofit(gsl_multifit_robust_bisquare, x, y, c, cov);

    for (j = 0; j < c->size; j++)
    {
        coefs[j] = gsl_vector_get(c, j);
    }

    /* Free the memories */
    gsl_matrix_free (x);
    gsl_vector_free (y);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return SUCCESS;
}

/******************************************************************************
MODULE:  auto_mask

PURPOSE:  Multitemporal cloud, cloud shadow, & snow masks (global version)

RETURN VALUE:
Type = int

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int auto_mask
(
    int *clrx,
    int16 **clry,
    int start,
    int end,
    float years,
    int t_b1,
    int t_b2,
    int n_t,
    int *bl_ids
)
{
    char FUNC_NAME[] = "auto_mask";
    int year;
    float w, w2;
    float *coefs, *coefs2;
    int i;
    float **x;
    int status;
    float *pred_b2, *pred_b5;
    int nums;

    nums = end - start + 1;
    /* Allocate memory */
    x = (float **)allocate_2d_array(4, nums, sizeof(float));
    if (x == NULL)
        RETURN_ERROR("ERROR allocating x memory", FUNC_NAME, FAILURE);
    coefs = (float *)malloc(5 * sizeof(float));
    if (coefs == NULL)
        RETURN_ERROR("ERROR allocating coefs memory", FUNC_NAME, FAILURE);
    coefs2 = (float *)malloc(5 * sizeof(float));
    if (coefs2 == NULL)
        RETURN_ERROR("ERROR allocating coefs2 memory", FUNC_NAME, FAILURE);
    pred_b2 = (float *)malloc(nums * sizeof(float));
    if (pred_b2 == NULL)
        RETURN_ERROR("ERROR allocating pred_b2 memory", FUNC_NAME, FAILURE);
    pred_b5 = (float *)malloc(nums * sizeof(float));
    if (pred_b5 == NULL)
        RETURN_ERROR("ERROR allocating pred_b5 memory", FUNC_NAME, FAILURE);

    year = ceil(years);
    w = TWO_PI / 365.25;
    w2 = w / (float)year;

    for (i = 0; i < nums; i++)
    {
        x[0][i] = cos(w * (float)clrx[i+start]);
        x[1][i] = sin(w * (float)clrx[i+start]);
        x[2][i] = cos(w2 * (float)clrx[i+start]);
        x[3][i] = sin(w2 * (float)clrx[i+start]);
    }

    /* Do robust fitting for band 2 */
    status = auto_robust_fit(x, clry, nums, start, 1, coefs);

    /* Do robust fitting for band 5 */
    status = auto_robust_fit(x, clry, nums, start, 4, coefs2);

    /* predict band 2 * band 5 refs, bl_ids value of 0 is clear and 
       1 otherwise */
    for (i = 0; i < nums; i++)
    {
        pred_b2[i] = coefs[0] + coefs[1] * cos((float)clrx[i+start] * w ) + 
                  coefs[2] * sin((float)clrx[i+start] * w ) + coefs[3] * 
                  cos((float)clrx[i+start] * w2 ) 
                  + coefs[4] * sin((float)clrx[i+start] * w2);
        pred_b5[i] = coefs2[0] + coefs2[1] * cos((float)clrx[i+start] * w ) + 
                  coefs2[2] * sin((float)clrx[i+start] * w ) + coefs2[3] * 
                  cos((float)clrx[i+start] * w2 ) 
                  + coefs2[4] * sin((float)clrx[i+start] * w2);
        if ((((float)clry[i+start][1]-pred_b2[i]) > (float)(n_t * t_b1)) || 
            (((float)clry[i+start][4]-pred_b5[i]) < -(float)(n_t * t_b2)))
            bl_ids[i] = 1;
        else
            bl_ids[i] = 0;
    }

    /* Free allocated memory */
    free(coefs);
    free(coefs2);
    free(pred_b2);
    free(pred_b5);
    if (free_2d_array ((void **) x) != SUCCESS)
        RETURN_ERROR ("Freeing memory: x\n", FUNC_NAME, FAILURE);

    return SUCCESS;
}

/******************************************************************************
MODULE:  auto_ts_predict

PURPOSE:  Using lasso regression fitting coefficients to predict new values

RETURN VALUE:
Type = int

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
void auto_ts_predict
(
    int *clrx,
    float **coefs,
    int band_index,
    int start,
    int end,
    float *pred_y
)
{
    int i;
    int nums = end - start + 1;
    float w, w2, w3;
    w = TWO_PI / 365.25;
    w2 = 2.0 * w;
    w3 = 3.0 * w;

    for (i = 0; i < nums; i++)
    { 
        pred_y[i]  = coefs[0][band_index] + coefs[1][band_index] * 
              (float)clrx[i+start] + coefs[2][band_index] * 
              cos((float)clrx[i+start] * w ) + coefs[3][band_index] * 
              sin((float)clrx[i+start] * w ) + coefs[4][band_index] * 
              cos((float)clrx[i+start] * w2 ) + coefs[5][band_index] * 
              sin((float)clrx[i+start] * w2) + coefs[6][band_index] * 
              cos((float)clrx[i+start] * w3 ) +coefs[7][band_index] * 
              sin((float)clrx[i+start] * w3);
    }
}

/******************************************************************************
MODULE:  auto_ts_fit

PURPOSE:  Lasso regression fitting with full outputs

RETURN VALUE:
Type = int

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int auto_ts_fit
(
    int *clrx,
    int16 **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse
)
{
    char FUNC_NAME[] = "auto_ts_fit";
    float w;
    int i;
    float **x;
    int status;
    float **v_dif;
    float *yhat;
    float v_dif_norm;
    int nums;
    FILE *fd;

    nums = end - start + 1;
    w = TWO_PI / 365.25;
    /* Allocate memory */
    if (df ==2)
    {
        x = (float **)allocate_2d_array(1, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 4)
    {
        x = (float **)allocate_2d_array(3, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 6)
    {
        x = (float **)allocate_2d_array(5, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 8)
    {
        x = (float **)allocate_2d_array(7, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else 
        RETURN_ERROR("Unsupported df value", FUNC_NAME, FAILURE);

    v_dif = (float **)allocate_2d_array(end - start + 1, TOTAL_BANDS - 1, 
             sizeof(float));
    if (v_dif == NULL)
        RETURN_ERROR("Allocating v_dif memory", FUNC_NAME, FAILURE);

    yhat = (float *)malloc(nums * sizeof(float));
    if (yhat == NULL)
        RETURN_ERROR("Allocating yhat memory", FUNC_NAME, FAILURE);

    for (i = 0; i < nums; i++)
    {
        if (df == 2)
        {
            x[0][i] = (float)clrx[i+start];
        }
        else if (df == 4)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
        }
        else if (df == 6)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
            x[3][i] = cos(2.0 * w * (float)clrx[i+start]);
            x[4][i] = sin(2.0 * w * (float)clrx[i+start]);
        }
        else if (df == 8)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
            x[3][i] = cos(2.0 * w * (float)clrx[i+start]);
            x[4][i] = sin(2.0 * w * (float)clrx[i+start]);
            x[5][i] = cos(3.0 * w * (float)clrx[i+start]);
            x[6][i] = sin(3.0 * w * (float)clrx[i+start]);
        }
        else
            RETURN_ERROR("Unsupported df value", FUNC_NAME, FAILURE);
    }

    /* Save the inputs for lasso fitting */

    if (df == 2)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%d\n", x[0][i], clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }

        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df2.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f", &coefs[0][band_index], &coefs[1][band_index]);
        fclose(fd);
    }
    else if (df == 4)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }
        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df4.r");
        if (status !=SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index]);
        fclose(fd);
    }
    else if (df == 6)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         x[3][i], x[4][i], clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }

        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df6.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index], 
               &coefs[4][band_index], &coefs[5][band_index]);
        fclose(fd);
    }
    else if (df == 8)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%f,%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         x[3][i], x[4][i], x[5][i], x[6][i],
                         clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }
        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df8.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f %f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index], &coefs[4][band_index], 
               &coefs[5][band_index], &coefs[6][band_index], &coefs[7][band_index]);
        fclose(fd);
    }
    else
         RETURN_ERROR("Incorrect df value", FUNC_NAME, FAILURE);

    /* predict lasso model results */
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        auto_ts_predict(clrx, coefs, band_index, start, end, yhat);
        for (i = 0; i < nums; i++)
        {
            v_dif[i][band_index] = (float)clry[i+start][band_index] - yhat[i];
        }
        matlab_2d_array_norm(v_dif, band_index, nums, &v_dif_norm);
        *rmse = v_dif_norm / sqrt((float)(nums - df));
    }

    /* Free allocated memory */
    free(yhat);
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        if (free_2d_array ((void **) x) != SUCCESS)
            RETURN_ERROR ("Freeing memory: x\n", FUNC_NAME, FAILURE);
    }
    if (free_2d_array ((void **) v_dif) != SUCCESS)
        RETURN_ERROR ("Freeing memory: v_dif\n", FUNC_NAME, FAILURE);

    /* Remove the temporary file */
    status = system("rm glmnet_fit_inputs.txt");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_inputs.txt file", FUNC_NAME, FAILURE);
    status = system("rm glmnet_fit_outputs.txt");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_outputs.txt file", FUNC_NAME, FAILURE);
    status = system("rm glmnet_fit_*.r.Rout");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_*.r.Rout file", FUNC_NAME, FAILURE);

    return SUCCESS;
}

/******************************************************************************
MODULE:  auto_ts_fit_full

PURPOSE:  Lasso regression fitting with full outputs

RETURN VALUE:
Type = int

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/5/2015   Song Guo         Original Development

NOTES:
******************************************************************************/
int auto_ts_fit_full
(
    int *clrx,
    int16 **clry,
    int band_index,
    int start,
    int end,
    int df,
    float **coefs,
    float *rmse,
    float **v_dif
)
{
    char FUNC_NAME[] = "auto_ts_fit_full";
    float w;
    int i;
    float **x;
    int status;
    float *yhat;
    float v_dif_norm;
    int nums;
    FILE *fd;

    nums = end - start + 1;
    w = TWO_PI / 365.25;

    /* Allocate memory */
    if (df ==2)
    {
        x = (float **)allocate_2d_array(1, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 4)
    {
        x = (float **)allocate_2d_array(3, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 6)
    {
        x = (float **)allocate_2d_array(5, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else if (df == 8)
    {
        x = (float **)allocate_2d_array(7, nums, sizeof(float));
        if (x == NULL)
            RETURN_ERROR("Allocating x memory", FUNC_NAME, FAILURE);
    }
    else 
        RETURN_ERROR("Unsupported df value", FUNC_NAME, FAILURE);

    yhat = (float *)malloc(nums * sizeof(float));
    if (yhat == NULL)
        RETURN_ERROR("Allocating yhat memory", FUNC_NAME, FAILURE);

    for (i = 0; i < nums; i++)
    {
        if (df == 2)
        {
            x[0][i] = (float)clrx[i+start];
        }
        else if (df == 4)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
        }
        else if (df == 6)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
            x[3][i] = cos(2.0 * w * (float)clrx[i+start]);
            x[4][i] = sin(2.0 * w * (float)clrx[i+start]);
        }
        else if (df == 8)
        {
            x[0][i] = (float)clrx[i+start];
            x[1][i] = cos(w * (float)clrx[i+start]);
            x[2][i] = sin(w * (float)clrx[i+start]);
            x[3][i] = cos(2.0 * w * (float)clrx[i+start]);
            x[4][i] = sin(2.0 * w * (float)clrx[i+start]);
            x[5][i] = cos(3.0 * w * (float)clrx[i+start]);
            x[6][i] = sin(3.0 * w * (float)clrx[i+start]);
        }
        else
            RETURN_ERROR("Unsupported df value", FUNC_NAME, FAILURE);
    }

    /* Save the inputs for lasso fitting */

    if (df == 2)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%d\n", x[0][i], clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }

        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df2.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f", &coefs[0][band_index], &coefs[1][band_index]);
        fclose(fd);
    }
    else if (df == 4)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }
#if 0
            printf("i,start,x[0][i], x[1][i], x[2][i],clry[i+start][band_index]=%d,%d,%f,%f,%f,%d\n",
                   i,start,x[0][i], x[1][i], x[2][i],clry[i+start][band_index]);
#endif
        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df4.r");
        if (status !=SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index]);
        fclose(fd);
    }
    else if (df == 6)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         x[3][i], x[4][i], clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }

        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df6.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index], 
               &coefs[4][band_index], &coefs[5][band_index]);
        fclose(fd);
    }
    else if (df == 8)
    {
        fd = fopen("glmnet_fit_inputs.txt", "w");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file1", FUNC_NAME, FAILURE);

        for (i = 0; i < nums; i++)
        {
            if (fprintf (fd, "%f,%f,%f,%f,%f,%f,%f,%d\n", x[0][i], x[1][i], x[2][i], 
                         x[3][i], x[4][i], x[5][i], x[6][i],
                         clry[i+start][band_index]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before nums"
                              " lines", FUNC_NAME, FAILURE);
            }
        }
        fclose(fd);

        /* Call R script to do lasso fitting */
        status = system("R CMD BATCH $BIN/glmnet_fit_df8.r");
        if (status != SUCCESS)
            RETURN_ERROR ("Running glmnet fit R scripts", FUNC_NAME, FAILURE);

        fd = fopen("glmnet_fit_outputs.txt", "r");
        if (fd == NULL)
            RETURN_ERROR("ERROR opening temporary file2", FUNC_NAME, FAILURE);

        /* Read out the lasso fit coefficients */
        fscanf(fd, "%f %f %f %f %f %f %f %f", &coefs[0][band_index], &coefs[1][band_index], 
               &coefs[2][band_index], &coefs[3][band_index], &coefs[4][band_index], 
               &coefs[5][band_index], &coefs[6][band_index], &coefs[7][band_index]);
        fclose(fd);
    }
    else
         RETURN_ERROR("Incorrect df value", FUNC_NAME, FAILURE);

    /* predict lasso model results */
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        auto_ts_predict(clrx, coefs, band_index, start, end, yhat);
        for (i = 0; i < nums; i++)
        {
            v_dif[i][band_index] = (float)clry[i+start][band_index] - yhat[i];
        }
        matlab_2d_array_norm(v_dif, band_index, nums, &v_dif_norm);
        *rmse = v_dif_norm / sqrt((float)(nums - df));
    }

    /* Free allocated memory */
    free(yhat);
    if (df == 2 || df == 4 || df == 6 || df == 8)
    {
        if (free_2d_array ((void **) x) != SUCCESS)
            RETURN_ERROR ("Freeing memory: x\n", FUNC_NAME, FAILURE);
    }

    /* Remove the temporary file */
    status = system("rm glmnet_fit_inputs.txt");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_inputs.txt file", FUNC_NAME, FAILURE);
    status = system("rm glmnet_fit_outputs.txt");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_outputs.txt file", FUNC_NAME, FAILURE);
    status = system("rm glmnet_fit_*.r.Rout");
    if (status != SUCCESS)
        RETURN_ERROR ("Deleting glmnet_fit_*.r.Rout file", FUNC_NAME, FAILURE);

    return SUCCESS;
}

