#include <string.h>
#include <stdarg.h>
#include <getopt.h>

#include "classification.h"
#include "utilities.h"

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
8/5/2015    Song Guo         Original Development
******************************************************************************/
int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    int *rows,                /* O: number of rows */
    int *cols,                /* O: number of columns */
    int *ref_rows,            /* O: number of rows for reference data */
    int *nclass,              /* O: number of classification types */
    bool *verbose             /* O: verbose flag */
)
{
    int c;                          /* current argument index */
    int option_index;               /* index for the command-line option */
    static int verbose_flag = 0;    /* verbose flag */
    static int cols_default = 71;   /* Default buffer for number of columns */
    static int nclass_default = 11; /* Default buffer for number of classes */
    char errmsg[MAX_STR_LEN];       /* error message */
    char FUNC_NAME[] = "get_args";  /* function name */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"rows", required_argument, 0, 'r'},
        {"cols", required_argument, 0, 'c'},
        {"ref_rows", required_argument, 0, 'f'},
        {"nclass", required_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Assign the default values */
    *cols = cols_default;
    *nclass = nclass_default;

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
                return FAILURE;
                break;

            case 'r':             
                *rows = atoi (optarg);
                break;

            case 'c':             
                *cols = atoi (optarg);
                break;

            case 'f':             
                *ref_rows = atoi (optarg);
                break;

            case 'n':             
                *nclass = atoi (optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
                break;
        }
    }

    /* Check the input values */
    if (*rows < 0)
    {
        sprintf (errmsg, "number of rows must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    if (*cols < 0)
    {
        sprintf (errmsg, "number of columns must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    if (*ref_rows < 0)
    {
        sprintf (errmsg, "number of reference rows must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    if (*nclass < 0)
    {
        sprintf (errmsg, "number of classes must be > 0");
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }


    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    if (*verbose)
    {
        printf ("rows = %d\n", *rows);
        printf ("cols = %d\n", *cols);
        printf ("ref_rows = %d\n", *ref_rows);
        printf ("nclass = %d\n", *nclass);
        printf ("verbose = %d\n", *verbose);
    }

    return SUCCESS;
}

