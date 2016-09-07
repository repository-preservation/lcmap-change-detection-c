/*****************************************************************************
!File: input.c
*****************************************************************************/

#include "input.h"
#include "utilities.h"
#include "defines.h"

//#define TOTAL_IMAGE_BANDS 7

const char raw_binary_format[][4] = {"rb", "wb", "rb+"};

FILE *open_raw_binary
(
    char *infile,        /* I: name of the input file to be opened */
    char *access_type    /* I: string for the access type for reading the
                               input file; use the raw_binary_format
                               array at the top of this file */
)
{
    FILE *rb_fptr = NULL;    /* pointer to the raw binary file */
    char FUNC_NAME[] = "open_raw_binary"; /* function name */

    /* Open the file with the specified access type */
    rb_fptr = fopen (infile, access_type);
    if (rb_fptr == NULL)
    {
         ERROR_MESSAGE("Opening raw binary", FUNC_NAME);
	 return NULL;
    }

    /* Return the file pointer */
    return rb_fptr;
}

void close_raw_binary
(
    FILE *fptr      /* I: pointer to raw binary file to be closed */
)
{
    fclose (fptr);
}

int write_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to write to the file */
    int nsamps,         /* I: number of samples to write to the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    void *img_array     /* I: array of nlines * nsamps * size to be written
                              to the raw binary file */
)
{
    int nvals;               /* number of values written to the file */
    char FUNC_NAME[] = "write_raw_binary"; /* function name */

    /* Write the data to the raw binary file */
    nvals = fwrite (img_array, size, nlines * nsamps, rb_fptr);
    if (nvals != nlines * nsamps)
    {
        RETURN_ERROR("Incorrect amount of data written", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}

int read_raw_binary
(
    FILE *rb_fptr,      /* I: pointer to the raw binary file */
    int nlines,         /* I: number of lines to read from the file */
    int nsamps,         /* I: number of samples to read from the file */
    int size,           /* I: number of bytes per pixel (ex. sizeof(uint8)) */
    void *img_array     /* O: array of nlines * nsamps * size to be read from
                              the raw binary file (sufficient space should
                              already have been allocated) */
)
{
    int nvals;               /* number of values read from the file */
    char FUNC_NAME[] = "read_raw_binary"; /* function name */

    /* Read the data from the raw binary file */
    nvals = fread (img_array, size, nlines * nsamps, rb_fptr);
    if (nvals != nlines * nsamps)
    {
        RETURN_ERROR("Incorrect amount of data read", FUNC_NAME, ERROR);
    }

    return (SUCCESS);
}

/******************************************************************************
MODULE: trimwhitespace

PURPOSE: Trim leading spaces of a sting
 
RETURN VALUE:
Type = string without trailing space

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2015    Song Guo         Modified from online code

NOTES:
*****************************************************************************/
char *trimwhitespace(char *str)
{
  char *end;

  /* Trim leading space */
  while(isspace(*str)) str++;

  if(*str == 0)  
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  /* Write new null terminator */
  *(end+1) = 0;

  return str;
}


/******************************************************************************
MODULE: read_envi_header

PURPOSE: Reads envi header info into input structure
 
RETURN VALUE:
Type = None

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
1/16/2015    Song Guo         Original development
5/03/2016    Brian Davis      extracted all of the name creation from main.c
                              to here.

NOTES:
*****************************************************************************/

int read_envi_header
(
    char *data_type,       /* I: input data type        */
    char *scene_name,      /* I: scene name             */
    Input_meta_t *meta     /* O: saved header file info */
)
{
    char  buffer[MAX_STR_LEN] = "\0"; /* for retrieving fields        */
    char  *label = NULL;              /* for retrieving string tokens */
    char  *tokenptr = NULL;           /* for retrieving string tokens */
    char  *tokenptr2 = NULL;          /* for retrieving string tokens */
    char  *seperator = "=";           /* for retrieving string tokens */
    char  *seperator2 = ",";          /* for retrieving string tokens */
    FILE *in;                         /* file ptr to input hdr file   */
    int ib;                           /* loop index                   */
    char map_info[10][MAX_STR_LEN];   /* projection information fields*/
    char FUNC_NAME[] = "read_envi_header"; /* function name           */
    char filename[MAX_STR_LEN];       /* scene name                   */
    int len;                          /* for strlen                   */
    char short_scene[MAX_STR_LEN];    /* char string for text manipulation */
    char directory[MAX_STR_LEN];      /* for constucting path/file names */
    char tmpstr[MAX_STR_LEN];         /* char string for text manipulation */
    char scene_list_name[MAX_STR_LEN];/* char string for text manipulation */
    int landsat_number;               /* mission number defines file name */


    /******************************************************************/
    /*                                                                */
    /* Determine the file name.                                       */
    /*                                                                */
    /******************************************************************/

    if (strcmp(data_type, "tifs") == 0)
    {
        len = strlen(scene_name);
        landsat_number = atoi(sub_string(scene_name,(len-19),1));
        if (landsat_number == 8)
            sprintf(filename, "%s_sr_band2.hdr", scene_name);
        else
            sprintf(filename, "%s_sr_band1.hdr", scene_name);
    }
    else if (strcmp(data_type, "bip") == 0)
    {
        len = strlen(scene_name);
        strncpy(short_scene, scene_name, len-5);
        split_directory_scenename(scene_name, directory, scene_list_name);
        if (strncmp(short_scene, ".", 1) == 0)
        {
            strncpy(tmpstr, short_scene + 2, len - 2);
            sprintf(filename, "%s/%s_MTLstack.hdr", tmpstr, scene_list_name);
        }
        else
            sprintf(filename, "%s/%s_MTLstack.hdr", short_scene, scene_list_name);
    }

    in=fopen(filename, "r");
    if (in == NULL)
    {
        RETURN_ERROR ("opening header file", FUNC_NAME, FAILURE);
    }

    /* process line by line */
    while(fgets(buffer, MAX_STR_LEN, in) != NULL) 
    {

        char *s;
        s = strchr(buffer, '=');
        if (s != NULL)
        {
            /* get string token */
            tokenptr = strtok(buffer, seperator);
            label=trimwhitespace(tokenptr);

            if (strcmp(label,"lines") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->lines = atoi(tokenptr);
            }

            if (strcmp(label,"data type") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->data_type = atoi(tokenptr);
            }

            if (strcmp(label,"byte order") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->byte_order = atoi(tokenptr);
            }

            if (strcmp(label,"samples") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                meta->samples = atoi(tokenptr);
            }

            if (strcmp(label,"interleave") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
                strcpy(meta->interleave, tokenptr);
            }

            if (strcmp(label,"UPPER_LEFT_CORNER") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
            }

            if (strcmp(label,"map info") == 0)
            {
                tokenptr = trimwhitespace(strtok(NULL, seperator));
            }

            if (strcmp(label,"map info") == 0)
            {
                tokenptr2 = strtok(tokenptr, seperator2);
                ib = 0;
                while(tokenptr2 != NULL)
                {
                    strcpy(map_info[ib], tokenptr2);
                    if (ib == 3)
                        meta->upper_left_x = atoi(map_info[ib]);
                    if (ib == 4)
                        meta->upper_left_y = atoi(map_info[ib]);
                    if (ib == 5)
                        meta->pixel_size = atoi(map_info[ib]);
                    if(ib == 7)
                        meta->utm_zone = atoi(map_info[ib]);
                    tokenptr2 = strtok(NULL, seperator2);
                    ib++;
                }
            }
        }
    }
    fclose(in);

    return (SUCCESS);
}


/*******************************************************************************
MODULE: read_stdin

PURPOSE: Reads pixels values from stdin and assings to array/pointer.
         Part of an on-going effort to remove all read I/O from main.c.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160421     Brian Davis      Original development

NOTES:
*******************************************************************************/

int read_stdin
(
    int           *updated_sdate_array, /* I/O: pointer to date values buffer */
    int           *buf,                 /* I/O: pointer to image bands buffer */
    unsigned char *updated_cfmask_buf,  /* I/O: pointer to cfmask pixel buffer*/
    int           num_bands,            /* I:   total number of bands         */
    int           *clear_sum,           /* O:   accumulator for clear pixels  */
    int           *water_sum,           /* O:   accumulator for clear pixels  */
    int           *shadow_sum,          /* O:   accumulator for clear pixels  */
    int           *snow_sum,            /* O:   accumulator for clear pixels  */
    int           *cloud_sum,           /* O:   accumulator for clear pixels  */
    int           *fill_sum,            /* O:   accumulator for clear pixels  */
    int           *all_sum,             /* O:   accumulator for clear pixels  */
    int           *valid_num_scenes,    /* O:   total scenes read             */
    bool          debug                 /* I:   flag for printing mesgs/info  */
)

{
    bool end_of_file = 0;               /* To identify end of stdin.          */
    int i, j;                           /* loop counters.                     */
    int status;                         /* function return status.            */

    /**************************************************************/
    /*                                                            */
    /* For stdin:                                                 */
    /* this assumes order of: julian date value, then 6 SR and 1  */
    /* thermal band values, then cfmask band values, each         */
    /* set/group together, each element of the group separated by */
    /* white space, culiminated with a newline per scene, for     */
    /* number of scenes. For example:                             */
    /* 2456445 94 156 164 758 807 492 2809 0                      */
    /*                                                            */
    /**************************************************************/

    i = 0;
    while (!end_of_file)
    {
        /**************************************************************/
        /*                                                            */
        /* Read Julian date/time value.                               */
        /*                                                            */
        /**************************************************************/

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

        /**************************************************************/
        /*                                                            */
        /* Read surface reflectance and thermal band values.          */
        /*                                                            */
        /**************************************************************/

        for (j = 0; j < num_bands; j++)
        {
            scanf("%d", &buf[i + j]);  // offsets 
            if (debug)
                //printf( "You entered: %d\n", buf[j][i]);
                printf( "You entered: %d\n", buf[i + j]);  // offsets 
        }

        /**************************************************************/
        /*                                                            */
        /* Read cfmask band value.                                    */
        /*                                                            */
        /**************************************************************/

        scanf("%hhu", &updated_cfmask_buf[i]);
        if (debug)
            printf( "You entered: %u\n", updated_cfmask_buf[i]);

        /**************************************************************/
        /*                                                            */
        /* Call the function with the case statement for cfmask       */
        /* values, because it is used for all input type optoins,     */
        /* and we put it in a function so we did not have to make     */
        /* multiple updates each time something changed.              */
        /*                                                            */
        /**************************************************************/

        status = assign_cfmask_values (updated_cfmask_buf[i], clear_sum,
                                       water_sum, shadow_sum, snow_sum,
                                       cloud_sum, fill_sum, all_sum);
        if (status != SUCCESS)
        {
            printf ("Error calling assign_cfmask_values.");
            return (FAILURE);
        }

        i++;
        (*all_sum)++;
    }

    *valid_num_scenes = i;
    if (debug)
        printf ("\n");

    return (SUCCESS);

}


/*******************************************************************************
MODULE: assign_cfmask_values

PURPOSE: updates pixel type accumulators for cfmask values.  This logic is
         needed in several places so separated out into its own function, to
         eliminate missing future updates in some places.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160421     Brian Davis      Original development

NOTES:
*******************************************************************************/

int assign_cfmask_values
(
    unsigned char cfmask_value,         /* I: current cfmask pixel value.   */
    int           *clear_sum,           /* O: accumulator for clear  pixels */
    int           *water_sum,           /* O: accumulator for water  pixels */
    int           *shadow_sum,          /* O: accumulator for shadow pixels */
    int           *snow_sum,            /* O: accumulator for snow   pixels */
    int           *cloud_sum,           /* O: accumulator for cloud  pixels */
    int           *fill_sum,            /* O: accumulator for full   pixels */
    int           *all_sum              /* O: accumulator for all    pixels */
)

{

    /******************************************************************/
    /*                                                                */
    /* Add the current cfmask pixel value to the correct category,    */
    /* for reporting totals at the end of reading all input.          */
    /*                                                                */
    /******************************************************************/

    switch (cfmask_value)
    {
        case CFMASK_CLEAR:
            (*clear_sum)++;
            break;
        case CFMASK_WATER:
            (*water_sum)++;
            (*clear_sum)++;
            break;
        case CFMASK_SHADOW:
            (*shadow_sum)++;
            break;
        case CFMASK_SNOW:
            (*snow_sum)++;
            break;
        case CFMASK_CLOUD:
            (*cloud_sum)++;
            break;
        case CFMASK_FILL:
            (*fill_sum)++;
            break;
        default:
            printf ("Unknown cfmask value %d", cfmask_value);
            return (FAILURE);
            break;
        }

    if (cfmask_value < CFMASK_FILL)
        (*all_sum)++;

    return (SUCCESS);
}


/*******************************************************************************
MODULE: read_tifs

PURPOSE: Creates the image band file names to open and read, and fills image
         buffer with values read. For data-type=tifs, all bands are in separate
         tif image files. Part of an on-going effort to remove all read I/O
         from main.c.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160425     Brian Davis      Original development
                              Fixed a bug from original code extracted
                              from main incorrectly checking for 7th
                              (cfmask) band in the image band loop.

NOTES:
*******************************************************************************/

int read_tifs
(
    char *sceneID_name,  /* I:   current file name in list of sceneIDs  */
    FILE ***fp_tifs,     /* I/O: file pointer array for band file names */
    int  curr_scene_num, /* I:   current num. in list of scenes to read */
    int  row,            /* I:   the row (Y) location within img/grid   */
    int  col,            /* I:   the col (X) location within img/grid   */
    int  nrows,          /* I:   number of rows to read in tile of data */
    int  ncols,          /* I:   number of colums to read (1 for bip)   */
    int  num_samples,    /* I:   number of image samples (X width)      */
    bool debug,          /* I:   flag for printing debug messages       */
    int  *image_buf      /* O:   pointer to 2-D image band values array */
)

{
    int  k;                     /* band loop counter.                   */
    int  row_inx, col_inx;      /* 2d array loop counters               */
    int  len;                   /* for string length call.              */
    int  landsat_number;        /* numeric mission number to make names */
    char filename[MAX_STR_LEN]; /* file name constructed from sceneID   */
    int  status;                /* return status of system call(s)      */
    int  offset;                /* pixel locatoffset into image buffer  */
    int  row_size;              /* size of a row to calculate offset    */
    int  col_size;              /* size of a col to calculate offset    */


    row_size = num_samples;
    col_size = TOTAL_BANDS;

    /******************************************************************/
    /*                                                                */
    /* Read the image bands for this scene. using the                 */
    /*                                                                */
    /******************************************************************/

    for (k = 0; k < TOTAL_IMAGE_BANDS; k++)
    {
        len = strlen(sceneID_name);
        landsat_number = atoi(sub_string(sceneID_name,(len-19),1));
        if (landsat_number != 8)
        {
            if (k == 5)
                sprintf(filename, "%s_sr_band%d.img", sceneID_name, k+2);
            else if (k == 6)
                sprintf(filename, "%s_toa_band6.img", sceneID_name);
            else
                sprintf(filename, "%s_sr_band%d.img", sceneID_name, k+1);
        }
        else
        {
            if (k == 6)
                sprintf(filename, "%s_toa_band10.img", sceneID_name);
            else 
                sprintf(filename, "%s_sr_band%d.img",  sceneID_name, k+2);
        }

        fp_tifs[k][curr_scene_num] = open_raw_binary(filename,"rb");
        if (fp_tifs[k][curr_scene_num] == NULL)
            printf("error open %d scene, %d bands files\n",curr_scene_num, k+1);

        status = fseek(fp_tifs[k][curr_scene_num], ((row * num_samples) + col) * sizeof(short int), SEEK_SET);
        if (status != 0)
            printf("error seeking %d scene, %d bands\n", curr_scene_num, (k + 1));

        for (row_inx = 0; row_inx < nrows; row_inx++)
        {
            for (col_inx = 0; col_inx < ncols; col_inx++)
            {
                offset = (row_inx * row_size) + (col_inx * col_size) + curr_scene_num + k;
                //if (read_raw_binary(fp_tifs[k][curr_scene_num], 1, 1,
                //                    sizeof(short int),
                //                    &image_buf[k][curr_scene_num]) != 0)
                if (read_raw_binary(fp_tifs[k][curr_scene_num], 1, 1, sizeof(short int), image_buf[offset]) != 0)
                    printf("error reading %d scene, %d bands\n", curr_scene_num, (k + 1));
                if (debug)
                {
                    printf("%d ", (short int)image_buf[offset]);
                }
            }
        }
    
        close_raw_binary(fp_tifs[k][curr_scene_num]);


    }

    return (SUCCESS);
}



/*******************************************************************************
MODULE: read_bip

PURPOSE: Creates the image band file names to open and read, and fills image
         buffer with values read. For data-type=bip, all image bands are in a
         single envi bip format image file. Part of an on-going effort to
         remove all read I/O from main.c.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160428     Brian Davis      Original development

NOTES:
*******************************************************************************/

int read_bip
(
    char *curr_scene_name,    /* I:   current file name in list of sceneIDs  */
    FILE **fp_bip,            /* I/O: file pointer array for BIP  file names */
    int  curr_scene_num,      /* I:   current num. in list of scenes to read */
    int  row,                 /* I:   the row (Y) location within img/grid   */
    int  col,                 /* I:   the col (X) location within img/grid   */
    int  nrows,               /* I:   number of rows to read in tile of data */
    int  ncols,               /* I:   number of colums to read (1 for bip)   */
    int  num_samples,         /* I:   number of image samples (X width)      */
    int  *image_buf           /* O:   pointer to 2-D image band values array */
)

{

    int  row_inx, col_inx, band_inx; /* for relative to start numbers.  */
    int  row_count, col_count; /* row, column, band loop counters.      */
    int  len;                   /* for string length call.              */
    char filename[MAX_STR_LEN]; /* file name constructed from sceneID   */
    char shorter_name[MAX_STR_LEN];/* file name constructed from sceneID*/
    char directory[MAX_STR_LEN];
    char scene_name[MAX_STR_LEN];
    char tmpstr[MAX_STR_LEN];   /* for string manipulation              */
    char errmsg[MAX_STR_LEN];   /* for printing error text to the log.  */
    bool debug = true;          /* for debug printing                   */
    int  in_offset;             /* pixel location offset into input buf */
    int  out_offset;            /* pixel location offset into output buf*/


    /******************************************************************/
    /*                                                                */
    /* Determine the BIP file name, open, fseek.                      */
    /*                                                                */
    /******************************************************************/

    len = strlen(curr_scene_name);
    strncpy(shorter_name, curr_scene_name, len-5);
    // somebody is trashing memory and I cannot find him...
    shorter_name[len-5] = '\0';

    split_directory_scenename(curr_scene_name, directory, scene_name);

    if (strncmp(shorter_name, ".", 1) == 0)
    {
        strncpy(tmpstr, shorter_name + 2, len - 2);
        sprintf(filename, "%s/%s_MTLstack", tmpstr, scene_name);
    }
    else
        sprintf(filename, "%s/%s_MTLstack", shorter_name, scene_name);

    fp_bip[curr_scene_num] = open_raw_binary(filename,"rb");
    if (fp_bip[curr_scene_num] == NULL)
    {
        sprintf(errmsg, "Opening %d scene files\n", curr_scene_num);
        printf(errmsg);
        return (FAILURE);
    }

    in_offset = ((row - 1) * num_samples + col - 1) * TOTAL_BANDS * sizeof(short int);
    fseek(fp_bip[curr_scene_num], in_offset, SEEK_SET);

    /******************************************************************/
    /*                                                                */
    /* Read the image bands for this scene.                           */
    /*                                                                */
    /******************************************************************/

    for (row_inx = (row - 1), row_count = 0; row_count < nrows; row_count++, row_inx++)
    {
        for (col_inx = (col - 1), col_count = 0; col_count < ncols; col_count++, col_inx++)
        {
            for (band_inx = 0; band_inx < TOTAL_IMAGE_BANDS; band_inx++)
            {
                out_offset = (row_count * num_samples) + (col_count * TOTAL_BANDS) + curr_scene_num + band_inx;
                if (read_raw_binary(fp_bip[curr_scene_num], 1, 1, sizeof(short int),
                                    &image_buf[out_offset]) != 0)
                    //sizeof(short int), &image_buf[(col_inx * TOTAL_BANDS) + band_inx][curr_scene_num]) != 0)
                {
        	    sprintf(errmsg, "error reading %d scene, row %d, col %d, %d bands\n",
                            curr_scene_num, row_inx + 1, col_inx + 1, band_inx + 1);
                    printf(errmsg);
                    return (FAILURE);
                }
                if (debug)
                    printf("%d ", (short int)image_buf[out_offset]);
            }
        }
    }

    close_raw_binary(fp_bip[curr_scene_num]);

    return (SUCCESS);
}


/*******************************************************************************
MODULE: read_cfmask

PURPOSE: Creates the cfmask file names to open and read, and fills cfmask
         buffer with values read. For data-type=bip, all image bands are in a
         single envi bip format image file. For tif files, all bands, including
         cfmask, are in separate files.  Part of an on-going effort to
         remove all read I/O from main.c.
 
RETURN VALUE:
Type = int

HISTORY:
Date         Programmer       Reason
----------   --------------   -------------------------------------
20160428     Brian Davis      Original development

NOTES:
*******************************************************************************/

int read_cfmask
(
    int  curr_scene_num, /* I:   current num. in list of scenes to read       */
    char *data_type,     /* I:   type of flies, tifs or single BIP            */
    char **scene_list,   /* I:   current scene name in list of sceneIDs       */
    int  row,            /* I:   the row (Y) location within img/grid         */
    int  col,            /* I:   the col (X) location within img/grid         */
    int  nrows,          /* I:   number of rows to read in a tile of some data*/
    int  ncols,          /* I:   number of cols to read in a line of bip data */
    int  num_samples,    /* I:   number of image samples (X width)            */
    FILE ***fp_tifs,     /* I/O: file ptr array for tif band file names       */
    FILE **fp_bip,       /* I/O: file pointer array for BIP file names        */
    unsigned char *fmask_buf,/* O:   pointer to cfmask band values            */
                         /* I/O: Worldwide Reference System path and row for  */
                         /* I/O: the current swath, this group of variables   */
                         /* I/O: is for filtering out swath overlap, and      */
    int *prev_wrs_path,  /* I/O: using the first of two scenes in a swath,    */
    int *prev_wrs_row,   /* I/O: , because it is recommended to use the meta  */
    int *prev_year,      /* I/O: data from the first for things like sun      */
    int *prev_jday,      /* I/O: angle, etc. However, always removing a       */
    unsigned char *prev_fmask_buf,/* I/O: redundant x/y location specified    */
    int *valid_scene_count,/* I/O: x/y is not always valid for gridded data,  */
    int *swath_overlap_count,/* I/O: it may/may not be in swap overlap area.  */
    char **valid_scene_list,/* I/O: 2-D array for list of filtered            */
    int *clear_sum,      /* I/O: Total number of clear cfmask pixels          */
    int *water_sum,      /* I/O: counter for cfmask water pixels.             */
    int *shadow_sum,     /* I/O: counter for cfmask shadow pixels.            */
    int *sn_sum,         /* I/O: Total number of snow cfmask pixels           */
    int *cloud_sum,      /* I/O: counter for cfmask cloud pixels.             */
    int *fill_sum,       /* I/O: counter for cfmask fill pixels.              */
    int *all_sum,        /* I/O: Total of all cfmask pixels                   */
    unsigned char *updated_fmask_buf, /* I/O: new entry in valid fmask values */
    int *updated_sdate_array, /* I/O: new buf of valid date values            */
    int *sdate,          /* I:   Original array of julian date values         */
    int *valid_num_scenes/* I/O: number of valid scenes after reading cfmask  */
)

{

    int len;             /* for strlen call                             */
    int landsat_number;  /* mission number for determining file names   */
    char scene_name[MAX_STR_LEN]; /* current scene id name              */
    char filename[MAX_STR_LEN];   /* temp for constructing file name    */
    int wrs_path;        /* Worldwide Reference System path             */
    int wrs_row = 0;     /* WRS row                                     */
    int year;            /* Year of acquisition date of current scene   */
    int jday;            /* Julian day since 0 of current scene date    */
    bool debug = 1;      /* for debug printing                          */
    int status;          /* for return status of function calls         */
    char short_scene[MAX_STR_LEN]; /* for parsing file names            */
    char directory[MAX_STR_LEN]; /* for parsing file names              */
    char tmpstr[MAX_STR_LEN]; /* for parsing file names                 */
    char errmsg[MAX_STR_LEN]; /* for printing errors before log/quit    */
    char FUNC_NAME[] = "read_cfmask"; /* for printing errors messages   */
    int int_buf;         /* for reading cfmask value then type cast     */
    int row_inx, col_inx;/* for reading 2d arrays of tiles              */

    if (strcmp(data_type, "tifs") == 0)
    {

        /**************************************************************/
        /*                                                            */
        /* Determine the cfmask file name to read.                    */
        /*                                                            */
        /**************************************************************/
    
        len = strlen(scene_list[curr_scene_num]);
        landsat_number = atoi(sub_string(scene_list[curr_scene_num],(len-19),1));
        wrs_path = atoi(sub_string(scene_list[curr_scene_num],(len-18),3));
        wrs_row =  atoi(sub_string(scene_list[curr_scene_num],(len-15),3));
        year = atoi(sub_string(scene_list[curr_scene_num],(len-12),4));
        jday = atoi(sub_string(scene_list[curr_scene_num],(len- 8),3));
        sprintf(filename, "%s_cfmask.img", scene_list[curr_scene_num]);
    
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
    
        fp_tifs[CFMASK_BAND][curr_scene_num] = open_raw_binary(filename,"rb");
        if (fp_tifs[CFMASK_BAND][curr_scene_num] == NULL)
            printf("error open %d scene, %d bands files\n", curr_scene_num, CFMASK_BAND+1);
    
        fseek(fp_tifs[CFMASK_BAND][curr_scene_num], ((row * num_samples) + col) *
              sizeof(unsigned char), SEEK_SET);
    
        for (row_inx = 0; row_inx < nrows; row_inx++)
        {
            for (col_inx = 0; col_inx < ncols; col_inx++)
            {
                if (read_raw_binary(fp_tifs[CFMASK_BAND][curr_scene_num], 1, 1,
                    sizeof(unsigned char), &fmask_buf[curr_scene_num]) != 0)
                    printf("error reading %d scene, %d bands\n", curr_scene_num, CFMASK_BAND+1);
            }
        }

        close_raw_binary(fp_tifs[CFMASK_BAND][curr_scene_num]);
    }

    else if (strcmp(data_type, "bip") == 0)

    {

        len = strlen(scene_list[curr_scene_num]);
        strncpy(short_scene, scene_list[curr_scene_num], len-5);
        split_directory_scenename(scene_list[curr_scene_num], directory, scene_name);

        len = strlen(scene_name);
        landsat_number = atoi(sub_string(scene_name,(len-19),1));
        wrs_path = atoi(sub_string(scene_name,(len-18),3));
        wrs_row =  atoi(sub_string(scene_name,(len-15),3));
        year = atoi(sub_string(scene_name,(len-12),4));
        jday = atoi(sub_string(scene_name,(len- 8),3));

        if (strncmp(short_scene, ".", 1) == 0)
        {
            strncpy(tmpstr, short_scene + 2, len - 2);
            sprintf(filename, "%s/%s_MTLstack", tmpstr, scene_name);
        }
        else
            sprintf(filename, "%s/%s_MTLstack", short_scene, scene_name);
        fp_bip[curr_scene_num] = open_raw_binary(filename,"rb");
        if (fp_bip[curr_scene_num] == NULL)
        {
            sprintf(errmsg, "Opening %d scene files\n", curr_scene_num);
            RETURN_ERROR (errmsg, FUNC_NAME, ERROR);
        }

        for (row_inx = 0; row_inx < nrows; row_inx++)
        {
            fseek(fp_bip[curr_scene_num], ((row - 1) * num_samples + col - 1) *
                  TOTAL_BANDS * sizeof(short int) + (TOTAL_IMAGE_BANDS * sizeof(short int)), SEEK_SET);
            for (col_inx = 0; col_inx < ncols; col_inx++)
            {
                if (read_raw_binary(fp_bip[curr_scene_num], 1, 1,
                        sizeof(short int), &int_buf) != 0)
                {
                    sprintf(errmsg, "error reading %d scene, %d bands\n",curr_scene_num, CFMASK_BAND+1);
                    RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
                }
                fmask_buf[curr_scene_num + ((row_inx * num_samples) + (col - 1) + col_inx)] = (unsigned char)int_buf;
            }
        }

        close_raw_binary(fp_bip[curr_scene_num]);

    }

    /******************************************************************/
    /*                                                                */
    /* Check for swath overlap pixels.  If consecutive temporal       */
    /* are in the same path, and in adjacent rows, are not fill, and  */
    /* have the same acquisition date, then they are essentially the  */
    /* same pixel, so use the first, because the metadata such as     */
    /* sun angle, etc. are more closely associated with the first.    */
    /*                                                                */
    /******************************************************************/

    if ((wrs_path == *prev_wrs_path) && (wrs_row == (*prev_wrs_row - 1)) &&
        (year == *prev_year) && (jday == *prev_jday) &&
        (fmask_buf[curr_scene_num] != CFMASK_FILL) && (*prev_fmask_buf != CFMASK_FILL))

    {
        (*swath_overlap_count)++;
        strcpy(valid_scene_list[(*valid_scene_count) - 1], scene_list[curr_scene_num]);
    }

    else
 
    {

        /**************************************************************/
        /*                                                            */
        /* Call the function with the case statement for totalling    */
        /* cfmask values, because it is used for all input type       */
        /* options, and we put it in a function so we did not have to */
        /* make multiple updates if something changed.                */
        /*                                                            */
        /**************************************************************/

        status = assign_cfmask_values (fmask_buf[curr_scene_num], clear_sum,
                                       water_sum, shadow_sum, sn_sum,
                                       cloud_sum, fill_sum, all_sum);
        if (status != SUCCESS)
        {
            RETURN_ERROR ("Calling assign_cfmask_values", 
                          "read_cfmask", FAILURE);
        }

        /**************************************************************/
        /*                                                            */
        /* After succesfully reading the cfmask data, update the list */
        /* of valid julian dates.                                     */
        /*                                                            */
        /**************************************************************/

        if (fmask_buf[curr_scene_num] < CFMASK_FILL)
        {
            updated_fmask_buf[*valid_scene_count] = fmask_buf[curr_scene_num];
            updated_sdate_array[*valid_scene_count] = sdate[curr_scene_num];
            strcpy(valid_scene_list[*valid_scene_count], scene_list[curr_scene_num]);
            (*valid_scene_count)++;
            (*valid_num_scenes)++;
        }
    }

    *prev_wrs_path = wrs_path;
    *prev_wrs_row  = wrs_row;
    *prev_year = year;
    *prev_jday = jday;
    *prev_fmask_buf = fmask_buf[curr_scene_num];

    return(SUCCESS);

}
