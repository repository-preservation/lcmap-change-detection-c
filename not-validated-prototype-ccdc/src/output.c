/*****************************************************************************
File: output.c

Description: Functions creating and writting data to the product output file.

Revision History:
 Song Guo
 Original Version - borrowed and modified from some of the LEDAPS libraries.

Team Unique Header:
  This software was developed by the Landsat Science, Research, and
  Development (LSRD) Team at the USGS EROS.

Design Notes:
 1. The following public functions handle the output files:
    CreateOutput - Create new output file.
    OutputFile - Setup 'output' data structure.
    CloseOutput - Close the output file.
    FreeOutput - Free the 'output' data structure memory.
    PutMetadata - Write the output product metadata.
    WriteOutput - Write a line of data to the output product file.
*****************************************************************************/

#include "utilities.h"
#include "output.h"


#define LST_BAND "lst"
#define LST_LONG_NAME "lst_band"


/******************************************************************************
Description: 'OutputFile' sets up the 'output' data structure and opens the
 output file for write access.

Input Parameters:this
 in_meta        input XML metadata structure (band-related info)
 input          input structure with input image metadata (nband, iband, size)

Output Parameters:
 (returns)      'output' data structure or NULL when an error occurs

HISTORY:
Date         Programmer       Reason
---------    ---------------  -------------------------------------
2/19/2014    Gail Schmidt     Modified to work with ESPA internal raw binary
                              file format

Design Notes:
1. MASK_INDEX "0 clear; 1 water; 2 cloud_shadow; 3 snow; 4 cloud"
*****************************************************************************/
Output_t *OpenOutput
(
    Espa_internal_meta_t *in_meta, /* I: input metadata structure */
    Input_t *input                 /* I: input reflectance band data */
)
{
    Output_t *this = NULL;
    char *mychar = NULL;        /* pointer to '_' */
    char scene_name[STR_SIZE];  /* scene name for the current scene */
    char file_name[STR_SIZE];   /* output filename */
    char production_date[MAX_DATE_LEN + 1]; /* current date/time for
                                               production */
    time_t tp;                  /* time structure */
    struct tm *tm = NULL;       /* time structure for UTC time */
    int indx = -1;              /* thermal band index in XML file */
    Espa_band_meta_t *bmeta = NULL;     /* pointer to the band metadata array
                                           within the output structure */

    /* Create the Output data structure */
    this = (Output_t *) malloc (sizeof (Output_t));
    if (this == NULL)
    {
        RETURN_ERROR ("allocating Output data structure", "OpenOutput", NULL);
    }

    /* Use thermal band for band metadata information */
    if (input->nband_th == 1)
        indx = 5;
    else
        indx = 9;

    /* Initialize the internal metadata for the output product. The global
       metadata won't be updated, however the band metadata will be updated
       and used later for appending to the original XML file. */
    init_metadata_struct (&this->metadata);

    /* Allocate memory for the output band */
    if (allocate_band_metadata (&this->metadata, 1) != SUCCESS)
    {
        RETURN_ERROR ("allocating band metadata", "OpenOutput", NULL);
    }

    bmeta = this->metadata.band;

    /* Determine the scene name */
    strcpy (scene_name, in_meta->band[indx].file_name);
    mychar = strchr (scene_name, '_');
    if (mychar != NULL)
        *mychar = '\0';

    /* Get the current date/time (UTC) for the production date of each band */
    if (time (&tp) == -1)
    {
        RETURN_ERROR ("unable to obtain current time", "OpenOutput", NULL);
    }

    tm = gmtime (&tp);
    if (tm == NULL)
    {
        RETURN_ERROR ("converting time to UTC", "OpenOutput", NULL);
    }

    if (strftime (production_date, MAX_DATE_LEN, "%Y-%m-%dT%H:%M:%SZ", tm)
        == 0)
    {
        RETURN_ERROR ("formatting the production date/time", "OpenOutput",
                      NULL);
    }

    /* Populate the data structure */
    this->open = false;
    this->fp_bin = NULL;
    this->nband = 1;
    this->size.l = input->size_th.l;
    this->size.s = input->size_th.s;

    strncpy (bmeta[0].short_name, in_meta->band[indx].short_name, 3);
    bmeta[0].short_name[3] = '\0';
    strcat (bmeta[0].short_name, "LST");
    strcpy (bmeta[0].product, "lst");
    strcpy (bmeta[0].source, "level1");
    strcpy (bmeta[0].category, "image");
    bmeta[0].nlines = this->size.l;
    bmeta[0].nsamps = this->size.s;
    bmeta[0].pixel_size[0] = input->meta.pixel_size[0];
    bmeta[0].pixel_size[1] = input->meta.pixel_size[1];
    strcpy (bmeta[0].pixel_units, "meters");
    sprintf (bmeta[0].app_version, "lst_%s", LST_VERSION);
    strcpy (bmeta[0].production_date, production_date);
    bmeta[0].data_type = ESPA_UINT16;
    bmeta[0].fill_value = FILL_VALUE;
    strcpy (bmeta[0].name, LST_BAND);
    strcpy (bmeta[0].long_name, LST_LONG_NAME);
    strcpy (bmeta[0].data_units, "quality/feature classification");
    bmeta[0].scale_factor = 10000.0;
    bmeta[0].add_offset = 0.0;

    /* Set up the filename with the scene name and band name and open the
       file for write access */
    sprintf (file_name, "%s_%s.img", scene_name, bmeta[0].name);
    strcpy (bmeta[0].file_name, file_name);
    this->fp_bin = open_raw_binary (file_name, "w");
    if (this->fp_bin == NULL)
    {
        RETURN_ERROR ("unable to open output file", "OpenOutput", NULL);
    }

    this->open = true;

    return this;
}


/******************************************************************************
Description: 'CloseOutput' closes the output files which are open.

Input Parameters:
 this           'output' data structure

Output Parameters:
 this           'output' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
                  'false' = error return
*****************************************************************************/
bool
CloseOutput (Output_t *this)
{
    if (!this->open)
    {
        RETURN_ERROR ("image files not open", "CloseOutput", false);
    }

    close_raw_binary (this->fp_bin);
    this->open = false;

    return true;
}


/******************************************************************************
Description: 'FreeOutput' frees the 'output' data structure memory.

Input Parameters:
 this           'output' data structure for which the fields are freed

Output Parameters:
 this           'output' data structure
 (returns)      status:
                  'true' = okay
                  'false' = error occurred
*****************************************************************************/
bool
FreeOutput (Output_t *this)
{
    if (this->open)
    {
        RETURN_ERROR ("file still open", "FreeOutput", false);
    }

    free (this);
    this = NULL;

    return true;
}


/******************************************************************************
Description: 'PutOutputLine' writes nlines of data to the output file.

Input Parameters:
 this           'output' data structure; the following fields are written:
                buf -- contains the line to be written
 final_mask     current nlines of data to be written (0-based)

Output Parameters:
 this           'output' data structure; the following fields are modified:
 (returns)      status:
                  'true' = okay
                  'false' = error return
*****************************************************************************/
bool
PutOutput (Output_t * this, unsigned char **final_mask)
{
    int il;
    void *buf = NULL;

    /* Check the parameters */
    if (this == (Output_t *) NULL)
    {
        RETURN_ERROR ("invalid input structure", "PutOutputLine", false);
    }

    if (!this->open)
    {
        RETURN_ERROR ("file not open", "PutOutputLine", false);
    }

    for (il = 0; il < this->size.l; il++)
    {
        buf = (void *) final_mask[il];
        if (write_raw_binary (this->fp_bin, 1, this->size.s,
                              sizeof (unsigned char), buf) != SUCCESS)
        {
            RETURN_ERROR ("writing output line", "PutOutput", false);
        }
    }

    return true;
}
