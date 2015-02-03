#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <stdio.h>

#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "ccdc.h"

#define NUM_LASSO_BANDS 5
#define TOTAL_BANDS 8
#define MAX_NUM_FC 10 /* Values change with number of pixels run */
int lasso_band_list[NUM_LASSO_BANDS] = {2, 3, 4, 5, 7};

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
--------    ---------------  -------------------------------------
1/15/2013   Song Guo         Original Development

NOTES: type ./ccdc --help for information to run the code
******************************************************************************/
int
main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";
    char msg_str[MAX_STR_LEN];  /* input data scene name */
    char filename[MAX_STR_LEN];         /* input binary filenames */
    char directory[MAX_STR_LEN];        /* input/output data directory */
    char extension[MAX_STR_LEN];        /* input TOA file extension */
    char data_dir[MAX_STR_LEN];         /* input/output data directory */
    char appendix[MAX_STR_LEN];         /* input TOA file extension */
    Input_t *input = NULL;      /* input data and meta data */
    char scene_name[MAX_STR_LEN];       /* input data scene name */
    char command[MAX_STR_LEN];
    float min_rmse;
    float t_cg;
    float t_max_cg;
    int conse;
    int status;                 /* return value from function call */
    Output_t *rec_cg[MAX_NUM_FC-1] = NULL;/* output structure and metadata */
    bool verbose;               /* verbose flag for printing messages */
    float alb = 0.1;
    int i, k;
    int num_points;
    char **scene_list = NULL;
    float **results = NULL;
    FILE *fd;
    int num_scenes;
    int i;
    int min_num_c = 4;
    int mid_num_c = 6;
    int max_num_c = 8;
    int num_c = 8;            /* max number of coefficients for the model */
    int n_times = 3;          /* number of clear observations/coefficients*/
    int num_fc = 0;           /* intialize NUM of Functional Curves */
    float num_yrs = 365.25;   /* number of days per year */
    int num_byte = 2;         /* number of bytes: int16 */
    int nbands = 8;           /* bands 1-7, cfmask */
    int num_b1 = 2;           /* Band for multitemporal cloud/snow detection 
                                 (green) */ 
    int num_b2 = 5;           /* Band for multitemporal shadow/snow shadow 
                                 detection (SWIR) */
    int t_const = 400;        /* Threshold for cloud, shadow, and snow detection */
    int mini_yrs = 1;         /* minimum year for model intialization */
    int num_detect = NUM_LASSO_BANDS;/* number of bands for change detection */
    float p_min = 0.1;        /* percent of ref for mini_rmse */
    float t_ws = 0.95;        /* no change detection for permanent water pixels */
    float t_sn = 0.6;         /* no change detection for permanent snow pixels */ 
    float t_cs = 0.6;         /* Fmask fails threshold */
    int *sdate;
    Input_meta_t *meta;
    int row, col;
    int landsat_number;
    int fmask_sum = 0;
    int claer_water_sum = 0;
    int clear_land_sum = 0;
    int snow_sum = 0;
    int cloud_and_shadow_sum = 0;
    float ws_pct;
    float sn_pct;
    float cs_pct;
    int n_ws = 0;
    int n_sn = 0;
    int n_clr = 0;
    int *clrx;
    int **clry;
    int i_start;
    int end;
    float **fit_cft;
    float **rmse;
    int i_span;
    int update_num_c;
    int v_qa;
    int8 *id_label;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    status = get_args (argc, argv, &row, &col, &min_rmse, &t_cg, &t_max_cg, 
                       &conse, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /* allocate memory for scene_list */
    scene_list = (char **) allocate_2d_array (MAX_SCENE_LIST, MAX_STR_LEN,
                                         sizeof (char));
    if (scene_list == NULL)
    {
        RETURN_ERROR ("Allocating scene_list memory", FUNC_NAME, FAILURE);
    }

    /* check if scene_list.txt file exists, if not, create the scene_list
       from existing files in the current data working directory */
    if (access("scene_list.txt", F_OK) != -1) /* File exists */
    {
        num_scenes = MAX_SCENE_LIST;
    }
    else /* File not exists */
    {
        stats = create_scene_list("L*_sr_band1.hdr", num_scenes, scene_list); 
    }

    fd = fopen("scene_list.txt", "r");
    if (fd == NULL)
    {
        RETURN_ERROR("Opening scene_list file", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_scenes; i++)
    {
        if (fscanf(fd, "%s", scene_list[i]) == EOF)
        {
            RETURN_ERROR("Reading scene_list file", FUNC_NAME, FAILURE);
            num_scenes = i;
            break;
        }
    }

    /* Allocate memory */
    sdate = malloc(num_scenes * sizeof(int));
    if (sdate == NULL)
        RETURN_ERROR("ERROR allocating sdate memory", FUNC_NAME, FAILURE);

    clrx = malloc(num_scenes * sizeof(int));
    if (clrx == NULL)
        RETURN_ERROR("ERROR allocating clrx memory", FUNC_NAME, FAILURE);

    id_label = malloc(num_scenes * sizeof(int8));
    if (clrx == NULL)
        RETURN_ERROR("ERROR allocating id_label memory", FUNC_NAME, FAILURE);

    clry = (int16 **) allocate_2d_array (num_scenes, TOTAL_BANDS - 1, 
                                         sizeof (int16));
    if (clry == NULL)
    {
        RETURN_ERROR ("Allocating clry memory", FUNC_NAME, FAILURE);
    }

    fit_cft = (float **) allocate_2d_array (max_num_c, TOTAL_BANDS - 1, 
                                         sizeof (float));
    if (fir_ctf == NULL)
    {
        RETURN_ERROR ("Allocating fit_cft memory", FUNC_NAME, FAILURE);
    }

    rmse = malloc((TOTAL_BANDS - 1) * sizeof(float));
    if (rmse == NULL)
    {
        RETURN_ERROR ("Allocating rmse memory", FUNC_NAME, FAILURE);
    }

    /* sort scene_list based on year & julian_day */
    status = sort_scene_based_on_year_doy(scene_list, num_scenes, sdate);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Create the Input metadata structure */
    meta = (Input_meta_t *)malloc(sizeof(Input_meta_t));
    if (meta == NULL) 
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, FAILURE);

    /* Get the metadata, all scene metadata are the same for stacked scenes */
    status = read_envi_header(scene_list[0], meta);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Calling sort_scene_based_on_year_jday", 
                      FUNC_NAME, EXIT_FAILURE);
    }

    if (verbose)
    {
        /* Print some info to show how the input metadata works */
        printf ("DEBUG: Number of input lines: %d\n", meta->lines);
        printf ("DEBUG: Number of input samples: %d\n", meta->samples);
        printf ("DEBUG: UL_MAP_CORNER: %f, %f\n", meta->upper_left_x,
                meta->upper_left_y);
        printf ("DEBUG: ENVI data type: %d\n", meta->data_type);
        printf ("DEBUG: ENVI byte order: %d\n", meta->byte_order);
        printf ("DEBUG: UTM zone number: %d\n", meta->utm_zone);
        printf ("DEBUG: Pixel size: %d\n", meta->pixel_size);
        printf ("DEBUG: Envi save format: %s\n", meta->interleave);
    }

    /* Open input files */
    FILE *fp_bin[num_scenes][TOTAL_BANDS];
    short int buf[num_scenes][TOTAL_BANDS-1];
    unsigned char fmask_buf[num_scenes];
    /* Open input files */
    for (i = 0; i < num_scenes; i++)
    {
        for (k = 0; k < TOTAL_BANDS; k++)
        {
            landsat_number = atoi(sub_string(scene_list[i],2,1));
            if (landsat_number != 8)
            {
                if (k == 5)
                    sprintf(filename, "%s_toa_band6.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+1);
            }
            else
            {
                if (k == 5)
                    sprintf(filename, "%s_toa_band10.img", scene_list[i]);
                else if (k == 7)
                    sprintf(filename, "%s_cfmask.img", scene_list[i]);
                else if (k == 6)
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+1);
                else 
                    sprintf(filename, "%s_sr_band%d.img", scene_list[i], k+2);
            }

            fp_bin[i][k] = open_raw_binary(filename,"rb");
            if (fp_bin[i][k] == NULL)
                printf("error open %d scene, %d bands files\n",i, k+1);
            if (k != TOTAL_BANDS-1)
            {
                fseek(fp_bin[i][k], (row * meta->samples + col)*sizeof(short int), SEEK_SET);
                if (read_raw_binary(fp_bin[i][k], meta->lines, meta->samples,
                    sizeof(short int), &buf[i][k]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
                printf("scene_number,band_number,buf[i][k] = %d, %d, %d\n",
                   i,k+1,buf[i][k]);
            }
            else
            {
                fseek(fp_bin[i][k], (row * meta->samples + col)*sizeof(unsigned char), 
                    SEEK_SET);
                if (read_raw_binary(fp_bin[i][k], meta->lines, meta->samples,
                    sizeof(unsigned char), &fmask_buf[i]) != 0)
                    printf("error reading %d scene, %d bands\n",i, k+1);
                printf("scene_number,band_number,buf[i][k] = %d, %d, %d\n",
                       i,k+1,fmask_buf[i]);
            }
        close_raw_binary(fp_bin[i][k]);
        }
    }

    /* Only run CCDC for places where more than 50% of images has data */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] < 255)
            fmask_sum++;
    }
    if (fmask_sum < (uint) 0.5 * num_scenes)
        RETURN_ERROR ("Not enough clear-sky pisels", FUNC_NAME, EXIT_FAILURE);
    else
        printf("Clear-sky pixel percentage = %f7.2\n", fmask_sum / num_scenes);

    /* pixel value ranges should follow physical rules */
    for (i = 0; i < num_scenes; i++)
    { 
        if ((buf[i][0] > 0) && (buf[i][0] < 10000) &&
            (buf[i][1] > 0) && (buf[i][1] < 10000) &&
            (buf[i][2] > 0) && (buf[i][2] < 10000) &&
            (buf[i][3] > 0) && (buf[i][3] < 10000) &&
            (buf[i][4] > 0) && (buf[i][4] < 10000) &&
            (buf[i][5] > -9320) && (buf[i][5] < 7070) &&
            (buf[i][6] > 0) && (buf[i][6] < 10000))
            id_label[i] = 1;
        else
            id_label[i] = 0;
    }

    /* Get each mask pixel totals */
    for (i = 0; i < num_scenes; i++)
    { 
        if (fmask_buf[i] == 0)
            clear_land_sum++;
        if (fmask_buf[i] == 1)
            clear_water_sum++;
        if (fmask_buf[i] == 2 || fmask_sum[i] ==4)
            cloud_and_shadow_sum++;
        if (fmask_buf[i] == 3)
            snow_sum++;
    }

    /* percent of water pixels */
    ws_pct = (float) clear_water_sum / (float) (clear_land_sum + clear_water_sum);

    /* percent of snow observations */
    sn_pct =  (float) snow_sum / (float)(clear_land_sum + 
               clear_water_sum + snow_sum); 

    /* percent of cloud and cloud shadow pixels */
    cs_pct = (float)cloud_and_shadow_sum / (float)fmask_sum;  

    /* fit only water pixels (permanet water) */
    if (ws_pct > t_ws)
    {
        for (i = 0; i < num_scenes; i++)
        { 
            if (fmask_buf[i] == 1 && id_label[i] == 1)
            {
                clrx[n_ws] = sdate[i];
                for (k = 0; k < TOTAL_BANDS - 1; k++)
                     clry[n_ws][k] = buf[i][k];
                n_ws++;
            }   
        }
        end = n_ws;

        printf("Permanent water pixel = %f\n", 100.0 * ws_pct);

        if (n_ws < ntimes * min_num_c ) /* not enough water pixels */
            RETURN_ERROR ("Not enough good water observations\n", 
                 FUNC_NAME, EXIT_FAILURE);
        else
        {
            /* start model fit for snow persistent pixels */
            printf ("Fit permanent water observations\n"); 

            /* the first observation for TSFit */
            i_start = 0; /* the first observation for TSFit */

            /* use all water and snow values to fit */
            i_span = n_ws;
            update_cft(i_span, n_times, min_num_c, mid_num_c, max_num_c, num_c,
                       &update_num_c);

            for (k = 0; k < TOTAL_BANDS - 1; k++)
            {
             // TO DO: [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),update_num_c);   
            }

            /* update information at each iteration */
            rec_cg[num_fc].t_start = clrx[i_start]; /* record time of curve start */
            rec_cg[num_fc].t_end = clrx[end]; /* record time of curve end */
            rec_cg[num_fc].t_break = 0; /* no break at the moment */
            rec_cg[num_fc].pos.row = row; /* record postion of the pixel */
            rec_cg[num_fc].pos.col = col; /* record postion of the pixel */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
            {
                for (k = 0; k < update_num_c; k++)
                    rec_cg[num_fc].coefs[i][k] = fit_cft[i][k]; /* record fitted coefficients */
                rec_cg[num_fc].rmse[i] = rmse[i]; /* record rmse of the pixel */
            }
            rec_cg[num_fc].change_prob = 0.0; /* record change probability */
            rec_cg[num_fc].num_obs = n_ws; /* record number of observations */
            rec_cg[num_fc].category = 20 + update_num_c; /* record fit category */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
                rec_cg[num_fc].magnitude = 0.0; /* record change magnitude */
            num_fc++;    /* NUM of Fitted Curves (num_fc) */
        }
    }
    else if (sn_pct > t_sn) /* fit only snow pixels (permanet snow) */
    {
        for (i = 0; i < num_scenes; i++)
        { 
            if (fmask_buf[i] == 3)
            {
                clrx[n_sn] = sdate[i];
                for (k = 0; k < TOTAL_BANDS - 1; k++)
                     clry[n_sn][k] = buf[i][k];
                n_sn++;
            }   
        }
        end = n_sn;

        printf("Permanent snow pixel = %f\n", 100.0 * sn_pct);

        if (n_sn < ntimes * min_num_c ) /* not enough snow pixels */
            RETURN_ERROR ("Not enough good snow observations\n", 
                 FUNC_NAME, EXIT_FAILURE);
        else
        {
            /* start model fit for snow persistent pixels */
            printf ("Fit permanent snow observations\n"); 

            /* the first observation for TSFit */
            i_start = 0; /* the first observation for TSFit */

            /* treat saturated and unsaturated pixels differently */
            for (k = 0; k < TOTAL_BANDS -1; k++)
            {
                i_span = 0;
                if (k != TOTAL_BANDS - 3) /* for optical bands */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > 0 && clry[i][k] < 10000)
                            i_span++;
                        if (i_span < min_num_c * n_times)
                            fit_cft[i][k] = 10000; /* fixed value for saturated
                                                       pixels */
                        else
                        {
                            update_cft(i_span, n_times, min_num_c, mid_num_c, 
                                       max_num_c, num_c, &update_num_c);

                            for (k = 0; k < TOTAL_BANDS - 1; k++)
                            {
                                // TO DO: [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),update_num_c);  
                            }
                        }
                    } 
                }
                else /* for thermal band */
                {
                    for (i = 0; i < num_scenes; i++)
                    {
                        if (clry[i][k] > -9300 && clry[i][k] < 7070)
                            i_span++;
                        update_cft(i_span, n_times, min_num_c, mid_num_c, 
                                   max_num_c, num_c, &update_num_c);

                        // TO DO: [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),update_num_c);  
                    }
                }


            /* update information at each iteration */
            rec_cg[num_fc].t_start = clrx[i_start]; /* record time of curve start */
            rec_cg[num_fc].t_end = clrx[end]; /* record time of curve end */
            rec_cg[num_fc].t_break = 0; /* no break at the moment */
            rec_cg[num_fc].pos.row = row; /* record postion of the pixel */
            rec_cg[num_fc].pos.col = col; /* record postion of the pixel */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
            {
                for (k = 0; k < update_num_c; k++)
                    rec_cg[num_fc].coefs[i][k] = fit_cft[i][k]; /* record fitted coefficients */
                rec_cg[num_fc].rmse[i] = rmse[i]; /* record rmse of the pixel */
            }
            rec_cg[num_fc].change_prob = 0.0; /* record change probability */
            rec_cg[num_fc].num_obs = n_ws; /* record number of observations */
            rec_cg[num_fc].category = 20 + update_num_c; /* record fit category */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
                rec_cg[num_fc].magnitude = 0.0; /* record change magnitude */
            num_fc++;    /* NUM of Fitted Curves (num_fc) */
        }
    }
    else /* clear land or water pixels */
    {
        printf("Land pixel (water/snow) = %f/%f\n", 100.0 * ws_pct, 100.0 * sn_pct);

        if (cs_pct < t_cs)  /* Fmask works */
        {
            for (i = 0; i < num_scenes; i++)
            { 
                if ((fmask_buf[i] == 0 || fmask_buf[i] == 1) && id_label[i] == 1)
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_BANDS - 1; k++)
                         clry[n_clr][k] = buf[i][k];
                    n_clr++;
                }   
            }
            end = n_clr;
            v_qa = 40; /* QA var for normal procedure */
        }
        else /* Fmask fails in persistent commission */
        {
            /* catch back all valid wthin range observations */
            for (i = 0; i < num_scenes; i++)
            { 
                if ((fmask_buf[i] == 0 || fmask_buf[i] == 1
                     fmask_buf[i] == 2 || fmask_buf[i] == 4) && id_label[i] == 1)
                {
                    clrx[n_clr] = sdate[i];
                    for (k = 0; k < TOTAL_BANDS - 1; k++)
                         clry[n_clr][k] = buf[i][k];
                    n_clr++;
                }   
            }
            end = n_clr;
            v_qa = 30; /* QA var for normal procedure */
        }

        if (n_ws < ntimes * min_num_c ) /* not enough snow pixels */
            RETURN_ERROR ("Not enough good water or snow observations\n", 
                 FUNC_NAME, EXIT_FAILURE);
        else
        {
            /* start model fit for snow persistent pixels */
            printf ("Fit permanent water or snow observations\n"); 

            /* the first observation for TSFit */
            i_start = 0; /* the first observation for TSFit */

            /* use all water and snow values to fit */
            i_span = n_ws;
            update_cft(i_span, n_times, min_num_c, mid_num_c, max_num_c, num_c,
                       &update_num_c);

            for (k = 0; k < TOTAL_BANDS - 1; k++)
            {
             // TO DO: [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),update_num_c);   
            }

            /* update information at each iteration */
            rec_cg[num_fc].t_start = clrx[i_start]; /* record time of curve start */
            rec_cg[num_fc].t_end = clrx[end]; /* record time of curve end */
            rec_cg[num_fc].t_break = 0; /* no break at the moment */
            rec_cg[num_fc].pos.row = row; /* record postion of the pixel */
            rec_cg[num_fc].pos.col = col; /* record postion of the pixel */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
            {
                for (k = 0; k < update_num_c; k++)
                    rec_cg[num_fc].coefs[i][k] = fit_cft[i][k]; /* record fitted coefficients */
                rec_cg[num_fc].rmse[i] = rmse[i]; /* record rmse of the pixel */
            }
            rec_cg[num_fc].change_prob = 0.0; /* record change probability */
            rec_cg[num_fc].num_obs = n_ws; /* record number of observations */
            rec_cg[num_fc].category = 20 + update_num_c; /* record fit category */
            for (i = 0; i < TOTAL_BANDS - 1; i++)
                rec_cg[num_fc].magnitude = 0.0; /* record change magnitude */
            num_fc++;    /* NUM of Fitted Curves (num_fc) */
        }
    }

#if 0
    /* Open the output file */
    output = OpenOutput (&xml_metadata, input);
    if (output == NULL)
    {                           /* error message already printed */
        RETURN_ERROR ("Opening output file", FUNC_NAME, EXIT_FAILURE);
    }

    if (!PutOutput (output, pixel_mask))
    {
        RETURN_ERROR ("Writing output LST in HDF files\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Close the output file */
    if (!CloseOutput (output))
    {
        RETURN_ERROR ("closing output file", FUNC_NAME, EXIT_FAILURE);
    }

    /* Create the ENVI header file this band */
    if (create_envi_struct (&output->metadata.band[0], &xml_metadata.global,
                            &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Creating ENVI header structure.", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Write the ENVI header */
    strcpy (envi_file, output->metadata.band[0].file_name);
    cptr = strchr (envi_file, '.');
    if (cptr == NULL)
    {
        RETURN_ERROR ("error in ENVI header filename", FUNC_NAME,
                      EXIT_FAILURE);
    }

    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Writing ENVI header file.", FUNC_NAME, EXIT_FAILURE);
    }

    /* Append the LST band to the XML file */
    if (append_metadata (output->nband, output->metadata.band, xml_name)
        != SUCCESS)
    {
        RETURN_ERROR ("Appending spectral index bands to XML file.",
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Free the structure */
    if (!FreeOutput (output))
    {
        RETURN_ERROR ("freeing output file structure", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the input file and free the structure */
    CloseInput (input);
    FreeInput (input);

    free (xml_name);
    printf ("Processing complete.\n");
#endif

    /* Free memory allocation */
    free(sdate);
    free(clrx);
    free(rmse);
    free(id_label);
    status = free_2d_array ((void **) scene_list);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: scene_list\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    status = free_2d_array ((void **) clry);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: clry\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    status = free_2d_array ((void **) fit_cft);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: fit_cft\n", FUNC_NAME,
                      EXIT_FAILURE);
    }



    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC end_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    return EXIT_SUCCESS;
}


/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/15/2013   Song Guo         Original Development
8/15/2013   Song Guo         Modified to use TOA reflectance file 
                             as input instead of metadata file
2/19/2014   Gail Schmidt     Modified to utilize the ESPA internal raw binary
                             file format

******************************************************************************/
void
usage ()
{
    printf ("Landsat Surface Temperature\n");
    printf ("\n");
    printf ("usage: scene_based_lst"
            " --xml=input_xml_filename"
            " --dem=input_dem_filename"
            " --emi=input_emissivity_filename" " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    -xml: name of the input XML file\n");
    printf ("\n");
    printf ("where the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("scene_based_lst --help will print the usage statement\n");
    printf ("\n");
    printf ("Example: scene_based_lst"
            " --xml=LE70390032010263EDC00.xml"
            " --dem=17_30_DEM.tif"
            " --emi=AG100B.v003.-20.122.0001.bin" " --verbose\n");
    printf ("Note: The scene_based_lst must run from the directory"
            " where the input data are located.\n\n");
}
