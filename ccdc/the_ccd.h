int the_ccd
(
    int valid_num_scenes,            /* Number of scenes after fmask filtering*/
    int row,                         /* The input row of the data frame.      */
    int col,                         /* The input column of the data frame.   */
    int nrows,                       /* The number of rows in "tile" of data  */
    int ncols,                       /* Number of columns in "tile" of data   */
    int num_samples,                 /* Number of samples in input grid       */
    int *updated_sdate_array,        /* sdate array after cfmask filtering    */
    unsigned char *updated_fmask_buf,/*sub-set of fmask buf, valid pixels only*/
    int *buf,                        /* This is the image bands buffer.       */
    int clr_sum,                     /* Total number of clear cfmask pixels   */
    int sn_sum,                      /* Total number of snow  cfmask pixels   */
    int water_sum,                   /* Total of cfmask water pixels.         */
    int all_sum,                     /* Total of all cfmask pixels            */
    bool verbose,                    /* Verbose flag for printing messages    */
    Output_t *rec_cg,                /* Output structure of vals and metadata */
    int *num_fc                      /* Number of fitted curves, cumulative   */
                                     /* for all rows/cols and models          */
);
