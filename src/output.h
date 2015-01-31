#ifndef OUTPUT_H
#define OUTPUT_H


#include <time.h>


#include "ccdc.h"
#include "input.h"
#include "date.h"

#define FILL_VALUE 255
#define NUM_COEFFS 5
#define NUM_BANDS 7

/* Structure for the 'output' data type */

typedef struct
{
    bool open;             /* Flag to indicate whether output file is open 
                              for access; 'true' = open, 'false' = not open */
    FILE *fp_bin;          /* File pointer for binary output file */
    Date_t t_start;        /* time when series model gets started */
    Date_t t_end;          /* time when series model gets ended */
    Date_t t_break;        /* time when the first break (change) is observed */
    float model_coefs[NUM_BANDS][NUM_COEFFS];
                           /*  coefficients for each time series model for each 
                               spectral band*/    
    float model_rmse[NUM_BANDS][NUM_COEFFS];
                           /*  RMSE for each time series model for each 
                               spectral band*/    
    int pos;               /* the location of each time series model */
    float change_prob;     /* the probability of a pixel that have undergone 
                              change (between 0 and 100) */
    int num_obs;           /* the number of "good" observations used for model 
                              estimation */
    float category;        /* the quality of the model estimation (what model 
                              is used, what process is used)*/
    float magnitude;       /* the magnitude of change (difference between model 
                              prediction and observation for each spectral band)*/
} Output_t;


/* Prototypes */
Output_t *OpenOutput (Espa_internal_meta_t *in_meta, Input_t *input);


bool PutOutput (Output_t *this, unsigned char **final_mask);


bool CloseOutput (Output_t *this);


bool FreeOutput (Output_t *this);


#endif
