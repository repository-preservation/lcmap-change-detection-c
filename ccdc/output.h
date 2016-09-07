#ifndef OUTPUT_H
#define OUTPUT_H

/* Aall defines in one file to avoid conflicts */
#include "defines.h"

typedef struct {
  int row;
  int col;
} Position_t;

/* Structure for the 'output' data type */
typedef struct
{
    int t_start;           /* time when series model gets started */
    int t_end;             /* time when series model gets ended */
    int t_break;           /* time when the first break (change) is observed */
    float coefs[NUM_BANDS][NUM_COEFFS];
                           /*  coefficients for each time series model for each 
                               spectral band*/    
    float rmse[NUM_BANDS];
                           /*  RMSE for each time series model for each 
                               spectral band*/    
    //Position_t pos;        /* the location of each time series model */
    int pos;               /* the pixel location of each time series model */
    float change_prob;     /* the probability of a pixel that have undergone 
                              change (between 0 and 100) */
    int num_obs;           /* the number of "good" observations used for model 
                              estimation */
    int category;          /* the quality of the model estimation (what model 
                              is used, what process is used) 
                              1x: persistent snow    2x: persistent water 
                              3x: Fmask fails        4x: normal precedure
                              x1: mean value (1)     x4: simple fit (4)
                              x6: basic fit (6)      x8: full fit (8) */
    float magnitude[NUM_BANDS];/* the magnitude of change (difference between model 
                                  prediction and observation for each spectral band)*/
} Output_t;

#endif
