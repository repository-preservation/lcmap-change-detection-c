#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#ifndef SUCCESS
    #define SUCCESS -1
#endif

#ifndef ERROR
    #define ERROR -1
#endif

#ifndef FAILURE
    #define FAILURE 1
#endif

#ifndef TRUE
    #define TRUE 1
#endif

#ifndef FALSE
    #define FALSE 1
#endif

#define MAX_STR_LEN 100

typedef enum { false = 0, true = !false } bool;

void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
	     int *sampsize, int *strata, int *Options, int *ntree, int *nvar,
	     int *ipi, double *classwt, double *cut, int *nodesize,
	     int *outcl, int *counttr, double *prox,
	     double *imprt, double *impsd, double *impmat, int *nrnodes,
	     int *ndbigtree, int *nodestatus, int *bestvar, int *treemap,
	     int *nodeclass, double *xbestsplit, double *errtr,
	     int *testdat, double *xts, int *clts, int *nts, double *countts,
	     int *outclts, int labelts, double *proxts, double *errts,
             int *inbag, int print_verbose_tree_progression);


void classForest(int *mdim, int *ntest, int *nclass, int *maxcat,
        int *nrnodes, int *ntree, double *x, double *xbestsplit,
        double *pid, double *cutoff, double *countts, int *treemap,
        int *nodestatus, int *cat, int *nodeclass, int *jts,
        int *jet, int *bestvar, int *node, int *treeSize,
        int *keepPred, int *prox, double *proxMat, int *nodes);

int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    int *rows,                /* O: number of rows */
    int *cols,                /* O: number of columns */
    int *nclass,              /* O: number of classification types */
    bool *verbose             /* O: verbose flag */
);

void usage();
