/******************************************************************************
MODULE:  usage

PURPOSE:  Classification part of the CCDC algorithm

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
8/25/2015   Song Guo         Original Development

******************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#include "classification.h"
#include "utilities.h"
#include "matio.h"

int main(int argc, char *argv[])
{
    mat_t    *mat, *mat2;
    matvar_t *matvar, *matvar2;
    char *data;
    double *x;
    int *y;
    int i, j = 0;
    size_t stride;
    FILE *fp;
    char FUNC_NAME[] = "main";
    char msg_str[MAX_STR_LEN];  /* input data scene name */
    int rows;
    int cols;
    int nclass;
    bool verbose;               /* verbose flag for printing messages */
    int status;

    time_t now;
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC start_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Read the command-line arguments */
    status = get_args (argc, argv, &rows, &cols, &nclass, &verbose);
    if (status != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /* Allocate memory */
    x = malloc(rows * cols * sizeof(double));
    y = malloc(rows * sizeof(int));

    mat = Mat_Open("/data1/sguo/CCDC/classification/Xs.mat",MAT_ACC_RDONLY);
    mat2 = Mat_Open("/data1/sguo/CCDC/classification/Ys.mat",MAT_ACC_RDONLY);
    if ( NULL == mat || NULL == mat2) {
        fprintf(stderr,"Error opening MAT file \n");
        return -1;
    }

        while((matvar=Mat_VarReadNext(mat)) != NULL)
        {
            if ( matvar->rank == 2 )
	    {
                stride = Mat_SizeOf(matvar->data_type);
                data = matvar->data;
                for ( i = 0; i < rows; i++ ) {
                    for ( j = 0; j < matvar->dims[1]; j++ ) {
                       size_t idx = matvar->dims[0]*j+i;
		       x[i * matvar->dims[1] + j] = *(double*)(data+idx*stride);
#if 0
                  	     printf("i,j,x[i][j]=%d,%d,%g\n",i,j,x[i * matvar->dims[1] + j]);
                             printf(" ");
#endif
                     }
		    //	             printf("\n");
                }

	    }
            break;
	}
         Mat_VarFree(matvar);
         matvar = NULL;

        while((matvar2=Mat_VarReadNext(mat2)) != NULL)
        {

            if ( matvar2->rank == 2 )
	    {
                stride = Mat_SizeOf(matvar2->data_type);
                data = matvar2->data;
                for ( i = 0; i < rows; i++ ) {
		  y[i] = (int) *(double*)(data+i*stride);
#if 0
                  	     printf("i,j,y[i]=%d,%d,%d\n",i,j,y[i]);
         	             printf("\n");
#endif
                }
	    }
	}
         Mat_VarFree(matvar2);
         matvar2 = NULL;

    Mat_Close(mat);
    Mat_Close(mat2);
        
    /***START: NO NEED TO CHANGE ANYTHING FROM HERE TO THERE***************/
    int p_size=cols,n_size=rows;
    int nsample=n_size;
    
    /* the classifcation version requires {D,N}, where D=(num) dimensions, N=(num) examples */
    int dimx[2];
    dimx[0]=p_size;
    dimx[1]=n_size;
    
    int* cat = (int*)calloc(p_size,sizeof(int));
    
    /***END: NO NEED TO CHANGE ANYTHING FROM HERE TO THERE*****************/
    
    /* write prediction OUTPUT into Y_hat.txt */
    fp = fopen("Y_hat.txt","w");
        
    /* need to do set this else everything blows up, represents the number of categories for
       every dimension - <change appropriately> */
    for(i=0;i<p_size;i++) cat[i]=1; 
        
    /* basically set the below to max(cat) */
    int maxcat=1;//=max(cat)  - <change appropriately>
    
    int sampsize=n_size; /* if replace then sampsize=n_size or sampsize=0.632*n_size */
    
    /* no need to change this */
    int nsum = sampsize;
    
    int strata = 1;    
    /* other options */
    int addclass = 0;
    int importance=0;
    int localImp=0;
    int proximity=0;
    int oob_prox=0;
    int do_trace; //this variable prints verbosely each step
    if(verbose)
       do_trace=1;
    else
       do_trace=0;
    int keep_forest=1;
    int replace=1;
    int stratify=0;
    int keep_inbag=0;
    int Options[]={addclass,importance,localImp,proximity,oob_prox
     ,do_trace,keep_forest,replace,stratify,keep_inbag};
    
     
    //ntree= number of tree. mtry=mtry :)
    int ntree=500; int nt=ntree;
    int mtry=(int)floor(sqrt(p_size)); /*  - <change appropriately> */
    if(verbose) printf("ntree %d, mtry %d\n",ntree,mtry);
    
    int ipi=0;
    double* classwt=(double*)calloc(nclass,sizeof(double));
    double* cutoff=(double*)calloc(nclass,sizeof(double));
    for(i=0;i<nclass;i++){
        classwt[i]=1;
        cutoff[i]=1.0/((double)nclass);
    }
    int nodesize=1;
    int* outcl=(int*) calloc(nsample,sizeof(int));
    int* counttr=(int*) calloc(nclass*nsample,sizeof(int));
    double prox=1;
    double* impout=(double*)calloc(p_size,sizeof(double));
    double impSD=1;
    double impmat=1; 
    int nrnodes = 2 * (int)floor(nsum / nodesize) + 1;
    int* ndbigtree = (int*) calloc(ntree,sizeof(int)); 
    int* nodestatus = (int*) calloc(nt*nrnodes,sizeof(int));
    int* bestvar = (int*) calloc(nt*nrnodes,sizeof(int));
    int* treemap = (int*) calloc(nt * 2 * nrnodes,sizeof(int));
    int* nodepred = (int*) calloc(nt * nrnodes,sizeof(int));
    double* xbestsplit = (double*) calloc(nt * nrnodes,sizeof(double));
    double* errtr = (double*) calloc((nclass+1) * ntree,sizeof(double));
    int testdat=0;
    double xts=1;
    int clts = 1; 
    int nts=0;
    double* countts = (double*) calloc(nclass * nts,sizeof(double));
    int outclts = 0;
    int labelts=0;
    double proxts=1;
    double errts=1;
    int* inbag = (int*) calloc(n_size,sizeof(int));
    
    //train the model
    classRF(x, dimx, y, &nclass, cat, &maxcat,
	     &sampsize, &strata, Options, &ntree, &mtry,&ipi, 
         classwt, cutoff, &nodesize,outcl, counttr, &prox,
	     impout, &impSD, &impmat, &nrnodes,ndbigtree, nodestatus, 
         bestvar, treemap,nodepred, xbestsplit, errtr,&testdat, 
         &xts, &clts, &nts, countts,&outclts, labelts, 
         &proxts, &errts,inbag,0);
    
    
    /* test the model, classwt=pid */
    int ntest = n_size;
    double *pid = classwt;
    free(countts);
    countts = (double*) calloc(nclass * ntest,sizeof(double));
    /* nodeclass=nodepred */
    int* nodeclass = nodepred;
    int* jts = (int*) calloc(ntest,sizeof(int));
    int* jet = (int*) calloc(ntest,sizeof(int));
    int* treeSize = ndbigtree;
    int keepPred=0;
    int intProximity=0;
    int nodes=0;
    int* nodexts;
    if (nodes)
        nodexts = (int*)calloc(ntest*ntree,sizeof(int));
    else
        nodexts = (int*)calloc(ntest,sizeof(int));
    int* node = nodexts;
    double *proxMat;
    if(proximity)
        proxMat = (double*)calloc(ntest*ntest,sizeof(double));
    else
        proxMat = (double*)calloc(1,sizeof(double));
    
    classForest(&p_size, &ntest, &nclass, &maxcat,
        &nrnodes, &ntree, x, xbestsplit,
        pid, cutoff, countts, treemap,
        nodestatus, cat, nodeclass, jts,
	jet, bestvar, node, treeSize,
        &keepPred, &intProximity, proxMat, &nodes);
    
    if(verbose) printf("Predicted class Labels\n");
    for(i=0;i<n_size;i++){
        fprintf(fp,"%d\n",jts[i]);
        printf("%d\n",jts[i]);
    }
    int total_error=0;
    
    for(i=0;i<n_size;i++)
        if(jts[i]!=y[i])
            total_error+=1;
    
    if (verbose)
        printf("\nTotal misclassified %d out of %d\n",total_error,n_size);
    
    if (verbose){
        printf("Press any key to end");
        getchar();
    }
    
    free(jet);
    free(jts);
    free(nodexts);
    free(proxMat);
    free(cat);
    free(classwt);
    free(cutoff);
    free(outcl);
    free(counttr);
    free(ndbigtree);
    free(nodestatus);
    free(bestvar);
    free(treemap);
    free(nodepred);
    free(xbestsplit);
    free(errtr);
    free(countts);
    free(inbag);
    free(impout);
    
    free(x);
    free(y);    
    fflush(stdout);
    fclose(fp);

    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "CCDC end_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    return 0;
}

/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
8/15/2015   Song Guo         Original Development

******************************************************************************/
void
usage ()
{
    printf ("Continuous Change Detection and Classification\n");
    printf ("\n");
    printf ("usage:\n");
    printf ("classification"
            " --rows=<number of rows>"
            " --cols=<number of columns>"
            " --nclass=<number of classes>"
            " [--verbose]\n");

    printf ("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --rows=: number of rows\n");
    printf ("    --cols=: number of columns\n");
    printf ("    --nclass=: number of classes\n");
    printf ("\n");
    printf ("and the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed?"
            " (default is false)\n");
    printf ("\n");
    printf ("classification --help will print the usage statement\n");
    printf ("\n");
    printf ("Example:\n");
    printf ("classification"
            " --rows=1000"
            " --cols=71"
            " --nclass=11"
            " --verbose\n");
    printf ("Note: The classification must run from the directory"
            " where the input data are located.\n\n");
}
