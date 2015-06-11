/**
 * This function is based on /extern/examples/mex/arrayProduct.c
 * 
 * To compile, type in the matlab shell: `mex helloWorld.c`
 * 
 * For more documentation go to:
 * http://fr.mathworks.com/help/matlab/mex-library.html 
 */

#include <stdio.h>
#include <string.h>
#include "mex.h"

/**
 * The function must be called mexFunction with inputs:
 *    - int nlhs = number of OUTPUTS  (l as left)
 *    - mxArray *plhs[] = OUTPUTS
 *    - int nrhs = number of INPUT (r as right)
 *    - mxArray *rlhs[] = INTPUTS
 * 
 * In matlab it will be called by
 *      [a,b] = helloWorld(c)
 * where    a = plhs[0]
 *          b = plhs[1]
 *          c = prhs[2]
 */

/* The computational routine */
void helloWorld(char * name)
{
    char text[30];
    strcpy(text,"Hello ");
    strcat(text,name);
    strcat(text," !\n");
    printf(text);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *name;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","One input required.");
    }
    if(nlhs!=0) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","No output required.");
    }
    
    /* make sure the first input argument is a string */
    if( !mxIsChar(prhs[0]) ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notString","Input name must be a string/double.");
    }
    
    /* get the value of the String input  */
    name = (char *) mxCalloc(mxGetN(prhs[0]), sizeof(char));
    mxGetString(prhs[0], name, mxGetN(prhs[0])+1);    
    
    /* call the computational routine */
    helloWorld(name);
}
