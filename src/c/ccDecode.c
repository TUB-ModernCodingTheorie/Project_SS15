#include "mex.h"

#define INFINITY     1.e10

/* Decoder */
int ccDecode(int **prior,
                  int ***backward,
                  int init_state,
                  int final_state,
                  const int Memory,
                  const int Kbit,
                  const int Mbit,
                  const int Nbit,
                  const int frame_length
                )
{
    if (Nbit % Mbit != 0) {
        printf("ERROR: Nbit is not a multiple of Mbit\n");
        exit(0);
    } 
    
    int asize,
        sdim,
        input_size,     /* */
        state_size,     /* Number of state elements */
        state,          /* index for the states */
        old_state,      /* index of a statets */
        t,x,i,i0,jj,d,
        survstate,      /* */
        survinput,
        survoutput,
        b,
        mask;
        
    double path,max,maxmax;
   
    double *path0 = calloc(state_size, sizeof(double));
    double *path1 = calloc(state_size, sizeof(double));
    int *detected_info_frame  = calloc(frame_length,sizeof(int));
    int *detected_code_frame  = calloc(frame_length,sizeof(int));
    int ***trace_back = calloc(state_size,sizeof(*trace_back));
    
    for (i = 0; i < state_size; ++i)
    {
        trace_back[i] = calloc(frame_length, sizeof(*trace_back[i]));
        if (trace_back[i] == NULL) {
            fprintf (stderr, "Memory allocation failure on trace_back");
        }
    }
    
    state_size = (1 << Memory);
    asize = (1 << Mbit);
    input_size = (1 << Kbit);
    
    sdim = Nbit/Mbit;
    
    /* initialize all states */
    for (state = 0 ; state < state_size ; state++)
        path0[state] = -INFINITY;
    
    /* begin */
    path0[init_state] = 0;

    mask = asize - 1;
	
    for (t = 0 ; t < frame_length ; t++) {
        maxmax = -INFINITY;
        i0  = t*sdim;
        
        for(state = 0; state < state_size; state++) {
            max = -INFINITY;
            for( b = 0 ; b < input_size ; b++) {
                x = backward[state][b][2];
                old_state = backward[state][b][0];
                path = path0[old_state];
                for (d = 0 ; d < sdim ; d++) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    path += prior[i][jj];
                }
                if(path>=max) {
                    survinput  = backward[state][b][1];
                    survoutput = x;
                    survstate  = old_state;
                    max = path;
                }
            }
            trace_back[state][t][0] = survstate;
            trace_back[state][t][1] = survinput;
            trace_back[state][t][2] = survoutput;
            path1[state] = max;
            
            if (max >= maxmax)
                maxmax = max;
        }
        for(state = 0 ; state < state_size ; state++)
            path0[state] = path1[state] - maxmax;
    }

    /* Final trace back (with trellis termination) */

    state = final_state;
    for (t = frame_length-1 ; t >= 0 ; t--) {
        detected_info_frame[t] = trace_back[state][t][1];
        detected_code_frame[t] = trace_back[state][t][2];
        state = trace_back[state][t][0];
    }
    
    free(path0);
    free(path1);
    free(detected_info_frame);
    free(detected_code_frame);
    for (i = 0; i < state_size; ++i)
    {
        free(trace_back[i]);
    }
    free(trace_back);
    
    return(state);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *c;	        /* output codewords */
    mxArray *sN;	/* final state */

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
    }
    
    /* make sure the first input argument is a string */
    if( !mxIsDouble(prhs[0]) ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input must be a double.");
    }
	if( !mxIsDouble(prhs[1]) ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input must be a double.");
    }
	if( !mxIsDouble(prhs[2]) ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input must be a double.");
    }

    if( mxGetNumberOfDimensions(prhs[0]) != 3){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:dimensionMismatch","The forward table must have three dimensions.");
    }
    if( mxGetM(prhs[1]) != 1){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
    if( mxGetM(prhs[2]) != 1){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }

    int **prior;
    int ***backward;
    int init_state;
    int final_state;
    const int Memory;
    const int Kbit;
    const int Mbit;
    const int Nbit;
    const int frame_length;
    
    ccDecode(prior,
                  backward,
                  init_state,
                  final_state,
                  Memory,
                  Kbit,
                  Mbit,
                  Nbit,
                  frame_length
                )

    /* Define output */
}
