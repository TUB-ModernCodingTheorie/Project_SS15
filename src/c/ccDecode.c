#include "mex.h"

#define INFINITY     1.e10

/**
 * Viterbi Decoder
 * 
 * @param *seq input encoded sequence
 * @param *bwd backward with
 *          bwd[current_state][idx][0] = previous_state
 *          bwd[current_state][idx][1] = input uncoded/plain bit
 *          bwd[current_state][idx][2] = output bit
 * @param initState state at which we begin
 * @param finalState state in which we are supposed to finish with
 * @param *decodedSeq decoded sequence output
 * 
 */
void ccDecode(mxArray *encodedFrameArray, mxArray *bwdArray, int initState, int finalState, mxArray *decodedFrameArray, mxArray *codedBitsArray)
{
    uint64_T *encodedFrame = (uint64_T*) mxGetData(encodedFrameArray);
    uint64_T *decodedFrame = (uint64_T*) mxGetData(decodedFrameArray);
    uint64_T *codedBits = (uint64_T*) mxGetData(codedBitsArray);
    uint64_T ***bwd = (uint64_T***) mxGetData(bwdArray);

    int stateSize = mxGetM(bwdArray);    /* Number of possible states: */
    int inputSize = mxGetN(bwdArray);    /* Number of bits per symbol  */
    int frameLength = mxGetN(encodedFrameArray)/*/??*/;  /* Number of unencoded input symbols per frame (ie: information bits) */
    
    int state,          /* index for the states */
        old_state,      /* index of a statets */
        t,x,i,i0,jj,d,
        survstate,      /* */
        survinput,
        survoutput,
        bit;
    
    double path,max,maxmax;
   
    double *path0 = calloc(stateSize, sizeof(double));
    double *path1 = calloc(stateSize, sizeof(double));
    int *decoded_frame  = calloc(frameLength,sizeof(int));
    int *detected_code_frame  = calloc(frameLength,sizeof(int));
    int ***trace_back = calloc(stateSize,sizeof(*trace_back));
    
    for (i = 0; i < stateSize; ++i)
    {
        trace_back[i] = calloc(frameLength, sizeof(*trace_back[i]));
        if (trace_back[i] == NULL) {
            fprintf (stderr, "Memory allocation failure on trace_back");
        }
    }
    
    /**
     * 1. Initialization
     * 
     * All metrics are -infty exept the initial state
     */
    for (state = 0 ; state < stateSize ; state++)
        path0[state] = -INFINITY;
    
    path0[initState] = 0;
    
    /**
     * 2. ACS
     *
     * For each information bit of the frame
     *  For each possible state transition
     *   For each possible parent state
     */
    
    for (t = frameLength-1 ; t >= 0 ; t--) {
        for (state = 0 ; state < stateSize ; state++) {
            max = -INFINITY;
            for (bit = 0 ; bit < inputSize ; bit++) {
                x = bwd[state][bit][2];
                old_state = bwd[state][bit][0];
                path = path0[old_state];
                
                if (path > max) {
                    survinput = bwd[state][bit][1];
                    survoutput = x;
                    max = path;
                }
                trace_back[state][t][0] = survstate;
                trace_back[state][t][1] = survinput;
                trace_back[state][t][2] = survoutput;
                path1[state] = max;
                
                if (max >= maxmax)
                    maxmax = max;
                
                for(state = 0 ; state < stateSize ; state++)
                    path0[state] = path1[state] - maxmax;
            }
        }
    }
    
    /* Final trace back (with trellis termination) */

    state = finalState;
    for (t = frameLength-1 ; t >= 0 ; t--) {
        decoded_frame[t] = trace_back[state][t][1];
        detected_code_frame[t] = trace_back[state][t][2];
        state = trace_back[state][t][0];
    }
    
    if (state != initState) {
        printf("- Error: trace-back does not recover initial state\n");
        exit(0);
    }
    
    /**
     * Free pointers
     */
    
  /*  free(path0);
    free(path1);
    for (i = 0; i < stateSize; ++i) {
        free(trace_back[i]);
    }
    free(trace_back);*/
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *encodedFrame, *bwd, *decodedFrame, *codedBits;
    int initState, finalState; 

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

    
    
    ccDecode(encodedFrame, bwd, initState, finalState, decodedFrame, codedBits);
}
