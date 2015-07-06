#include "mex.h"

#define INFINITY     1.e10


/**
 * Viterbi Decoder
 * 
 * @param *encodedFrame input encoded sequence
 * @param *bwdArray bwd with
 *          bwd[current_state][idx][0] = previous_state
 *          bwd[current_state][idx][1] = input uncoded/plain bit
 *          bwd[current_state][idx][2] = output bit
 * @param initState state at which we begin
 * @param finalState state in which we are supposed to finish with
 * @param *decodedSeqArray decoded sequence (output)
 * @param *recodedBitsArray encoded sequence going with the decoded one (output)
 * 
 */
void ccDecode(  mxArray *encodedFrame,
                int initState,
                int finalState,
                mxArray *priorArray,
                mxArray *bwdArray,
                int stateSize,
                int inputSize,
                int outputSize,
                mxArray **decodedFrameArray,
                mxArray **recodedFrameArray
            )
{
    int max, maxmax, survoutput, survinput,survstate,old_state,mask,Mbit,i,i0,sdim,b,jj,d,x,t, state;
    double M; /* Metric */
   
    uint64_T *encodedFrame = (uint64_T*) mxGetData(encodedFrameArray);
    uint64_T *bwd = (uint64_T*) mxGetData(bwdArray);
    double *prior = (double*) mxGetData(priorArray);
    
    printf("\ninput: %d\nstate: %d\n", inputSize, stateSize);
    /**
     * Number of unencoded input symbols per frame
     * For example: inputSize = 2
     *              encodedFrame = {1 2 0 3} = {01100011}
     *              frameLength = 4 = 8/2
     */
    int frameLength = (int) mxGetM(encodedFrameArray);

    double M = calloc(stateSize, sizeof(double));
    double *path1 = calloc(stateSize, sizeof(double));
    uint64_T *decodedFrame  = calloc(frameLength,sizeof(uint64_T));
    uint64_T *recodedFrame  = calloc(frameLength,sizeof(uint64_T));
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
    for (state = 0 ; state < stateSize ; state++) {
        M[state] = -INFINITY;
    }
    M[init_state] = 0;


    /**
     * 2. ACS
     *
     * For each symbol of the frame
     *  For each possible state transition
     *   For each possible symbol
     */
    for (t = frameLength-1 ; t >= 0 ; t--) {
        for (state = 0 ; state < stateSize ; state++) {
            max = -INFINITY;
            i0 = t*outputSize;
            for (b = 0 ; b < inputSize ; b++) {
                x = bwd[state][b][2];
                old_state = bwd[state][b][0];
                M_tmp = M[old_state];
                for (d = 0 ; d < outputSize ; d++) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    M_tmp += prior[i + jj*outputSize] /* P(i | jj) */
                }
                if(M_path>=max) {
                    survinput  = bwd[state][b][1];
                    survoutput = x;
                    survstate  = old_state;
                    max = M_tmp;
                }
            }
            trace_back[state][t][0] = survstate;
            trace_back[state][t][1] = survinput;
            trace_back[state][t][2] = survoutput;
            path1[state] = max;
            
            if (max >= maxmax)
                maxmax = max;
        }
        for(state = 0 ; state < stateSize ; state++)
            M[state] = path1[state] - maxmax;
    }

    /* Final trace back (with trellis termination) */
    state = finalState;
    for (t = frameLength-1 ; t >= 0 ; t--) {
        decodedFrame[t] = trace_back[state][t][1];
        recodedFrame[t] = trace_back[state][t][2];
        state = trace_back[state][t][0];
    }
    
    if (state != initState) {
        printf("- Error: trace-back does not recover initial state\n");
        exit(0);
    }

    /**
     * Free pointers
     */
    
    mxFree(path0);
    mxFree(path1);
    for (i = 0; i < stateSize; ++i)
    {
        free(trace_back[i]);
    }
    free(trace_back);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *encodedFrame;
    int initState, finalState;
    
    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:nrhs","Four inputs required.");
    }
    
    /* make sure the first input argument is a string 
    if( !mxIsClass(prhs[0], "uint64") ){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:notUint64","Encoded frame must be uint64.");
    }
	if( !mxIsClass(prhs[1], "uint64") ){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:notUint64","Backward table must be uint64.");
    }
	if( !mxIsScalar(prhs[2]) ){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:notScalar","Initial state must be an integer.");
    }
        if( !mxIsScalar(prhs[3]) ){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:notScalar","Final state must be an integer.");
    }

    if( mxGetNumberOfDimensions(prhs[1]) != 3){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:dimensionMismatch","The bwd table must have three dimensions.");
    }
    if( mxGetM(prhs[0]) != 1){
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:notRowVector","Encoded frame must be a row vector.");
    }*/
    
    encodedFrame = mxGetPr(prhs[0]);
    
    mxArray *bwd = prhs[1];
    
    initState = (int) mxGetScalar(prhs[1]);
    finalState = (int) mxGetScalar(prhs[2]);
    
    mxArray *decodedFrame;  
    mxArray *recodedFrame;
    
    printf("0");
    ccDecode(encodedFrame,
             initState,
             finalState,
             &bwd,
             stateSize,
             inputSize,
             &decodedFrame, 
             &recodedFrame
            );
    printf("0");
    
    decodedFrame = mxCreateDoubleScalar(3);
    plhs[0] = decodedFrame;
    plhs[1] = recodedFrame;
    printf("FIN");
}

