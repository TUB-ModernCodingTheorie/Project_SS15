#include "mex.h"

#define INFINITY     1.e10

typedef struct
{
    void *next;
    void *previous;
    double *stateMetrics;
} stateMetrics_t;

/**
 * Viterbi Decoder
 * 
 * @param *bwdArray bwd with
 *          bwd[current_state][output][0] = previous_state
 *          bwd[current_state][output][1] = input uncoded/plain
 *          bwd[current_state][output][2] = output
 * @param *encodedFrame input encoded sequence
 * @param initState state at which we begin
 * @param finalState state in which we are supposed to finish with
 * @param m number of bits per state
 * @param k number of bits of the input symbols
 * @param n number of bits of the output symbols
 * @param **decodedSeqArray decoded sequence (output)
 * @param **recodedBitsArray encoded sequence going with the decoded one (output)
 */
void ccDecode(  const mxArray *bwdArray,
                const mxArray *encodedFrameArray,
                const mxArray *metricArray,
                const uint64_T initState,
                const uint64_T finalState,                
                const int m,
                const int k,
                const int n,
                mxArray **decodedFrameArray,
                mxArray **recodedFrameArray
            )
{
    int frameSize = (int) mxGetN(encodedFrameArray);
    int stateSize = 1 << m;
    int outputSize = 1 << n;
    int inputSize = 1 << k;
    
    int t,d,i0,i;
    uint64_T s,b,x,old_s,jj,survinput,survoutput,survstate;
    double path,max,maxmax;
    
    double *encodedFrame = mxGetPr(encodedFrameArray);
    double *metric = mxGetPr(metricArray);
    double *decodedFrame,*recodedFrame;

    uint64_T *bwd = (uint64_T*) mxGetData(bwdArray);

    /**
     * Allocate pointers memory
     */
    uint64_T *traceBack = mxCalloc(stateSize*frameSize*3,sizeof(uint64_T)); 
    double *path0 = mxCalloc(stateSize, sizeof(double));
    double *path1 = mxCalloc(stateSize, sizeof(double));
    
    if (path0 == NULL || path1 == NULL)
            mexErrMsgIdAndTxt("CodingLibrary:ccDecode:noMemory",
                              "Could not allocate array. Not enough free memory.");
    if (traceBack == NULL)
            mexErrMsgIdAndTxt("CodingLibrary:ccDecode:noMemory",
                              "Could not allocate array. Not enough free memory.");
    
    (*decodedFrameArray) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
    (*recodedFrameArray) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
    decodedFrame = mxGetPr(*decodedFrameArray);
    recodedFrame = mxGetPr(*recodedFrameArray);
    
    /**
     * 1. Initialization
     * 
     * All metrics are -infty exept the initial state
     */
		
    for (s = 0 ; s < stateSize ; s++)
        path0[s] = -INFINITY;
    path0[initState] = 0;

    /**
     * 2. ACS
     *
     * For each symbol of the frame
     *  For each possible state transition
     *   For each possible symbol
     */
    for (t = 0; t < frameSize; t++) {
        maxmax = -INFINITY;
        i0  = t*n;
        for (s = 0 ; s < stateSize ; s++) {
            max = -INFINITY;
            for (b = 0 ; b < inputSize ; b++) {
                x     = bwd[s + b*stateSize + 2*stateSize*inputSize];
                old_s = bwd[s + b*stateSize];
                path  = path0[old_s];
                for (d = 0 ; d < n ; d++) {
                    i = i0 + d;
                    jj = (x >> d)&1;
                    path += metric[i + jj*frameSize*n];
                }
                if (path >= max) {
                   survinput  = bwd[s + b*stateSize + 1*stateSize*inputSize];
                   survoutput = x;
                   survstate  = old_s;
                   max = path;
                }
            }
            traceBack[s + t*stateSize] = survstate;
            traceBack[s + t*stateSize + 1*stateSize*frameSize] = survinput;
            traceBack[s + t*stateSize + 2*stateSize*frameSize] = survoutput;
            path1[s] = max;
            if (max >= maxmax)
                maxmax = max;
        }
        for (s = 0 ; s < stateSize ; s++)
            path0[s] = path1[s] - maxmax;
    }

    /**
     * Final trace back (with trellis termination)
     */
    s = finalState;
    for (t = frameSize-1 ; t >= 0 ; t--) {
        decodedFrame[t] = (double)traceBack[s + t*stateSize + 1*stateSize*frameSize];
        recodedFrame[t] = (double)traceBack[s + t*stateSize + 2*stateSize*frameSize];
        s = traceBack[s + t*stateSize];
    }
    
    if (s != initState) {
	mexErrMsgIdAndTxt("CodingLibrary:ccDecode:decodingError","Decoding error! The initial state doesn't match the traced back state.\n");
    }

    /**
     * Free pointers
     */
    mxFree(path0);
    mxFree(path1);
    mxFree(traceBack);
}

/**
 * The gateway function
 * [m, c] = ccDecode(bwdTrellis, encodedFrame, metric, initState, finalState)
**/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    uint64_T initState, finalState;
    int n,m,k;
    mxArray *bwd;
    mxArray *decodedFrame;
    mxArray *recodedFrame;
    
    /** 
     * Check that the arguments are corrects
     */
    if (nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:nrhs",
                          "Five inputs required.");
    }
    if (nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:nlhs",
                          "Two outputs required.");
    }
    if (!mxIsScalar(prhs[3])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccDecode:notScalar",
                          "4th argument (initial state) must be a scalar.");
    }
    if (!mxIsScalar(prhs[4])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccDecode:notScalar",
                          "5th argument (final state) must be a scalar.");
    }
    if (!mxIsClass(prhs[2], "double")) {
        mexErrMsgIdAndTxt("CodingLibrary:ccDecode:notDouble",
                          "3rd argmument (metric) must be a matrix of double.");
    }
    if (!mxIsClass(prhs[0], "struct") || mxGetNumberOfFields(prhs[0]) != 4) {
        mexErrMsgIdAndTxt("CodingLibrary:ccDecode:notStruct",
                          "1st argument (backward structure) must be a structure with 4 fields");
    }
    if (mxGetM(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("CodingLibrary:ccDescode:notRowVector",
                          "2rd argument (encoded frame) must be a row vector.");
    }
    
    /**
     * Define variables and call the encoding function
     */
    bwd = mxGetFieldByNumber(prhs[0],0,0);
    m = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,1));
    n = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,2));
    k = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,3));
    
    initState = (uint64_T) mxGetScalar(prhs[3]);
    finalState = (uint64_T) mxGetScalar(prhs[4]);
    
    ccDecode(bwd, prhs[1], prhs[2], initState, finalState, m, k, n,
             &decodedFrame, &recodedFrame);
    
    /**
     * Define the outputs
     */
    plhs[0] = decodedFrame;
    plhs[1] = recodedFrame;
}
