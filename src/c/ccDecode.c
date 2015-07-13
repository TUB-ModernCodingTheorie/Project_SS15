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
void ccDecode(  mxArray *encodedFrameArray,
                uint64_T initState,
                uint64_T finalState,
                mxArray *priorArray,
                mxArray *bwdArray,
                int m,
                int k,
                int n,
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
    uint64_T *bwd = (uint64_T*) mxGetData(bwdArray);
    double *prior = mxGetPr(priorArray);
	double *decodedFrame,*recodedFrame;

    uint64_T *traceBack = mxCalloc(stateSize*frameSize*3,sizeof(uint64_T)); 
	double *path0 = mxCalloc(stateSize, sizeof(double));
	double *path1 = mxCalloc(stateSize, sizeof(double));
	
	if (path0 == NULL || path1 == NULL)
		mexErrMsgIdAndTxt("CodingLibrary:ccDecode:noMemory","Could not allocate array. Not enough free memory.\n");
	
	if (traceBack == NULL)
		mexErrMsgIdAndTxt("CodingLibrary:ccDecode:noMemory","Could not allocate array. Not enough free memory.\n");
	
	(*decodedFrameArray) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
	(*recodedFrameArray) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
	decodedFrame = mxGetPr(*decodedFrameArray);
	recodedFrame = mxGetPr(*recodedFrameArray);
    
	
	/**
     * Number of unencoded input symbols per frame
     * For example: inputSize = 2
     *              encodedFrame = {1 2 0 3} = {01100011}
     *              frameLength = 4 = 8/2
     */

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
			for(b=0;b<inputSize;b++) {
				x     = bwd[s + b*stateSize + 2*stateSize*inputSize];
                old_s = bwd[s + b*stateSize];
                path  = path0[old_s];
				for(d=0;d<n;d++) {
                    i = i0 + d;
                    jj = (x >> d)&1;
					//printf("Prior: %f i: %d jj: %d\n",prior[i + jj*frameSize*n],i+1,jj+1);
                    path += prior[i + jj*frameSize*n];
                }
				if(path>=max) {
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
            if(max>=maxmax) maxmax = max;
		}
		for(s=0;s<stateSize;s++) path0[s] = path1[s] - maxmax;
	}

    /* Final trace back (with trellis termination) */
    s = finalState;
    for (t = frameSize-1 ; t >= 0 ; t--) {
        decodedFrame[t] = (double)traceBack[s + t*stateSize + 1*stateSize*frameSize];
        recodedFrame[t] = (double)traceBack[s + t*stateSize + 2*stateSize*frameSize];
        s = traceBack[s + t*stateSize];
    }
    
    if (s != initState) {
        printf("s doesn't match initial state\n");
		//mexErrMsgIdAndTxt("CodingLibrary:ccDecode:decodingError","Decoding error! The initial state doesn't match the traced back state.\n");
    }

    /**
     * Free pointers
     */
    mxFree(path0);
    mxFree(path1);
	mxFree(traceBack);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *encodedFrame;
    uint64_T initState, finalState;
	int n,m,k;
	mxArray *bwd;
	mxArray *Bmetric;
	mxArray *decodedFrame;  
    mxArray *recodedFrame;
    
    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:ccDecode:nrhs","Four inputs required.");
    }
	
	if(nlhs!=2) {
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
    
    
    bwd = mxGetFieldByNumber(prhs[0],0,0);
	m = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,1));
    n = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,2));
    k = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,3));
	
	encodedFrame = prhs[1];
	Bmetric = prhs[2];
    
    initState = (uint64_T) mxGetScalar(prhs[3]);
    finalState = (uint64_T) mxGetScalar(prhs[4]);
    
    ccDecode(encodedFrame, initState, finalState, Bmetric, bwd, m, k, n,
             &decodedFrame, &recodedFrame);
    
    plhs[0] = decodedFrame;
    plhs[1] = recodedFrame;
}

