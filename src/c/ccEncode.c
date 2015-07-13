#include "mex.h"
#include "tools.h"

/**
 * Actual encoders
 * 
 * @param *state current state (output: new state)
 * @param *codeword output codeword based on current state (output)
 * @param info input code
 * @param *fwd forward trellis
 * @param stateSize number of states
 * @param inputSize size of inputs symbols
 **/
void BinEnc( uint64_T *state,
             uint64_T *codeWord,
             uint64_T info,
             uint64_T *fwd,
             mwSize stateSize,
             mwSize inputSize
           )
{
    /**
     * codeword = fwd[currentState][input][1]
     * newState = fwd[currentState][input][0]
     */
    *codeWord = fwd[*state + info*stateSize + 1*(stateSize*inputSize)];
    *state = fwd[*state + info*stateSize + 0*(stateSize*inputSize)];
}

/**
 * Encoding function for convolutional codes
 * 
 * @param *fwd forward trellis
 * @param s0 initial state of the encoding
 * @param *seq input sequence to encode
 * @param m number of bits per state
 * @param k number of bits of the input symbols
 * @param n number of bits of the output symbols
 * @param **c codeword (output)
 * @param *sN final state (output)
 **/
void ccEncode(  const mxArray *fwd,
                const uint64_T s0,
                const mxArray *seq,
                const int m,
                const int k,
                const int n,
                mxArray **c,
                uint64_T *sN
             )
{
    int i;
    uint64_T currState = 0;
    uint64_T codeword;
    int frameSize = mxGetN(seq);
    double *cData;        /* pointer to the codeword *c */
    double *uData;        /* pointer to the inout sequence */
    uint64_T *fwdTrellis;
    int inputSize = (1 << k);
    int stateSize = (1 << m);
    
    (*c) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
    cData = mxGetPr(*c);
    uData = mxGetPr(seq);
    
    fwdTrellis = (uint64_T*)mxGetData(fwd);
    
    currState = s0;
    /**
     * For each symbol of the frame
     *  encode the current symbol with the trellis
     */
    for(i=0; i < frameSize; i++) {
        BinEnc(&currState, &codeword, (int)uData[i], fwdTrellis, stateSize, inputSize);
        cData[i] = (double)codeword;
    }
    
    (*sN) = currState; /* final state */
}

/**
 * The gateway function
 * [c, sN] = ccEncode(fwdStruct, s0, seq)
 **/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *c;     /* output codewords */
    uint64_T sN;    /* final state */
    mxArray *fwd;   /* forward trellis */
    int m,k,n;
    
    /** 
     * Check that the arguments are corrects
     */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:nrhs",
                          "Three inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:nlhs",
                          "Two outputs required.");
    }
    if (!mxIsClass(prhs[1], "double")
        || mxGetM(prhs[1]) != 1 
        || mxGetN(prhs[1]) != 1
    ) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notScalar",
                          "2nd argument (initial state) must be a scalar.");
    }
    if (!mxIsClass(prhs[2], "double")) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notDouble",
                          "3rd argmument (input sequence) must be double.");
    }
    if (!mxIsClass(prhs[0], "struct") || mxGetNumberOfFields(prhs[0]) != 4) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notStruct",
                          "1st argument (forward structure) must be a structure with 4 fields");
    }
    if (mxGetM(prhs[2]) != 1){
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notRowVector",
                          "3rd argument (input sequence) must be a row vector.");
    }
    
    /**
     * Define variables and call the encoding function
     */
    fwd = mxGetFieldByNumber(prhs[0],0,0);
    m = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,1));
    n = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,2));
    k = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,3));
    
    ccEncode(fwd, (uint64_T) mxGetScalar(prhs[1]), prhs[2], m, k, n, &c, &sN);
    
    /**
     * Define the outputs
     */
    plhs[0] = c;
    plhs[1] = mxCreateDoubleScalar((double)sN);
}
