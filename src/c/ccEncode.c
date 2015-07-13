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
     * Structure of the fwd table:
     *  codeword = fwd[currentState][input][1]
     *  newState = fwd[currentState][input][0]
     */
    *codeWord = fwd[*state + info*stateSize + 1*(stateSize*inputSize)];
    *state = fwd[*state + info*stateSize + 0*(stateSize*inputSize)];
}

/**
 * Encoding function for convolutional codes
 * 
 * @param *fwd forward trellis
 *          fwd[current_state][input][0] = next_state
 *          fwd[current_state][input][1] = output
 * @param *seq input sequence to encode
 * @param initialState initial state of the encoding
 * @param m number of bits per state
 * @param k number of bits of the input symbols
 * @param n number of bits of the output symbols
 * @param **codeword codeword (output)
 * @param *finalState final state (output)
 **/
void ccEncode(  const mxArray *fwd,
                const mxArray *seq,
                const uint64_T initialState,
                const int m,
                const int k,
                const int n,
                mxArray **codeword,
                uint64_T *finalState
             )
{
    int i;
    uint64_T currState = 0;
    uint64_T codeword_tmp;
    int frameSize = mxGetN(seq);
    double *codewordData;
    double *seqData;
    uint64_T *fwdTrellis;
    int inputSize = (1 << k);
    int stateSize = (1 << m);
    
    (*codeword) = mxCreateNumericMatrix(1, frameSize, mxDOUBLE_CLASS, 0);
    codewordData = mxGetPr(*codeword);
    seqData = mxGetPr(seq);
    
    fwdTrellis = (uint64_T*)mxGetData(fwd);
    
    currState = initialState;
    /**
     * For each symbol of the frame
     *  encode the current symbol with the trellis
     */
    for(i=0; i < frameSize; i++) {
        BinEnc(&currState, &codeword_tmp, (int)seqData[i], fwdTrellis, stateSize, inputSize);
        codewordData[i] = (double)codeword_tmp;
    }
    
    (*finalState) = currState;
}

/**
 * The gateway function
 * [c, finalState] = ccEncode(fwdStruct, seq, initialState)
 **/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    uint64_T initState;
    int m,k,n;
    mxArray *fwd;           /* forward trellis */
    mxArray *codeword;      /* output codewords */
    uint64_T finalState;    /* final state */
    
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
    if (!mxIsScalar(prhs[2])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notScalar",
                          "3rd argument (initial state) must be a scalar.");
    }
    if (!mxIsClass(prhs[1], "double")) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notDouble",
                          "2nd argmument (input sequence) must be double.");
    }
    if (!mxIsClass(prhs[0], "struct") || mxGetNumberOfFields(prhs[0]) != 4) {
        mexErrMsgIdAndTxt("CodingLibrary:ccEncode:notStruct",
                          "1st argument (forward structure) must be a structure with 4 fields");
    }
    if (mxGetM(prhs[1]) != 1){
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
    
    initState = (uint64_T) mxGetScalar(prhs[2]);
    
    ccEncode(fwd, prhs[1], initState, m, k, n, &codeword, &finalState);
    
    /**
     * Define the outputs
     */
    plhs[0] = codeword;
    plhs[1] = mxCreateDoubleScalar((double)finalState);
}
