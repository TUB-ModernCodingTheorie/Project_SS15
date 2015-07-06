/* Include header files */

#include "mex.h"

/* This yields the i-th binary digit of an integer (LSB is in position i = 0)*/
#define Int2bin(x_,i_) ((x_ >> i_)&1)

/* This is just a way of computing efficiently a*2^i, where a in {0,1} */
 #define Bin2int(x_,i_) ((x_) << (i_))

/* Encoder */

void BinEnc( uint64_T *state,
    uint64_T *codeWord,
    uint64_T info,
    uint64_T *fwd,
    mwSize fwdDim1,
    mwSize fwdDim2
    )
{
    *codeWord = fwd[*state + info*fwdDim1 + 1*(fwdDim1*fwdDim2)];
    *state = fwd[*state + info*fwdDim1 + 0*(fwdDim1*fwdDim2)];
}

void ccEncode(mxArray *fwd, const mxArray *s0, const mxArray *seq, int inputSize, int outputSize, mxArray **c, mxArray **sN)
{
    int i, j, rc = 0;
    uint64_T currState = 0;
    uint64_T codeWord, inputSeq;
    
    int stateSize = mxGetN(s0);
    int frameSize = mxGetN(seq);
    double *cData;
    double *uData;
    uint64_T *fwdTrellis;
    double *s0_ptr;
    mwSize *fwdSize = mxGetDimensions(fwd);
    
    (*sN) = mxDuplicateArray(s0);
    (*c) = mxCreateNumericMatrix(1, frameSize*outputSize, mxDOUBLE_CLASS, 0);
    
    cData = mxGetPr(*c);
    uData = mxGetPr(seq);
    fwdTrellis = (uint64_T*)mxGetData(fwd);

    s0_ptr = mxGetPr(s0);
    for(i=0; i < stateSize; i++)
        currState ^= Bin2int((int)s0_ptr[i],i);

    for(i=0; i < frameSize; i++) {
        inputSeq = 0;
        for(j=0; j<inputSize; j++)
            inputSeq ^= Bin2int((int)uData[i*inputSize + j],j);
        
        BinEnc(&currState, &codeWord, inputSeq, fwdTrellis, fwdSize[0], fwdSize[1]);
        
        for(j=0; j < outputSize ; j++)
            cData[i*outputSize + j] = Int2bin(codeWord,j);
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *c;     /* output codewords */
    mxArray *sN;    /* final state */
    mxArray *fwd;
    int k,n;
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
    }
    
    /*if (!mxIsClass(prhs[1], "double")) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Initial state must be double.");
    }
    if (!mxIsClass(prhs[2], "double")) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notUint64","Input sequence must be double.");
    }
    if (!mxIsClass(prhs[3], "double")) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input size must be double.");
    }
    if (!mxIsClass(prhs[4], "double")) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Output size must be double.");
    }
    if (mxGetNumberOfDimensions(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:dimensionMismatch","The forward table must have three dimensions.");
    }
    if (mxGetM(prhs[1]) != 1){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
    if (mxGetM(prhs[2]) != 1){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
    if (mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input size must be a single double.");
    }
    if (mxGetM(prhs[4]) != 1 || mxGetN(prhs[4]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Output size must be a single double.");
    }
*/
    fwd = mxGetFieldByNumber(prhs[0],0,0);
    n = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,2));
    k = (int) mxGetScalar(mxGetFieldByNumber(prhs[0],0,3));
    
    ccEncode(fwd, prhs[1], prhs[2], k, n, &c, &sN);

    plhs[0] = c;
    plhs[1] = sN;
}
