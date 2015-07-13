#include "mex.h"
#include "tools.h"

/**
 * Build the trellis
 * 
 * @param *A
 * @param *B
 * @param *C
 * @param *D
 * @param m number of rows of A = number of bits per states
 * @param n number of rows of B = number of bits per output symbol
 * @param k number of rows of C = number of bits per input symbol
 * @param **fwd forward trellis
 *          fwd[current_state][input][0] = next_state
 *          fwd[current_state][input][1] = output
 * @param **bwd backward trellis
 *          bwd[current_state][output][0] = previous_state
 *          bwd[current_state][output][1] = input uncoded/plain
 *          bwd[current_state][output][2] = output
 * 
 * @comment:
 *      In mex a matrix is read [1 4 7
 *                               2 5 8
 *                               3 6 9]
 **/
void makeTrellis(double *A, double *B, double *C, double *D,
                 int m, int n, int k, mxArray **fwd, mxArray **bwd)
{
    int i,j,l;
    uint64_T *fwdData, *bwdData;
    uint64_T state, inp;
    uint64_T codeword, nextState;
    uint64_T numberOfStates = 2 << (m-1);
    uint64_T numberOfInputs = 2 << (k-1);
    uint64_T numberOfOutputs = 2 << (n-1);
    const mwSize dimsFwd[3] = {numberOfStates, numberOfInputs, 2};
    const mwSize dimsBwd[3] = {numberOfStates, numberOfInputs, 3};
    
    (*fwd) = mxCreateNumericArray(3, dimsFwd, mxUINT64_CLASS, 0);
    (*bwd) = mxCreateNumericArray(3, dimsBwd, mxUINT64_CLASS, 0);
    fwdData = (uint64_T*)mxGetData(*fwd);
    bwdData = (uint64_T*)mxGetData(*bwd);
    
    for (state = 0; state < numberOfStates; state++) {
        for (inp = 0; inp < numberOfInputs ; inp++) {
            codeword = 0;
            nextState = 0;
            
            for(i = 0; i < m; i++) {
                for(j = 0; j < m; j++) {
                    nextState ^= Bin2int(Int2bin(state, j)*(int)A[j + i*m], i);
                }
                for(j = 0; j < k; j++) {
                    nextState ^= Bin2int(Int2bin(inp, j)*(int)B[j + i*k], i);
                }
            }
            for(i = 0; i < n; i++) {
                for(j = 0; j < m; j++) {
                    codeword ^= Bin2int(Int2bin(state, j)*(int)C[j + i*m], i);
                }
                for(j = 0; j < k; j++) {
                    codeword ^= Bin2int(Int2bin(inp, j)*(int)D[j + i*k], i);
                }
            }
            
            fwdData[state + inp*numberOfStates + 1*(numberOfStates*numberOfInputs)] = codeword;
            fwdData[state + inp*numberOfStates + 0*(numberOfStates*numberOfInputs)] = nextState;
        }
    }
    
    for (nextState = 0; nextState < numberOfStates; nextState++) {
        i = 0;
        for (inp = 0; inp < numberOfInputs ; inp++) {
            for (state = 0; state < numberOfStates; state++) {
                if (nextState == fwdData[state + inp*numberOfStates]) {
                    bwdData[nextState + i*numberOfStates] = state;
					bwdData[nextState + i*numberOfStates + 1*(numberOfStates*numberOfInputs)] = inp;
                    bwdData[nextState + i*numberOfStates + 2*(numberOfStates*numberOfInputs)] = fwdData[state + inp*numberOfStates + 1*(numberOfStates*numberOfInputs)];
					i++;
                }
            }
        }
    }
}

/**
 * The gateway function
 * [forwardStruct, backwardStruct] = ccInitialize(A,B,C,D)
 **/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    const int nbFields = 4;
    const char *field_names[] = {"trellis", "ldStates", "ldOutputs", "ldInputs"};
    mwSize structDims[2] = {1,1};
    
    mxArray *fwdArray, *bwdArray;
    uint64_T *fwd, *bwd;
    double *A, *B, *C, *D;
    int m,n,k;
    int i,j;
    
    /** 
     * Check that the arguments are corrects
     * 
     * A should be m x m
     * B should be k x m
     * C should be n x n
     * D should be k x n 
     **/
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:nrhs",
                          "Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:nlhs",
                          "Two outputs required.");
    }
    if (!mxIsClass(prhs[0], "double")
        || !mxIsClass(prhs[1], "double")
        || !mxIsClass(prhs[2], "double")
        || !mxIsClass(prhs[3], "double")
    ) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:notDouble",
                          "All input must be double.");
    }
    if (mxGetM(prhs[0]) != mxGetN(prhs[0])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:wrongSize",
                          "1st argument (A) must be mxm matrix.");
    }
    if (mxGetM(prhs[2]) != mxGetN(prhs[2])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:wrongSize",
                          "1st argument (C) must be nxn matrix.");
    }
    if (mxGetM(prhs[0]) != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:wrongSize",
                          "A should be of size mxm and C kxm. Check the m value");
    }
    if (mxGetM(prhs[1]) != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:wrongSize",
                          "B should be kxm and D kxn. Check the k value.");
    }
    if (mxGetM(prhs[2]) != mxGetN(prhs[3])) {
        mexErrMsgIdAndTxt("CodingLibrary:ccInitialize:wrongSize",
                          "C should be nxn and D kxn. Check the n value.");
    }
    
    /**
     * Define variables and call the encoding function
     */
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    D = mxGetPr(prhs[3]);
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[2]);
    k = mxGetM(prhs[1]);
    
    makeTrellis(A, B, C, D, m, n, k, &fwdArray, &bwdArray);
    
    /**
     * Define the outputs
     */
    plhs[0] = mxCreateStructArray(2, structDims, nbFields, field_names);
    plhs[1] = mxCreateStructArray(2, structDims, nbFields, field_names);
    
    mxSetFieldByNumber(plhs[0],0,0,fwdArray);                /* forward */
    mxSetFieldByNumber(plhs[0],0,1,mxCreateDoubleScalar(m)); /* ldStates */
    mxSetFieldByNumber(plhs[0],0,2,mxCreateDoubleScalar(n)); /* ldOutputs */
    mxSetFieldByNumber(plhs[0],0,3,mxCreateDoubleScalar(k)); /* ldInputs */
    
    mxSetFieldByNumber(plhs[1],0,0,bwdArray);                /* backward */
    mxSetFieldByNumber(plhs[1],0,1,mxCreateDoubleScalar(m)); /* ldStates */
    mxSetFieldByNumber(plhs[1],0,2,mxCreateDoubleScalar(n)); /* ldOutputs */
    mxSetFieldByNumber(plhs[1],0,3,mxCreateDoubleScalar(k)); /* ldInputs */
}
