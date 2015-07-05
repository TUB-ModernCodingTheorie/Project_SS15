#include "mex.h"

/* This yields the i-th binary digit of an integer (LSB is in position i = 0)*/
#define Int2bin(x_,i_) ((x_ >> i_)&1)

/* This is just a way of computing efficiently a*2^i, where a in {0,1} */
#define Bin2int(x_,i_) ((x_) << (i_))

/*struct treillis_t
{
    int ldStates;
    int ldOutputs;
    int ldInputs;
    mxArray **forward;
    mxArray **backward;
}*/

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
    const mwSize dimsBwd[3] = {numberOfStates, numberOfOutputs, 2};
    
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
					bwdData[nextState + i*numberOfStates + 1*(numberOfStates*numberOfOutputs)] = fwdData[state + inp*numberOfStates + 1*(numberOfStates*numberOfInputs)];
				}
			}
		}
	}
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

 /*   struct trellis_t trellisData;*/
    const char *field_names[] = {"forward", "backward", "ldStates", "ldOutputs", "ldInputs"};
    mwSize structDims[2] = {1,1};
    int nbFields = (sizeof(field_names)/sizeof(*field_names));
    
    mxArray *fwdArray, *bwdArray;
    uint64_T *fwd, *bwd;
    double *A, *B, *C, *D;
    int m,n,k;
    int i,j;
    
    /**
     * @TODO:  check that everything is ok 
     * 
     * A should be m x m
     * B sould be k x m
     * C shoudl be n x n
     * D should be k x n
     * 
     * In mex a matrix is read [1 4 7
     *                          2 5 8
     *                          3 6 9]
     **/
    
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    D = mxGetPr(prhs[3]);
    
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[2]);
    k = mxGetM(prhs[1]);
    
    makeTrellis(A, B, C, D, m, n, k, &fwdArray, &bwdArray);

/*    trellis.forward = fwdArray;
    trellis.backward = bwdArray;
    trellis.ldStates = m;
    trellis.ldOutputs = n;
    trellis.ldInputs = k;*/
    
    plhs[0] = mxCreateStructArray(2, structDims, nbFields, field_names);
    mxSetFieldByNumber(plhs[0],0,0,fwdArray);
    mxSetFieldByNumber(plhs[0],0,1,bwdArray);
    mxSetFieldByNumber(plhs[0],0,2,mxCreateDoubleScalar(m));
    mxSetFieldByNumber(plhs[0],0,3,mxCreateDoubleScalar(n));
    mxSetFieldByNumber(plhs[0],0,4,mxCreateDoubleScalar(k));
    
/*    plhs[0] = fwdArray;
    plhs[1] = bwdArray;*/
    
}
