#include "mex.h"

#define INFINITY     1.e10

/* This yields the i-th binary digit of an integer (LSB is in position i = 0)*/
#define Int2bin(x_,i_) ((x_ >> i_)&1)

/* This is just a way of computing efficiently a*2^i, where a in {0,1} */
 #define Bin2int(x_,i_) ((x_) << (i_))

 
void makeTrellis(const int *A, const int *B, const int *C, const int *D,
                 int m, int n, int k, uint64_T ***fwd, uint64_T ***bwd)
{

    printf("%d\n", A[1]);
    printf("%d\n", A[2]);
    printf("%d\n", A[3]);
    printf("%d\n", A[4]);
    
    int i,j,l;
    int numberOfState = 2 << (m-1);
    int numberOfInput = 2 << (k-1);
    int numberOfOutput = 2 << (n-1);
    int si;
    int tmp1, tmp2;

    fwd = mxCalloc(numberOfState, sizeof *fwd);
    for (i = 0; i < numberOfState; ++i) {
        fwd[i] = mxCalloc(numberOfInput, sizeof *fwd[i]);
        if (fwd[i] == NULL) {
            fprintf (stderr, "Memory allocation failure on fwd");
        }
        for (j = 0; j < numberOfInput ; j++) {
           fwd[i][j] = mxCalloc(2, sizeof *fwd[i][j]);
            if (fwd[i][j] == NULL) {
                fprintf (stderr, "Memory allocation failure on fwd");
            }
        }
    }
    bwd = mxCalloc(numberOfState, sizeof *bwd);
    for (i = 0; i < numberOfState; ++i) {
        bwd[i] = mxCalloc(numberOfOutput, sizeof *bwd[i]);
        if (bwd[i] == NULL) {
            fprintf (stderr, "Memory allocation failure on bwd");
        }
        for (j = 0; j < numberOfOutput ; j++) {
           bwd[i][j] = mxCalloc(2, sizeof *bwd[i][j]);
            if (bwd[i][j] == NULL) {
                fprintf (stderr, "Memory allocation failure on bwd");
            }
        }
    }
    
    for(l = 0 ; l < 2 ; l++){
        for(i = 0 ; i < numberOfState ; i++) {
            for(j = 0 ; j < numberOfOutput ; j++) {
                bwd[i][j][l] = -1;
            }
        }
    }

    printf("bwd is ok\n");
    for (si = 0 ; si < numberOfState ; si++) {
        int ui;
        for (ui = 0; ui < numberOfInput ; ui++) {
            int tmp1,tmp2;
            int next_state = 0;
            int output = 0;

            for(i = 0; i < m; i++){ /* number of columns */
                tmp1 = 0;
                tmp2 = 0;
                for(j = 0; j < m; j++){
                    tmp1 = (Int2bin(si,i)*A[m-1-j+i*m]+tmp1)%2;
                }
                for(j = 0; j < k; j++){
                    tmp2 = (Int2bin(ui,j)*B[k-1-j+i*m]+tmp2)%2;
                }

                next_state += Bin2int((tmp1+tmp2)%2,m-i-1);
            }
            for(i = 0; i < n; i++){
                tmp1 = 0;
                tmp2 = 0;
                for(j = 0; j < n; j++){
                    tmp1 = (Int2bin(si,j)*C[n-1-j+i*n]+tmp1)%2;
                }
                for(j = 0; j < k;j++){
                    tmp2 = (Int2bin(ui,j)*D[k-1-j+i*n]+tmp2)%2;
                }
                output += Bin2int((tmp1+tmp2)%2,n-i-1);
            }
            printf("%d, %d, %d ,%d\n",si,ui, output, next_state);
            fwd[si][ui][0] = next_state;
            fwd[si][ui][1] = output;
            bwd[next_state][output][0] = si;
            bwd[next_state][output][1] = ui;
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mxArray *fwdArray, *bwdArray;
    uint64_T ***fwd, ***bwd;
    int *A, *B, *C, *D;
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
    
    A = (int*) mxGetPr(prhs[0]);
    B = (int*) mxGetPr(prhs[1]);
    C = (int*) mxGetPr(prhs[2]);
    D = (int*) mxGetPr(prhs[3]);
    
    m = mxGetScalar(prhs[4]);
    n = mxGetScalar(prhs[5]);
    k = mxGetScalar(prhs[6]);
    
    
    
    makeTrellis(A, B, C, D, m, n, k, fwd, bwd);
    
    mxSetData(plhs[0], fwd);
    mxSetData(plhs[1], bwd);
    
}
