/* Include header files */

#include "mex.h"

/* This yields the i-th binary digit of an integer (LSB is in position i = 0)*/
#define Int2bin(x_,i_) ((x_ >> i_)&1)

/* This is just a way of computing efficiently a*2^i, where a in {0,1} */
 #define Bin2int(x_,i_) ((x_) << (i_))

/* Encoder */

int BinEnc( int *state,
            int *x,
            const int info,
            const double ***fwd,
            const int state_size,
            const int sizeK
        )
{
    if (*state >= state_size || info >= sizeK)
        return 0;
	
    *x = fwd[*state][info][1];
    *state = fwd[*state][info][0];

    return 1;	
}

/* Frame encoder */

/*
int FrameEncoder(int state)
{
    int i,ii,j,x,index,mask;
    
    mask = Asize - 1;

    for(i=0;i<Frame_length;i++) {
        x = BinEnc(&state,Info_frame[i]);
        Coded_frame[i] = x;
        for(j=0;j < Sdim;j++) {
            index = (x >> j*Mbit)&mask;
            ii = i*Sdim + j;
            XX[ii][0] = Signal_set[index][0];
            XX[ii][1] = Signal_set[index][1];
        }
    }
    return(state);
}
*/

void ccEncode(const mxArray *fwd, const mxArray *s0, const mxArray *seq, mxArray *c, mxArray *sN)
{
    int i, rc, x, s = 0;
    double *s0_ptr = mxGetPr(s0);
    double *seq_ptr = mxGetPr(seq);

    for(i=0; i < mxGetN(s0); i++)
        s ^= Bin2int(s0_ptr[i]==1,i);

    for(i=0; i < mxGetN(seq); i++) {
            rc = BinEnc(&s, &x, seq_ptr[i], mxGetPr(fwd), mxGetN(fwd), mxGetM(fwd));
            if(!rc) {
                mexErrMsgIdAndTxt("MyToolbox:ccEncode:encodeError","An encode error occured.");
                return;
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mxArray *c;	        /* output codewords */
    mxArray *sN;	/* final state */

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

    ccEncode(prhs[0], prhs[1], prhs[2], c, sN);

    plhs[0] = c;
    plhs[1] = sN;
}