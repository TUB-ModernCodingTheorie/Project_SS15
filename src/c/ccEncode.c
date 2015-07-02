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
            uint64_T ***fwd
        )
{
    *codeWord = fwd[*state][info][1];
    *state = fwd[*state][info][0];
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

void ccEncode(mxArray *fwd, mxArray *s0, mxArray *seq, mxArray **c, mxArray **sN)
{
    int i, rc = 0;

    uint64_T currState = 0;
    int stateSize = mxGetN(s0);
    int frameSize = mxGetN(seq);
    uint64_T *cData;
    uint64_T *uData;
    uint64_T *fwdTrellis;
    double *s0_ptr;

    (*sN) = mxDuplicateArray(s0);
    (*c) = mxCreateNumericMatrix(1, frameSize, mxUINT64_CLASS, 0);

    cData = (uint64_T*)mxGetData(*c);
    uData = (uint64_T*)mxGetData(seq);
    fwdTrellis = (uint64_T*)mxGetData(fwd);

    s0_ptr = mxGetPr(s0);
    for(i=0; i < stateSize; i++)
        currState ^= Bin2int((int)s0_ptr[i],i);

    for(i=0; i < frameSize; i++) {
		printf("state: %d data: %d\n",currState, uData[i]);
        BinEnc(&currState, &cData[i], uData[i], fwdTrellis);
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
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
    }
    
    /* make sure the first input argument is a string */
    if( !mxIsClass(prhs[0], "uint64") ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notUint64","Forward table must be uint64.");
    }
	if( !mxIsClass(prhs[1], "double") ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Initial state must be double.");
    }
	if( !mxIsClass(prhs[2], "uint64") ){
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notUint64","Input sequence must be uint64.");
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

    ccEncode(prhs[0], prhs[1], prhs[2], &c, &sN);

    plhs[0] = c;
    plhs[1] = sN;
}