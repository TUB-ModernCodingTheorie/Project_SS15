// Include header files 

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "mex.h"

// Numerical constants 

#define INFINITY     1.e10
#define LN2          0.69314718
#define LN10         2.3025851
#define SQRT2        1.41421356
#define PIG          3.14159265
#define TINY         1.e-10
#define BIG          36.736
#define SMALL        1.67e-16

// Memory allocation functions

int   chk_bytes = 0;
size_t Alloc_array_size[20];

void *Check_calloc (size_t nitems, size_t size)
{
  void *tmp;
  tmp = calloc (nitems, size);
  if (tmp == NULL) {
    fprintf (stderr, "Memory allocation failure (%d items of size %d bytes)\n"
                   , nitems, size);
    fprintf (stderr, "Allocated memory before this call: %d bytes\n", chk_bytes);
    exit (0);
  }
  chk_bytes += nitems*size;
  return tmp;
}

void *Alloc_array1 (size_t elem_size, int dim)
{
  if (dim == 1) {
    return Check_calloc (Alloc_array_size[0], elem_size);
  } else {
    unsigned int i; void **ret;

    assert (dim < sizeof (Alloc_array_size));
    dim--;
    ret = Check_calloc (Alloc_array_size[dim], sizeof (void **));
    for (i = 0; i < Alloc_array_size[dim]; i++) {
      ret[i] = Alloc_array1 (elem_size, dim);
    }
    return ret;
  }
}

// NOTICE: THIS IS THE ONE TO USE, BUT IT NEEDS THE ABOVE TWO FUNCTIONS IN ORDER TO WORK

void *Alloc_array (size_t elem_size, int dim, ...)
{
  int i;
  va_list dimpt;

  va_start (dimpt, dim);
  for (i = 0; i < dim; i++) {
    assert ((Alloc_array_size[dim-1-i] = va_arg (dimpt, size_t)) > 0);
  }
  va_end (dimpt);
  assert (dim > 0);
  return Alloc_array1 (elem_size, dim);
}

// ****************************************************************************
// Input from file functions

// Read a vector of integers from file, this is given as a sequence of integers, on the same line, separated by spaces

int GetLongVector (FILE *ff, int n, int *vector)
{
  int i;
  char line[1024], *pt;

  if (n <= 0) return 1;
  fgets (line, sizeof (line)-1, ff);
  pt = strtok (line, " ,\t");
  for (i = 0; i < n; i++) {
    if (!pt || sscanf (pt, "%d", &vector[i]) != 1) return 1;
    pt = strtok (NULL, " ,\t");
  }
  return 0;
}

// Read a matrix of integers from file, this is given as a multi-row table, each row is a vector as defined before

int GetLongMatrix (FILE *ff, int nrow, int ncol, int **matrix)
{
  int i, j;
  char line[1024], *pt;

  if (nrow <= 0 || ncol <= 0) return 1;
  for (i = 0; i < nrow; i++) {
    fgets (line, sizeof (line)-1, ff);
    pt = strtok (line, " ,\t");
    for (j = 0; j < ncol; j++) {
      if (!pt || sscanf (pt, "%d", &matrix[i][j]) != 1) return 1;
      pt = strtok (NULL, " ,\t");
    }
  }
  return 0;
}

// Read a vector of doubles, same as before, for doubles

int GetDoubleVector (FILE *ff, int n, double *vector)
{
  int i;
  char line[1024], *pt;

  if (n <= 0) return 1;
  fgets (line, sizeof (line)-1, ff);
  pt = strtok (line, " ,\t");
  for (i = 0; i < n; i++) {
    if (!pt || sscanf (pt, "%lf", &vector[i]) != 1) return 1;
    pt = strtok (NULL, " ,\t");
  }
  return 0;
}

// Read a matrix of doubles, same as before, for doubles

int GetDoubleMatrix (FILE *ff, int nrow, int ncol, double **matrix)
{
  int i, j;
  char line[1024], *pt;

  if (nrow <= 0 || ncol <= 0) return 1;
  for (i = 0; i < nrow; i++) {
    fgets (line, sizeof (line)-1, ff);
    pt = strtok (line, " ,\t");
    for (j = 0; j < ncol; j++) {
      if (!pt || sscanf (pt, "%lf", &matrix[i][j]) != 1) return 1;
      pt = strtok (NULL, " ,\t");
    }
  }
  return 0;
}

// ****************************************************************************
// Random generators

#define MODULE       2147483647
#define FACTOR       16807
#define LASTXN       127773
#define UPTOMOD     -2836

// Uniform [0,1] 

static int Random=0; 
// Notice: Random is a **global variable**, if it is set to zero, it initializes by reading the time, otherwise, it contains the initial seed of the generator

double Uniform (void)
{
  static int  first_call=1;
  static int seed;
  static int times, rest, prod1, prod2;

  if(first_call) {
    if(Random<=0) {
      seed = time(NULL)%100000000;
    } else {
      seed = Random;
    }

    /* printf("RANDOM GENERATOR\n");
    printf("Seed initial value   : %d\n",seed);
    printf("-------------------------------------------------\n"); */

    first_call = 0;
  }
  times    = seed / LASTXN;
  rest     = seed - times * LASTXN;
  prod1    = times * UPTOMOD;
  prod2    = rest * FACTOR;
  seed     = prod1 + prod2;
  if (seed < 0) seed = seed + MODULE;
  return (double) seed / MODULE;
}

// Uniform on the integers {0,q-1}

int Bitq (int q)
{
  int a;

  do {
    a = (int) (q * Uniform());
  } while (a>=q);

  return(a);
}

// Normal bivariate zero-mean with covariance diag(1/2,1/2)

void Gaussian2 (double *x1, double *x2)
{
  double u1, u2, s;

  do {
    u1 = Uniform ();
  } while (u1<=0);
  u2  = PIG*(2*Uniform () - 1);
  s   = sqrt(-log(u1));
  *x1 = s*cos(u2);
  *x2 = s*sin(u2);
}

//*****************************************************************************************
// Miscellanea functions of eneral utility

double Minimum(
       double a,
       double b
       )
{
  if(a<=b) return(a);
  else     return(b);
}

double Maximum(
       double a,
       double b
       )
{
  if(a>=b) return(a);
  else     return(b);
}

// Notice: this is an efficient and numerically stable way of computing log-sum

double LogSum(
       double a,
       double b
       )
{
  double d;

  d = b - a;
  if(d<=-4.0)      return(a);
  else if (d>=4.0) return(b);
  else if(d<=0)    return(a + log(1 + exp(d)));
  else             return(b + log(1 + exp(-d)));
}

// This yields the i-th binary digit of an integer (LSB is in position i = 0)
//Baris: more efficient way
#define Int2bin(x_,i_) ((x_ >> i_)&1)
/*int Int2bin(
    int a,
    int i
    )
{
  return((a >> i)&1);
}*/

// This is just a way of computing efficiently a*2^i, where a in {0,1}
//Baris: more efficient way
 #define Bin2int(x_,i_) (x_ << i_)
/*int Bin2int(
    int a,
    int i
    )
{
  return((a << i));
}*/

// This function computes the bit-reversal of an integer (i.e., exchanging MSB with LSB and so on)

//Baris: this function is unused
/*int Reverse(a,n)
{
	int i,s;

	s = 0;
	for(i=0;i<n;i++) s += Bin2int(Int2bin(a,i),n-i-1);
	return(s);
}*/

//****************************************************************************************************
// Convolutional encoder over GF(2)

static int   Kbit,
             Nbit,
			 Sdim,
			 Mbit,
             Memory,
             State_size,
             Input_size,
             Output_size,
			 Asize,
           **Generator,
          ***Forward,
          ***Backward;

static int Code_length,
           Frame_length,
		  *Info_frame,
          *Coded_frame;

static double **XX,
              **Prior,
              **Signal_set;

double BinEnc_init(char *ccfilename, char *ssfilename)
{

    int s,next_s,a,b,x,i,j,mask;

	double Rate,A;

	char line[256];

	FILE *file;

	// Read CC encoder definition file

    file = fopen(ccfilename,"r");
    if(file == NULL) {
       printf("- ERROR: File not found -\n");
       exit(-1);
	}

    fgets(line,sizeof(line)-1,file);
    sscanf(line,"%d %d %d",&Nbit,&Kbit,&Memory);

    Generator = (int **) Alloc_array(sizeof(int),2,Nbit+Memory,Kbit+Memory);

    i = GetLongMatrix(file,Nbit+Memory,Kbit+Memory,Generator);
    if(i) {
       printf("- ERROR: error reading the generator matrix\n");
	   exit(-1);
	}

	fclose(file);

    Input_size  = (1 << Kbit);
    Output_size = (1 << Nbit);
    State_size  = (1 << Memory);

// Trellis matrices

    Forward   = (int ***) Alloc_array(sizeof(int),3,State_size,Input_size,3);
    Backward  = (int ***) Alloc_array(sizeof(int),3,State_size,Input_size,3);

    for(s=0;s<State_size;s++) {
      for(b=0;b<Input_size;b++) {
        x = 0;
        for(i=0;i<Nbit;i++) {
          a = 0;
          for(j=0;j<Kbit;j++)    a ^= Generator[i][j]&Int2bin(b,j);
          for(j=0;j<Memory;j++)  a ^= Generator[i][j+Kbit]&Int2bin(s,j);
          x += Bin2int(a,i);
        }
        next_s = 0;
        for(i=0;i<Memory;i++) {
          a = 0;
          for(j=0;j<Kbit;j++)    a ^= Generator[i+Nbit][j]&Int2bin(b,j);
          for(j=0;j<Memory;j++)  a ^= Generator[i+Nbit][j+Kbit]&Int2bin(s,j);
          next_s += Bin2int(a,i);
        }
        Forward[s][b][0] = next_s;
        Forward[s][b][1] = b;
        Forward[s][b][2] = x;
      }
    }

    for(next_s=0;next_s<State_size;next_s++) {
      i = 0;
      for(b=0;b<Input_size;b++) {
        for(s=0;s<State_size;s++) {
          if(Forward[s][b][0]==next_s) {
            Backward[next_s][i][0] = s;
            Backward[next_s][i][1] = b;
            Backward[next_s][i][2] = Forward[s][b][2];
            i += 1;
          }
        }
      }
    }
	
    // This part is a check .. uncomment if you want to output the trellis matrices
	
    /* printf("Forward trellis array :\n");
    for(i=0;i<State_size;i++) {
      for(j=0;j<Input_size;j++) printf("%d,%d,%d ",Forward[i][j][0],Forward[i][j][1],Forward[i][j][2]);
      printf("\n");
    }
    printf("Backward trellis array :\n");
    for(i=0;i<State_size;i++) {
      for(j=0;j<Input_size;j++) printf("%d,%d,%d ",Backward[i][j][0],Backward[i][j][1],Backward[i][j][2]);
      printf("\n");
    }
    exit(0);*/ 
	
  // Read the Signal set definition file

  file = fopen(ssfilename,"r");
  if(file == NULL) {
     printf("- ERROR: Signal file not found -\n");
     exit(0);
  }
  fgets(line,sizeof(line)-1,file);
  sscanf(line,"%d",&Mbit);

  Asize = (1 << Mbit);
  Signal_set = (double **) Alloc_array(sizeof(double),2,Asize,2);

  i = GetDoubleMatrix(file,Asize,2,Signal_set);
  if(i) {
     printf("-- WARNING: error reading the signal set\n");
  }
  mask = Asize - 1;
  fclose(file);

  A = 0;
  for(i=0;i<Asize;i++) A += Signal_set[i][0]*Signal_set[i][0] + Signal_set[i][1]*Signal_set[i][1];
  A = sqrt(A/Asize);
  for(i=0;i<Asize;i++) {
      Signal_set[i][0] = Signal_set[i][0]/A;
      Signal_set[i][1] = Signal_set[i][1]/A;
  }

  Sdim = Nbit/Mbit;
  if(Sdim*Mbit!=Nbit) {
        printf("ERROR: Nbit is not a multiple of Mbit\n");
        exit(0);
  }
  
  Code_length = Frame_length*Sdim;
  
  Info_frame = (int *) Alloc_array(sizeof(int),1,Frame_length);
  Coded_frame = (int *) Alloc_array(sizeof(int),1,Frame_length);
  XX = (double **) Alloc_array(sizeof(double),2,Code_length,2);
  Prior = (double **) Alloc_array(sizeof(double),2,Code_length,Asize);
    
  Rate = (double) Kbit/ (double) (Nbit/Mbit);
  return(Rate);
  
}

// Encoder

int BinEnc(
     int  *state,
     int   info
     )
{
    int x;
	
    x = Forward[*state][info][2];
    *state = Forward[*state][info][0];
    return(x);	
}

// Frame encoder 

int FrameEncoder(
     int state
     )
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

// AWGN channel and receiver front-end 

void AWGNChannel(
     double amp
     )
{
  double x,y,w1,w2;

  int i,j;

      for(j=0;j<Code_length;j++) {
          Gaussian2(&w1,&w2);
          x = amp*XX[j][0] + w1;
          y = amp*XX[j][1] + w2;
          for(i=0;i<Asize;i++) Prior[j][i] = -(x - amp*Signal_set[i][0])*(x - amp*Signal_set[i][0]) - (y - amp*Signal_set[i][1])*(y - amp*Signal_set[i][1]);
      }
}


// Viterbi algorithm 

static int  ***Trace_back,
              *Detected_info_frame,
              *Detected_code_frame;

static double    *Path0,
                 *Path1;


void  Viterbi_init(void)
{

    Path0      = (double *) Alloc_array(sizeof(double),1,State_size);
    Path1      = (double *) Alloc_array(sizeof(double),1,State_size);
    Trace_back = (int ***) Alloc_array(sizeof(int),3,State_size,Frame_length,3);
    Detected_info_frame  = (int *) Alloc_array(sizeof(int),1,Frame_length);
    Detected_code_frame  = (int *) Alloc_array(sizeof(int),1,Frame_length);
}


int  Viterbi(
     int init_state,
     int final_state
     )
{
  int    t,s,x,i,i0,jj,d,survstate,survinput,survoutput,odd,even,old_s,b,mask;
  double path,max,maxmax;


    for(s=0;s<State_size;s++) Path0[s] = -INFINITY;
    Path0[init_state] = 0;

    mask = Asize - 1;
	
    for(t=0;t<Frame_length;t++) {
        maxmax = -INFINITY;
        odd = t%2;
        even = odd^1;
        i0  = t*Sdim;
        for(s=0;s<State_size;s++) {
            max = -INFINITY;
            for(b=0;b<Input_size;b++) {
                x     = Backward[s][b][2];
                old_s = Backward[s][b][0];
                path  = Path0[old_s];
                for(d=0;d<Sdim;d++) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    path += Prior[i][jj];
                }
                if(path>=max) {
                   survinput  = Backward[s][b][1];
                   survoutput = x;
                   survstate  = old_s;
                   max = path;
                }
            }
            Trace_back[s][t][0] = survstate;
            Trace_back[s][t][1] = survinput;
			Trace_back[s][t][2] = survoutput;
            Path1[s] = max;
            if(max>=maxmax) maxmax = max;
        }
        for(s=0;s<State_size;s++) Path0[s] = Path1[s] - maxmax;
    }

    // Final trace back (with trellis termination)

    s = final_state;
    for(t=Frame_length-1;t>=0;t--) {
        Detected_info_frame[t] = Trace_back[s][t][1];
        Detected_code_frame[t] = Trace_back[s][t][2];
        s = Trace_back[s][t][0];
    }
    return(s);
	
}

// BCJR algorithm 

static double    *Beta0,
                 *Beta1,
                **Alpha,
                **APPout,
                 *APPin;

void BCJR_init(void)
{

    Beta0      = (double *) Alloc_array(sizeof(double),1,State_size);
    Beta1      = (double *) Alloc_array(sizeof(double),1,State_size);
    Alpha      = (double **) Alloc_array(sizeof(double),2,State_size,Frame_length+1);
    APPout     = (double **) Alloc_array(sizeof(double),2,Code_length,Asize);
    APPin      = (double *) Alloc_array(sizeof(double),1,Input_size);
    Detected_info_frame  = (int *) Alloc_array(sizeof(int),1,Frame_length);
    Detected_code_frame  = (int *) Alloc_array(sizeof(int),1,Frame_length);
	
}

void BCJR(
     int init_state,
     int final_state
     )
{
  int    t,s,x,i,i0,jj,d,old_s,b,in,mask;
  double path,max,maxmax,branch;


    for(s=0;s<State_size;s++) Alpha[s][Frame_length] = -INFINITY;
    Alpha[final_state][Frame_length] = 0;

    mask = Asize - 1;
	
    for(t=Frame_length-1;t>=0;t--) {
        maxmax = -INFINITY;
        i0  = t*Sdim;
        for(s=0;s<State_size;s++) {
            max = -INFINITY;
            for(b=0;b<Input_size;b++) {
                x     = Forward[s][b][2];
                old_s = Forward[s][b][0];
                path  = Alpha[old_s][t+1];
                for(d=Sdim-1;d>=0;d--) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    path += Prior[i][jj];
                }
                max = LogSum(max,path);
            }
            Alpha[s][t] = max;
            maxmax = LogSum(maxmax,max);
        }
        for(s=0;s<State_size;s++) Alpha[s][t] -= maxmax;
    }

    for(s=0;s<State_size;s++) Beta0[s] = -INFINITY;
    Beta0[init_state] = 0;

    for(t=0;t<Frame_length;t++) {
        maxmax = -INFINITY;
        i0  = t*Sdim;
        for(d=0;d<Sdim;d++) {
            for(i=0;i<Asize;i++) APPout[i0+d][i] = -INFINITY;
        }
        for(i=0;i<Input_size;i++) APPin[i] = -INFINITY;

        for(s=0;s<State_size;s++) {
            max = -INFINITY;
            for(b=0;b<Input_size;b++) {
                x     = Backward[s][b][2];
                in    = Backward[s][b][1];
                old_s = Backward[s][b][0];
                path  = 0;
                for(d=0;d<Sdim;d++) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    path += Prior[i][jj];
                }
                max = LogSum(max,path+Beta0[old_s]);
                branch = Beta0[old_s] + path + Alpha[s][t+1];
                for(d=0;d<Sdim;d++) {
                    i = i0 + d;
                    jj = (x >> d*Mbit)&mask;
                    // APPout[i][jj] = LogSum(APPout[i][jj],branch - Prior[i][jj]); 
                   APPout[i][jj] = LogSum(APPout[i][jj],branch);
                }
                APPin[in] = LogSum(APPin[in],branch);
            }
            Beta1[s] = max;
            maxmax = LogSum(maxmax,max);
        }
        
		for(s=0;s<State_size;s++) Beta0[s] = Beta1[s] - maxmax;
        
		x = 0;
        for(d=0;d<Sdim;d++) {
            max = -INFINITY;
            maxmax = -INFINITY;
            for(i=0;i<Asize;i++) {
                max = LogSum(max,APPout[i0+d][i]);
                if(APPout[i0+d][i] > maxmax) {
                   maxmax = APPout[i0+d][i];
                   jj = i;
                }
            }
            x += (jj << d*Mbit);
            for(i=0;i<Asize;i++) APPout[i0+d][i] -= max;
        }
        Detected_code_frame[t] = x;
        max = -INFINITY;
        for(b=0;b<Input_size;b++) {
            if(APPin[b]>max) {
               max = APPin[b];
               jj = b;
            }
        }
        Detected_info_frame[t] = jj;
    }
}


// -----------------------------------------------------------------
// -----------------------------------------------------------------
// -----------------------------------------------------------------
// -----------------------------------------------------------------
// CONVOLUTIONAL CODE SIMULATOR on BI-AWGN
// -----------------------------------------------------------------
// -----------------------------------------------------------------
// -----------------------------------------------------------------
// -----------------------------------------------------------------

#define MIN_BIT_ERRORS 500
#define MIN_FRAME_ERRORS 50

main (int argc, char *argv[])
{
  int max_frames,frame_error_cnt,symb_error_cnt,bit_error_cnt,
      frame_cnt,frame,init_state,final_state,check_state,i,j,a,b,out,mask,verbose;

  double EbN0_init,EbN0,SNR,EbN0_step,EbN0_final,Rate,A,ps,pb,pf;

  char  output_file_name[256],encoder_file_name[256],
        signal_file_name[256],line[256];

  FILE *file;

// Get input parameters

    if(argc < 4) {
        printf ("CC-basic simulator: v1.0\n");
        printf ("Usage: %s <output file> <cc file> <ss file> [S <range>] [N <frame length>] [F <max number of frames>] [V]\n\n", argv[0]);
        exit (0);
	}
    strcpy (output_file_name  , argv[1]);
    strcpy (encoder_file_name  , argv[2]);
    strcpy (signal_file_name  , argv[3]);

    if(argc >= 4) {

	    for(i=4;i<argc;i++) {
            if (*argv[i] == 'F') {
                sscanf(argv[i+1],"%d",&max_frames);
			    i++;
			} else if (*argv[i] == 'V') {
				verbose = 1;
			} else if (*argv[i] == 'N') {
                sscanf(argv[i+1],"%d",&Frame_length);
			    i++;
			} else if (*argv[i] == 'S') {
                sscanf(argv[i+1],"%lf",&EbN0_init);
                sscanf(argv[i+2],"%lf",&EbN0_step);
                sscanf(argv[i+3],"%lf",&EbN0_final);
			    i += 3;
			} else {
			    printf("-ERROR: input option not supported-\n");
			    exit(0);
			}
		}
	}


// CC encoder initialization

    Rate = BinEnc_init(encoder_file_name,signal_file_name);
	
// CC decoder initialization

    Viterbi_init();
	//BCJR_init();

// Simulation parameter summary

	if(verbose) {
    printf("\n*************************************************\n");
    printf("OUTPUT       FILE    : %s\n",output_file_name);
    printf("-------------------------------------------------\n");
    printf("BINARY CC ENCODER\n");
    printf("k                    : %d\n",Kbit);
    printf("n                    : %d\n",Nbit);
    printf("Memory               : %d\n",Memory);
    printf("Generators           :\n");
    for(i=0;i<Nbit+Memory;i++) {
        for(j=0;j<Kbit+Memory;j++) printf("%d ",Generator[i][j]);
        printf("\n");
	}
    printf("Coding rate          : %lf\n",Rate);
    printf("-------------------------------------------------\n");
    printf("SIGNAL SET           : %s\n",signal_file_name);
    printf("Coordinates          :\n");
    for(i=0;i<Asize;i++) printf("%lf %lf\n",Signal_set[i][0],Signal_set[i][1]);
	printf("-------------------------------------------------\n");
    printf("CHANNEL              : AWGN\n");
    printf("EbN0 range (dB)      : %lf : %lf : %lf\n",EbN0_init,EbN0_step,EbN0_final);
    printf("MAX FRAMES           : %d\n",max_frames);
    printf("-------------------------------------------------\n");
    printf("FRAME LENGHT         : %d\n",Frame_length);
    printf("-------------------------------------------------\n");
    //printf("DECODER              : Viterbi Algorithm\n");
    printf("DECODER              : BCJR\n");
    printf("-------------------------------------------------\n");

    fflush(stdout);
	}

    file = fopen(output_file_name,"w");
    if(file == NULL) {
        printf("- ERROR: File not found -\n");
        exit(0);
	}

// EbN0 loop
   
	mask = Asize - 1;
	
    for(EbN0=EbN0_init;EbN0<=EbN0_final + 1e-7;EbN0 += EbN0_step) {

// SNR calculation 

    A = sqrt(Rate*pow(10.,EbN0/10));
    SNR = EbN0 + 10*log(Rate)/LN10;

    symb_error_cnt  = 0;
    bit_error_cnt   = 0;
	frame_error_cnt = 0;
    frame_cnt       = 0;

// Simulation loop: this is the main Monte Carlo simulation loop

    for(frame=0;frame<max_frames;frame++) {

       // Information bits
	   for(i=0;i<Frame_length;i++) Info_frame[i] = Bitq(Input_size);

       // Encoding
	   init_state = 0;
	   final_state = FrameEncoder(init_state);
  
	   // Transmission and demodulation
	   
       AWGNChannel(A);
	   
	   // Decoding

       check_state = Viterbi(init_state,final_state);
       if(check_state != init_state) {
	       printf("- Error: trace-back does not recover initial state\n");
		   exit(0);
		}	   
	   
        /* BCJR(init_state,final_state); */

// Counting errors

      out = 0;
      for(i=0;i<Frame_length;i++) {
          b = Info_frame[i]^Detected_info_frame[i];
          for(j=0;j<Kbit;j++) {
		      a = (b >> j)&1;
			  if(a > 0) out = 1;
		      bit_error_cnt += a;
		  }
          a = Coded_frame[i]^Detected_code_frame[i];
          for(j=0;j<Sdim;j++) symb_error_cnt += (((a >> j*Mbit)&mask) != 0);
      }
	  frame_error_cnt += out;
      frame_cnt += 1;

      if((bit_error_cnt >= MIN_BIT_ERRORS) && (frame_error_cnt >= MIN_FRAME_ERRORS)) break;

    }

/* ----------------------------------------------------------------- */
/* Output section */

    if(verbose) {
    printf("Number of simulated frames : %d\n",frame_cnt);
    printf("Number of simulated bits   : %d\n",frame_cnt*Frame_length*Kbit);
    printf("Number of symbol errors    : %d\n",symb_error_cnt);
    printf("Number of bit    errors    : %d\n",bit_error_cnt);
    if(bit_error_cnt< MIN_BIT_ERRORS) {
       printf("-- WARNING: simulation terminated with too few bit errors\n");
    }
    if(frame_error_cnt < MIN_FRAME_ERRORS) {
       printf("-- WARNING: simulation terminated with too few frame errors\n");
    }
    ps = (double) symb_error_cnt / (double) frame_cnt / (double) Frame_length / (double) Sdim;
    pb = (double) bit_error_cnt / (double) frame_cnt / (double) Frame_length / (double) Kbit;
    pf = (double) frame_error_cnt / (double) frame_cnt;
    printf("*************************************************\n");
    printf("SER = %.4le\n",ps);
    printf("BER = %.4le\n",pb);
    printf("FER = %.4le\n",pf);
    fflush(stdout);
    }
 
    fprintf(file,"%.3le %.3le %.3le %.3le %.3le\n",EbN0,SNR,ps,pb,pf);
    fflush(file);

  }

  fclose(file);
  return(0);

}
