#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This yields the i-th binary digit of an integer (LSB is in position i = 0)*/
#define Int2bin(x_,i_) ((x_ >> i_)&1)
#define Bin2int(x_,i_) ((x_) << (i_))
#define K 1
#define N 2


void matrixShow(int a[][2][2],int m,int n);
//void matrixMulti(int A[][N], int B[][N], int C[][N], int m, int p, int n);

int main()
{
//    int a = pow(2,7);
    int A[N][N] = {0,1,0,0};
    int B[K][N] = {1,0};
    int C[N][N] = {0,1,1,1};
    int D[K][N] = {1,1};
//    matrixMulti(B,A,C,K,N,N);
//    printf("%d",);

    int Fwd[4][2][2];//[][][0]:next_state;[][][1]:output
    int si;
    for (si = 0;si < pow(2,N);si++){
        int ui;
        for (ui = 0;ui < pow(2,K);ui++){
            int i1;
            int next_state = 0;
            int output = 0;

            for (i1 = 0; i1 < N; i1++){
                int i21,i22;
                int tmp11 = 0;
                int tmp12 = 0;
                int tmp21 = 0;
                int tmp22 = 0;

                for (i21 = 0; i21 < N; i21++){
                    tmp11 = (Int2bin(si,i21)*A[N-1-i21][i1]+tmp11)%2;
                    tmp12 = (Int2bin(si,i21)*C[N-1-i21][i1]+tmp12)%2;
                }

                for (i22 = 0; i22 < K; i22++){
                    tmp21 = (Int2bin(ui,i22)*B[K-1-i22][i1]+tmp21)%2;
                    tmp22 = (Int2bin(ui,i22)*D[K-1-i22][i1]+tmp22)%2;
                }

                next_state += Bin2int((tmp11+tmp21)%2,N-i1-1);
                output += Bin2int((tmp12+tmp22)%2,N-i1-1);
            }
            Fwd[si][ui][0] = next_state;
            Fwd[si][ui][1] = output;
        }

    }




    printf("fertig!\n");
    matrixShow(Fwd,4,2);




    return 0;
}

void matrixShow(int a[][2][2],int m,int n)
{
    int i,j;
    for (i = 0; i < m; i++)
    {
        for(j = 0;j < n; j++ )
        {
            printf(" %d",a[i][j][0]);
        }
        printf("\n");
    }
}

//void matrixMulti(int A[][N], int B[][N], int C[][N], int m, int p, int n)
//{
//    int i,j,k;
//    for (i = 0; i < m; i++)
//    {
//        for (j = 0; j < p; j++)
//        {
//            C[i][j] = 0;
//            for (k = 0; k < n; k++)
//            {
//                C[i][j] += A[i][k] * B[k][j];
//            }
//         }
//     }
//}
