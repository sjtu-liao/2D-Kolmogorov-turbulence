#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include"mpi.h"
#include"gmp.h"
#include"mpfr.h"
#include"mpf2mpfr.h"
#define _GMP_H
#define _MPFR_H
#include "mpi_gmp.h"
#include "mpi_mpfr.h"

#define M 1024
#define N 1024
#define p 10    /*M=N=2^p*/
#define L 140    /*Taylor*/
#define prec 864    /*Ns*/

//#define HE 5
//#define np 32    /* pow(2,HE) */
//#define NP 1024    /* np*np */
//#define evel 32    /* N/np */
//#define evel2 1024    /* evel*evel */

#define HE 6
#define np 64    /* pow(2,HE) */
#define NP 4096    /* np*np */
#define evel 16    /* N/np */
#define evel2 256    /* evel*evel */

mpfr_t psi[evel][evel];
mpfr_t pX[L][evel][evel];
mpfr_t pZ[L][evel][evel];
mpfr_t H1[L][evel][evel];
mpfr_t H2[L][evel][evel];
mpfr_t Coe_psi_real[evel][evel];
mpfr_t Coe_psi_vir[evel][evel];
mpfr_t Coe_psi0_real[evel][evel];
mpfr_t Coe_psi0_vir[evel][evel];
mpfr_t MN[evel][evel];
mpfr_t F[M];
mpfr_t hl[L+1];

mpfr_t CH1[evel][evel];
mpfr_t CH2[evel][evel];
mpfr_t CH1N[evel][evel];
mpfr_t CH2N[evel][evel];
mpfr_t CH[evel][evel];
mpfr_t MM[evel][evel];
mpfr_t NN[evel][evel];
mpfr_t MMN[evel][evel];
mpfr_t NNN[evel][evel];

mpfr_t Spectrum[501];

int myid, numprocs;
long int l, m, n, i, j, k;
long int tempZ, tempi, tempii;
long int vj, IG, IU, CC, ss, SS, vv, I1, I2, vr, vg, vk;
long int d[p+1], d2[p+1];
long int DAO[N][p], RD[N], rd[N];
long int NOWid[N][N], NOWm[N][N], NOWn[N][N];
long int Mm[evel][evel], Nn[evel][evel];

mpfr_t tempreal, tempvir, tempreal2, tempvir2, temp, temp1, temp2;
mpfr_t TEMPreal[evel][evel], TEMPvir[evel][evel], VIR[evel][evel];
mpfr_t Wnzreal[N][p+1], Wnzvir[N][p+1], Wnreal[N][p+1], Wnvir[N][p+1];

void *packed_REAL;
void *packed_VIR;
void *packed_REAL2;
void *packed_VIR2;

mpfr_t Coe[evel][evel*2];
mpfr_t CoeALL[NP][evel][evel*2];
void *packed_Coe;
void *packed_CoeALL[NP];

MPI_Status status, status1;

double Gauss[M];
int i_gauss;

void gauss(double ex, double dx) //ex:均值；dx:方差
{
    for(i_gauss=0;i_gauss<N;i_gauss++)
    {
        Gauss[i_gauss]=(sqrt(-2*log((double)rand()/RAND_MAX))*cos((double)rand()/RAND_MAX*2*3.1415926))*sqrt(dx)+ex;
        if(Gauss[i_gauss]>100.0||Gauss[i_gauss]<(-1.0)*100.0)
            Gauss[i_gauss]=0.0;
    }
}

void fft(int ZHE, int NN, int ID, int IP, mpfr_t aareal[evel][evel], mpfr_t aavir[evel][evel], mpfr_t bbreal[evel][evel], mpfr_t bbvir[evel][evel])
{
tempZ=ID*evel;
for(vj=1; vj<=p; vj++)
{
    if(vj<=HE)
    {
    IG=NN/d2[vj];
    IU=NN/d[vj];
    SS=tempZ%IG;
    ss=tempZ/IG;
    tempi=(IU/evel)*ZHE;
    if(SS<IU)
    {
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Send(packed_REAL, evel2, MPI_MPF, myid+tempi, 1, MPI_COMM_WORLD);
        MPI_Send(packed_VIR, evel2, MPI_MPF, myid+tempi, 2, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(packed_REAL2, evel2, MPI_MPF, myid-tempi, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(packed_VIR2, evel2, MPI_MPF, myid-tempi, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    if(SS>=IU)
    {
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Send(packed_REAL, evel2, MPI_MPF, myid-tempi, 1, MPI_COMM_WORLD);
        MPI_Send(packed_VIR, evel2, MPI_MPF, myid-tempi, 2, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(packed_REAL2, evel2, MPI_MPF, myid+tempi, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(packed_VIR2, evel2, MPI_MPF, myid+tempi, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    tempi=tempZ-ss*IG-IU;
    for(vr=0; vr<evel; vr++)
    {
    tempii=tempi+vr;
    for(vk=0; vk<evel; vk++)
    {
    if(SS<IU)
    {
        mpfr_add(aareal[vr][vk], aareal[vr][vk], bbreal[vr][vk], GMP_RNDN);
        mpfr_add(aavir[vr][vk], aavir[vr][vk], bbvir[vr][vk], GMP_RNDN);
    }
    else
    {
        mpfr_sub(tempreal, bbreal[vr][vk], aareal[vr][vk], GMP_RNDN);
        mpfr_sub(tempvir, bbvir[vr][vk], aavir[vr][vk], GMP_RNDN);
        mpfr_mul(temp1, tempreal, Wnreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, tempvir, Wnvir[tempii][vj], GMP_RNDN);
        mpfr_sub(aareal[vr][vk], temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, tempvir, Wnreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, tempreal, Wnvir[tempii][vj], GMP_RNDN);
        mpfr_add(aavir[vr][vk], temp1, temp2, GMP_RNDN);
    }
    }
    }
    }
    else
    {
    IG=NN/d2[vj];
    IU=NN/d[vj];
    CC=NN/d[vj]-1;
    vv=evel/IG;
    for(vg=0; vg<=vv-1; vg++)
    {
    I1=vg*IG;
    I2=vg*IG+IU;
    for(vr=0; vr<=CC; vr++)
    {
    tempi=vr+I1;
    tempii=vr+I2;
    for(vk=0; vk<evel; vk++)
    {
        mpfr_set(tempreal, aareal[tempi][vk], GMP_RNDN);
        mpfr_set(tempvir, aavir[tempi][vk], GMP_RNDN);
        mpfr_add(aareal[tempi][vk], tempreal, aareal[tempii][vk], GMP_RNDN);
        mpfr_add(aavir[tempi][vk], tempvir, aavir[tempii][vk], GMP_RNDN);
        mpfr_sub(aareal[tempii][vk], tempreal, aareal[tempii][vk], GMP_RNDN);
        mpfr_sub(aavir[tempii][vk], tempvir, aavir[tempii][vk], GMP_RNDN);
        mpfr_set(tempreal, aareal[tempii][vk], GMP_RNDN);
        mpfr_set(tempvir, aavir[tempii][vk], GMP_RNDN);
        mpfr_mul(temp1, tempreal, Wnreal[vr][vj], GMP_RNDN);
        mpfr_mul(temp2, tempvir, Wnvir[vr][vj], GMP_RNDN);
        mpfr_sub(aareal[tempii][vk], temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, tempvir, Wnreal[vr][vj], GMP_RNDN);
        mpfr_mul(temp2, tempreal, Wnvir[vr][vj], GMP_RNDN);
        mpfr_add(aavir[tempii][vk], temp1, temp2, GMP_RNDN);
    }
    }
    }
    }
}
tempZ=IP*evel;
for(vj=1; vj<=p; vj++)
{
    if(vj<=HE)
    {
    IG=NN/d2[vj];
    IU=NN/d[vj];
    SS=tempZ%IG;
    ss=tempZ/IG;
    tempi=IU/evel;
    if(SS<IU)
    {
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Send(packed_REAL, evel2, MPI_MPF, myid+tempi, 1, MPI_COMM_WORLD);
        MPI_Send(packed_VIR, evel2, MPI_MPF, myid+tempi, 2, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(packed_REAL2, evel2, MPI_MPF, myid-tempi, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(packed_VIR2, evel2, MPI_MPF, myid-tempi, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    if(SS>=IU)
    {
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Send(packed_REAL, evel2, MPI_MPF, myid-tempi, 1, MPI_COMM_WORLD);
        MPI_Send(packed_VIR, evel2, MPI_MPF, myid-tempi, 2, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(packed_REAL2, evel2, MPI_MPF, myid+tempi, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(packed_VIR2, evel2, MPI_MPF, myid+tempi, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    tempi=tempZ-ss*IG-IU;
    for(vr=0; vr<evel; vr++)
    {
    tempii=tempi+vr;
    for(vk=0; vk<evel; vk++)
    {
    if(SS<IU)
    {
        mpfr_add(aareal[vk][vr], aareal[vk][vr], bbreal[vk][vr], GMP_RNDN);
        mpfr_add(aavir[vk][vr], aavir[vk][vr], bbvir[vk][vr], GMP_RNDN);
    }
    else
    {
        mpfr_sub(tempreal, bbreal[vk][vr], aareal[vk][vr], GMP_RNDN);
        mpfr_sub(tempvir, bbvir[vk][vr], aavir[vk][vr], GMP_RNDN);
        mpfr_mul(temp1, tempreal, Wnreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, tempvir, Wnvir[tempii][vj], GMP_RNDN);
        mpfr_sub(aareal[vk][vr], temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, tempvir, Wnreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, tempreal, Wnvir[tempii][vj], GMP_RNDN);
        mpfr_add(aavir[vk][vr], temp1, temp2, GMP_RNDN);
    }
    }
    }
    }
    else
    {
    IG=NN/d2[vj];
    IU=NN/d[vj];
    CC=NN/d[vj]-1;
    vv=evel/IG;
    for(vg=0; vg<=vv-1; vg++)
    {
    I1=vg*IG;
    I2=vg*IG+IU;
    for(vr=0; vr<=CC; vr++)
    {
    tempi=vr+I1;
    tempii=vr+I2;
    for(vk=0; vk<evel; vk++)
    {
        mpfr_set(tempreal, aareal[vk][tempi], GMP_RNDN);
        mpfr_set(tempvir, aavir[vk][tempi], GMP_RNDN);
        mpfr_add(aareal[vk][tempi], tempreal, aareal[vk][tempii], GMP_RNDN);
        mpfr_add(aavir[vk][tempi], tempvir, aavir[vk][tempii], GMP_RNDN);
        mpfr_sub(aareal[vk][tempii], tempreal, aareal[vk][tempii], GMP_RNDN);
        mpfr_sub(aavir[vk][tempii], tempvir, aavir[vk][tempii], GMP_RNDN);
        mpfr_set(tempreal, aareal[vk][tempii], GMP_RNDN);
        mpfr_set(tempvir, aavir[vk][tempii], GMP_RNDN);
        mpfr_mul(temp1, tempreal, Wnreal[vr][vj], GMP_RNDN);
        mpfr_mul(temp2, tempvir, Wnvir[vr][vj], GMP_RNDN);
        mpfr_sub(aareal[vk][tempii], temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, tempvir, Wnreal[vr][vj], GMP_RNDN);
        mpfr_mul(temp2, tempreal, Wnvir[vr][vj], GMP_RNDN);
        mpfr_add(aavir[vk][tempii], temp1, temp2, GMP_RNDN);
    }
    }
    }
    }
}
}

void ifft(int ZHE, int NN, int ID, int IP, mpfr_t aareal[evel][evel], mpfr_t aavir[evel][evel], mpfr_t bbreal[evel][evel], mpfr_t bbvir[evel][evel])
{
tempZ=IP*evel;
for(vj=1; vj<=p; vj++)
{
    I1=d[vj]/2;
    if(vj<=p-HE)
    {
    for(vg=0; vg<=I1-1; vg++)
    {
    for(vr=vg; vr<evel; vr+=d[vj])
    {
    tempi=vr+I1;
    for(vk=0; vk<evel; vk++)
    {
        mpfr_mul(temp1, aareal[vk][tempi], Wnzreal[vg][vj], GMP_RNDN);
        mpfr_mul(temp2, aavir[vk][tempi], Wnzvir[vg][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, aavir[vk][tempi], Wnzreal[vg][vj], GMP_RNDN);
        mpfr_mul(temp2, aareal[vk][tempi], Wnzvir[vg][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_set(tempreal2, aareal[vk][vr], GMP_RNDN);
        mpfr_set(tempvir2, aavir[vk][vr], GMP_RNDN);
        mpfr_add(aareal[vk][vr], tempreal2, tempreal, GMP_RNDN);
        mpfr_add(aavir[vk][vr], tempvir2, tempvir, GMP_RNDN);
        mpfr_sub(aareal[vk][tempi], tempreal2, tempreal, GMP_RNDN);
        mpfr_sub(aavir[vk][tempi], tempvir2, tempvir, GMP_RNDN);
    }
    }
    }
    }
    else
    {
    SS=tempZ%d[vj];
    ss=tempZ/d[vj];
    tempi=I1/evel;
    if(SS<I1)
    {
        tempii=myid+tempi;
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Sendrecv(packed_REAL, evel2, MPI_MPF, tempii, 1, packed_REAL2, evel2, MPI_MPF, tempii, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(packed_VIR, evel2, MPI_MPF, tempii, 2, packed_VIR2, evel2, MPI_MPF, tempii, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    else
    {
        tempii=myid-tempi;
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Sendrecv(packed_REAL, evel2, MPI_MPF, tempii, 1, packed_REAL2, evel2, MPI_MPF, tempii, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(packed_VIR, evel2, MPI_MPF, tempii, 2, packed_VIR2, evel2, MPI_MPF, tempii, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    tempi=tempZ-ss*d[vj];
    for(vg=0; vg<evel; vg++)
    {
    tempii=tempi+vg;
    I2=tempii-I1;
    for(vk=0; vk<evel; vk++)
    {
    if(SS<I1)
    {
        mpfr_mul(temp1, bbreal[vk][vg], Wnzreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, bbvir[vk][vg], Wnzvir[tempii][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, bbvir[vk][vg], Wnzreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, bbreal[vk][vg], Wnzvir[tempii][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_add(aareal[vk][vg], aareal[vk][vg], tempreal, GMP_RNDN);
        mpfr_add(aavir[vk][vg], aavir[vk][vg], tempvir, GMP_RNDN);
    }
    else
    {
        mpfr_mul(temp1, aareal[vk][vg], Wnzreal[I2][vj], GMP_RNDN);
        mpfr_mul(temp2, aavir[vk][vg], Wnzvir[I2][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, aavir[vk][vg], Wnzreal[I2][vj], GMP_RNDN);
        mpfr_mul(temp2, aareal[vk][vg], Wnzvir[I2][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_sub(aareal[vk][vg], bbreal[vk][vg], tempreal, GMP_RNDN);
        mpfr_sub(aavir[vk][vg], bbvir[vk][vg], tempvir, GMP_RNDN);
    }
    }
    }
    }
}
tempZ=ID*evel;
for(vj=1; vj<=p; vj++)
{
    I1=d[vj]/2;
    if(vj<=p-HE)
    {
    for(vg=0; vg<=I1-1; vg++)
    {
    for(vr=vg; vr<evel; vr+=d[vj])
    {
    tempi=vr+I1;
    for(vk=0; vk<evel; vk++)
    {
        mpfr_mul(temp1, aareal[tempi][vk], Wnzreal[vg][vj], GMP_RNDN);
        mpfr_mul(temp2, aavir[tempi][vk], Wnzvir[vg][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, aavir[tempi][vk], Wnzreal[vg][vj], GMP_RNDN);
        mpfr_mul(temp2, aareal[tempi][vk], Wnzvir[vg][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_set(tempreal2, aareal[vr][vk], GMP_RNDN);
        mpfr_set(tempvir2, aavir[vr][vk], GMP_RNDN);
        mpfr_add(aareal[vr][vk], tempreal, tempreal2, GMP_RNDN);
        mpfr_add(aavir[vr][vk], tempvir, tempvir2, GMP_RNDN);
        mpfr_sub(aareal[tempi][vk], tempreal2, tempreal, GMP_RNDN);
        mpfr_sub(aavir[tempi][vk], tempvir2, tempvir, GMP_RNDN);
    }
    }
    }
    }
    else
    {
    SS=tempZ%d[vj];
    ss=tempZ/d[vj];
    tempi=I1/evel*ZHE;
    if(SS<I1)
    {
        tempii=myid+tempi;
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Sendrecv(packed_REAL, evel2, MPI_MPF, tempii, 1, packed_REAL2, evel2, MPI_MPF, tempii, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(packed_VIR, evel2, MPI_MPF, tempii, 2, packed_VIR2, evel2, MPI_MPF, tempii, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    else
    {
        tempii=myid-tempi;
        pack_mpf(aareal[0][0], evel2, packed_REAL);
        pack_mpf(aavir[0][0], evel2, packed_VIR);
        MPI_Sendrecv(packed_REAL, evel2, MPI_MPF, tempii, 1, packed_REAL2, evel2, MPI_MPF, tempii, 1, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(packed_VIR, evel2, MPI_MPF, tempii, 2, packed_VIR2, evel2, MPI_MPF, tempii, 2, MPI_COMM_WORLD, &status1);
        unpack_mpf(packed_REAL2, bbreal[0][0], evel2);
        unpack_mpf(packed_VIR2, bbvir[0][0], evel2);
    }
    tempi=tempZ-ss*d[vj];
    for(vg=0; vg<evel; vg++)
    {
    tempii=tempi+vg;
    I2=tempii-I1;
    for(vk=0; vk<evel; vk++)
    {
    if(SS<I1)
    {
        mpfr_mul(temp1, bbreal[vg][vk], Wnzreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, bbvir[vg][vk], Wnzvir[tempii][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, bbvir[vg][vk], Wnzreal[tempii][vj], GMP_RNDN);
        mpfr_mul(temp2, bbreal[vg][vk], Wnzvir[tempii][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_add(aareal[vg][vk], aareal[vg][vk], tempreal, GMP_RNDN);
        mpfr_add(aavir[vg][vk], aavir[vg][vk], tempvir, GMP_RNDN);
    }
    else
    {
        mpfr_mul(temp1, aareal[vg][vk], Wnzreal[I2][vj], GMP_RNDN);
        mpfr_mul(temp2, aavir[vg][vk], Wnzvir[I2][vj], GMP_RNDN);
        mpfr_sub(tempreal, temp1, temp2, GMP_RNDN);
        mpfr_mul(temp1, aavir[vg][vk], Wnzreal[I2][vj], GMP_RNDN);
        mpfr_mul(temp2, aareal[vg][vk], Wnzvir[I2][vj], GMP_RNDN);
        mpfr_add(tempvir, temp1, temp2, GMP_RNDN);
        mpfr_sub(aareal[vg][vk], bbreal[vg][vk], tempreal, GMP_RNDN);
        mpfr_sub(aavir[vg][vk], bbvir[vg][vk], tempvir, GMP_RNDN);
    }
    }
    }
    }
}
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    mpfr_set_default_prec(prec);
    commit_mpf(&(MPI_MPF),prec,MPI_COMM_WORLD);
    create_mpf_op(&(MPI_MPF_SUM), _mpi_mpf_add, MPI_COMM_WORLD);
    double st, et;
    long int k_2D;
    double K_2D;
    long int Nint, Nint2, Ss;
    Nint=N;
    Nint2=N*N;
    Ss=0;
    int Na=N/3;
    int ID, IP, ZHE;
    ZHE=np;
    ID=(myid%NP)/ZHE;
    IP=(myid%NP)%ZHE;
    int T_compare, Compare;
    FILE *fp;
    char filename[64];
    char Char[1000];
    char Char1[1000];
    char Char2[1000];
    mpfr_t PI, ZERO, Small, POINTONE, HALF, NEGONE, ONE, TWO, FOUR, NK, Re;
    mpfr_inits2(prec, tempreal, tempvir, tempreal2, tempvir2, temp, temp1, temp2, PI, ZERO, Small, POINTONE, HALF, NEGONE, ONE, TWO, FOUR, NK, Re, (mpfr_ptr) 0);
    mpfr_set_str(PI, "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132000568127145263560827785771342757789609173637178721468440901224953430146549585371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303598253490428755468731159562863882353787593751957781857780532171226806613001927876611195909216420199", 10, GMP_RNDN);
    mpfr_set_str(ZERO, "0.0", 10, GMP_RNDN);
    mpfr_set_str(Small, "0.0001", 10, GMP_RNDN);
    mpfr_set_str(POINTONE, "0.1", 10, GMP_RNDN);
    mpfr_set_str(HALF, "0.5", 10, GMP_RNDN);
    mpfr_set_str(NEGONE, "-1.0", 10, GMP_RNDN);
    mpfr_set_str(ONE, "1.0", 10, GMP_RNDN);
    mpfr_set_str(TWO, "2.0", 10, GMP_RNDN);
    mpfr_set_str(FOUR, "4.0", 10, GMP_RNDN);
    mpfr_set_str(NK, "16.0", 10, GMP_RNDN);    /*n_K*/
    mpfr_set_str(Re, "2000.0", 10, GMP_RNDN);    /*Re*/
    packed_REAL=allocbuf_mpf(prec, evel2);
    packed_VIR=allocbuf_mpf(prec, evel2);
    packed_REAL2=allocbuf_mpf(prec, evel2);
    packed_VIR2=allocbuf_mpf(prec, evel2);
    packed_Coe=allocbuf_mpf(prec, evel2*2);
    for(i=0; i<NP; i++)
    {
        packed_CoeALL[i]=allocbuf_mpf(prec, evel2*2);
    }
    mpfr_t h, T, t, TS;
    mpfr_inits2(prec, h, T, t, TS, (mpfr_ptr) 0);
    mpfr_set_str(h, "0.001", 10, GMP_RNDN);    /*h*/
    mpfr_set_str(T, "300.0001", 10, GMP_RNDN);    /*T*/
    mpfr_set_str(TS, "0.1", 10, GMP_RNDN);    /*T-out*/
    mpfr_set(t, ZERO, GMP_RNDN);
    T_compare=0;
    for(k=0; k<evel; k++)
    {
        for(j=0; j<evel; j++)
        {
            mpfr_inits2(prec, psi[j][k], TEMPreal[j][k], TEMPvir[j][k], VIR[j][k], Coe_psi_real[j][k], Coe_psi_vir[j][k], Coe_psi0_real[j][k], Coe_psi0_vir[j][k], MN[j][k], CH1[j][k], CH2[j][k], CH1N[j][k], CH2N[j][k], CH[j][k], MM[j][k], NN[j][k], MMN[j][k], NNN[j][k], Coe[j][k], Coe[j][k+evel], (mpfr_ptr) 0);
            for(i=0; i<NP; i++)
            {
                mpfr_inits2(prec, CoeALL[i][j][k], CoeALL[i][j][k+evel], (mpfr_ptr) 0);
            }
            for(l=0; l<L; l++)
            {
                mpfr_inits2(prec, pX[l][j][k], pZ[l][j][k], H1[l][j][k], H2[l][j][k], (mpfr_ptr) 0);
            }
        }
    }
    for(i=0; i<=500; i++)
    {
        mpfr_inits2(prec, Spectrum[i], (mpfr_ptr) 0);
    }
    for(j=1; j<=p; j++)
    {
        d[j]=pow(2, j);
        d2[j]=pow(2, j-1);
    }
    for(k=0; k<N; k++)
    {
        mpfr_mul(temp, TWO, PI, GMP_RNDN);
        mpfr_mul_si(temp, temp, k, GMP_RNDN);
        for(j=1; j<=p; j++)
        {
            mpfr_inits2(prec, Wnzreal[k][j], Wnzvir[k][j], Wnreal[k][j], Wnvir[k][j], (mpfr_ptr) 0);
            mpfr_div_si(temp1, temp, d[j], GMP_RNDN);
            mpfr_cos(Wnzreal[k][j], temp1, GMP_RNDN);
            mpfr_sin(Wnzvir[k][j], temp1, GMP_RNDN);
            tempi=N/d2[j];
            mpfr_div_si(temp2, temp, tempi, GMP_RNDN);
            mpfr_cos(Wnreal[k][j], temp2, GMP_RNDN);
            mpfr_sin(Wnvir[k][j], temp2, GMP_RNDN);
            mpfr_mul(Wnvir[k][j], Wnvir[k][j], NEGONE, GMP_RNDN);
        }
    }
    for(i=0; i<N; i++)
    {
        tempi=i;
        for(j=0; j<p; j++)
        {
            DAO[i][j]=tempi%2;
            tempi=tempi/2;
        }
    }
    for(i=0; i<N; i++)
    {
        RD[i]=0;
        for(j=0; j<p; j++)
        {
            RD[i]=RD[i]+DAO[i][j]*(pow(2, p-j-1));
        }
    }
    for(i=0; i<N; i++)
    {
        if(RD[i]<N/2)
            rd[i]=RD[i];
        else
            rd[i]=RD[i]-N;
    }
    for(n=0; n<evel; n++)
    {
        for(m=0; m<evel; m++)
        {
            j=ID*evel+m;
            k=IP*evel+n;
            Mm[m][n]=RD[j];
            Nn[m][n]=RD[k];
            if(Mm[m][n]>=N/2)
                Mm[m][n]-=N;
            if(Nn[m][n]>=N/2)
                Nn[m][n]-=N;
        }
    }
    
    for(j=0; j<N; j++)
    {
        for(k=0; k<N; k++)
        {
            tempi=0;
            for(m=0; m<N; m++)
            {
                for(n=0; n<N; n++)
                {
                    if(RD[m]==j&&RD[n]==k)
                    {
                        tempi=1;
                        break;
                    }
                }
                if(tempi==1)
                    break;
            }
            NOWid[j][k]=(m/evel)*np+(n/evel);
            NOWm[j][k]=m%evel;
            NOWn[j][k]=n%evel;
        }
    }
    
    for(n=0; n<evel; n++)
    {
        for(m=0; m<evel; m++)
        {
            mpfr_mul_si(temp1, ONE, Nn[m][n], GMP_RNDN);
            mpfr_mul(temp1, temp1, temp1, GMP_RNDN);
            mpfr_mul_si(temp2, ONE, Mm[m][n], GMP_RNDN);
            mpfr_mul(temp2, temp2, temp2, GMP_RNDN);
            mpfr_add(MN[m][n], temp1, temp2, GMP_RNDN);
            mpfr_mul(MN[m][n], MN[m][n], NEGONE, GMP_RNDN);
        }
    }
    
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul_si(temp1, ONE, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp1, temp1, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp1, temp1, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, ONE, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, temp2, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, temp2, Mm[j][k], GMP_RNDN);
                    mpfr_add(CH1[j][k], temp1, temp2, GMP_RNDN);
                    mpfr_mul(CH1N[j][k], CH1[j][k], NEGONE, GMP_RNDN);
                }
            }
            
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul_si(temp1, ONE, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp1, temp1, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp1, temp1, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, ONE, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, temp2, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, temp2, Mm[j][k], GMP_RNDN);
                    mpfr_add(CH2[j][k], temp1, temp2, GMP_RNDN);
                    mpfr_mul(CH2N[j][k], CH2[j][k], NEGONE, GMP_RNDN);
                }
            }
            
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul_si(temp1, ONE, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp1, temp1, Nn[j][k], GMP_RNDN);
                    mpfr_mul(temp1, temp1, temp1, GMP_RNDN);
                    mpfr_mul_si(temp2, ONE, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp2, temp2, Mm[j][k], GMP_RNDN);
                    mpfr_mul(temp2, temp2, temp2, GMP_RNDN);
                    mpfr_mul_si(temp, TWO, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp, temp, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(temp, temp, Nn[j][k], GMP_RNDN);
                    mpfr_mul_si(temp, temp, Nn[j][k], GMP_RNDN);
                    mpfr_add(temp, temp, temp1, GMP_RNDN);
                    mpfr_add(temp, temp, temp2, GMP_RNDN);
                    mpfr_div(CH[j][k], temp, Re, GMP_RNDN);
                }
            }
            
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul_si(MM[j][k], ONE, Mm[j][k], GMP_RNDN);
                    mpfr_mul_si(NN[j][k], ONE, Nn[j][k], GMP_RNDN);
                    mpfr_mul(MMN[j][k], MM[j][k], NEGONE, GMP_RNDN);
                    mpfr_mul(NNN[j][k], NN[j][k], NEGONE, GMP_RNDN);
                }
            }
    
    /*initial*/
    for(n=0; n<evel; n++)
    {
        for(m=0; m<evel; m++)
        {
            j=ID*evel+m;
            k=IP*evel+n;
            mpfr_mul(temp, TWO, PI, GMP_RNDN);
            mpfr_div_si(temp, temp, Nint, GMP_RNDN);
            mpfr_mul_si(temp1, temp, j, GMP_RNDN);
            mpfr_mul_si(temp2, temp, k, GMP_RNDN);
            mpfr_cos(temp1, temp1, GMP_RNDN);
            mpfr_cos(temp2, temp2, GMP_RNDN);
            mpfr_mul(temp1, temp1, temp2, GMP_RNDN);
            mpfr_mul(psi[m][n], temp1, NEGONE, GMP_RNDN);
        }
    }
            
                    for(k=0; k<evel; k++)    /*Coe_p*/
                    {
                        for(j=0; j<evel; j++)
                        {
                            mpfr_set(Coe_psi_real[j][k], psi[j][k], GMP_RNDN);
                            mpfr_set(Coe_psi_vir[j][k], ZERO, GMP_RNDN);
                        }
                    }
                    fft(ZHE, N, ID, IP, Coe_psi_real, Coe_psi_vir, TEMPreal, TEMPvir);
                    for(k=0; k<evel; k++)
                    {
                        for(j=0; j<evel; j++)
                        {
                            mpfr_div_si(Coe_psi_real[j][k], Coe_psi_real[j][k], Nint2, GMP_RNDN);
                            mpfr_div_si(Coe_psi_vir[j][k], Coe_psi_vir[j][k], Nint2, GMP_RNDN);
                            mpfr_set(Coe_psi0_real[j][k], Coe_psi_real[j][k], GMP_RNDN);
                            mpfr_set(Coe_psi0_vir[j][k], Coe_psi_vir[j][k], GMP_RNDN);
                        }
                    }
                    
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_set(Coe[j][k], Coe_psi_real[j][k], GMP_RNDN);
                    mpfr_set(Coe[j][k+evel], Coe_psi_vir[j][k], GMP_RNDN);
                }
            }
            
            pack_mpf(Coe[0][0], evel2*2, packed_Coe);    /*Gather*/
            if(myid!=0)
            {
                MPI_Send(packed_Coe, evel2*2, MPI_MPF, 0, myid, MPI_COMM_WORLD);
            }
            if(myid==0)
            {
                for(i=1; i<NP; i++)
                {
                    MPI_Recv(packed_CoeALL[i], evel2*2, MPI_MPF, i, i, MPI_COMM_WORLD, &status);
                }
            }
            if(myid==0)
            {
                unpack_mpf(packed_Coe, CoeALL[0][0][0], evel2*2);
                for(i=1; i<NP; i++)
                {
                    unpack_mpf(packed_CoeALL[i], CoeALL[i][0][0], evel2*2);
                }
            }
            
            if(myid==0)
            {
                {
                    sprintf(filename, "ALL_t%d.dat", Ss);
                    fp=fopen(filename, "w+");
                    for(k=0; k<N; k++)
                    {
                        for(j=0; j<N; j++)
                        {
                            mpfr_sprintf(Char1, "%.262Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]]);
                            mpfr_sprintf(Char2, "%.262Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]+evel]);
                            fprintf(fp, "%d\t%d\t%.272s\t%.272s\n", k, j, Char1, Char2);
                        }
                    }
                    fclose(fp);
                }
            }
            
    for(l=0; l<=L; l++)
    {
        mpfr_inits2(prec, hl[l], (mpfr_ptr) 0);
        if(l==0)
            mpfr_set(hl[l], ONE, GMP_RNDN);
        else
            mpfr_mul(hl[l], hl[l-1], h, GMP_RNDN);
    }
    
    for(j=0; j<M; j++)
    {
        mpfr_inits2(prec, F[j], (mpfr_ptr) 0);
        mpfr_mul(temp, TWO, PI, GMP_RNDN);
        mpfr_div_si(temp, temp, Nint, GMP_RNDN);
        mpfr_mul_si(temp, temp, j, GMP_RNDN);
        mpfr_mul(temp, NK, temp, GMP_RNDN);
        mpfr_cos(temp, temp, GMP_RNDN);
        mpfr_mul(F[j], NK, temp, GMP_RNDN);
    }

    st=MPI_Wtime();
            
    for(; T_compare<=0; )
    {
        for(l=1; l<=L; l++)
        {
            
            for(k=0; k<evel; k++)    /*pX*/
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(pX[l-1][j][k], Coe_psi_vir[j][k], NNN[j][k], GMP_RNDN);
                    mpfr_mul(VIR[j][k], Coe_psi_real[j][k], NN[j][k], GMP_RNDN);
                }
            }
            ifft(ZHE, N, ID, IP, pX[l-1], VIR, TEMPreal, TEMPvir);
            
            for(k=0; k<evel; k++)    /*pZ*/
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(pZ[l-1][j][k], Coe_psi_vir[j][k], MMN[j][k], GMP_RNDN);
                    mpfr_mul(VIR[j][k], Coe_psi_real[j][k], MM[j][k], GMP_RNDN);
                }
            }
            ifft(ZHE, N, ID, IP, pZ[l-1], VIR, TEMPreal, TEMPvir);
            
            for(k=0; k<evel; k++)    /*H1*/
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(H1[l-1][j][k], Coe_psi_vir[j][k], CH1[j][k], GMP_RNDN);
                    mpfr_mul(VIR[j][k], Coe_psi_real[j][k], CH1N[j][k], GMP_RNDN);
                }
            }
            ifft(ZHE, N, ID, IP, H1[l-1], VIR, TEMPreal, TEMPvir);
            
            for(k=0; k<evel; k++)    /*H2*/
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(H2[l-1][j][k], Coe_psi_vir[j][k], CH2[j][k], GMP_RNDN);
                    mpfr_mul(VIR[j][k], Coe_psi_real[j][k], CH2N[j][k], GMP_RNDN);
                }
            }
            ifft(ZHE, N, ID, IP, H2[l-1], VIR, TEMPreal, TEMPvir);
            
            for(k=0; k<evel; k++)    /*pXXXX+pZZZZ+pXXZZ*/
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(Coe_psi_real[j][k], Coe_psi_real[j][k], CH[j][k], GMP_RNDN);
                    mpfr_mul(Coe_psi_vir[j][k], Coe_psi_vir[j][k], CH[j][k], GMP_RNDN);
                }
            }
            ifft(ZHE, N, ID, IP, Coe_psi_real, Coe_psi_vir, TEMPreal, TEMPvir);
            
            for(k=0; k<evel; k++)    /*f*/
            {
                for(j=0; j<evel; j++)
                {
                    m=ID*evel+j;
                    mpfr_set(temp1, ZERO, GMP_RNDN);
                    mpfr_set(temp2, ZERO, GMP_RNDN);
                    for(i=0; i<=l-1; i++)
                    {
                        mpfr_mul(temp, pZ[i][j][k], H1[l-1-i][j][k], GMP_RNDN);
                        mpfr_add(temp1, temp1, temp, GMP_RNDN);
                        mpfr_mul(temp, pX[i][j][k], H2[l-1-i][j][k], GMP_RNDN);
                        mpfr_add(temp2, temp2, temp, GMP_RNDN);
                    }
                    mpfr_add(Coe_psi_real[j][k], Coe_psi_real[j][k], temp1, GMP_RNDN);
                    mpfr_sub(Coe_psi_real[j][k], Coe_psi_real[j][k], temp2, GMP_RNDN);
                    if(l==1)
                    {
                        mpfr_sub(Coe_psi_real[j][k], Coe_psi_real[j][k], F[m], GMP_RNDN);
                    }
                    mpfr_div_si(Coe_psi_real[j][k], Coe_psi_real[j][k], l, GMP_RNDN);
                    mpfr_set(Coe_psi_vir[j][k], ZERO, GMP_RNDN);
                }
            }
            
            fft(ZHE, N, ID, IP, Coe_psi_real, Coe_psi_vir, TEMPreal, TEMPvir);    /*Coe_p*/
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    if(myid==0&&j==0&&k==0)
                    {
                        mpfr_set(Coe_psi_real[j][k], ZERO, GMP_RNDN);
                        mpfr_set(Coe_psi_vir[j][k], ZERO, GMP_RNDN);
                    }
                    else
                    {
                        mpfr_div_si(Coe_psi_real[j][k], Coe_psi_real[j][k], Nint2, GMP_RNDN);
                        mpfr_div_si(Coe_psi_vir[j][k], Coe_psi_vir[j][k], Nint2, GMP_RNDN);
                        mpfr_div(Coe_psi_real[j][k], Coe_psi_real[j][k], MN[j][k], GMP_RNDN);
                        mpfr_div(Coe_psi_vir[j][k], Coe_psi_vir[j][k], MN[j][k], GMP_RNDN);
                    }
                    if(Mm[j][k]>Na||Mm[j][k]<-Na||Nn[j][k]>Na||Nn[j][k]<-Na)
                    {
                        mpfr_set(Coe_psi_real[j][k], ZERO, GMP_RNDN);
                        mpfr_set(Coe_psi_vir[j][k], ZERO, GMP_RNDN);
                    }
                }
            }
            
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_mul(temp1, Coe_psi_real[j][k], hl[l], GMP_RNDN);
                    mpfr_mul(temp2, Coe_psi_vir[j][k], hl[l], GMP_RNDN);
                    mpfr_add(Coe_psi0_real[j][k], Coe_psi0_real[j][k], temp1, GMP_RNDN);
                    mpfr_add(Coe_psi0_vir[j][k], Coe_psi0_vir[j][k], temp2, GMP_RNDN);
                }
            }
            
        }

            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_set(Coe_psi_real[j][k], Coe_psi0_real[j][k], GMP_RNDN);
                    mpfr_set(Coe_psi_vir[j][k], Coe_psi0_vir[j][k], GMP_RNDN);
                }
            }
        
        mpfr_add(t, t, h, GMP_RNDN);
        mpfr_sub(temp, t, TS, GMP_RNDN);
        mpfr_abs(temp, temp, GMP_RNDN);
        Compare=mpfr_cmp(temp, Small);
        if(Compare<0)
        {
            mpfr_add(TS, TS, POINTONE, GMP_RNDN);
            Ss++;
            
            for(k=0; k<evel; k++)
            {
                for(j=0; j<evel; j++)
                {
                    mpfr_set(Coe[j][k], Coe_psi_real[j][k], GMP_RNDN);
                    mpfr_set(Coe[j][k+evel], Coe_psi_vir[j][k], GMP_RNDN);
                }
            }
            
            pack_mpf(Coe[0][0], evel2*2, packed_Coe);    /*Gather*/
            if(myid!=0)
            {
                MPI_Send(packed_Coe, evel2*2, MPI_MPF, 0, myid, MPI_COMM_WORLD);
            }
            if(myid==0)
            {
                for(i=1; i<NP; i++)
                {
                    MPI_Recv(packed_CoeALL[i], evel2*2, MPI_MPF, i, i, MPI_COMM_WORLD, &status);
                }
            }
            if(myid==0)
            {
                unpack_mpf(packed_Coe, CoeALL[0][0][0], evel2*2);
                for(i=1; i<NP; i++)
                {
                    unpack_mpf(packed_CoeALL[i], CoeALL[i][0][0], evel2*2);
                }
            }
            
            if(myid==0)
            {
                if(Ss%10==0)
                {
                    sprintf(filename, "ALL_t%d.dat", Ss);
                    fp=fopen(filename, "w+");
                    for(k=0; k<N; k++)
                    {
                        for(j=0; j<N; j++)
                        {
                            mpfr_sprintf(Char1, "%.262Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]]);
                            mpfr_sprintf(Char2, "%.262Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]+evel]);
                            fprintf(fp, "%d\t%d\t%.272s\t%.272s\n", k, j, Char1, Char2);
                        }
                    }
                    fclose(fp);
                }
                else
                {
                    sprintf(filename, "ALL_t%d.dat", Ss);
                    fp=fopen(filename, "w+");
                    for(k=0; k<N; k++)
                    {
                        for(j=0; j<N; j++)
                        {
                            mpfr_sprintf(Char1, "%.5Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]]);
                            mpfr_sprintf(Char2, "%.5Re", CoeALL[NOWid[j][k]][NOWm[j][k]][NOWn[j][k]+evel]);
                            fprintf(fp, "%d\t%d\t%.15s\t%.15s\n", k, j, Char1, Char2);
                        }
                    }
                    fclose(fp);
                }
            }
        }
        
        T_compare=mpfr_cmp(t, T);
    }
    
    et=MPI_Wtime();
    if(myid==0)
    {
        printf("\nCPU's time is %fs\n", et-st);
    }
    
    free_mpf_op(&(MPI_MPF_SUM));
    free_mpf(&(MPI_MPF));
    MPI_Finalize();
    return 0;
}
