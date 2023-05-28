#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <windows.h>
#include <omp.h>
#include <string.h>

#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_poly.h"

#include "Constantfinal.h"

char *name[LAYER];
int N[2],layer,smt[LAYER],divnum[LAYER],NX[LAYER],valley[LAYER],cond[LAYER];
double d[LAYER],die[LAYER],bar[LAYER],mass[LAYER],Ef[LAYER],Q[LAYER],Work[LAYER],massxy[LAYER];
double dx,temp;

int main(int argc, char *argv[])
{
	if(argc < 2){
		printf("aaaaaaa");
		return 0;
	}
	int number = atoi(argv[1]);
	time_t t = time(NULL);
	struct tm *local = localtime(&t);
	char days[50];
	char filename_potential[50];
	char filename_wave[50];
	strftime(days, sizeof(days), "%Y_%m_%d",local);
	sprintf(filename_potential, "%s_%s_%d.csv", "potential", days, number);
	sprintf(filename_wave, "%s_%s_%d.csv", "wave", days, number);

	
	//int p=2; 									/* プログラムの番号				*/
	int j, n, loop; //齋藤が変数loop(ループを回す為の変数)を追加 (2016.10.13)
	double QW,VRTD;
	double *v = calloc(N[0]+1, sizeof(double));
	
	double *vRTD = calloc(N[1]+1, sizeof(double));				/* ポテンシャルの格納用			*/
	QW=0*ELEC*1e4*5e12;
	
	setmaterial(-100);
	VRTD=0.2;
	potential(v,vRTD,VRTD,QW);
	FILE *fp;
	fp = fopen(filename_potential, "w");
	for(j=0; j<N[0]+1; j++){
		fprintf(fp, "%d\t%e\t%e\n", j, j*dx*1e9, v[j]);
	}
	fclose(fp);

	
	//double V = v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);
	n=60;								/* 計算したい準位の数	*/
	double En[n];						/* 準位の格納			*/
	double wavestore[n][DIV*N[0]];  //波動関数格納
	getconfinedstates(n,v,0,En);
	for(j=0;j<n;j++){
		wavefunction(v,En[j],0, wavestore[j]); //齋藤が引数wavestore[j]を追加 (2016.10.13)
	}

	fp = fopen(filename_wave, "w");
	for(j=0; j<DIV*N[0]; j++){
		fprintf(fp, "%g\t", (double)j*dx*1e9/DIV);
		for(loop=0; loop<n; loop++)
			fprintf(fp, "%e\t", wavestore[loop][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	//Beep( 440, 1200 );
	return 0;
}

int mt(int n)
{
	int x;
	for(x=0; x < layer; x++)
	{
		if (n <= NX[x])
			return x;
	}
	return x;
}

int nRTD(int np)
{
	if(np==0)
		return 0;
	else if(np==1)
		return NX[1];
	else
	{
		printf("想定していない値が引き数npに使われています。\n");
		exit(1);
	}
}


void potential(double v[], double vRTD[], double VRTD, double QW)
{
	int n,output;
	double DR,DL;
	selfpotential(QW,VRTD,vRTD);
	for(n=0;n<N[1]+1;n++)
	{
		v[n+NX[1]]=vRTD[n];
	}
	DL=die[LBAR]*(vRTD[1]-vRTD[0])/dx;
	DR=die[RBAR]*(vRTD[N[1]]-vRTD[N[1]-1])/dx;
	for(n=LBAR-1;n>=0;n--)
		DL=calVRL(-1,n,DL,v);
	for(n=RBAR+1;n<=layer;n++)
		DR=calVRL(1,n,DR,v);
	output=OUTPUT;
	setpotential(v,0,output);
	return;
}

void selfpotential(double Q, double VRTD, double vRTD[])
{
	int np=RTD;								/*npはRTD構造か、全体かを表している。np=1はRTD構造。=0は素子全体 */
	int aa=-1;								/*計算終了の判断のために使用*/
	double E;
	potential0(Q,VRTD,vRTD,np);				/*近似式 初期値として使用*/
	setpotential(vRTD,np,0);				/*階段近似適用*/
	gsl_complex A[N[np]],B[N[np]],k[N[np]];	/*Aは前進波、Bは後進波の係数。kは波数 */
	A[0]=gsl_complex_rect(0,0);				/* 初期値の設定                    */
	B[0]=gsl_complex_rect(1,0);				/* 初期値の設定                    */
	while(aa==-1)							/*自己無撞着計算* 計算が収束しない可能性あり。*/
	{
		E=confinedstate(1,vRTD,np);
		cal2(E,vRTD,A,B,np);				/*エネルギーEにおける波動関数の係数AとBを計算*/
		setk(E,vRTD,k,np);					/*波数計算*/
		aa=makewave2(Q,E,VRTD,vRTD,k,A,B,np);/*電荷Qを考慮したポテンシャル計算 出力は、前のポテンシャルと計算後のポテンシャルの差が1e-6以下なら-1　それ以外が1*/
	}
}

void potential0 (double Q, double VRTD, double v[], int np)
{
	int n,nrtd;
	nrtd=NX[LBAR-1];
	double x,DL,DR;
	DL = (VRTD-(d[WELL]/2/die[WELL]+d[RBAR]/die[RBAR])*Q) / (d[LBAR]/die[LBAR]+d[WELL]/die[WELL]+d[RBAR]/die[RBAR]);
	DR = DL+Q;
	for(n=0; n<nrtd+1; n++)
	{
		x=n*dx;
		switch( mt(n+nrtd) )
		{
		case 0: case 1:	v[n] = VRTD; break;
		case 2: v[n] = -DL*x/die[LBAR]+VRTD; break;
		case 3: v[n] = -Q*( pow((x-d[LBAR]),2) + d[WELL]*d[WELL]*(cos(2*PI*(x-d[LBAR])/d[WELL])-1)/2/PIS)/2/die[WELL]/d[WELL] - DL*(x-d[LBAR])/die[WELL]-DL*d[LBAR]/die[LBAR]+VRTD; break;
		case 4: v[n] = -DR*(x-d[LBAR]-d[WELL]-d[RBAR])/die[RBAR]; break;
		case 5: case 6: case 7: v[n] = 0; break;
		}
	}
	return;
}


void setpotential (double v[], int np, int output)
{
	int n,nrtd;								/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
	nrtd=nRTD(np);
	double vr[N[np]+1],vl[N[np]+1];			/* ポテンシャルの左側からの極限vrと右側からの極限			*/
	for(n=0; n<N[np]+1; n++){
		vl[n]=v[n]+bar[mt(n+nrtd)];
		vr[n]=v[n]+bar[mt(n+nrtd+1)];
	}
	for(n=0; n<N[np]; n++){
		v[n]=(vr[n]+vl[n+1])/2;				/* 階段近似適用	*/
	}
	return;
}


void setk(double E, double v[],gsl_complex k[],int np)
{
	int n,nrtd;					/* ポテンシャルから波数kを計算して格納*/
	nrtd=nRTD(np);
	for(n=0; n<N[np]; n++)
	{
		k[n]=gsl_complex_div_real(gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n])*MSTAR*ELEC),HBAR);
	}
}

double cal(double E, double v[], int np, int f) // 出力は、透過率|T|の2乗=1/|T(1,1)|の2乗 npは0が全体、1がRTD内。fは、0のとき、バンドギャップ内も計算する。
{
	int n,nrtd;
	double m,t;
	
	nrtd=nRTD(np);

	if(f==1 && (E-v[N[np]-1]<=0 || E-v[0]<=0))	//バンドギャップ内と判断。計算せず。
		return 0;

	gsl_complex kk,kn,kn1,pp,pm,zp,zm;
	gsl_matrix_complex *temp  = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(temp);
	gsl_matrix_complex *dummy = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(dummy);
	gsl_matrix_complex *trans = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(trans);

	for(n=0; n<N[np]-1; n++)
	{
		if((E-v[n])*(E-v[n+1])!=0){
			m=sqrt(mass[mt(n+nrtd+2)]/mass[mt(n+nrtd+1)]);
			kk=gsl_complex_div(gsl_complex_sqrt_real(E-v[n]),gsl_complex_sqrt_real(E-v[n+1]));

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]を計算
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]を計算

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //エネルギーE、領域n+1における波数を計算 HBARで割っていないことに注意

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)を計算
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)を計算
			
			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0));                         //全ての要素に0.5をかける
		}

		else if(E-v[n+1]==0 && E-v[n]!=0)
		{
			kn=gsl_complex_sqrt_real(2*mass[mt(n+nrtd)]*(E-v[n])*MSTAR*ELEC/(HBAR*HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];

			zp=gsl_complex_mul_imag(kn, m);
			zm=gsl_complex_mul_imag(kn,-m);

			gsl_matrix_complex_set(temp,0,0,zp);
			gsl_matrix_complex_set(temp,0,1,zm);
			
			zp=gsl_complex_mul_imag(kn, m*dx);
			zm=gsl_complex_mul_imag(kn,-m*dx);

			gsl_matrix_complex_set(temp,1,0,gsl_complex_add_real(zp,1));
			gsl_matrix_complex_set(temp,1,1,gsl_complex_add_real(zm,1));

		}

		else if(E-v[n]==0 && E-v[n+1]!=0)
		{
			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC);
			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];
			kn1=gsl_complex_sqrt_real(HBAR*HBAR/(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC));
			pp=gsl_complex_mul_imag(kn1, m);
			pm=gsl_complex_mul_imag(kn1,-m);

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //全ての要素に0.5をかける
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,dummy,gsl_complex_rect(0,0),trans); //行列の掛け算

		gsl_matrix_complex_memcpy(dummy, trans);
	}
	if(f==1)
	{
		if(E-v[N[np]-1]>0 && E-v[0]>0)
			t=sqrt ( (mass[layer]/mass[0]) * ( (E-v[0]) / (E-v[N[np]-1]) ) );
		else
			t=0;
	}
	else
		t=1;
	t*=1/gsl_complex_abs2(gsl_matrix_complex_get(trans,1,1));
	
	gsl_matrix_complex_free(temp);
	gsl_matrix_complex_free(dummy);
	gsl_matrix_complex_free(trans);

	return t;
}

void cal2(double E, double v[], gsl_complex A[], gsl_complex B[], int np)
{
	int n,nrtd;
	double m;
	
	nrtd=nRTD(np);

	gsl_complex kk,kn,kn1,pp,pm,zp,zm;
	gsl_matrix_complex *temp  = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(temp);
	gsl_matrix_complex *FN = gsl_matrix_complex_calloc(2,1);
	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);

	for(n=0; n<N[np]-1; n++)
	{
		gsl_matrix_complex_set(FN,0,0,A[n]);
		gsl_matrix_complex_set(FN,1,0,B[n]);

		if((E-v[n])*(E-v[n+1])!=0){
			m=sqrt(mass[mt(n+nrtd+2)]/mass[mt(n+nrtd+1)]);
			kk=gsl_complex_div(gsl_complex_sqrt_real(E-v[n]),gsl_complex_sqrt_real(E-v[n+1]));

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]を計算
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]を計算

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //エネルギーE、領域n+1における波数を計算 HBARで割っていないことに注意

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)を計算
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)を計算
			
			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0));                         //全ての要素に0.5をかける
		}

		else if(E-v[n+1]==0 && E-v[n]!=0)
		{
			kn=gsl_complex_sqrt_real(2*mass[mt(n+nrtd)]*(E-v[n])*MSTAR*ELEC/(HBAR*HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];

			zp=gsl_complex_mul_imag(kn, m);
			zm=gsl_complex_mul_imag(kn,-m);

			gsl_matrix_complex_set(temp,0,0,zp);
			gsl_matrix_complex_set(temp,0,1,zm);
			
			zp=gsl_complex_mul_imag(kn, m*dx);
			zm=gsl_complex_mul_imag(kn,-m*dx);

			gsl_matrix_complex_set(temp,1,0,gsl_complex_add_real(zp,1));
			gsl_matrix_complex_set(temp,1,1,gsl_complex_add_real(zm,1));

		}

		else if(E-v[n]==0 && E-v[n+1]!=0)
		{
			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC);
			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];
			kn1=gsl_complex_sqrt_real(HBAR*HBAR/(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC));
			pp=gsl_complex_mul_imag(kn1, m);
			pm=gsl_complex_mul_imag(kn1,-m);

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //全ての要素に0.5をかける
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納
		}

		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,FN,gsl_complex_rect(0,0),FNN); //行列の掛け算 FNN=temp*FN

		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);

		gsl_matrix_complex_memcpy(FN, FNN);
	}
	
	gsl_matrix_complex_free(FN);
	gsl_matrix_complex_free(FNN);
	gsl_matrix_complex_free(temp);

	return;
}

void cal3(double E, double v[], gsl_complex A[], gsl_complex B[], int np)
{
	int i,n,nrtd;
	double m;

	nrtd=nRTD(np);	
	gsl_complex kk,kn1,pp,pm,zp,zm;
	gsl_matrix_complex *temp[N[np]];
	for(n=0; n<N[np]; n++){
		temp[n]=gsl_matrix_complex_alloc(2,2);
	}
	gsl_matrix_complex *dummy = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(dummy);
	gsl_matrix_complex *trans = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(trans);
	for(n=0; n<N[np]-1; n++)
	{
		if((E-v[n])*(E-v[n+1])!=0){
			m=sqrt(mass[mt(n+nrtd+2)]/mass[mt(n+nrtd+1)]);
			kk=gsl_complex_div(gsl_complex_sqrt_real(E-v[n]),gsl_complex_sqrt_real(E-v[n+1]));

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]を計算
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]を計算

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //エネルギーE、領域n+1における波数を計算 HBARで割っていないことに注意

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)を計算
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)を計算
			
			gsl_matrix_complex_set(temp[n],0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp[n],0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp[n],1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp[n],1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp[n], gsl_complex_rect(0.5,0));                         //全ての要素に0.5をかける
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],dummy,gsl_complex_rect(0,0),trans); //行列の掛け算

		gsl_matrix_complex_memcpy(dummy, trans);
	}
	if(v[0]<E)
	{
		if(v[N[np]]<E)
			i=0;
		else
			i=0;
	}
	else if(v[0]>E)
		i=1;
	switch(i)
	{
	case 0:	A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_negative( gsl_complex_div( gsl_matrix_complex_get(trans,1,0), gsl_matrix_complex_get(trans,1,1) ) );
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
	}
	gsl_matrix_complex *FN  = gsl_matrix_complex_calloc(2,1);
	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);
	
	for(n=0; n<N[np]-1; n++)
	{
		gsl_matrix_complex_set(FN,0,0,A[n]);
		gsl_matrix_complex_set(FN,1,0,B[n]);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],FN,gsl_complex_rect(0,0),FNN); //行列の掛け算 FNN=temp*FN
		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);
	}
	for(n=0; n<N[np]; n++){
		gsl_matrix_complex_free(temp[n]);
	}
	gsl_matrix_complex_free(dummy);
	gsl_matrix_complex_free(trans);
}

double confinedstate(int num, double v[], int np)		/* 引数numは、量子数nのこと。nは別のところで使われているのでnumにした。 */
{														/*バグあり。恐らく、左右の一番端のエネルギーが準位として検出されてしまう。*/
	int i,j,n;
	double E,T,S,dE;
	double emax=5;                                    /* whileの無限ループを避けるため */
	double min=5;

	for(n=0;n<N[np]+1;n++)		
	{
		if(min>v[n])
			min=v[n];                                   /* 分割されたポテンシャルの最小値を求める */
	}
	E=min;
	for(i=1;i<num+1;i++)
	{
		dE=1e-3;                                      /* 刻み幅(粗い)の設定                                              */
		T=j=0;
		S=cal(E,v,np,0);
		while(E<emax)                                 /* 透過率の減少区間の検索                                          */
		{                                             /* 以下の準位を求めるプログラムは                                  */
			E+=dE;                                    /* エネルギーが上昇していくとともに、透過率が上昇していくと仮定し、*/
			T=cal(E,v,np,0);                          /* その最大値をとる箇所が準位であるとしている                      */
			if(S>T)                                   /* よって、透過率が減少していく区間は不要                          */
				S=T;
			else
				break;                                /* 透過率の減少区間の終了 */
		}
		S=T=0;                                        /* S、Tの初期化                          */
		while(E<emax)                                 /* 基底準位の探索 刻み幅：粗い           */
		{
			E+=dE;                                    /* エネルギーの設定                      */
			T=cal(E,v,np,0);                              /* 透過率の計算                          */
			if(S<T)                                   /* 増減の判断　増加関数ならばwhile文続行 */
				S=T;
			else
			{
				E-=2*dE;                              /* エネルギーを少し戻して、刻み幅の細かい方へ移行 */
				dE/=10;
				//j++;
				S=0;
				if(dE<1e-10)
					break;
			}
		}
		dE/=10;                                      /* 刻み幅(細かい)の設定                           */
		S=0;                                          /* Sの初期化                                      */
		while(E<emax)                                 /* 基底準位の探索 刻み幅：細かい                  */
		{
			E+=dE;                                    /* エネルギーの設定                               */
			T=cal(E,v,np,0);                              /* 透過率の計算                                   */
			if(S<T)
				S=T;                                  /* 増減の判断　増加関数ならばwhile文続行          */
			else{
				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
				{										/*cal関数の影響で、E=v[0]のところが準位に見えてしまうので*/
					E+=10*dE;							/*それを避けるために導入*/
					i--;								/*従って、準位がv[0]±1e-10[eV]の近辺に存在した場合、*/
					break;								/*検知できない可能性あり*/
				}
				if(i==num)
					return E-dE;
				E+=10*dE;
				break;                                /* 透過率の増加区間の終了 その一個前の計算結果を  */
			}
		} /* 簡略化したい */
	}
	E-=dE;
	return E;                                         /* num番目の準位を返す                           */
}

void makewave(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np, double wavestore[]) //齋藤が引数wavestoreを追加 (2016.10.13)
{
	int i,j,n;
	double max=0;                     // 規格化用
	double sum=0;
	double xn,px,py, pyy[N[np]*DIV]; //齋藤が配列追加(pyy, 2016.10.13)
	gsl_complex tempA,tempB;
	j=max=0;//p[0]=-ML;q[0]=E;
	for(n=0;n<N[np];n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=py*dx*1e9/DIV;
			j++;
		}
	}j=0;//2015/5/20追加　思考停止 とりあえず波動関数の面積を求める(↑のmaxに格納)。波動関数の2乗の全区間で積分すると1になるように規格化。
	 //量子井戸の場合、膜厚が厚すぎるとexp(±ikz)の掛け算でオーバーフロー(？)を起こし、減衰項ではなく増幅項が支配的になるので注意。
	for(n=0;n<N[np];n++)
	{

	//ここから先，出力の仕方を大幅に変更(旧QCL.cとほぼ同じにした) by齋藤 (2016.10.13)

		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			pyy[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			if(pyy[j]>max)
				max=pyy[j];
			j++;
		}
	}
	for(i=0; i<N[np]*DIV; i++)
		wavestore[i] = E+pyy[i]/max;
	return;
}

int makewave2(double Q, double E, double VRTD, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[],int np)
{
	int i,j,n;
	double max=0;                     // 規格化用
	double qmax=0;
	double xn,p[DIV*N[np]],q[DIV*N[np]]; // pがx座標、qがy座標 qは波動関数の絶対値の2乗を出力予定
	double qtemp[N[np]+1], vnew[N[np]];
	gsl_complex tempA,tempB;
	j=0;//	p[0]=-ML;
	for(n=0; n<N[np]; n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			p[j]=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],p[j]-xn)));     /* tempA= A[n] × exp( ik[n](x-x[n])) */
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-p[j])));     /* tempB= B[n] × exp(-ik[n](x-x[n])) */
			q[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));                                /* 波動関数の絶対値の2乗を計算        */
			qmax+=q[j]/DIV;
			j++;
		}
	}
	j=0;
	for(n=1;n<N[np];n++)
	{
		qtemp[n]=-Q*q[DIV*n]/qmax;
	}
	qtemp[0]=qtemp[N[np]]=0;
	potential0(0,VRTD,vnew,np);
	calpotential(vnew,qtemp,np);
	double opv[N[np]+1];
	for(n=0;n<N[np]+1;n++)
		opv[n]=vnew[n];
	setpotential(vnew,np,0);
	for(n=1; n<N[np]-1; n++)
	{
		if(fabs(v[n]-vnew[n])>max)
			max=fabs(v[n]-vnew[n]);
	}
	for(n=0; n<N[np]; n++)
			v[n]=vnew[n];
	if(max<1e-6){			/* 重要　自己無撞着計算の収束判定*/
		for(n=0;n<N[np];n++)
			v[n]=opv[n];
		v[N[np]]=0;
		return 1;
	}
	else
		return -1;
	//return max;
}

void calpotential(double vnew[], double q[],int np)
{
	int n,nrtd;
	int div;
	double temp,dietemp;
	
	nrtd=nRTD(np);
	div=N[np]-1;

	gsl_vector *V = gsl_vector_alloc (div);
	gsl_vector *X = gsl_vector_alloc (div);
	gsl_matrix *S = gsl_matrix_calloc(div,div);
	gsl_permutation *P = gsl_permutation_calloc(div);
	
	for(n=0;n<div;n++)
	{
		temp=q[n+1]*dx;
		if(n+nrtd+1==NX[2])
			dietemp=(die[LBAR]+die[WELL])/2;
		else if(n+nrtd+1==NX[3])
			dietemp=(die[WELL]+die[RBAR])/2;
		else
			dietemp=die[mt(n+nrtd+1)];

		if(n==0){
			gsl_matrix_set(S, 0, 0, -2*dietemp);
			gsl_matrix_set(S, 1, 0,  1*dietemp);
			gsl_vector_set(V, 0   , temp);
		}else if(n==div-1){
			gsl_matrix_set(S, n-1, n,  1*dietemp);
			gsl_matrix_set(S, n  , n, -2*dietemp);
			gsl_vector_set(V, n  , temp);
		}else{
			gsl_matrix_set(S, n-1, n,  1*dietemp);
			gsl_matrix_set(S, n  , n, -2*dietemp);
			gsl_matrix_set(S, n+1, n,  1*dietemp);
			gsl_vector_set(V, n  , temp);
		}
	}
	gsl_linalg_LU_decomp(S, P, &n);
	gsl_linalg_LU_solve(S, P, V, X);

	for(n=0;n<div;n++)
	{
		vnew[n+1]+=gsl_vector_get(X,n);
	}
	gsl_vector_free(V);
	gsl_vector_free(X);
	gsl_matrix_free(S);
	gsl_permutation_free(P);

	return;
}


double calVRL(int direction, int n, double D, double v[])
{
	int sm,la,i,j;
	double A,E,tempV,fai,q,w,wmax,x;
	
	la=NX[n];
	if(n-1<0)
		sm=0;
	else
		sm=NX[n-1];
	if(smt[n]==7 || smt[n]==5 || smt[n]==23 || smt[n]==24 || smt[n]==22 )//i-Si smt[n]==24を追加 07/28大野 8/19ていがn=22を追加
	{
		E=D/die[n];
		if(direction==-1)
		{
			for(i=la;i>=sm;i--)
				v[i-1]=v[i]-E*dx;
		}
		else
		{
			for(i=sm;i<=la;i++)
				v[i+1]=v[i]+E*dx;
		}
		return D;
	}
	if(cond[n]==0)//金属
	{
		if(direction==-1)
		{
			for(i=la;i>=sm;i--)
				v[i-1]=v[n];
		}
		else
		{
			for(i=sm;i<=la;i++)
				v[i+1]=v[NX[n-1]];
		}
		return 0;
	}
	if(cond[n]==1)//n-SiSub.
	{
		if(direction<0)
				j=la;
			else
				j=sm;
		if( (direction<0 && D>=0) || (direction>0 && D<=0)) //空乏層
		{
			A=ELEC*Q[n]*1e6/2/(die[n]);
			w=fabs(D/ELEC/(Q[n]*1e6));
			wmax=sqrt(4*die[n]*KB*temp*log(Q[n]*1e6/NI)/(ELEC*ELEC*Q[n]*1e6));
			x=0;
			if(1)//w<d[n] || wmax<d[n])
			{
					if(fabs(w)>wmax)
					{
						w=wmax;
					}
					if(direction<0)
					{
						for(i=la;i>=sm;i--)
						{
							if(x<w)
								v[i]=A*x*(x-2*w)+v[j];
							else
								v[i]=-A*w*w+v[j];
							x+=dx;
						}
					}
					else
					{
						for(i=sm;i<=la;i++)
						{
							if(x<w)
								v[i]=A*x*(x-2*w)+v[j];
							else
								v[i]=-A*w*w+v[j];
							x+=dx;
						}
					}
			}
			else //エラー処理
			{
				printf("\n%d層目:空乏層幅が膜厚以上です。w=%e\n膜厚を増やして下さい。\n空乏層幅Max=%e",n,w,wmax);
				int wML=wmax/ML;
				printf("=%6.4gML<%dML \n現在の膜厚d[%d]=%e\n",wmax/ML,wML+1,n,d[n]);
				exit(1);
			}
			return 0;
		}
		else //蓄積層
		{
			A=sqrt(2*KB*temp*Q[n]*1e6/die[n]);
			q=ELEC/(KB*temp);
			fai=-vacc(D,Q[n]);
			E=A*sqrt(exp(-q*fai)+q*fai-1);
			tempV=0;
			if(direction>0)
			{
				for(i=sm;i<=la;i++)
				{
					v[i]=tempV+v[j];
					fai+=E*dx;
					tempV+=direction*E*dx;
					E=A*sqrt(exp(-q*fai)+q*fai-1);
				}
			}
			else
			{
				for(i=la;i>=sm;i--)
				{
					v[i]=tempV+v[j];
					fai+=E*dx;
					tempV-=direction*E*dx;
					E=A*sqrt(exp(-q*fai)+q*fai-1);
				}
			}
		}
	}
	if(cond[n]==2) //p-Si Sub.
	{
		if(direction<0)
				j=la;
			else
				j=sm;
		if( (direction<0 && D<=0) || (direction>0 && D>=0) ) //空乏層
		{
			A=-ELEC*Q[n]*1e6/2/die[n];
			w=fabs(D/ELEC/(Q[n]*1e6));
			wmax=sqrt(4*die[n]*KB*temp*log(Q[n]*1e6/NI)/(ELEC*ELEC*Q[n]*1e6));
			x=0;
			if(w<d[n] || wmax<d[n])
			{
					if(fabs(w)>wmax)
					{
						w=wmax;
					}
					if(direction<0)
					{
						for(i=la;i>=sm;i--)
						{
							if(x<w)
								v[i]=A*x*(x-2*w)+v[j];
							else
								v[i]=-A*w*w+v[j];
							x+=dx;
						}
					}
					else
					{
						for(i=sm;i<=la;i++)
						{
							if(x<w)
								v[i]=A*x*(x-2*w)+v[j];
							else
								v[i]=-A*w*w+v[j];
							x+=dx;
						}
					}
			}
			else //エラー処理
			{
				printf("\n%d層目:空乏層幅が膜厚以上です。w=%e\n膜厚を増やして下さい。\n空乏層幅Max=%e",n,w,wmax);
				int wML=wmax/ML;
				printf("=%6.4gML<%dML\n",wmax/ML,wML+1);
				exit(1);
			}
			return 0;
		}
		else //蓄積層
		{
			A=sqrt(2*KB*temp*Q[n]*1e6/die[n]);
			q=ELEC/(KB*temp);
			fai=-vacc(D,Q[n]);
			E=A*sqrt(exp(q*fai)-q*fai-1);
			tempV=0;
			if(direction>0)
			{
				for(i=sm;i<=la;i++)
				{
					v[i]=tempV+v[j];
					fai+=E*dx;
					tempV-=direction*E*dx;
					E=A*sqrt(exp(q*fai)-q*fai-1);
				}
			}
			else
			{
				for(i=la;i>=sm;i--)
				{
					v[i]=tempV+v[j];
					fai+=E*dx;
					tempV+=direction*E*dx;
					E=A*sqrt(exp(-q*fai)+q*fai-1);
				}
			}
		}
	}
	return 0;
}

double vacc(double D, double nd)
{
	int loop=1;
	double a,q,xk;
	double Vacc=0.1;
	double sub=1;
	double f1,f2;

	nd*=1e6*ELEC;
	q=ELEC/KB/temp;
	a=q*D*D/2/die[1]/nd;
	if(D==0)
		return 0;
	while(sub>1e-6)
	{
		if(Vacc>20 && loop<=5)
		{
			loop++;
			Vacc=0.1*loop;
		}
		else if(loop>5)
		{
			printf("蓄積層のポテンシャル計算が収束しませんでした。\nドーピング濃度や極端に高い電圧をかけていないか等確認して下さい。");
			exit(1);
		}
		xk=Vacc;
		f1=exp(q*xk)-q*xk-1-a;
		f2=q*(exp(q*xk)-1);
		Vacc=xk-f1/f2;
		sub=fabs(Vacc-xk);
	}
	return Vacc;
}

void wavefunction(double v[], double E,int np, double wavestore[]) //齋藤が引数wavestoreを追加 (2016.10.13)
{
	gsl_complex A[N[np]],B[N[np]],k[N[np]];
	cal3(E,v,A,B,np);                               /* 波動関数の計算                  */
	setk(E,v,k,np);                                 /* 波数の計算                      */		
	makewave(E,v,k,A,B,np, wavestore);                         /* 波動関数の計算と表示            */
}

double func(double E, double v, int n)
{
	if(fabs(v-E)<1e-15)
		return 0;
	double t;
	t=1.0+exp(ELEC*(v+Ef[n]-E)/KB/temp);
	return log(t);
}

double fermia(int n, double t, int flag)			/* フェルミエネルギーの計算　要検討							*/
{													/* 本来は、温度を変数とした計算を行いたかったが、現在凍結中	*/
	double Ed=-0.05; //ドナーorアクセプター準位。本来ならドープ材料によって変化?今は50meVとしている。

	if(flag==1 || flag==2 || flag==19 || flag==20)// n-Si用
	{
		double ef;
		ef=-EG_SI/2+KB*t*log((Q[n]*1e6)/NI)/ELEC;
		return ef;
	}
	if(flag==3 || flag==4 || flag==10)//p-Si用
	{
		double nv,ef;
		nv=2*pow((MSTAR*KB*t/HBAR/HBAR/PI/2.0),1.5)*(pow(0.16,1.5)+pow(0.49,1.5));
		ef=Ed+KB*t*log(-1.0/4.0+pow(1+8*Q[n]*1e6*exp(-Ed/KB/t)/nv,0.5)/4);
		return ef/ELEC;
	}
	return -0.56;
}

double calconfinedstate(double v[], int np, double E1, double E2)
{
	/*	E1からE2のエネルギー範囲において、最も低い準位を探す	*/
	/*	ない場合には、EMAXを返す。								*/
	/*	変数	v	ポテンシャル								*/
	/*			np	0なら空間全体、1ならRTD部分のみ				*/
	double E,T,S,dE;
	E=E1;											/* 探索範囲の下限													*/
	while(E<E2)
	{											
		dE=1e-3;									/* 刻み幅(粗い)の設定                                              */
		T=0;
		S=cal(E,v,np,0);
		while(E<E2)									/* 透過率の減少区間の検索                                          */
		{											/* 以下の準位を求めるプログラムは                                  */
			E+=dE;									/* エネルギーが上昇していくとともに、透過率が上昇していくと仮定し、*/
			T=cal(E,v,np,0);						/* その最大値をとる箇所が準位であるとしている                      */
			if(S>T)									/* よって、透過率が減少していく区間は不要                          */
				S=T;
			else
				break;								/* 透過率の減少区間の終了 */
		}
		S=T=0;										/* S、Tの初期化                          */
		while(E<E2)									/* 基底準位の探索 刻み幅：粗い           */
		{
			E+=dE;									/* エネルギーの設定                      */
			T=cal(E,v,np,0);						/* 透過率の計算                          */
			if(S<T)									/* 増減の判断　増加関数ならばwhile文続行 */
				S=T;
			else
			{
				E-=2*dE;							/* エネルギーを少し戻して、刻み幅を細かくしていく */
				dE/=10;
				S=0;
				if(dE<1e-10)						/*　*/
					break;
			}
		}
		dE/=10;										/* 刻み幅(細かい)の設定                           */
		S=0;										/* Sの初期化                                      */
		while(E<E2)									/* 基底準位の探索 刻み幅：細かい                  */
		{
			E+=dE;									/* エネルギーの設定                               */
			T=cal(E,v,np,0);						/* 透過率の計算                                   */
			if(S<T)
				S=T;								/* 増減の判断　増加関数ならばwhile文続行          */
			else
			{
				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
				{									/*cal関数の影響で、E=v[0]のところが準位に見えてしまうので*/
					E+=10*dE;						/*それを避けるために導入*/
					break;							/*検知できない可能性あり*/
				}
				E-=dE;
				return E;
			}
		}
	}
	E2=EMAX;
	return E2;
}

double getconfinedstate(int n, double v[], int np)
{											/* 第n準位を出力												*/
	int i;
	double E,E1,E2;
	E1=E2=EMAX;
	for(i=0;i<N[np]+1;i++)										
	{
		if(E1>v[i])							/* エネルギーの探索範囲の決定									*/
			E1=v[i];						/* 分割されたポテンシャルの最小値を求める						*/
	}
	for(i=0;i<n;i++)
	{
		E=calconfinedstate(v,np,E1,E2);
		E1=E+1e-9;
	}
	return E;
}

void getconfinedstates(int n, double v[], int np, double En[])
{											/* 複数の準位を取得し、配列に出力								*/
	int i;
	double E1,E2;
	E1=E2=EMAX;
	for(i=0;i<N[np]+1;i++)										
	{
		if(E1>v[i])							/* エネルギーの探索範囲の決定									*/
			E1=v[i];						/* 分割されたポテンシャルの最小値を求める						*/
	}
	for(i=0;i<n;i++)
	{
		En[i]=calconfinedstate(v,np,E1,E2);
		E1=En[i]+1e-9;
	}
}

double cal4(double E, double v[], int np)
{
	int i,j,n,nrtd;
	double m;
	gsl_complex A[N[np]],B[N[np]],k[N[np]];

	nrtd=nRTD(np);	
	gsl_complex kk,kn1,pp,pm,zp,zm;
	gsl_matrix_complex *temp[N[np]];
	for(n=0; n<N[np]; n++){
		temp[n]=gsl_matrix_complex_alloc(2,2);
	}
	gsl_matrix_complex *dummy = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(dummy);
	gsl_matrix_complex *trans = gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex_set_identity(trans);
	for(n=0; n<N[np]-1; n++)
	{
		if((E-v[n])*(E-v[n+1])!=0){
			m=sqrt(mass[mt(n+nrtd+2)]/mass[mt(n+nrtd+1)]);
			kk=gsl_complex_div(gsl_complex_sqrt_real(E-v[n]),gsl_complex_sqrt_real(E-v[n+1]));

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]を計算
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]を計算

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //エネルギーE、領域n+1における波数を計算 HBARで割っていないことに注意

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)を計算
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)を計算
			
			gsl_matrix_complex_set(temp[n],0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp[n],0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp[n],1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp[n],1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp[n], gsl_complex_rect(0.5,0));                         //全ての要素に0.5をかける
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],dummy,gsl_complex_rect(0,0),trans); //行列の掛け算

		gsl_matrix_complex_memcpy(dummy, trans);
	}
	if(v[0]<E)
	{
		if(v[N[np]]<E)
			i=0;
		else
			i=0;
	}
	else if(v[0]>E)
		i=1;
	switch(i)
	{
	case 0:	A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_negative( gsl_complex_div( gsl_matrix_complex_get(trans,1,0), gsl_matrix_complex_get(trans,1,1) ) );
			break;//左側から右側に透過する場合
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
			break;//右側から左側に透過する場合
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
			break;//要検討。v[0]<v[N+1]のとき適用可能？2015/3/25
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
			break;//簡易版。透過率が1の場合適用できる？
	}
	gsl_matrix_complex *FN  = gsl_matrix_complex_calloc(2,1);
	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);
	
	for(n=0; n<N[np]-1; n++)
	{
		gsl_matrix_complex_set(FN,0,0,A[n]);
		gsl_matrix_complex_set(FN,1,0,B[n]);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],FN,gsl_complex_rect(0,0),FNN); //行列の掛け算 FNN=temp*FN
		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);
	}
	for(n=0; n<N[np]; n++){
		gsl_matrix_complex_free(temp[n]);
	}
	gsl_matrix_complex_free(dummy);
	gsl_matrix_complex_free(trans);

	setk(E,v,k,np);

	//double max=0;                     // 規格化用
	double sum=0;
	double xn,px,py;
	//double pr,pi;
	gsl_complex tempA,tempB;
	j=0;//p[0]=-ML;q[0]=E;
	for(n=0;n<N[np];n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=(double)py*dx/DIV;
			j++;
		}
	}
	return sum*func(E,v[0],0)*sqrt(1/(E-v[0]));
}
void setmaterial(int m)//m=-100だと、出力されない。実行はされる。引き数mは、物性値を変更したい場合に使用。
{
	struct material{
		double die;		/* 誘電率					*/
		double bar;		/* バンド不連続				*/
		double mass;	/* 有効質量					*/
		double massxy;	/* 横方向の有効質量			*/
		double ef;		/* フェルミエネルギー		*/
		int valley;		/* 伝導帯の谷の数			*/
		char *name;		/* 物質の名前				*/
		int cond;
	};

	int i,j;
	double base;		/* 伝導帯エネルギーの基準。大抵はSiの電子親和力	*/
	struct material mtconst[100];/* 物質の設定。配列を十分確保しないとエラーが起きるので注意。*/
	temp=TEMP;			/* 温度						*/
	int flag=-100;		/* 引き数がこの値と同じなら、実効はされるが出力はされなくなる。基本的には出力した方がいいと思う。*/
	/* 各物質における物性値の設定。name:名前、die：誘電率、bar：バンド不連続、mass:有効質量、valley:谷の数。valley:基本1。cond:導電性。0が金属、1がn-Si、2がp-Si、3が絶縁物。*/
	i=0; mtconst[i].name="Al"; mtconst[i].die=DIE_AL   * DIEELECSTAR; mtconst[i].bar=BAR_AL ;  mtconst[i].mass=MASS_AL;   mtconst[i].massxy=MASS_AL;		mtconst[i].valley=1;			mtconst[i].cond=0;
	i=1; mtconst[i].name="n-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=1;
	i=2; mtconst[i].name="n-Si"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=1;
	i=3; mtconst[i].name="p-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=2;
	i=4; mtconst[i].name="p-Si"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=2;
	i=5; mtconst[i].name="CaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=BAR_CAF2; mtconst[i].mass=MASS_CAF2; mtconst[i].massxy=MASS_CAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=6; mtconst[i].name="CdF2"; mtconst[i].die=DIE_CDF2 * DIEELECSTAR; mtconst[i].bar=BAR_CDF2; mtconst[i].mass=MASS_CDF2; mtconst[i].massxy=MASS_CDF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=7; mtconst[i].name="i-Si"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=3;
	i=8; mtconst[i].name="Au"; mtconst[i].die=DIE_AU   * DIEELECSTAR; mtconst[i].bar=BAR_AU ;  mtconst[i].mass=MASS_AU;   mtconst[i].massxy=MASS_AU;		mtconst[i].valley=1;			mtconst[i].cond=0;
	i=9; mtconst[i].name="nSub.CaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=BAR_CAF2; mtconst[i].mass=MASS_CAF2; mtconst[i].massxy=MASS_CAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=10;mtconst[i].name="p-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_pSI;  mtconst[i].mass=MASS_pSI;  mtconst[i].massxy=MASS_pSI;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=11;mtconst[i].name="pCaF2"; mtconst[i].die=DIE_pCAF2* DIEELECSTAR; mtconst[i].bar=BAR_pCAF2;mtconst[i].mass=MASS_pCAF2;mtconst[i].massxy=MASS_pCAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=12;mtconst[i].name="pAl"; mtconst[i].die=DIE_pAL  * DIEELECSTAR; mtconst[i].bar=BAR_pAL;  mtconst[i].mass=MASS_pAL;  mtconst[i].massxy=MASS_pAL;	mtconst[i].valley=1;			mtconst[i].cond=0;
	i=13;mtconst[i].name="5MLCaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.7;      mtconst[i].mass=0.7;       mtconst[i].massxy=0.7;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=14;mtconst[i].name="仮想AL"; mtconst[i].die=DIE_AL   * DIEELECSTAR; mtconst[i].bar=-5.17  ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=15;mtconst[i].name="仮想CaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=7      ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=16;mtconst[i].name="仮想p-Si"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;      mtconst[i].massxy=0.55;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=17;mtconst[i].name="SiO2"; mtconst[i].die=DIE_SIO2 * DIEELECSTAR; mtconst[i].bar=BAR_SIO2; mtconst[i].mass=MASS_SIO2; mtconst[i].massxy=MASS_SIO2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=18;mtconst[i].name="nc-Si."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=19;mtconst[i].name="nSi100Z"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_Z; mtconst[i].massxy=MASS_SI_XY;	mtconst[i].valley=NUMVALLY_SI_Z;	mtconst[i].cond=1;
	i=20;mtconst[i].name="nSi100XY"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_XY;mtconst[i].massxy=MASS_SI_Z;	mtconst[i].valley=NUMVALLY_SI_XY;	mtconst[i].cond=1;
	i=21;mtconst[i].name="n-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;
	i=22;mtconst[i].name="CaF2(4ML)"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.6;	  mtconst[i].mass=0.85;	     mtconst[i].massxy=0.85;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=23;mtconst[i].name="p-i-Si"; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;	     mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;	
	i=24;mtconst[i].name="3MLCaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.5;      mtconst[i].mass=1.0;       mtconst[i].massxy=1.0;		mtconst[i].valley=1;			mtconst[i].cond=3;
	/* 物質を追加したい場合は、この行の上をコピーして貼り付け。iを変更するのとmtconstの配列に収まるように注意。*/
	if(DX>0)
		dx=ML/DX;	/* DXは分割数	MLはこのプログラムのz軸の基本単位。そのMLをさらにDX分割する。*/
	else			/* ML=1e-9にすれば、z[nm]。ML=0.31e-9ならz[ML]	*/
	{
		printf("分割数DXが0以下です。適切な値を設定して下さい。\n");
		exit(1);
	}
	N[0]=0;			/* RTD全体の層数*/
	FILE *fd1;
	if((fd1 = fopen("set.dat", "r")))
	{
		if(m!=flag)
			printf("層番号\t材料番号\t物質名\tML数\t層厚[nm]\t比誘電率\t有効質量\t障壁の高さ\t谷\t分割数\tNX\tEF\t仕事関数\t電荷量\n");
		for(i=0; fscanf(fd1, "%d %d %lf\n", &smt[i], &j, &Q[i])!=EOF; i++)
		{
			if( modf((double)j*DX,&d[i])!=0)
			{
				printf("分割エラー　分割間隔DXを変更して下さい\n");
				exit(1);
			}
			divnum[i] = j*DX;					/* 分割数					*/
			N[ALL]        += divnum[i];			/* 総分割数					*/
			NX[i]     = N[0];
			d[i]      = j*ML;					/* 膜厚[m]					*/
			name[i]   = mtconst[smt[i]].name;	/* 材料名					*/
			die[i]    = mtconst[smt[i]].die;	/* 各層の誘電率[F/m]		*/
			bar[i]    = mtconst[smt[i]].bar;	/* 各層のバンド不連続ΔEc	*/
			mass[i]   = mtconst[smt[i]].mass;	/* 各層の有効質量			*/
			massxy[i] = mtconst[smt[i]].massxy;	/* 各層の横方向の有効質量	*/
			valley[i] = mtconst[smt[i]].valley;	/* 伝導帯の谷の数			*/
			cond[i]	  = mtconst[smt[i]].cond;
			layer=i;							/* 層番号					*/
			switch(smt[i]){/*フェルミエネルギーの設定。金属は物性値、半導体は教科書の式から計算、絶縁物には対応していない。*/
				case 0 : Ef[i]=-0.17+EF_AL;base=BASE;break;
				case 1 : Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;//-EG_SI/2+KT*log(Q[i]*1e6/NI)/ELEC;base=BASE;break;
				case 2 : Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 3 : Ef[i]=-EG_SI-fermia(i,temp,smt[i]);base=BASE;break;
				case 4 : Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;//-EG_SI/2-KT*log(Q[i]*1e6/NI)/ELEC;base=BASE;break;
				case 7 : Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 8 : Ef[i]=EF_AU;base=BASE;break;
				case 9 : Ef[i]=-0.3;base=BASE;break;
				case 10: Ef[i]=fermia(i,temp,smt[i]);base=pBASE;break;
				case 12: Ef[i]=EF_pAL;base=pBASE;break;
				case 14: Ef[i]=4.37;base=pBASE;break;
				case 16: Ef[i]=-EG_SI/2+KB*temp*log(Q[i]*1e6/NI)/ELEC;base=pBASE;break;
				case 19: Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 20: Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 23: Ef[i]=fermia(i,temp,smt[i]);base=pBASE;break;
				default: Ef[i]=0;base=pBASE;break;
			}
			switch(mtconst[smt[i]].cond){ //2016.09.30追加．絶縁体は仕事関数を計算させないようにした．
				case 0 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 1 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 2 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 3 : break;
			}
	//			Work[i]   = -Ef[i]-bar[i]+base;		/* 仕事関数。当然、絶縁物には対応していない。*/
			if(m!=flag)
				printf("%d\t%d\t%s\t%d\t%.4g\t%.3g\t%.2g\t%.4g\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\n", i, smt[i], name[i], j, j*ML*1e9, die[i]/DIEELECSTAR, mass[i], bar[i], valley[i], divnum[i], NX[i],Ef[i],Work[i],Q[i]);//出力。
		}
		fclose(fd1);
	}
	else 
	{
		printf("ファイルが開けていません\n");
		exit(1);
	}
	N[1]=divnum[LBAR]+divnum[WELL]+divnum[RBAR];	/* RTD構造のみの層数*/
	if(m!=flag){
		printf("\n1ML[nm]\t分割数\t積分範囲\tΔE\tピーク時ΔE\tN\tNRTD\tTEMP\tII*temp\tbar[%d]\tbar[%d]\n%8.3g\t%6d\t%8.2g\t%8.2g\t%8.2g\t%5d\t%4d\t%4.3g\t%8.2g",layer,LBAR-2,ML*1e9,DX,VI,DELTA,DELTAE,N[0],N[1],temp,II*temp);
		printf("なし\tなし\n");
	}
	return;
}