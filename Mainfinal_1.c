#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include <windows.h>
#include <omp.h>

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

char filename_potential[50];
char filename_wave[50];

int main(int argc, char *argv[])
{
	clock_t start,end;							/* �v���O�����̌v�Z���Ԃ̑���	*/
	start = clock();							/* �J�n�����擾					*/
	int p=2; printf("�v���O����No.%d\n",p);		/* �v���O�����̔ԍ�				*/
	setmaterial(-100);							/* �z�肵�Ă���\����ޗ��̎擾�@��������-100���Əo�͂���Ȃ� */
	int a,j,m,n, loop; //�V�����ϐ�loop(���[�v���񂷈ׂ̕ϐ�)��ǉ� (2016.10.13)
	double I,QW,VRTD,V, dv;
	double v[N[0]+1],vRTD[N[1]+1];				/* �|�e���V�����̊i�[�p			*/
	for(m=0;m<1;m++)
	{
		if(p==0)								/* �d���v�Z�p*/
		{
			if(OUTPUT==1)
			{
				printf("OUTPUT=1�ł��B�����I�����܂��B\n�|�e���V�����v�Z�Əo�͂̓v���O����No.1�ł���Ă��������B\n����p=0\n");
				exit(1);
			}
			else
			{
				setmaterial(0);					/* �ޗ����ēx�ݒ�@��ɕ����l��ύX���邽�߂ɓ���*/
				QW=0*ELEC*1e4*5e12;				/* ��˓��ɒ~�ς���Ă���d�ׂ̐ݒ�@set.dat�t�@�C�����ł��ݒ�ł��邪�A�ϐ��Ƃ������ꍇ�͂����ŉ��߂Đݒ�*/
				printf("\nNQW[/cm2],Q[C/m2]\n%e,%e\n\nVRTD,V,I[A/cm2]\n",QW*1e-4/ELEC,QW);
				dv=0.0025;
				for(j=-3;j<400;j++)				/* �d���v�Z�J�n */
				{	
					VRTD=dv*j;					/* RTD�Ɉ������Ă���d��������@���̌㑼�̑w�̓d�����v�Z���A�f�q�ɂ������Ă���d�����v�Z����B�@�S�̂ɂ������Ă���d������|�e���V�������v�Z���邱�Ƃ͍��� */
					potential(v,vRTD,VRTD,QW);	/* �|�e���V�����v�Z	���Ȗ������Ȍv�Z�𓱓����Ă���@�v�Z���������Ȃ��\������	*/
					V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);	/* �d���v�Z�@���T�C�h�̃t�F���~�G�l���M�[�̈ʒu����d�����v�Z			*/
					//printf("%e %e\n",VRTD,V);
					if(V<0)
						continue;
					else if(V<5)//+(4.7-4.05+Ef[0]))
					{
						I=current(v);			/* �d���v�Z	�d�����v�Z���邽�߂ɂ́A�|�e���V�������������Ƃ��ĕK�v �������A�萔���������Ă��Ȃ��̂Œ���*/
					}
					else
						break;
					//printf("%e %e %e\n",VRTD,V,III*I/1e4);	/* �o��*/
					printf("%e,%e,%e\n",VRTD,V,II*temp*I/1e4);	/* �o��*/
				}
			}
		}
		else if(p==1)							/* �|�e���V�����Ɠ��ߗ��v�Z�p*/
		{
			double E,T,dE;
			setmaterial(0);						/* �ޗ����ēx�ݒ�@��ɕ����l��ύX���邽�߂ɓ����@���ɕύX����K�v���Ȃ���Εs�v							*/
			QW=0*ELEC*1e4*1e13;					/* ��˓��ɒ~�ς���Ă���d�ׂ̐ݒ�@set.dat�t�@�C�����ł��ݒ�ł��邪�A�ϐ��Ƃ������ꍇ�͂����ŉ��߂Đݒ�	*/
			VRTD=2.855;
			potential(v,vRTD,VRTD,QW);
			//setpotential(vRTD,1,1);
			V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);
			printf("\nVRTD,V,Q[/cm2],Q[C/m2]\n%e,%e,%e,%e\n\n\nE,T\n",VRTD,V,QW*1e-4/ELEC,QW);
			dE=1e-6;							/* �G�l���M�[�̕���	���̒l��K�؂ɂƂ邱�Ɓ@��ʓI�ɂׂ͍����Ƃ�Ηǂ��@���̕��v�Z�Ɏ��Ԃ�������*/
			for(j=3200000;j<3500000;j++)
			{
				E=dE*j;					/* �G�l���M�[[eV]�̐ݒ�	*/
				T=cal(E,v,0,1);					/* ���ߗ��v�Z			*/
				printf("%e,%e\n",E,T);
			}
		}
		else if(p==2)							/* �g���֐��v�Z�p*/
		{
			if(argc < 2){
				printf("aaaaaaa");
				return 0;
			}
			int number = atoi(argv[1]);
			time_t t = time(NULL);
			struct tm *local = localtime(&t);
			char days[50];
			strftime(days, sizeof(days), "%Y_%m_%d",local);
			sprintf(filename_potential, "%s_%s_%d.csv", "potential", days, number);
			sprintf(filename_wave, "%s_%s_%d.csv", "wave", days, number);
			QW=0*ELEC*1e4*5e12;
			setmaterial(0);
			printf("NQW[/cm2],Q[C/m2]\n%e,%e\n",QW*1e-4/ELEC,QW);
			//VRTD=0.3375;
			VRTD=atof(argv[2]);
			
			potential(v,vRTD,VRTD,QW);
			V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);printf("\nVRTD,V\n%e %e\n\n",VRTD,V);
			n=atoi(argv[3]);								/* �v�Z���������ʂ̐�	*/
			double En[n];						/* ���ʂ̊i�[			*/
			double wavestore[n][DIV*N[0]];  //�g���֐��i�[
			//for(j=0;j<n;j++)
				//En[j]=confinedstate(j+1,v,0);
			confinedstates(v,0);
			getconfinedstates(n,v,0,En);
			for(a=0;a<4;a++)
			{
				if(a==1)
					printf("%e,",NX[LBAR]*dx*1e9);
				else if(a==2)
					printf("%e,",(NX[LBAR]+divnum[WELL]/2)*dx*1e9);
				else if(a==3)
					printf("%e,",(NX[LBAR]+divnum[WELL])*dx*1e9);
				for(j=0;j<n;j++)
				{
					if(a==0)
						printf("��%d����,",j+1);
					else
						printf("%e,",En[j]);
				}
				printf("\n");
			}
			for(j=0;j<n;j++)
				wavefunction(v,En[j],0, wavestore[j]); //�V��������wavestore[j]��ǉ� (2016.10.13)

				
			FILE *fp;
			fp = fopen(filename_wave, "w");
			for(j=0; j<DIV*N[0]; j++){
				fprintf(fp, "%g,", (double)j*dx*1e9/DIV);
				for(loop=0; loop<n; loop++)
					fprintf(fp, "%e,", wavestore[loop][j]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
		else if(p==3)
		{
			int a;
			int i=0;
			int in=15;
			double VRTQWS[in];
			double vout[in][N[0]+1];
			double Eout[in];
			double Vout[in];

			setmaterial(0);						/* �ޗ����ēx�ݒ�@��ɕ����l��ύX���邽�߂ɓ����@���ɕύX����K�v���Ȃ���Εs�v							*/
			QW=m*ELEC*1e4*5e12;					/* ��˓��ɒ~�ς���Ă���d�ׂ̐ݒ�@set.dat�t�@�C�����ł��ݒ�ł��邪�A�ϐ��Ƃ������ꍇ�͂����ŉ��߂Đݒ�	*/
			for(i=0;i<in;i++)
			{
				VRTQWS[i]=-0.1*i;
				potential(vout[i],vRTD,VRTQWS[i],QW);
				Vout[i]=vout[i][0]+Ef[0]-(vout[i][N[0]-1]+Ef[layer]);
				Eout[i]=getconfinedstate(1,vout[i],0);
			}
			printf("Q[/cm2],%e,Q[C/m2],%e\nz,",QW*1e-4/ELEC,QW);
			for(i=0;i<in;i++)
			{
				printf("%e,",VRTQWS[i]);
			}printf("\n");
			for(i=0;i<N[0];i++)
			{
				printf("%d,",i);
				for(j=0;j<in;j++)
				{
					printf("%e,",vout[j][i]);
				}
				printf("\n");
			}
			for(a=0;a<3;a++)
			{
				printf("\n%d,",NX[LBAR+a]);
				for(i=0;i<in;i++)
					printf("%e,",Eout[i]);
			}
			printf("\n\nV,");
			for(i=0;i<in;i++)
				printf("%e,",Vout[i]);
			printf("\n");
		}
		else if(p==4)//�e�X�g���[�h
		{
			int i,iE;
			double C;
			//double E;
			double Enum;
			double dE;
			double NQ;
			//gsl_complex A[N[np]],B[N[np]],k[N[np]];
			setmaterial(m);					/* �ޗ����ēx�ݒ�@��ɕ����l��ύX���邽�߂ɓ���*/
			QW=0*ELEC*5e4*1e12;				/* ��˓��ɒ~�ς���Ă���d�ׂ̐ݒ�@set.dat�t�@�C�����ł��ݒ�ł��邪�A�ϐ��Ƃ������ꍇ�͂����ŉ��߂Đݒ�*/
			C=massxy[0]*MSTAR*KB*temp/HH/HH/1e4*sqrt(2*mass[0]*MSTAR/ELEC)/HBAR*ELEC;
			iE=1e5;
			dE=VI/iE;
			printf("\nNQW[/cm2],Q[C/m2]\n%e,%e\n\nVRTD,V,NQ[A/cm2],iE,dE\n",QW*1e-4/ELEC,QW);
			for(j=-4;j<20;j++){
				NQ=0;
				VRTD=0.05*j;					/* RTD�Ɉ������Ă���d��������@���̌㑼�̑w�̓d�����v�Z���A�f�q�ɂ������Ă���d�����v�Z����B�@�S�̂ɂ������Ă���d������|�e���V�������v�Z���邱�Ƃ͍��� */
				potential(v,vRTD,VRTD,QW);	/* �|�e���V�����v�Z	���Ȗ������Ȍv�Z�𓱓����Ă���@�v�Z���������Ȃ��\������	*/
				V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);	/* �d���v�Z�@���T�C�h�̃t�F���~�G�l���M�[�̈ʒu����d�����v�Z			*/
				//E=calconfinedstate(v,0,v[0],EMAX);
				for(a=0;a<1;a++)
				{
					NQ=0;
					// #pragma omp parallel for private(Enum) reduction(+:NQ) OpenMP���g���Ȃ��ׁA�R�����g�A�E�g�ɂ����B(2016/09/25)
					for(i=1;i<iE+1;i++)
					{
						Enum=v[0]+dE*i;
						NQ+=cal4(Enum,v,0)*dE;
						//printf("%d %e %e\n",i,Enum,cal4(Enum,v,0));
					}
					printf("%e,%e,%e,%d,%e\n",VRTD,V,valley[0]*C*NQ,iE,dE);	/* �o��*/
					//iE*=3;
					//dE=VI/iE;
				}
			}
		
		}
		else
		{
			
		}
		printf("\n");
	}
	end = clock(); printf("\n���Z���� %f �b\n",(double)(end-start)/CLOCKS_PER_SEC);
	//Beep( 440, 1200 );
	return 0;
}

void setmaterial(int m)//m=-100���ƁA�o�͂���Ȃ��B���s�͂����B������m�́A�����l��ύX�������ꍇ�Ɏg�p�B
{
	struct material{
		double die;		/* �U�d��					*/
		double bar;		/* �o���h�s�A��				*/
		double mass;	/* �L������					*/
		double massxy;	/* �������̗L������			*/
		double ef;		/* �t�F���~�G�l���M�[		*/
		int valley;		/* �`���т̒J�̐�			*/
		char *name;		/* �����̖��O				*/
		int cond;
	};

	int i,j;
	double base;		/* �`���уG�l���M�[�̊�B����Si�̓d�q�e�a��	*/
	struct material mtconst[100];/* �����̐ݒ�B�z����\���m�ۂ��Ȃ��ƃG���[���N����̂Œ��ӁB*/
	temp=TEMP;			/* ���x						*/
	int flag=-100;		/* �����������̒l�Ɠ����Ȃ�A�����͂���邪�o�͂͂���Ȃ��Ȃ�B��{�I�ɂ͏o�͂������������Ǝv���B*/
	/* �e�����ɂ����镨���l�̐ݒ�Bname:���O�Adie�F�U�d���Abar�F�o���h�s�A���Amass:�L�����ʁAvalley:�J�̐��Bvalley:��{1�Bcond:���d���B0�������A1��n-Si�A2��p-Si�A3���≏���B*/
	i=0; mtconst[i].name="Al       "; mtconst[i].die=DIE_AL   * DIEELECSTAR; mtconst[i].bar=BAR_AL ;  mtconst[i].mass=MASS_AL;   mtconst[i].massxy=MASS_AL;		mtconst[i].valley=1;			mtconst[i].cond=0;
	i=1; mtconst[i].name="n-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=1;
	i=2; mtconst[i].name="n-Si     "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=1;
	i=3; mtconst[i].name="p-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=2;
	i=4; mtconst[i].name="p-Si     "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=2;
	i=5; mtconst[i].name="CaF2     "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=BAR_CAF2; mtconst[i].mass=MASS_CAF2; mtconst[i].massxy=MASS_CAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=6; mtconst[i].name="CdF2     "; mtconst[i].die=DIE_CDF2 * DIEELECSTAR; mtconst[i].bar=BAR_CDF2; mtconst[i].mass=MASS_CDF2; mtconst[i].massxy=MASS_CDF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=7; mtconst[i].name="i-Si     "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=3;
	i=8; mtconst[i].name="Au       "; mtconst[i].die=DIE_AU   * DIEELECSTAR; mtconst[i].bar=BAR_AU ;  mtconst[i].mass=MASS_AU;   mtconst[i].massxy=MASS_AU;		mtconst[i].valley=1;			mtconst[i].cond=0;
	i=9; mtconst[i].name="nSub.CaF2"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=BAR_CAF2; mtconst[i].mass=MASS_CAF2; mtconst[i].massxy=MASS_CAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=10;mtconst[i].name="p-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_pSI;  mtconst[i].mass=MASS_pSI;  mtconst[i].massxy=MASS_pSI;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=11;mtconst[i].name="pCaF2    "; mtconst[i].die=DIE_pCAF2* DIEELECSTAR; mtconst[i].bar=BAR_pCAF2;mtconst[i].mass=MASS_pCAF2;mtconst[i].massxy=MASS_pCAF2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=12;mtconst[i].name="pAl      "; mtconst[i].die=DIE_pAL  * DIEELECSTAR; mtconst[i].bar=BAR_pAL;  mtconst[i].mass=MASS_pAL;  mtconst[i].massxy=MASS_pAL;	mtconst[i].valley=1;			mtconst[i].cond=0;
	i=13;mtconst[i].name="5MLCaF2  "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.7;      mtconst[i].mass=0.7;       mtconst[i].massxy=0.7;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=14;mtconst[i].name="���zAL   "; mtconst[i].die=DIE_AL   * DIEELECSTAR; mtconst[i].bar=-5.17  ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=15;mtconst[i].name="���zCaF2 "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=7      ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=16;mtconst[i].name="���zp-Si "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;      mtconst[i].massxy=0.55;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=17;mtconst[i].name="SiO2     "; mtconst[i].die=DIE_SIO2 * DIEELECSTAR; mtconst[i].bar=BAR_SIO2; mtconst[i].mass=MASS_SIO2; mtconst[i].massxy=MASS_SIO2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=18;mtconst[i].name="nc-Si.   "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=19;mtconst[i].name="nSi100Z  "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_Z; mtconst[i].massxy=MASS_SI_XY;	mtconst[i].valley=NUMVALLY_SI_Z;	mtconst[i].cond=1;
	i=20;mtconst[i].name="nSi100XY "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_XY;mtconst[i].massxy=MASS_SI_Z;	mtconst[i].valley=NUMVALLY_SI_XY;	mtconst[i].cond=1;
	i=21;mtconst[i].name="n-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;
	i=22;mtconst[i].name="CaF2(4ML)"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.6;	  mtconst[i].mass=0.85;	     mtconst[i].massxy=0.85;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=23;mtconst[i].name="p-i-Si   "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;	     mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;	
	i=24;mtconst[i].name="3MLCaF2  "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.5;      mtconst[i].mass=1.0;       mtconst[i].massxy=1.0;		mtconst[i].valley=1;			mtconst[i].cond=3;
	/* ������ǉ��������ꍇ�́A���̍s�̏���R�s�[���ē\��t���Bi��ύX����̂�mtconst�̔z��Ɏ��܂�悤�ɒ��ӁB*/
	if(DX>0)
		dx=ML/DX;	/* DX�͕�����	ML�͂��̃v���O������z���̊�{�P�ʁB����ML�������DX��������B*/
	else			/* ML=1e-9�ɂ���΁Az[nm]�BML=0.31e-9�Ȃ�z[ML]	*/
	{
		printf("������DX��0�ȉ��ł��B�K�؂Ȓl��ݒ肵�ĉ������B\n");
		exit(1);
	}
	N[0]=0;			/* RTD�S�̂̑w��*/
	FILE *fd1;
	if((fd1 = fopen("set.dat", "r")))
	{
		if(m!=flag)
			printf("�w�ԍ�,�ޗ��ԍ�,������,ML��,�w��[nm],��U�d��,�L������,��ǂ̍���,�J,������,NX,EF,�d���֐�,�d�ח�\n");
		for(i=0; fscanf(fd1, "%d %d %lf\n", &smt[i], &j, &Q[i])!=EOF; i++)
		{
			if( modf((double)j*DX,&d[i])!=0)
			{
				printf("�����G���[�@�����ԊuDX��ύX���ĉ�����\n");
				exit(1);
			}
			divnum[i] = j*DX;					/* ������					*/
			N[0]        += divnum[i];			/* ��������					*/
			NX[i]     = N[0];
			d[i]      = j*ML;					/* ����[m]					*/
			name[i]   = mtconst[smt[i]].name;	/* �ޗ���					*/
			die[i]    = mtconst[smt[i]].die;	/* �e�w�̗U�d��[F/m]		*/
			bar[i]    = mtconst[smt[i]].bar;	/* �e�w�̃o���h�s�A����Ec	*/
			mass[i]   = mtconst[smt[i]].mass;	/* �e�w�̗L������			*/
			massxy[i] = mtconst[smt[i]].massxy;	/* �e�w�̉������̗L������	*/
			valley[i] = mtconst[smt[i]].valley;	/* �`���т̒J�̐�			*/
			cond[i]	  = mtconst[smt[i]].cond;
			layer=i;							/* �w�ԍ�					*/
			switch(smt[i]){/*�t�F���~�G�l���M�[�̐ݒ�B�����͕����l�A�����̂͋��ȏ��̎�����v�Z�A�≏���ɂ͑Ή����Ă��Ȃ��B*/
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
				case 14: Ef[i]=4.37;break;
				case 16: Ef[i]=-EG_SI/2+KB*temp*log(Q[i]*1e6/NI)/ELEC;base=BASE;break;
				case 19: Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 20: Ef[i]=fermia(i,temp,smt[i]);base=BASE;break;
				case 23: Ef[i]=fermia(i,temp,smt[i]);base=pBASE;break;
				default: Ef[i]=0;base=pBASE;break;
			}
			switch(mtconst[smt[i]].cond){ //2016.09.30�ǉ��D�≏�͎̂d���֐����v�Z�����Ȃ��悤�ɂ����D
				case 0 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 1 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 2 : Work[i]   = -Ef[i]-bar[i]+base; break;
				case 3 : break;
			}
//			Work[i]   = -Ef[i]-bar[i]+base;		/* �d���֐��B���R�A�≏���ɂ͑Ή����Ă��Ȃ��B*/
			if(m!=flag)
				printf("%d,%d,%s,%d,%.4g,%.3g,%.2g,%.4g,%d,%d,%d,%.4g,%.4g,%.4g\n", i, smt[i], name[i], j, j*ML*1e9, die[i]/DIEELECSTAR, mass[i], bar[i], valley[i], divnum[i], NX[i],Ef[i],Work[i],Q[i]);//�o�́B
		}
		fclose(fd1);
	}
	else 
	{
		printf("�t�@�C�����J���Ă��܂���\n");
		exit(1);
	}
	N[1]=divnum[LBAR]+divnum[WELL]+divnum[RBAR];	/* RTD�\���݂̂̑w��*/
	if(m!=flag)
		printf("\n 1ML[nm],������,�ϕ��͈�,��E,�s�[�N����E,N,NRTD TEMP,II*temp,bar[%d],bar[%d]\n%.3g,%d,%.2g,%.2g,%.2g,%d,%d,%.3g,%.2g ",layer,LBAR-2,ML*1e9,DX,VI,DELTA,DELTAE,N[0],N[1],temp,II*temp);
	if(AAA>0){									/* �t�F���~�G�l���M�[�s���j���O�̌��ۂ��Č��B�v�����B��{�g��Ȃ������ǂ�?*/
		if(d[RBAR+1]!=0)
		{
			bar[layer]-=Work[layer-1]-Work[layer];
			if(m!=flag)
				printf("%.3g,",bar[layer]);
		}
		else
		{
			if(m!=flag)
				printf("�Ȃ�,");
		}
		if(d[LBAR-1]!=0)//����K�v
		{
			bar[LBAR-2]-=Work[LBAR-1]-Work[LBAR-2];
			if(m!=flag)
				printf("%.3g\n",bar[LBAR-2]);
		}
		else
		{
			if(m!=flag)
				printf("�Ȃ�\n");
		}
	}
	else
		if(m!=flag)
			printf("�Ȃ�,�Ȃ�\n");

	return;
}

void potential(double v[], double vRTD[], double VRTD, double QW)
{
	int n,output;
	double DR,DL;
	selfpotential(QW,VRTD,vRTD);
	for(n=0;n<N[1]+1;n++)
	{
		v[n+NX[1]]=vRTD[n];
		//printf("%d %e\n",n,vRTD[n]);
	}
	DL=die[LBAR]*(vRTD[1]-vRTD[0])/dx;
	DR=die[RBAR]*(vRTD[N[1]]-vRTD[N[1]-1])/dx;
	for(n=LBAR-1;n>=0;n--)
		DL=calVRL(-1,n,DL,v);
	for(n=RBAR+1;n<=layer;n++)
		DR=calVRL(1,n,DR,v);
	output=OUTPUT;
	setpotential(v,0,output);
	//outputpotential(v,0);
	//setpotential(vRTD,1,output);
	return;
}

double calcurrent(double E1, double E2, double v[], double delta)
{
	int i,n;										/* E1����E2�܂ł̓d���l���v�Z�B���ݕ��͈���delta�B											*/
	double dE,E,T,S1,S2,I;							/* delta�����ߗ��̔��l�������傫���ꍇ�A���l�v�Z�Őϕ������s����ƌ덷���傫���Ȃ�̂ŁA	*/
	I=0;											/* delta�͏\���������Ƃ邱�ƁB�������A����������ƌv�Z�Ɏ��Ԃ�������B						*/
	n=(E2-E1)/delta+1;								/* ���ݕ��̌���																				*/
	dE=(E2-E1)/(n-1);								/* �덷����B���܂�傫���Ȃ��̂ŕ��u�B														*/
	//printf("E1 %e E2 %e n %d dE %e\n",E1,E2,n,dE);
// #pragma omp parallel for private(E,T,S1,S2) reduction(+:I) OpenMP���g���Ȃ��ׁA�R�����g�A�E�g�ɂ���(2016/09/25)	/* openmp���g���Ȃ��ꍇ�́A���̍s�͖��������̂Ŗ��Ȃ��B		*/
	for(i=0;i<n;i++)
	{
		E=E1+i*dE;
		T=cal(E,v,0,1);													/* ���ߗ��̌v�Z                             */
		S1=massxy[0]*valley[0]*func(E,v[0],0);								/* Supply function�̌v�Z                    */
		S2=massxy[mt(N[0]-1)]*valley[mt(N[0]-1)]*func(E,v[N[0]-1],layer);	/* Supply function�̌v�Z                    */
		if(i==0 || i==n-1)
			I+=0.5*T*(S1-S2);
		else
			I+=T*(S1-S2);
	}
	I*=dE;
	return I;
}

double current(double v[])						/* �d���l�̌v�Z�B�����ł́A�ϕ��͈͂����肵�Ă���B*/
{
	int i;
	int n=1;
	int num=10;									/* ���݂��鏀�ʂ̐�					*/
	double Emin,Emax,En[num],E1,E2,max1,max2;
	double I=0;

	if(DELTA<DELTAE)
	{
		printf("DELTA < DELTAE�̂��ߌv�Z�ł��܂���ł����B\n�����I�����܂��B\nDELTAE��K�؂Ȓl�ɂ��ĉ�����\n");
		exit(0);
	}
	max1=max(v[0],v[0]+Ef[0]);					//printf("v[0] %e v[0]+Ef[0] %e\n",v[0],v[0]+Ef[0]);
	max2=max(v[N[0]-1],	v[N[0]-1]+Ef[layer]);	//printf("v[N] %e v[N]+Ef[N] %e\n",v[N[0]-1],v[N[0]-1]+Ef[layer]);
	Emax=max(max1,max2)+VI;						//printf("Emax %e max1 %e max2 %e\n",Emax,max1,max2);
	if(v[0]>v[N[0]-1])							/* �ϕ��͈͂̍ŏ��l�̐ݒ�			*/
		Emin=v[0];
	else if(v[0]<v[N[0]-1])
		Emin=v[N[0]-1];
	else
		Emin=v[0];
	getconfinedstates(num,v,0,En);				/* ���ʂ̒T���B�z��Ɋi�[�B			*/
	E1=Emin;
	for(i=0;i<num;i++)
	{
		if(En[i]<E1)							/*�ϕ���ԓ���En[i]�̏��ʂȂ�									*/
			continue;
		else if(Emax<En[i])
		{
			I+=calcurrent(E1,Emax,v,DELTA);		/*�ϕ���Ԃ̍Ō�܂Őϕ�										*/
			break;
		}
		else if(E1<En[i] && En[i]<Emax)			/*�ϕ���ԓ��ɏ��ʂ���@�ׂ����G�l���M�[�𕪊����Čv�Z����K�v����*/
		{
			n=1;
			while(E2<Emin+VI)
			{
				E2=E1+(double)(n*DELTA);
				if(En[i]<=E2)
					break;
				n++;
			}
			if(n<3)
			{
				I+=calcurrent(E1,E2,v,DELTAE);
				E1=E2;
			}
			else
			{
				E2=Emin+(double)((n-2)*DELTA);
				I+=calcurrent(E1,E2,v,DELTA);
				E1=E2;
				E2+=3*DELTA;
				I+=calcurrent(E1,E2,v,DELTAE);
				E1=E2;
			}
		}
	}
	return I;
}

int mt(int n)
{
	int x;
	for(x=0;x<layer;x++)
	{
		if(n<=NX[x])
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
		printf("�z�肵�Ă��Ȃ��l��������np�Ɏg���Ă��܂��B\n");
		exit(1);
	}
}

void potential0 (double Q, double VRTD, double v[], int np)
{
	int n,nrtd;
	nrtd=nRTD(np);
	double x,DL,DR;
	DL = (VRTD-(d[WELL]/2/die[WELL]+d[RBAR]/die[RBAR])*Q) / (d[LBAR]/die[LBAR]+d[WELL]/die[WELL]+d[RBAR]/die[RBAR]);
	DR = DL+Q;
	for(n=0; n<N[np]+1; n++)
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
	int n,nrtd;								/* �d�ʕ��zv�ɓ`���o���h�s�A���𓱓����A�K�i�ߎ���K�p����B*/
	nrtd=nRTD(np);
	double vr[N[np]+1],vl[N[np]+1];			/* �|�e���V�����̍�������̋Ɍ�vr�ƉE������̋Ɍ�			*/
	FILE *fp;
	if(output==1){
		fp = fopen(filename_potential, "w");
	}
	for(n=0; n<N[np]+1; n++){
		vl[n]=v[n]+bar[mt(n+nrtd)];
		vr[n]=v[n]+bar[mt(n+nrtd+1)];
		if(output==1)
		{
			fprintf(fp, "%e,%e\n",n*dx*1e9,vl[n]);
		}
	}
	//if(output==1)
		//printf("\nn z[nm] Potential[eV]\n");
	for(n=0; n<N[np]; n++){
		//printf("%d %d %d %e\n",n,n+nrtd,mt(n+nrtd),bar[mt(n+nrtd)]);
		v[n]=(vr[n]+vl[n+1])/2;				/* �K�i�ߎ��K�p	*/
		//if(np==0 && output==1)
			//printf("%e %e %e\n",n+0.5,(n+0.5)*dx*1e9,v[n]);
	}
	if(output==1)
	{
		fclose(fp);
	}
	return;
}

//void outputpotential (double v[],int np)
//{
//	int n,nrtd;
//	nrtd=nRTD(np);
//	double vr[N[np]+1],vl[N[np]+1];
//	for(n=0; n<N[np]+1; n++){
//		vl[n]=v[n]+bar[mt(n+nrtd)];
//		vr[n]=v[n]+bar[mt(n+nrtd+1)];
//		printf("%d %e\n%d %e\n",n,vl[n],n,vr[n]);
//	}
//	printf("\n");
//	for(n=0; n<N[np]; n++){
//		v[n]=(vr[n]+vl[n+1])/2;
//		//printf("%e %e\n",n+0.5,v[n]);
//	}
//	return;
//}

void setk(double E, double v[],gsl_complex k[],int np)
{
	int n,nrtd;					/* �|�e���V��������g��k���v�Z���Ċi�[*/
	nrtd=nRTD(np);
	for(n=0; n<N[np]; n++)
	{
		k[n]=gsl_complex_div_real(gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n])*MSTAR*ELEC),HBAR);
		//printf("name[%d] %s,k[%d]=%g +i %g\n",n,name[mt(n+1)],n,GSL_REAL(k[n]),GSL_IMAG(k[n]));
	}
}

double cal(double E, double v[], int np, int f) // �o�͂́A���ߗ�|T|��2��=1/|T(1,1)|��2�� np��0���S�́A1��RTD���Bf�́A0�̂Ƃ��A�o���h�M���b�v�����v�Z����B
{
	int n,nrtd;
	double m,t;
	
	nrtd=nRTD(np);

	if(f==1 && (E-v[N[np]-1]<=0 || E-v[0]<=0))	//�o���h�M���b�v���Ɣ��f�B�v�Z�����B
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

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]���v�Z
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]���v�Z

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //�G�l���M�[E�A�̈�n+1�ɂ�����g�����v�Z HBAR�Ŋ����Ă��Ȃ����Ƃɒ���

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)���v�Z
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)���v�Z
			
			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0));                         //�S�Ă̗v�f��0.5��������
		}

		else if(E-v[n+1]==0 && E-v[n]!=0)
		{
			//printf("kn!=0,kn1=0\n");
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
			//printf("kn=0,kn1!=0\n");
			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC);
			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];
			kn1=gsl_complex_sqrt_real(HBAR*HBAR/(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC));
			pp=gsl_complex_mul_imag(kn1, m);
			pm=gsl_complex_mul_imag(kn1,-m);

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //�S�Ă̗v�f��0.5��������
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			//printf("kn=0,kn1=0\n");
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,dummy,gsl_complex_rect(0,0),trans); //�s��̊|���Z

		gsl_matrix_complex_memcpy(dummy, trans);

		/*printf("Transfer Matrix n=%d material %s\nT(0,0)=%e +i %e T(0,1)=%e +i%e\nT(1,0)=%e +i %e T(1,1)=%e +i%e \n\n",n,name[mt(n+nrtd+1)],
				GSL_REAL(gsl_matrix_complex_get(trans,0,0)), GSL_IMAG(gsl_matrix_complex_get(trans,0,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,0,1)), GSL_IMAG(gsl_matrix_complex_get(trans,0,1)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,0)), GSL_IMAG(gsl_matrix_complex_get(trans,1,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,1)), GSL_IMAG(gsl_matrix_complex_get(trans,1,1)));*/
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

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]���v�Z
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]���v�Z

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //�G�l���M�[E�A�̈�n+1�ɂ�����g�����v�Z HBAR�Ŋ����Ă��Ȃ����Ƃɒ���

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)���v�Z
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)���v�Z
			
			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0));                         //�S�Ă̗v�f��0.5��������
		}

		else if(E-v[n+1]==0 && E-v[n]!=0)
		{
			//printf("kn!=0,kn1=0\n");
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
			//printf("kn=0,kn1!=0\n");
			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC);
			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));

			m=mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)];
			kn1=gsl_complex_sqrt_real(HBAR*HBAR/(2*mass[mt(n+nrtd+1)]*(E-v[n+1])*MSTAR*ELEC));
			pp=gsl_complex_mul_imag(kn1, m);
			pm=gsl_complex_mul_imag(kn1,-m);

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //�S�Ă̗v�f��0.5��������
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			//printf("kn=0,kn1=0\n");
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[
		}

		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,FN,gsl_complex_rect(0,0),FNN); //�s��̊|���Z FNN=temp*FN

		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);

		gsl_matrix_complex_memcpy(FN, FNN);
		/*printf("\nTransfer Matrix\nn   %d material name %s mass[mt(%d)])=%e\nn+1 %d material name %s mass[mt(%d)])=%e\nT(00) %e i %e T(01) %e i %e\nT(10) %e i %e T(11) %e i %e \n\n",n,name[mt(n+nrtd+1)],n,mass[mt(n+nrtd+1)],n+1,name[mt(n+nrtd+2)],n+1,mass[mt(n+nrtd+2)],
				GSL_REAL(gsl_matrix_complex_get(temp,0,0)), GSL_IMAG(gsl_matrix_complex_get(temp,0,0)),
				GSL_REAL(gsl_matrix_complex_get(temp,0,1)), GSL_IMAG(gsl_matrix_complex_get(temp,0,1)),
				GSL_REAL(gsl_matrix_complex_get(temp,1,0)), GSL_IMAG(gsl_matrix_complex_get(temp,1,0)),
				GSL_REAL(gsl_matrix_complex_get(temp,1,1)), GSL_IMAG(gsl_matrix_complex_get(temp,1,1)));*/

		//printf("%d %e\n",n,gsl_complex_abs2( gsl_complex_add(A[n],B[n]) ));
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

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]���v�Z
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]���v�Z

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //�G�l���M�[E�A�̈�n+1�ɂ�����g�����v�Z HBAR�Ŋ����Ă��Ȃ����Ƃɒ���

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)���v�Z
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)���v�Z
			
			gsl_matrix_complex_set(temp[n],0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp[n],0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp[n],1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp[n],1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp[n], gsl_complex_rect(0.5,0));                         //�S�Ă̗v�f��0.5��������
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],dummy,gsl_complex_rect(0,0),trans); //�s��̊|���Z

		gsl_matrix_complex_memcpy(dummy, trans);
		/*printf("Transfer Matrix n=%d material %s\nT(0,0)=%e +i %e T(0,1)=%e +i %e\nT(1,0)=%e +i %e T(1,1)=%e +i %e \n\n",n,name[mt(n+1)],
				GSL_REAL(gsl_matrix_complex_get(temp[n],0,0)), GSL_IMAG(gsl_matrix_complex_get(temp[n],0,0)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],0,1)), GSL_IMAG(gsl_matrix_complex_get(temp[n],0,1)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],1,0)), GSL_IMAG(gsl_matrix_complex_get(temp[n],1,0)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],1,1)), GSL_IMAG(gsl_matrix_complex_get(temp[n],1,1)));

		printf("Tran n=%d material %s\nT(0,0)=%e +i %e T(0,1)=%e +i %e\nT(1,0)= %e +i %e T(1,1)=%e +i %e \n\n",n,name[mt(n+1)],
				GSL_REAL(gsl_matrix_complex_get(trans,0,0)), GSL_IMAG(gsl_matrix_complex_get(trans,0,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,0,1)), GSL_IMAG(gsl_matrix_complex_get(trans,0,1)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,0)), GSL_IMAG(gsl_matrix_complex_get(trans,1,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,1)), GSL_IMAG(gsl_matrix_complex_get(trans,1,1)));*/
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
			printf("\nA[0],1+0i,B[0],1-r\n");break;//��������E���ɓ��߂���ꍇ
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
			printf("\nA[0],0+0i,B[0],1+0i\n");break;//�E�����獶���ɓ��߂���ꍇ
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
			printf("\nA[0],1+0i,B[0],1+0i\n");break;//�v�����Bv[0]<v[N+1]�̂Ƃ��K�p�\�H2015/3/25
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
			printf("\nA[0],1+0i,B[0],0+0i\n");break;//�ȈՔŁB���ߗ���1�̏ꍇ�K�p�ł���H
	}
	gsl_matrix_complex *FN  = gsl_matrix_complex_calloc(2,1);
	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);
	
	for(n=0; n<N[np]-1; n++)
	{
		gsl_matrix_complex_set(FN,0,0,A[n]);
		gsl_matrix_complex_set(FN,1,0,B[n]);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],FN,gsl_complex_rect(0,0),FNN); //�s��̊|���Z FNN=temp*FN
		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);
	}
	//printf("T %e R %e\n", sqrt(1-v[0]/E)/gsl_complex_abs2(gsl_matrix_complex_get(trans,1,1)),gsl_complex_abs2(B[0]) );
	for(n=0; n<N[np]; n++){
		gsl_matrix_complex_free(temp[n]);
	}
	gsl_matrix_complex_free(dummy);
	gsl_matrix_complex_free(trans);
}

void confinedstates(double v[], int np)
{
	double dE = 1e-4;
	double minV,maxV;
	minV = v[0];
	maxV = v[0];
	for(int n = 1; n < N[np]; n++){
		if( minV > v[n] )	minV = v[n];
		if( maxV < v[n] )	maxV = v[n];
	}
	//int count_E = (maxV - minV)/dE;
	FILE *fp = fopen("confinedstates.csv", "w");
	for(double e = minV; e <= maxV; e+=dE){
		fprintf(fp, "%e,%e\n", e, cal(e, v, np, 0));
	}
	fclose(fp);
}

double confinedstate(int num, double v[], int np)		/* ����num�́A�ʎq��n�̂��ƁBn�͕ʂ̂Ƃ���Ŏg���Ă���̂�num�ɂ����B */
{														/*�o�O����B���炭�A���E�̈�Ԓ[�̃G�l���M�[�����ʂƂ��Č��o����Ă��܂��B*/
	int i,j,n;
	double E,T,S,dE;
	double emax=5;                                    /* while�̖������[�v������邽�� */
	double min=5;

	for(n=0;n<N[np]+1;n++)		
	{
		if(min>v[n])
			min=v[n];                                   /* �������ꂽ�|�e���V�����̍ŏ��l�����߂� */
	}
	E=min;
	//printf("min=%e\n",E);
	for(i=1;i<num+1;i++)
	{
		dE=1e-3;                                      /* ���ݕ�(�e��)�̐ݒ�                                              */
		T=j=0;
		S=cal(E,v,np,0);
		while(E<emax)                                 /* ���ߗ��̌�����Ԃ̌���                                          */
		{                                             /* �ȉ��̏��ʂ����߂�v���O������                                  */
			E+=dE;                                    /* �G�l���M�[���㏸���Ă����ƂƂ��ɁA���ߗ����㏸���Ă����Ɖ��肵�A*/
			T=cal(E,v,np,0);                          /* ���̍ő�l���Ƃ�ӏ������ʂł���Ƃ��Ă���                      */
			if(S>T)                                   /* ����āA���ߗ����������Ă�����Ԃ͕s�v                          */
				S=T;
			else
				break;                                /* ���ߗ��̌�����Ԃ̏I�� */
		}
		S=T=0;                                        /* S�AT�̏�����                          */
		while(E<emax)                                 /* ��ꏀ�ʂ̒T�� ���ݕ��F�e��           */
		{
			E+=dE;                                    /* �G�l���M�[�̐ݒ�                      */
			T=cal(E,v,np,0);                              /* ���ߗ��̌v�Z                          */
			if(S<T)                                   /* �����̔��f�@�����֐��Ȃ��while�����s */
				S=T;
			else
			{
				E-=2*dE;                              /* �G�l���M�[�������߂��āA���ݕ��ׂ̍������ֈڍs */
				dE/=10;
				//j++;
				S=0;
				if(dE<1e-10)
					break;
			}
		}
		//printf("j=%d\n",j);
		dE/=10;                                      /* ���ݕ�(�ׂ���)�̐ݒ�                           */
		S=0;                                          /* S�̏�����                                      */
		while(E<emax)                                 /* ��ꏀ�ʂ̒T�� ���ݕ��F�ׂ���                  */
		{
			E+=dE;                                    /* �G�l���M�[�̐ݒ�                               */
			T=cal(E,v,np,0);                              /* ���ߗ��̌v�Z                                   */
			if(S<T)
				S=T;                                  /* �����̔��f�@�����֐��Ȃ��while�����s          */
			else{
				//printf("����%d %+20.18lf T %e\n",i,E,S);
				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
				{										/*cal�֐��̉e���ŁAE=v[0]�̂Ƃ��낪���ʂɌ����Ă��܂��̂�*/
					E+=10*dE;							/*���������邽�߂ɓ���*/
					i--;								/*�]���āA���ʂ�v[0]�}1e-10[eV]�̋ߕӂɑ��݂����ꍇ�A*/
					//printf("aS %e T %e E %e E-dE-v[0] %e E-dE-v[N[np]-1]%e \n",S,T,E,E-dE-v[0],E-dE-v[N[np]-1]);
					break;								/*���m�ł��Ȃ��\������*/
				}
				if(i==num)
					return E-dE;
				E+=10*dE;
				break;                                /* ���ߗ��̑�����Ԃ̏I�� ���̈�O�̌v�Z���ʂ�  */
			}
		} /* �ȗ��������� */
	}
	E-=dE;
	//printf("����%d %e T %e\n",i,E,S);                 /* ���ʂ̏o��                                    */
	return E;                                         /* num�Ԗڂ̏��ʂ�Ԃ�                           */
}

void makewave(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np, double wavestore[]) //�V��������wavestore��ǉ� (2016.10.13)
{
	int i,j,n;
	double max=0;                     // �K�i���p
	double sum=0;
	double xn,px,py, pyy[N[np]*DIV]; //�V�����z��ǉ�(pyy, 2016.10.13)
	//double pr,pi;
	gsl_complex tempA,tempB;
	//printf("\nn z[nm] Ar Ai Br Bi Yr Yi Y*Y\n");
	j=max=0;//p[0]=-ML;q[0]=E;
	for(n=0;n<N[np];n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			//pr=GSL_REAL(gsl_complex_add(tempA,tempB));
			//pi=GSL_IMAG(gsl_complex_add(tempA,tempB));
			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			//if(py>max)
				//max=py;                                                                       /* �g���֐���2��̍ő�l�����߂�      */
			//printf("%g %e %e %e %e %e %e %e %e\n",(double)j/DIV,px*1e9,GSL_REAL(tempA),GSL_IMAG(tempA),GSL_REAL(tempB),GSL_IMAG(tempB),pr,pi,py);
			max+=py*dx*1e9/DIV;
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=py*dx*1e9/DIV;
			j++;
		}
	}j=0;//2015/5/20�ǉ��@�v�l��~ �Ƃ肠�����g���֐��̖ʐς����߂�(����max�Ɋi�[)�B�g���֐���2��̑S��ԂŐϕ������1�ɂȂ�悤�ɋK�i���B
	//printf("\nz[nm] Y*Y Y*Y\n"); //�ʎq��˂̏ꍇ�A���������������exp(�}ikz)�̊|���Z�ŃI�[�o�[�t���[(�H)���N�����A�������ł͂Ȃ����������x�z�I�ɂȂ�̂Œ��ӁB
	for(n=0;n<N[np];n++)
	{

	//���������C�o�͂̎d����啝�ɕύX(��QCL.c�Ƃقړ����ɂ���) by�V�� (2016.10.13)

		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			pyy[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			if(pyy[j]>max)
				max=pyy[j];
			//printf("%e %e %e\n",px*1e9, py/max+E, py/max);
			j++;
		}
	}
	for(i=0; i<N[np]*DIV; i++)
		wavestore[i] = E+pyy[i]/max;
	//printf("\nmax %e sum %e\n",max,sum);
	return;
}

int makewave2(double Q, double E, double VRTD, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[],int np)
{
	int i,j,n;
	double max=0;                     // �K�i���p
	double qmax=0;
	double xn,p[DIV*N[np]],q[DIV*N[np]]; // p��x���W�Aq��y���W q�͔g���֐��̐�Βl��2����o�͗\��
	double qtemp[N[np]+1], vnew[N[np]];
	gsl_complex tempA,tempB;
	j=0;//	p[0]=-ML;
	//printf("dx=%e\n",dx);
	for(n=0; n<N[np]; n++)
	{
		xn=(n+1)*dx;
		//if(E-v[n]!=0)
		//{
			for(i=0; i<DIV; i++)
			{
				p[j]=(double)j*dx/DIV;
				tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],p[j]-xn)));     /* tempA= A[n] �~ exp( ik[n](x-x[n])) */
				tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-p[j])));     /* tempB= B[n] �~ exp(-ik[n](x-x[n])) */
				q[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));                                /* �g���֐��̐�Βl��2����v�Z        */
				qmax+=q[j]/DIV;
				//printf("%d %e %e\n",j,p[j],q[j]);
				j++;
			}
		//}
	}
	//printf("\nmax %e\nX q\n",max);
	j=0;
	for(n=1;n<N[np];n++)
	{
		//for(i=0;i<DIV;i++)
		//{
		//	q[j]*=Q/qmax;
		//	//printf("%e %e\n",p[j],q[j]);
		//	j++;
		//}
		qtemp[n]=-Q*q[DIV*n]/qmax;
		//printf("%d %e\n",n,qtemp[n]);
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
	//printf(" max %e \n",max);
	for(n=0; n<N[np]; n++)
			v[n]=vnew[n];
	if(max<1e-6){			/* �d�v�@���Ȗ������v�Z�̎�������*/
		/*printf("\ni X q\n");
		for(j=0;j<DIV*N[np];j++)
			printf("%d %e %e\n",j,p[j],Q*q[j]/qmax);*/
		for(n=0;n<N[np];n++)
			v[n]=opv[n];
		v[N[np]]=0;
		return 1;
	}
	else
		return -1;
	//return max;
}

double calelecnuminwell(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np)
{
	int i,j,n;
	//double max=0;                     // �K�i���p
	double sum=0;
	double xn,px,py;
	//double pr,pi;
	gsl_complex tempA,tempB;
	//printf("\nn z[nm] Ar Ai Br Bi Yr Yi Y*Y\n");
	j=0;//p[0]=-ML;q[0]=E;
	for(n=0;n<N[np];n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			//pr=GSL_REAL(gsl_complex_add(tempA,tempB));
			//pi=GSL_IMAG(gsl_complex_add(tempA,tempB));
			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			//if(py>max)
				//max=py;                                                                       /* �g���֐���2��̍ő�l�����߂�      */
			//printf("%g %e %e %e %e %e %e %e %e\n",(double)j/DIV,px*1e9,GSL_REAL(tempA),GSL_IMAG(tempA),GSL_REAL(tempB),GSL_IMAG(tempB),pr,pi,py);
			//max+=py*dx*1e9/DIV;
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=py*dx/DIV;
			j++;
		}
	}
	//printf("\nmax %e sum %e\n",max,sum);
	return sum;
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
		//temp=q[n+i]*dx/die[mt(n+i)];
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
		//printf("%d %e %d %e %e\n",n+1,q[n+1],n,gsl_vector_get(V,n),dietemp);
	}
	gsl_linalg_LU_decomp(S, P, &n);
	gsl_linalg_LU_solve(S, P, V, X);

	for(n=0;n<div;n++)
	{
		vnew[n+1]+=gsl_vector_get(X,n);
		//printf("%d %e\n",n+1,vnew[n+1]);
	}
	gsl_vector_free(V);
	gsl_vector_free(X);
	gsl_matrix_free(S);
	gsl_permutation_free(P);

	return;
}

void selfpotential(double Q, double VRTD, double vRTD[])
{
	int np=1;								/*np��RTD�\�����A�S�̂���\���Ă���Bnp=1��RTD�\���B=0�͑f�q�S�� */
	int aa=-1;								/*�v�Z�I���̔��f�̂��߂Ɏg�p*/
	double E;
	potential0(Q,VRTD,vRTD,np);				/*�ߎ��� �����l�Ƃ��Ďg�p*/
	setpotential(vRTD,np,0);				/*�K�i�ߎ��K�p*/
	gsl_complex A[N[np]],B[N[np]],k[N[np]];	/*A�͑O�i�g�AB�͌�i�g�̌W���Bk�͔g�� */
	A[0]=gsl_complex_rect(0,0);				/* �����l�̐ݒ�                    */
	B[0]=gsl_complex_rect(1,0);				/* �����l�̐ݒ�                    */
	while(aa==-1)							/*���Ȗ������v�Z* �v�Z���������Ȃ��\������B*/
	{
		E=confinedstate(1,vRTD,np);
		cal2(E,vRTD,A,B,np);				/*�G�l���M�[E�ɂ�����g���֐��̌W��A��B���v�Z*/
		setk(E,vRTD,k,np);					/*�g���v�Z*/
		aa=makewave2(Q,E,VRTD,vRTD,k,A,B,np);/*�d��Q���l�������|�e���V�����v�Z �o�͂́A�O�̃|�e���V�����ƌv�Z��̃|�e���V�����̍���1e-6�ȉ��Ȃ�-1�@����ȊO��1*/
	}
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
	//printf("small %d large %d\n",sm,la);	
	if(smt[n]==7 || smt[n]==5 || smt[n]==23 || smt[n]==24 || smt[n]==22 || smt[n]==13 || smt[n]==6 )//i-Si smt[n]==24��ǉ� 07/28��� 8/19�Ă���n=22��ǉ�
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
	if(cond[n]==0)//����
	{
		if(direction==-1)
		{
			for(i=la;i>=sm;i--)
				v[i-1]=v[NX[n]];
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
		if( (direction<0 && D>=0) || (direction>0 && D<=0)) //��R�w
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
						//printf("��R�wMax\n");
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
			else //�G���[����
			{
				printf("\n%d�w��:��R�w���������ȏ�ł��Bw=%e\n�����𑝂₵�ĉ������B\n��R�w��Max=%e",n,w,wmax);
				int wML=wmax/ML;
				printf("=%6.4gML<%dML \n���݂̖���d[%d]=%e\n",wmax/ML,wML+1,n,d[n]);
				exit(1);
			}
			//printf("D %e w %e\n",D,w);
			return 0;
		}
		else //�~�ϑw
		{
			A=sqrt(2*KB*temp*Q[n]*1e6/die[n]);
			q=ELEC/(KB*temp);
			fai=-vacc(D,Q[n]);
			E=A*sqrt(exp(-q*fai)+q*fai-1);
			//printf("A %e fai %e E %e ��E %e D %e\n",A,fai,E,E*die[1],D);
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
					//printf("fai %e E %e tempV %e\n",fai,E,tempV);
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
		if( (direction<0 && D<=0) || (direction>0 && D>=0) ) //��R�w
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
						//printf("��R�wMax\n");
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
			else //�G���[����
			{
				printf("\n%d�w��:��R�w���������ȏ�ł��Bw=%e\n�����𑝂₵�ĉ������B\n��R�w��Max=%e",n,w,wmax);
				int wML=wmax/ML;
				printf("=%6.4gML<%dML\n",wmax/ML,wML+1);
				exit(1);
			}
			//printf("D %e w %e\n",D,w);
			return 0;
		}
		else //�~�ϑw
		{
			A=sqrt(2*KB*temp*Q[n]*1e6/die[n]);
			q=ELEC/(KB*temp);
			fai=-vacc(D,Q[n]);
			E=A*sqrt(exp(q*fai)-q*fai-1);
			//printf("A %e fai %e E %e\n",A,fai,E);
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
	//printf("D=%e nd=%e q=%e a=%e\n",D,nd,q,a);
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
			printf("�~�ϑw�̃|�e���V�����v�Z���������܂���ł����B\n�h�[�s���O�Z�x��ɒ[�ɍ����d���������Ă��Ȃ������m�F���ĉ������B");
			exit(1);
		}
		xk=Vacc;
		f1=exp(q*xk)-q*xk-1-a;
		f2=q*(exp(q*xk)-1);
		Vacc=xk-f1/f2;
		sub=fabs(Vacc-xk);
		//printf("xk=%e f1=%e f2=%e Vacc=%e sub=%e\n",xk,f1,f2,Vacc,sub);
	}
	//printf("Vacc %e\n",Vacc);
	return Vacc;
}

void wavefunction(double v[], double E,int np, double wavestore[]) //�V��������wavestore��ǉ� (2016.10.13)
{
	gsl_complex A[N[np]],B[N[np]],k[N[np]];
	cal3(E,v,A,B,np);                               /* �g���֐��̌v�Z                  */
	setk(E,v,k,np);                                 /* �g���̌v�Z                      */		
	makewave(E,v,k,A,B,np, wavestore);                         /* �g���֐��̌v�Z�ƕ\��            */
}

double func(double E, double v, int n)
{
	if(fabs(v-E)<1e-15)
		return 0;
	double t;
	t=1.0+exp(ELEC*(v+Ef[n]-E)/KB/temp);
	//printf("n %d v %e Ef[%d] %e E %e exp() %e t %e temp %e ",n,v,n,Ef[n],v+Ef[n]-E,ELEC*(v+Ef[n]-E)/KB/temp,t,log(t)); 
	return log(t);
}

double fermia(int n, double t, int flag)			/* �t�F���~�G�l���M�[�̌v�Z�@�v����							*/
{													/* �{���́A���x��ϐ��Ƃ����v�Z���s�������������A���ݓ�����	*/
	double Ed=-0.05; //�h�i�[or�A�N�Z�v�^�[���ʁB�{���Ȃ�h�[�v�ޗ��ɂ���ĕω�?����50meV�Ƃ��Ă���B

	if(flag==1 || flag==2 || flag==19 || flag==20)// n-Si�p
	{
		double ef;
		//double nc;
		//nc=2*pow((MSTAR*KB*t/HBAR/HBAR/PI/2.0),1.5)*pow(0.328,1.5)*valley[n];//*pow(0.98*0.19*0.19,0.5),1.5); //�u�����̕����vp139�`�Q�ƁB[/m3]
		//printf("nc %e\n",nc);
		//ef=KB*t*log((Q[n]*1e6)/nc)/ELEC;//printf("ef2 %e\n",ef);
		ef=-EG_SI/2+KB*t*log((Q[n]*1e6)/NI)/ELEC;//printf("ef3 %e\n",ef);
		//printf("nc %e\nQNd %e\nef %e\n",nc,Q[n]*1e6,ef);
		return ef;
	}
	if(flag==3 || flag==4 || flag==10)//p-Si�p
	{
		double nv,ef;
		nv=2*pow((MSTAR*KB*t/HBAR/HBAR/PI/2.0),1.5)*(pow(0.16,1.5)+pow(0.49,1.5));
		ef=Ed+KB*t*log(-1.0/4.0+pow(1+8*Q[n]*1e6*exp(-Ed/KB/t)/nv,0.5)/4);
		//efhh=KB*t*log(nv/(Q[n]*1e6));
		return ef/ELEC;
	}
	return -0.56;
}

double calconfinedstate(double v[], int np, double E1, double E2)
{
	/*	E1����E2�̃G�l���M�[�͈͂ɂ����āA�ł��Ⴂ���ʂ�T��	*/
	/*	�Ȃ��ꍇ�ɂ́AEMAX��Ԃ��B								*/
	/*	�ϐ�	v	�|�e���V����								*/
	/*			np	0�Ȃ��ԑS�́A1�Ȃ�RTD�����̂�				*/
	double E,T,S,dE;
	E=E1;											/* �T���͈͂̉���													*/
	while(E<E2)
	{											
		dE=1e-3;									/* ���ݕ�(�e��)�̐ݒ�                                              */
		T=0;
		S=cal(E,v,np,0);
		while(E<E2)									/* ���ߗ��̌�����Ԃ̌���                                          */
		{											/* �ȉ��̏��ʂ����߂�v���O������                                  */
			E+=dE;									/* �G�l���M�[���㏸���Ă����ƂƂ��ɁA���ߗ����㏸���Ă����Ɖ��肵�A*/
			T=cal(E,v,np,0);						/* ���̍ő�l���Ƃ�ӏ������ʂł���Ƃ��Ă���                      */
			if(S>T)									/* ����āA���ߗ����������Ă�����Ԃ͕s�v                          */
				S=T;
			else
				break;								/* ���ߗ��̌�����Ԃ̏I�� */
		}
		S=T=0;										/* S�AT�̏�����                          */
		while(E<E2)									/* ��ꏀ�ʂ̒T�� ���ݕ��F�e��           */
		{
			E+=dE;									/* �G�l���M�[�̐ݒ�                      */
			T=cal(E,v,np,0);						/* ���ߗ��̌v�Z                          */
			if(S<T)									/* �����̔��f�@�����֐��Ȃ��while�����s */
				S=T;
			else
			{
				E-=2*dE;							/* �G�l���M�[�������߂��āA���ݕ����ׂ������Ă��� */
				dE/=10;
				S=0;
				if(dE<1e-10)						/*�@*/
					break;
			}
		}
		dE/=10;										/* ���ݕ�(�ׂ���)�̐ݒ�                           */
		S=0;										/* S�̏�����                                      */
		while(E<E2)									/* ��ꏀ�ʂ̒T�� ���ݕ��F�ׂ���                  */
		{
			E+=dE;									/* �G�l���M�[�̐ݒ�                               */
			T=cal(E,v,np,0);						/* ���ߗ��̌v�Z                                   */
			if(S<T)
				S=T;								/* �����̔��f�@�����֐��Ȃ��while�����s          */
			else
			{
				//printf("����%d %+20.18lf T %e\n",i,E,S);
				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
				{									/*cal�֐��̉e���ŁAE=v[0]�̂Ƃ��낪���ʂɌ����Ă��܂��̂�*/
					E+=10*dE;						/*���������邽�߂ɓ���*/
					break;							/*���m�ł��Ȃ��\������*/
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
{											/* ��n���ʂ��o��												*/
	int i;
	double E,E1,E2;
	E1=E2=EMAX;
	for(i=0;i<N[np]+1;i++)										
	{
		if(E1>v[i])							/* �G�l���M�[�̒T���͈͂̌���									*/
			E1=v[i];						/* �������ꂽ�|�e���V�����̍ŏ��l�����߂�						*/
	}
	for(i=0;i<n;i++)
	{
		E=calconfinedstate(v,np,E1,E2); //printf("E %e \n",E);
		E1=E+1e-9;
	}
	return E;
}

void getconfinedstates(int n, double v[], int np, double En[])
{											/* �����̏��ʂ��擾���A�z��ɏo��								*/
	int i;
	double E1,E2;
	E1=E2=EMAX;
	for(i=0;i<N[np]+1;i++)										
	{
		if(E1>v[i])							/* �G�l���M�[�̒T���͈͂̌���									*/
			E1=v[i];						/* �������ꂽ�|�e���V�����̍ŏ��l�����߂�						*/
	}
	for(i=0;i<n;i++)
	{
		En[i]=calconfinedstate(v,np,E1,E2);//printf("En[%d] %e\n",i,En[i]);
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

			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]���v�Z
			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]���v�Z

			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //�G�l���M�[E�A�̈�n+1�ɂ�����g�����v�Z HBAR�Ŋ����Ă��Ȃ����Ƃɒ���

			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)���v�Z
			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)���v�Z
			
			gsl_matrix_complex_set(temp[n],0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)���v�Z �s��temp(1,1)�Ɋi�[
			gsl_matrix_complex_set(temp[n],0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)���v�Z �s��temp(1,2)�Ɋi�[
			gsl_matrix_complex_set(temp[n],1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)���v�Z �s��temp(2,1)�Ɋi�[
			gsl_matrix_complex_set(temp[n],1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)���v�Z �s��temp(2,2)�Ɋi�[

			gsl_matrix_complex_scale(temp[n], gsl_complex_rect(0.5,0));                         //�S�Ă̗v�f��0.5��������
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],dummy,gsl_complex_rect(0,0),trans); //�s��̊|���Z

		gsl_matrix_complex_memcpy(dummy, trans);
		/*printf("Transfer Matrix n=%d material %s\nT(0,0)=%e +i %e T(0,1)=%e +i %e\nT(1,0)=%e +i %e T(1,1)=%e +i %e \n\n",n,name[mt(n+1)],
				GSL_REAL(gsl_matrix_complex_get(temp[n],0,0)), GSL_IMAG(gsl_matrix_complex_get(temp[n],0,0)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],0,1)), GSL_IMAG(gsl_matrix_complex_get(temp[n],0,1)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],1,0)), GSL_IMAG(gsl_matrix_complex_get(temp[n],1,0)),
				GSL_REAL(gsl_matrix_complex_get(temp[n],1,1)), GSL_IMAG(gsl_matrix_complex_get(temp[n],1,1)));

		printf("Tran n=%d material %s\nT(0,0)=%e +i %e T(0,1)=%e +i %e\nT(1,0)= %e +i %e T(1,1)=%e +i %e \n\n",n,name[mt(n+1)],
				GSL_REAL(gsl_matrix_complex_get(trans,0,0)), GSL_IMAG(gsl_matrix_complex_get(trans,0,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,0,1)), GSL_IMAG(gsl_matrix_complex_get(trans,0,1)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,0)), GSL_IMAG(gsl_matrix_complex_get(trans,1,0)),
				GSL_REAL(gsl_matrix_complex_get(trans,1,1)), GSL_IMAG(gsl_matrix_complex_get(trans,1,1)));*/
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
			//printf("\nA[0] 1+0i B[0] 1-r\n");
			break;//��������E���ɓ��߂���ꍇ
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
			//printf("\nA[0] 0+0i B[0] 1+0i\n");
			break;//�E�����獶���ɓ��߂���ꍇ
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
			//printf("\nA[0] 1+0i B[0] 1+0i\n");
			break;//�v�����Bv[0]<v[N+1]�̂Ƃ��K�p�\�H2015/3/25
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
			//printf("\nA[0] 1+0i B[0] 0+0i\n");
			break;//�ȈՔŁB���ߗ���1�̏ꍇ�K�p�ł���H
	}
	gsl_matrix_complex *FN  = gsl_matrix_complex_calloc(2,1);
	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);
	
	for(n=0; n<N[np]-1; n++)
	{
		gsl_matrix_complex_set(FN,0,0,A[n]);
		gsl_matrix_complex_set(FN,1,0,B[n]);
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],FN,gsl_complex_rect(0,0),FNN); //�s��̊|���Z FNN=temp*FN
		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
		B[n+1]=gsl_matrix_complex_get(FNN,1,0);
	}
	//printf("T %e R %e\n", sqrt(1-v[0]/E)/gsl_complex_abs2(gsl_matrix_complex_get(trans,1,1)),gsl_complex_abs2(B[0]) );
	for(n=0; n<N[np]; n++){
		gsl_matrix_complex_free(temp[n]);
	}
	gsl_matrix_complex_free(dummy);
	gsl_matrix_complex_free(trans);

	setk(E,v,k,np);

	//double max=0;                     // �K�i���p
	double sum=0;
	double xn,px,py;
	//double pr,pi;
	gsl_complex tempA,tempB;
	//printf("\nn z[nm] Ar Ai Br Bi Yr Yi Y*Y\n");
	j=0;//p[0]=-ML;q[0]=E;
	for(n=0;n<N[np];n++)
	{
		xn=(n+1)*dx;
		for(i=0; i<DIV; i++)
		{
			px=(double)j*dx/DIV;
			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
			//pr=GSL_REAL(gsl_complex_add(tempA,tempB));
			//pi=GSL_IMAG(gsl_complex_add(tempA,tempB));
			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
			//if(py>max)
				//max=py;                                                                       /* �g���֐���2��̍ő�l�����߂�      */
			//printf("%g %e %e %e %e %e %e %e %e\n",(double)j/DIV,px*1e9,GSL_REAL(tempA),GSL_IMAG(tempA),GSL_REAL(tempB),GSL_IMAG(tempB),pr,pi,py);
			//max+=py*dx*1e9/DIV;
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=(double)py*dx/DIV;
			j++;
		}
	}
	//printf("\nmax %e sum %e\n",max,sum);
	return sum*func(E,v[0],0)*sqrt(1/(E-v[0]));
}

int max(int a, int b) {
    return (a > b) ? a : b;
}
int min(int a, int b) {
    return (a < b) ? a : b;
}
