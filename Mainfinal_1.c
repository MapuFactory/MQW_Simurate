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
	clock_t start,end;							/* プログラムの計算時間の測定	*/
	start = clock();							/* 開始時刻取得					*/
	int p=2; printf("プログラムNo.%d\n",p);		/* プログラムの番号				*/
	setmaterial(-100);							/* 想定している構造や材料の取得　引き数が-100だと出力されない */
	int a,j,m,n, loop; //齋藤が変数loop(ループを回す為の変数)を追加 (2016.10.13)
	double I,QW,VRTD,V, dv;
	double v[N[0]+1],vRTD[N[1]+1];				/* ポテンシャルの格納用			*/
	for(m=0;m<1;m++)
	{
		if(p==0)								/* 電流計算用*/
		{
			if(OUTPUT==1)
			{
				printf("OUTPUT=1です。強制終了します。\nポテンシャル計算と出力はプログラムNo.1でやってください。\n今はp=0\n");
				exit(1);
			}
			else
			{
				setmaterial(0);					/* 材料を再度設定　主に物性値を変更するために導入*/
				QW=0*ELEC*1e4*5e12;				/* 井戸内に蓄積されている電荷の設定　set.datファイル中でも設定できるが、変数としたい場合はここで改めて設定*/
				printf("\nNQW[/cm2],Q[C/m2]\n%e,%e\n\nVRTD,V,I[A/cm2]\n",QW*1e-4/ELEC,QW);
				dv=0.0025;
				for(j=-3;j<400;j++)				/* 電流計算開始 */
				{	
					VRTD=dv*j;					/* RTDに印加されている電圧を決定　その後他の層の電圧を計算し、素子にかかっている電圧を計算する。　全体にかかっている電圧からポテンシャルを計算することは困難 */
					potential(v,vRTD,VRTD,QW);	/* ポテンシャル計算	自己無撞着な計算を導入している　計算が収束しない可能性あり	*/
					V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);	/* 電圧計算　両サイドのフェルミエネルギーの位置から電圧を計算			*/
					//printf("%e %e\n",VRTD,V);
					if(V<0)
						continue;
					else if(V<5)//+(4.7-4.05+Ef[0]))
					{
						I=current(v);			/* 電流計算	電流を計算するためには、ポテンシャルが引き数として必要 ただし、定数がかかっていないので注意*/
					}
					else
						break;
					//printf("%e %e %e\n",VRTD,V,III*I/1e4);	/* 出力*/
					printf("%e,%e,%e\n",VRTD,V,II*temp*I/1e4);	/* 出力*/
				}
			}
		}
		else if(p==1)							/* ポテンシャルと透過率計算用*/
		{
			double E,T,dE;
			setmaterial(0);						/* 材料を再度設定　主に物性値を変更するために導入　特に変更する必要がなければ不要							*/
			QW=0*ELEC*1e4*1e13;					/* 井戸内に蓄積されている電荷の設定　set.datファイル中でも設定できるが、変数としたい場合はここで改めて設定	*/
			VRTD=2.855;
			potential(v,vRTD,VRTD,QW);
			//setpotential(vRTD,1,1);
			V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);
			printf("\nVRTD,V,Q[/cm2],Q[C/m2]\n%e,%e,%e,%e\n\n\nE,T\n",VRTD,V,QW*1e-4/ELEC,QW);
			dE=1e-6;							/* エネルギーの分割	この値を適切にとること　一般的には細かくとれば良い　その分計算に時間がかかる*/
			for(j=3200000;j<3500000;j++)
			{
				E=dE*j;					/* エネルギー[eV]の設定	*/
				T=cal(E,v,0,1);					/* 透過率計算			*/
				printf("%e,%e\n",E,T);
			}
		}
		else if(p==2)							/* 波動関数計算用*/
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
			n=atoi(argv[3]);								/* 計算したい準位の数	*/
			double En[n];						/* 準位の格納			*/
			double wavestore[n][DIV*N[0]];  //波動関数格納
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
						printf("第%d準位,",j+1);
					else
						printf("%e,",En[j]);
				}
				printf("\n");
			}
			for(j=0;j<n;j++)
				wavefunction(v,En[j],0, wavestore[j]); //齋藤が引数wavestore[j]を追加 (2016.10.13)

				
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

			setmaterial(0);						/* 材料を再度設定　主に物性値を変更するために導入　特に変更する必要がなければ不要							*/
			QW=m*ELEC*1e4*5e12;					/* 井戸内に蓄積されている電荷の設定　set.datファイル中でも設定できるが、変数としたい場合はここで改めて設定	*/
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
		else if(p==4)//テストモード
		{
			int i,iE;
			double C;
			//double E;
			double Enum;
			double dE;
			double NQ;
			//gsl_complex A[N[np]],B[N[np]],k[N[np]];
			setmaterial(m);					/* 材料を再度設定　主に物性値を変更するために導入*/
			QW=0*ELEC*5e4*1e12;				/* 井戸内に蓄積されている電荷の設定　set.datファイル中でも設定できるが、変数としたい場合はここで改めて設定*/
			C=massxy[0]*MSTAR*KB*temp/HH/HH/1e4*sqrt(2*mass[0]*MSTAR/ELEC)/HBAR*ELEC;
			iE=1e5;
			dE=VI/iE;
			printf("\nNQW[/cm2],Q[C/m2]\n%e,%e\n\nVRTD,V,NQ[A/cm2],iE,dE\n",QW*1e-4/ELEC,QW);
			for(j=-4;j<20;j++){
				NQ=0;
				VRTD=0.05*j;					/* RTDに印加されている電圧を決定　その後他の層の電圧を計算し、素子にかかっている電圧を計算する。　全体にかかっている電圧からポテンシャルを計算することは困難 */
				potential(v,vRTD,VRTD,QW);	/* ポテンシャル計算	自己無撞着な計算を導入している　計算が収束しない可能性あり	*/
				V=v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);	/* 電圧計算　両サイドのフェルミエネルギーの位置から電圧を計算			*/
				//E=calconfinedstate(v,0,v[0],EMAX);
				for(a=0;a<1;a++)
				{
					NQ=0;
					// #pragma omp parallel for private(Enum) reduction(+:NQ) OpenMPが使えない為、コメントアウトにした。(2016/09/25)
					for(i=1;i<iE+1;i++)
					{
						Enum=v[0]+dE*i;
						NQ+=cal4(Enum,v,0)*dE;
						//printf("%d %e %e\n",i,Enum,cal4(Enum,v,0));
					}
					printf("%e,%e,%e,%d,%e\n",VRTD,V,valley[0]*C*NQ,iE,dE);	/* 出力*/
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
	end = clock(); printf("\n演算時間 %f 秒\n",(double)(end-start)/CLOCKS_PER_SEC);
	//Beep( 440, 1200 );
	return 0;
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
	i=14;mtconst[i].name="仮想AL   "; mtconst[i].die=DIE_AL   * DIEELECSTAR; mtconst[i].bar=-5.17  ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=15;mtconst[i].name="仮想CaF2 "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=7      ;  mtconst[i].mass=1;         mtconst[i].massxy=1;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=16;mtconst[i].name="仮想p-Si "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;      mtconst[i].massxy=0.55;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=17;mtconst[i].name="SiO2     "; mtconst[i].die=DIE_SIO2 * DIEELECSTAR; mtconst[i].bar=BAR_SIO2; mtconst[i].mass=MASS_SIO2; mtconst[i].massxy=MASS_SIO2;	mtconst[i].valley=1;			mtconst[i].cond=3;
	i=18;mtconst[i].name="nc-Si.   "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=19;mtconst[i].name="nSi100Z  "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_Z; mtconst[i].massxy=MASS_SI_XY;	mtconst[i].valley=NUMVALLY_SI_Z;	mtconst[i].cond=1;
	i=20;mtconst[i].name="nSi100XY "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_SI_XY;mtconst[i].massxy=MASS_SI_Z;	mtconst[i].valley=NUMVALLY_SI_XY;	mtconst[i].cond=1;
	i=21;mtconst[i].name="n-Si.Sub."; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=MASS_iSI;  mtconst[i].massxy=MASS_AL;		mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;
	i=22;mtconst[i].name="CaF2(4ML)"; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.6;	  mtconst[i].mass=0.85;	     mtconst[i].massxy=0.85;		mtconst[i].valley=1;			mtconst[i].cond=3;
	i=23;mtconst[i].name="p-i-Si   "; mtconst[i].die=DIE_iSI  * DIEELECSTAR; mtconst[i].bar=BAR_iSI;  mtconst[i].mass=0.55;	     mtconst[i].massxy=MASS_iSI;	mtconst[i].valley=NUMVALLY_SI;		mtconst[i].cond=0;	
	i=24;mtconst[i].name="3MLCaF2  "; mtconst[i].die=DIE_CAF2 * DIEELECSTAR; mtconst[i].bar=1.5;      mtconst[i].mass=1.0;       mtconst[i].massxy=1.0;		mtconst[i].valley=1;			mtconst[i].cond=3;
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
			printf("層番号,材料番号,物質名,ML数,層厚[nm],比誘電率,有効質量,障壁の高さ,谷,分割数,NX,EF,仕事関数,電荷量\n");
		for(i=0; fscanf(fd1, "%d %d %lf\n", &smt[i], &j, &Q[i])!=EOF; i++)
		{
			if( modf((double)j*DX,&d[i])!=0)
			{
				printf("分割エラー　分割間隔DXを変更して下さい\n");
				exit(1);
			}
			divnum[i] = j*DX;					/* 分割数					*/
			N[0]        += divnum[i];			/* 総分割数					*/
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
				case 14: Ef[i]=4.37;break;
				case 16: Ef[i]=-EG_SI/2+KB*temp*log(Q[i]*1e6/NI)/ELEC;base=BASE;break;
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
				printf("%d,%d,%s,%d,%.4g,%.3g,%.2g,%.4g,%d,%d,%d,%.4g,%.4g,%.4g\n", i, smt[i], name[i], j, j*ML*1e9, die[i]/DIEELECSTAR, mass[i], bar[i], valley[i], divnum[i], NX[i],Ef[i],Work[i],Q[i]);//出力。
		}
		fclose(fd1);
	}
	else 
	{
		printf("ファイルが開けていません\n");
		exit(1);
	}
	N[1]=divnum[LBAR]+divnum[WELL]+divnum[RBAR];	/* RTD構造のみの層数*/
	if(m!=flag)
		printf("\n 1ML[nm],分割数,積分範囲,ΔE,ピーク時ΔE,N,NRTD TEMP,II*temp,bar[%d],bar[%d]\n%.3g,%d,%.2g,%.2g,%.2g,%d,%d,%.3g,%.2g ",layer,LBAR-2,ML*1e9,DX,VI,DELTA,DELTAE,N[0],N[1],temp,II*temp);
	if(AAA>0){									/* フェルミエネルギーピンニングの現象を再現。要検討。基本使わない方が良い?*/
		if(d[RBAR+1]!=0)
		{
			bar[layer]-=Work[layer-1]-Work[layer];
			if(m!=flag)
				printf("%.3g,",bar[layer]);
		}
		else
		{
			if(m!=flag)
				printf("なし,");
		}
		if(d[LBAR-1]!=0)//分岐必要
		{
			bar[LBAR-2]-=Work[LBAR-1]-Work[LBAR-2];
			if(m!=flag)
				printf("%.3g\n",bar[LBAR-2]);
		}
		else
		{
			if(m!=flag)
				printf("なし\n");
		}
	}
	else
		if(m!=flag)
			printf("なし,なし\n");

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
	int i,n;										/* E1からE2までの電流値を計算。刻み幅は引数delta。											*/
	double dE,E,T,S1,S2,I;							/* deltaが透過率の半値幅よりも大きい場合、数値計算で積分を実行すると誤差が大きくなるので、	*/
	I=0;											/* deltaは十分小さくとること。ただし、小さすぎると計算に時間がかかる。						*/
	n=(E2-E1)/delta+1;								/* 刻み幅の決定																				*/
	dE=(E2-E1)/(n-1);								/* 誤差あり。あまり大きくないので放置。														*/
	//printf("E1 %e E2 %e n %d dE %e\n",E1,E2,n,dE);
// #pragma omp parallel for private(E,T,S1,S2) reduction(+:I) OpenMPが使えない為、コメントアウトにした(2016/09/25)	/* openmpが使えない場合は、この行は無視されるので問題ない。		*/
	for(i=0;i<n;i++)
	{
		E=E1+i*dE;
		T=cal(E,v,0,1);													/* 透過率の計算                             */
		S1=massxy[0]*valley[0]*func(E,v[0],0);								/* Supply functionの計算                    */
		S2=massxy[mt(N[0]-1)]*valley[mt(N[0]-1)]*func(E,v[N[0]-1],layer);	/* Supply functionの計算                    */
		if(i==0 || i==n-1)
			I+=0.5*T*(S1-S2);
		else
			I+=T*(S1-S2);
	}
	I*=dE;
	return I;
}

double current(double v[])						/* 電流値の計算。ここでは、積分範囲を決定している。*/
{
	int i;
	int n=1;
	int num=10;									/* 存在する準位の数					*/
	double Emin,Emax,En[num],E1,E2,max1,max2;
	double I=0;

	if(DELTA<DELTAE)
	{
		printf("DELTA < DELTAEのため計算できませんでした。\n強制終了します。\nDELTAEを適切な値にして下さい\n");
		exit(0);
	}
	max1=max(v[0],v[0]+Ef[0]);					//printf("v[0] %e v[0]+Ef[0] %e\n",v[0],v[0]+Ef[0]);
	max2=max(v[N[0]-1],	v[N[0]-1]+Ef[layer]);	//printf("v[N] %e v[N]+Ef[N] %e\n",v[N[0]-1],v[N[0]-1]+Ef[layer]);
	Emax=max(max1,max2)+VI;						//printf("Emax %e max1 %e max2 %e\n",Emax,max1,max2);
	if(v[0]>v[N[0]-1])							/* 積分範囲の最小値の設定			*/
		Emin=v[0];
	else if(v[0]<v[N[0]-1])
		Emin=v[N[0]-1];
	else
		Emin=v[0];
	getconfinedstates(num,v,0,En);				/* 準位の探索。配列に格納。			*/
	E1=Emin;
	for(i=0;i<num;i++)
	{
		if(En[i]<E1)							/*積分区間内にEn[i]の準位なし									*/
			continue;
		else if(Emax<En[i])
		{
			I+=calcurrent(E1,Emax,v,DELTA);		/*積分区間の最後まで積分										*/
			break;
		}
		else if(E1<En[i] && En[i]<Emax)			/*積分区間内に準位あり　細かくエネルギーを分割して計算する必要あり*/
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
		printf("想定していない値が引き数npに使われています。\n");
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
	int n,nrtd;								/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
	nrtd=nRTD(np);
	double vr[N[np]+1],vl[N[np]+1];			/* ポテンシャルの左側からの極限vrと右側からの極限			*/
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
		v[n]=(vr[n]+vl[n+1])/2;				/* 階段近似適用	*/
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
	int n,nrtd;					/* ポテンシャルから波数kを計算して格納*/
	nrtd=nRTD(np);
	for(n=0; n<N[np]; n++)
	{
		k[n]=gsl_complex_div_real(gsl_complex_sqrt_real(2*mass[mt(n+nrtd+1)]*(E-v[n])*MSTAR*ELEC),HBAR);
		//printf("name[%d] %s,k[%d]=%g +i %g\n",n,name[mt(n+1)],n,GSL_REAL(k[n]),GSL_IMAG(k[n]));
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

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //全ての要素に0.5をかける
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			//printf("kn=0,kn1=0\n");
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納
		}
		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,dummy,gsl_complex_rect(0,0),trans); //行列の掛け算

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

			gsl_matrix_complex_set(temp,0,0,gsl_complex_mul(pm,zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,zp);                     //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_mul(pp,zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,zm);                     //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

			gsl_matrix_complex_scale(temp, gsl_complex_rect(0.5,0)); //全ての要素に0.5をかける
		}

		else if(E-v[n]==0 && E-v[n+1]==0)
		{
			//printf("kn=0,kn1=0\n");
			gsl_matrix_complex_set(temp,0,0,gsl_complex_rect(mass[mt(n+nrtd+1)]/mass[mt(n+nrtd)],0));    //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
			gsl_matrix_complex_set(temp,0,1,gsl_complex_rect(0,0));                            //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
			gsl_matrix_complex_set(temp,1,0,gsl_complex_rect(mass[mt(n+nrtd+1)]*dx/mass[mt(n+nrtd)],0)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
			gsl_matrix_complex_set(temp,1,1,gsl_complex_rect(1,0));                            //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納
		}

		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp,FN,gsl_complex_rect(0,0),FNN); //行列の掛け算 FNN=temp*FN

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
			printf("\nA[0],1+0i,B[0],1-r\n");break;//左側から右側に透過する場合
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
			printf("\nA[0],0+0i,B[0],1+0i\n");break;//右側から左側に透過する場合
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
			printf("\nA[0],1+0i,B[0],1+0i\n");break;//要検討。v[0]<v[N+1]のとき適用可能？2015/3/25
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
			printf("\nA[0],1+0i,B[0],0+0i\n");break;//簡易版。透過率が1の場合適用できる？
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
	//printf("min=%e\n",E);
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
		//printf("j=%d\n",j);
		dE/=10;                                      /* 刻み幅(細かい)の設定                           */
		S=0;                                          /* Sの初期化                                      */
		while(E<emax)                                 /* 基底準位の探索 刻み幅：細かい                  */
		{
			E+=dE;                                    /* エネルギーの設定                               */
			T=cal(E,v,np,0);                              /* 透過率の計算                                   */
			if(S<T)
				S=T;                                  /* 増減の判断　増加関数ならばwhile文続行          */
			else{
				//printf("準位%d %+20.18lf T %e\n",i,E,S);
				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
				{										/*cal関数の影響で、E=v[0]のところが準位に見えてしまうので*/
					E+=10*dE;							/*それを避けるために導入*/
					i--;								/*従って、準位がv[0]±1e-10[eV]の近辺に存在した場合、*/
					//printf("aS %e T %e E %e E-dE-v[0] %e E-dE-v[N[np]-1]%e \n",S,T,E,E-dE-v[0],E-dE-v[N[np]-1]);
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
	//printf("準位%d %e T %e\n",i,E,S);                 /* 準位の出力                                    */
	return E;                                         /* num番目の準位を返す                           */
}

void makewave(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np, double wavestore[]) //齋藤が引数wavestoreを追加 (2016.10.13)
{
	int i,j,n;
	double max=0;                     // 規格化用
	double sum=0;
	double xn,px,py, pyy[N[np]*DIV]; //齋藤が配列追加(pyy, 2016.10.13)
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
				//max=py;                                                                       /* 波動関数の2乗の最大値を求める      */
			//printf("%g %e %e %e %e %e %e %e %e\n",(double)j/DIV,px*1e9,GSL_REAL(tempA),GSL_IMAG(tempA),GSL_REAL(tempB),GSL_IMAG(tempB),pr,pi,py);
			max+=py*dx*1e9/DIV;
			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
				sum+=py*dx*1e9/DIV;
			j++;
		}
	}j=0;//2015/5/20追加　思考停止 とりあえず波動関数の面積を求める(↑のmaxに格納)。波動関数の2乗の全区間で積分すると1になるように規格化。
	//printf("\nz[nm] Y*Y Y*Y\n"); //量子井戸の場合、膜厚が厚すぎるとexp(±ikz)の掛け算でオーバーフロー(？)を起こし、減衰項ではなく増幅項が支配的になるので注意。
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
	double max=0;                     // 規格化用
	double qmax=0;
	double xn,p[DIV*N[np]],q[DIV*N[np]]; // pがx座標、qがy座標 qは波動関数の絶対値の2乗を出力予定
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
				tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],p[j]-xn)));     /* tempA= A[n] × exp( ik[n](x-x[n])) */
				tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-p[j])));     /* tempB= B[n] × exp(-ik[n](x-x[n])) */
				q[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));                                /* 波動関数の絶対値の2乗を計算        */
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
	if(max<1e-6){			/* 重要　自己無撞着計算の収束判定*/
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
	//double max=0;                     // 規格化用
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
				//max=py;                                                                       /* 波動関数の2乗の最大値を求める      */
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
	int np=1;								/*npはRTD構造か、全体かを表している。np=1はRTD構造。=0は素子全体 */
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
	if(smt[n]==7 || smt[n]==5 || smt[n]==23 || smt[n]==24 || smt[n]==22 || smt[n]==13 || smt[n]==6 )//i-Si smt[n]==24を追加 07/28大野 8/19ていがn=22を追加
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
						//printf("空乏層Max\n");
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
			//printf("D %e w %e\n",D,w);
			return 0;
		}
		else //蓄積層
		{
			A=sqrt(2*KB*temp*Q[n]*1e6/die[n]);
			q=ELEC/(KB*temp);
			fai=-vacc(D,Q[n]);
			E=A*sqrt(exp(-q*fai)+q*fai-1);
			//printf("A %e fai %e E %e εE %e D %e\n",A,fai,E,E*die[1],D);
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
						//printf("空乏層Max\n");
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
			//printf("D %e w %e\n",D,w);
			return 0;
		}
		else //蓄積層
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
			printf("蓄積層のポテンシャル計算が収束しませんでした。\nドーピング濃度や極端に高い電圧をかけていないか等確認して下さい。");
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
	//printf("n %d v %e Ef[%d] %e E %e exp() %e t %e temp %e ",n,v,n,Ef[n],v+Ef[n]-E,ELEC*(v+Ef[n]-E)/KB/temp,t,log(t)); 
	return log(t);
}

double fermia(int n, double t, int flag)			/* フェルミエネルギーの計算　要検討							*/
{													/* 本来は、温度を変数とした計算を行いたかったが、現在凍結中	*/
	double Ed=-0.05; //ドナーorアクセプター準位。本来ならドープ材料によって変化?今は50meVとしている。

	if(flag==1 || flag==2 || flag==19 || flag==20)// n-Si用
	{
		double ef;
		//double nc;
		//nc=2*pow((MSTAR*KB*t/HBAR/HBAR/PI/2.0),1.5)*pow(0.328,1.5)*valley[n];//*pow(0.98*0.19*0.19,0.5),1.5); //「半導体物性」p139～参照。[/m3]
		//printf("nc %e\n",nc);
		//ef=KB*t*log((Q[n]*1e6)/nc)/ELEC;//printf("ef2 %e\n",ef);
		ef=-EG_SI/2+KB*t*log((Q[n]*1e6)/NI)/ELEC;//printf("ef3 %e\n",ef);
		//printf("nc %e\nQNd %e\nef %e\n",nc,Q[n]*1e6,ef);
		return ef;
	}
	if(flag==3 || flag==4 || flag==10)//p-Si用
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
				//printf("準位%d %+20.18lf T %e\n",i,E,S);
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
		E=calconfinedstate(v,np,E1,E2); //printf("E %e \n",E);
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
			break;//左側から右側に透過する場合
	case 1: A[0]=gsl_complex_rect(0,0);
			B[0]=gsl_complex_rect(1,0);
			//printf("\nA[0] 0+0i B[0] 1+0i\n");
			break;//右側から左側に透過する場合
	case 2: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(1,0);
			//printf("\nA[0] 1+0i B[0] 1+0i\n");
			break;//要検討。v[0]<v[N+1]のとき適用可能？2015/3/25
	case 3: A[0]=gsl_complex_rect(1,0);
			B[0]=gsl_complex_rect(0,0);
			//printf("\nA[0] 1+0i B[0] 0+0i\n");
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
	//printf("T %e R %e\n", sqrt(1-v[0]/E)/gsl_complex_abs2(gsl_matrix_complex_get(trans,1,1)),gsl_complex_abs2(B[0]) );
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
				//max=py;                                                                       /* 波動関数の2乗の最大値を求める      */
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
