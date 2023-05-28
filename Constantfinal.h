#define SETFILE "set.dat"

#define PI          3.1415926      /* PI                              */
#define PIS         9.8696044      /* PI^2                            */
//#define ML 1e-10
#define ML          0.31e-9        /* [m]   [ML]→[m]                 */
//#define ML          0.155e-9        /* [m]   [ML]→[m]                 */
//#define ML          0.2715e-9        /* [m]   [ML]→[m]                 */
#define MSTAR       9.109e-31      /* [kg]                            */
#define ELEC        1.602e-19      /* [C]   Electron charge           */
#define DIEELECSTAR 8.854e-12      /* [F/m] Vacum Dielectric Constant */
#define HBAR        1.054e-34      /* [Js]                            */
#define HH	6.626e-34		/*[Js]								*/
#define KB          1.38e-23       /* [J/K]                           */
#define TEMP        3.0e+2            /* [K]                             */
#define KT          (KB*TEMP)
//#define III (ELEC*ELEC*MSTAR*KB*TEMP/2/PIS/HBAR/HBAR/HBAR)
/* 1.08872e12*/
#define II (ELEC*ELEC*MSTAR*KB/2/PIS/HBAR/HBAR/HBAR)
#define VI          0.6            /* 電流計算時の積分範囲 */
#define DELTA       1e-5
#define DELTAE      1e-9
#define DX          10				/* 分割間隔 = ML/DX   DXが大きいと分割間隔が細かい */
#define OUTPUT      1             /* ポテンシャルの出力の有無                        */

#define NI  1.45e16                   /* [/m3] */
#define EG_SI  1.12                     /* [eV]   */

#define BASE 4.05
#define pBASE (-5.17)

#define DIE_CAF2  6.76
#define EAFF_CAF2 (4.05-1.0) /* 2ML物性値 Delta Ec = 1.0 [eV] 適用 */
#define BAR_CAF2  -(EAFF_CAF2-BASE)
#define MASS_CAF2 1

#define DIE_CDF2  8.83
#define EAFF_CDF2 (4.05+0.6)
#define BAR_CDF2  -(EAFF_CDF2-BASE)
#define MASS_CDF2 0.4

#define DIE_iSI   (11.8)
#define EAFF_iSI  4.05
#define BAR_iSI   -(EAFF_iSI-BASE)
#define MASS_iSI  0.26
#define NUMVALLY_SI 6
#define MASS_SI_Z  0.98
#define MASS_SI_XY 0.19

#define DIE_AL    0
#define EAFF_AL  15.83
#define BAR_AL   -(EAFF_AL-BASE)
#define MASS_AL  1
#define EF_AL    (11.63)

#define DIE_AU    0
#define EAFF_AU  10.33
#define BAR_AU   -(EAFF_AU-BASE)
#define MASS_AU  1
#define EF_AU    5.51

#define DIE_pSI   (11.8)
#define EAFF_pSI  -(4.05+1.12)
#define BAR_pSI   -(EAFF_pSI-pBASE)
#define MASS_pSI  0.55

#define DIE_pAL    0
#define EAFF_pAL   0
#define BAR_pAL   -(EAFF_pAL-pBASE)
#define MASS_pAL  1
#define EF_pAL    4.37

#define DIE_pCAF2  6.76
#define EAFF_pCAF2 -(4.05+1.12+1)
#define BAR_pCAF2  -(EAFF_pCAF2-pBASE)
#define MASS_pCAF2 1

#define DIE_SIO2   3.9
#define EAFF_SIO2  (4.05-2.9)
#define BAR_SIO2   -(EAFF_SIO2-BASE)
#define MASS_SIO2  0.26
#define NUMVALLY_SI 6
#define NUMVALLY_SI_Z 2
#define NUMVALLY_SI_XY 4

//#define DIE_CoSi2   8.5
//#define EAFF_CoSi2  17.7
//#define BAR_CoSi2   -(EAFF_CoSi2-BASE)
//#define MASS_CoSi2  1
//#define EF_CoSi2    (13.0)

#define LBAR	2
#define WELL	(LBAR+1)
#define RBAR	(LBAR+2)
#define DIM		2
#define LAYER	100 //利用する材料の層数の最大数．(ML数ではない)
#define DIV     10
#define EMAX	10

#define AAA	0

#define RTD 1
#define ALL 0

void setmaterial(int m);
int mt(int n);
int nRTD(int np);
void potential0 (double Q, double VRTD, double v[], int np);
void setpotential (double v[], int np, int output);
//void outputpotential (double v[],int np);
void setk(double E, double v[],gsl_complex k[],int np);
double cal(double E, double v[], int np, int f);
void cal2(double E, double v[], gsl_complex A[], gsl_complex B[], int np);
void cal3(double E, double v[], gsl_complex A[], gsl_complex B[], int np);
double confinedstate(int num, double v[], int np);
void makewave(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np, double wavestore[]); //齋藤が引数wavestore[]を追加 (2016.10.13)
int makewave2(double Q, double E, double VRTD, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[],int np);
double calelecnuminwell(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np);
void calpotential(double vnew[], double q[], int np);
void selfpotential(double Q,double VRTD,double vRTD[]);
double calVRL(int direction,int n, double D, double v[]);
double vacc(double D, double nd);
void wavefunction(double v[], double E,int np, double wavestore[]); //齋藤が引数jを追加 (2016.10.13)
double func(double E, double v, int n);
double calcurrent(double a, double b, double v[], double delta);
double current(double v[]);
void potential(double v[], double vRTD[], double VRTD, double QW);
double fermia(int n, double t, int flag);
double calconfinedstate(double v[], int np, double E1, double E2);
double getconfinedstate(int n, double v[], int np);
void getconfinedstates(int n, double v[], int np, double En[]);
double cal4(double E, double v[], int np);