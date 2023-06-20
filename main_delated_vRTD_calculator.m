SETFILE = "set.dat";

count = 0 %csv保存用のカウンタ　手動で変える


filename = datestr(now, 'yyyy_mm_dd_') + num2str(count);

j = 0;
n = 0;
loop = 0; %齋藤が変数loop(ループを回す為の変数)を追加 (2016.10.13)
global RTD_Designs;
global layer;
global const;
const = Constant();

setmaterial();

				%/* ML=1e-9にすれば、z[nm]。ML=0.31e-9ならz[ML]	*/
global N;
N = RTD_Designs(layer).NX;


VRTD=1.35;

v = potential(VRTD);
zn = (1 : N)*const.dx*1e9;
plot(zn, v);
hold on

mass = zeros(1,N);
mass(1:RTD_Designs(1).NX) = RTD_Designs(1).mass*const.MSTAR;
for i = 2:layer
    mass(RTD_Designs(i-1).NX+1:RTD_Designs(i).NX) = RTD_Designs(i).mass*const.MSTAR; 
end
n = 20;								%/* 計算したい準位の数	*/
En = getconfinedstates(n, v, mass);
n = length(En)                      %見つかった準位の数

wavestore = zeros(n, const.DIV * N);
zn = (1 : const.DIV * N)*const.dx*1e9/const.DIV;
for j = 1 : n
 	wavestore(j,:) = wavefunction(v, En(j), const, mass, RTD_Designs, N); %//齋藤が引数wavestore[j]を追加 (2016.10.13)
    plot(zn, wavestore(j,:)+En(j));
end


function v = potential(VRTD)

	v = potential0(VRTD);			%/*近似式 初期値として使用*/
	v = setpotential(v);				%/*階段近似適用*/
	
end

function v = potential0 (VRTD)
    global RTD_Designs;
    global layer;
    global const;
    global N;

    D = VRTD / sum([RTD_Designs(2:layer-1).d] ./ [RTD_Designs(2:layer-1).die]);%全体の電束密度を計算
    % 先に全ての要素を0で初期化
    die = zeros(1,N);
    % 繰り返し処理でaの要素を更新
    v = zeros(1,N);
    die(1 : RTD_Designs(1).NX) = RTD_Designs(1).die;
    for i = 2:layer
        die(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).die;
    end
    die = die(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX);
    v(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX) = VRTD - cumsum(D./die .* const.dx); %1つ前のvに対して、-D/ε_n*dxを計算することで蓄電状態にない系のポテンシャルを計算
    
    A = sqrt(2*const.KB*const.TEMP*RTD_Designs(1).Q*1e6/RTD_Designs(1).die);
    q = const.ELEC / (const.KB * const.TEMP);
    fai = - vacc(D, RTD_Designs(1).Q);
    E = A * sqrt(exp(-q*fai) + q*fai - 1);
    tempV = 0;
    for i = RTD_Designs(1).NX : -1 :1
        v(i) = tempV + v(RTD_Designs(1).NX+1);
        fai = fai + E*const.dx;
        tempV= tempV + E*const.dx;
        E = A*sqrt(exp(-q*fai) + q*fai - 1);
    end
    
    A = const.ELEC*RTD_Designs(layer).Q*1e6/2/RTD_Designs(layer).die;
    w = abs(D / const.ELEC / ( RTD_Designs(layer).Q*1e6 ));
    wmax = sqrt(4 * RTD_Designs(layer).die * const.KB * const.TEMP * log( RTD_Designs(layer).Q*1e6/const.NI ) / (const.ELEC^2 * RTD_Designs(layer).Q * 1e6));
    x=0;
    if w > wmax
        w=wmax;
    end
    for i = RTD_Designs(layer-1).NX+1 : N
        if x < w
            v(i) = A * x * (x - 2*w) + v(RTD_Designs(layer-1).NX);
        else
            v(i) = -A * w * w + v(RTD_Designs(layer-1).NX);
        end
        x = x + const.dx;
    end
end
function v = setpotential (v)							%/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
    global N;
    global RTD_Designs;
    global layer
	vr = zeros(N+1);
	vl = zeros(N+1);									%/* ポテンシャルの左側からの極限vrと右側からの極限			*/
    bar = zeros(1,N);
    % 繰り返し処理でaの要素を更新

    bar(1 : RTD_Designs(1).NX) = RTD_Designs(1).bar;
    for i = 2:layer
        bar(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).bar;
    end
	for n = 1 : N
		vl(n+1) = v(n) + bar(n);
		vr(n) = v(n) + bar(n);
	end

	for n = 1 : N
		v(n) = ( vr(n) + vl(n+1) )/2;						%/* 階段近似適用	*/
	end
end
function Vacc =  vacc(D, nd)
    global const
    global RTD_Designs
	loop = 1;
	Vacc = 0.1;
	sub = 1;

	nd = nd*1e6*const.ELEC;
	q = const.ELEC/const.KB/const.TEMP;
	a = q*D*D/2/RTD_Designs(1).die/nd;
	if D == 0
		Vacc = 0;
		return;
	end
	while sub > 1e-6
		if Vacc > 20 && loop <= 5
			loop = loop + 1;
			Vacc = 0.1*loop;
		elseif loop > 5
			print("蓄積層のポテンシャル計算が収束しませんでした。\nドーピング濃度や極端に高い電圧をかけていないか等確認して下さい。");
			exit(1);
		end
		xk = Vacc;
		f1 = exp(q*xk) - q*xk - 1 - a;
		f2 = q*(exp(q*xk) - 1);
		Vacc=  xk - f1/f2;
		sub = abs(Vacc - xk);
	end
end

function En = getconfinedstates(n, v, mass)
    global N;
    global const;
    E1 = min(v);
    E2 = max(v);%const.EMAX;
    dE = 5e-4;
    E = E1+1e-9 :dE: E2-1e-9;
    t11 = zeros(1,length(N));
    t12 = zeros(1,length(N));
    t21 = zeros(1,length(N));
    t22 = zeros(1,length(N));
    t22s = zeros(1,length(E));
    T = zeros(1,length(E));
	for i = 1 : length(E)
		[T(i), t11, t12, t21, t22] = TransMatrix(v, E(i), const, mass);
        t22s(i) = t22(N);
    end    
    [p, Es] = findpeaks(log(abs(T)), E);
    if n < length(Es)
        En = Es(1:n);
    else
        En = Es;
    end
end

function [T, t11, t12, t21, t22] = TransMatrix(v, E, const, mass)
    N = length(v);
    kn = sqrt(mass.*(E-v).*const.ELEC) / const.HBAR;
    P = zeros(2);
    ex = zeros(2);
    %A = zeros(1, N);
    %B = zeros(1, N);
    %A(N) = AN;
    %B(N) = BN;
    dx = const.dx;
    %AB = [AN;BN];
    trans = diag([1 1]);
    for n = 2:N
        P(1,1) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
        P(1,2) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
        P(2,1) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
        P(2,2) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
        ex(1,1) = exp(1i * kn(n-1) * dx);
        ex(2,2) = exp(-1i * kn(n-1) * dx);
        trans = 0.5*P*ex*trans;
        t11(n) = trans(1,1);
        t12(n) = trans(1,2);
        t21(n) = trans(2,1);
        t22(n) = trans(2,2);
    end
    T = mass(n)*kn(1)/mass(1)/kn(n)./abs(t11(N)).^2;
end
    % double calconfinedstate(double v[], int np, double E1, double E2)

function wavestore = wavefunction(v, E, const, mass, RTD_Designs, N) %//齋藤が引数wavestoreを追加 (2016.10.13)
	[T, t11, t12, t21, t22] = TransMatrix(v, E, const, mass);                               %/* 波動関数の計算                  */
    A = zeros(1,N);
    B = zeros(1,N);
    if v(1) <= E
        A(1) = 1;
        B(1) = -conj(t12(N))/conj(t11(N))*A(1);
    else
        A(1) = 1;
        B(1) = 0;
    end
    for n = 2:N
        A(n) = t11(n)*A(1) + t12(n)*B(1);
        B(n) = t21(n)*A(1) + t22(n)*B(1);
    end
    k = sqrt( 2*mass.*(E-v)*const.MSTAR*const.ELEC) / const.HBAR; %   /* 波数の計算                      */		
	wavestore = makewave(k, A, B, const, N);                                          %/* 波動関数の計算と表示            */
end


function wavestore = makewave(k, A, B, const, N) %//齋藤が引数wavestoreを追加 (2016.10.13)
	pyy = zeros(1,N*const.DIV); %//齋藤が配列追加(pyy, 2016.10.13)
    j = 1;
    %//2015/5/20追加　思考停止 とりあえず波動関数の面積を求める(↑のmaxに格納)。波動関数の2乗の全区間で積分すると1になるように規格化。
	%//量子井戸の場合、膜厚が厚すぎるとexp(±ikz)の掛け算でオーバーフロー(？)を起こし、減衰項ではなく増幅項が支配的になるので注意。
	for n = 1 : N
	%//ここから先，出力の仕方を大幅に変更(旧QCL.cとほぼ同じにした) by齋藤 (2016.10.13)
		xn = (n + 1)*const.dx;
		for i = 0 : const.DIV-1
			px = (j-1) * const.dx/const.DIV;
			tempA = A(n) * exp( 1i * k(n) * (px - xn) );
			tempB = B(n) * exp( 1i * k(n) * (xn - px) );
			pyy(j) = abs(tempA + tempB)^2;
			j =j+1;
        end
    end
	wavestore = pyy./max(pyy);
end
%void cal3(double E, double v[], gsl_complex A[], gsl_complex B[], int np)
%{
%	int i,n,nrtd;
%	double m;

%	nrtd=nRTD(np);	
%	gsl_complex kk,kn1,pp,pm,zp,zm;
%	gsl_matrix_complex *temp[N[np]];
%	for(n=0; n<N[np]; n++){
%		temp[n]=gsl_matrix_complex_alloc(2,2);
%	}
%	gsl_matrix_complex *dummy = gsl_matrix_complex_alloc(2,2);
%	gsl_matrix_complex_set_identity(dummy);
%	gsl_matrix_complex *trans = gsl_matrix_complex_alloc(2,2);
%	gsl_matrix_complex_set_identity(trans);
%	for(n=0; n<N[np]-1; n++)
%	{
%		if((E-v[n])*(E-v[n+1])!=0){
%			m=sqrt(mass[mt(n+nrtd+2)]/mass[mt(n+nrtd+1)]);
%			kk=gsl_complex_div(gsl_complex_sqrt_real(E-v[n]),gsl_complex_sqrt_real(E-v[n+1]));

%			pp=gsl_complex_mul_real(kk, m);                                                  //pm=+p=+ m[n+1]*k[n] / m[n]*k[n+1]を計算
%			pm=gsl_complex_mul_real(kk,-m);                                                  //pm=-p=- m[n+1]*k[n] / m[n]*k[n+1]を計算

%			kn1=gsl_complex_sqrt_real(2*mass[mt(n+nrtd+2)]*(E-v[n+1])*MSTAR*ELEC);                //エネルギーE、領域n+1における波数を計算 HBARで割っていないことに注意

%			zp=gsl_complex_exp(gsl_complex_mul_imag(kn1, dx/HBAR));                          //exp( ik[n+1]h)を計算
%			zm=gsl_complex_exp(gsl_complex_mul_imag(kn1,-dx/HBAR));                          //exp(-ik[n+1]h)を計算
			
%			gsl_matrix_complex_set(temp[n],0,0,gsl_complex_mul(gsl_complex_add_real(pp,1),zp)); //(1+p)exp( ik[n+1]h)を計算 行列temp(1,1)に格納
%			gsl_matrix_complex_set(temp[n],0,1,gsl_complex_mul(gsl_complex_add_real(pm,1),zp)); //(1-p)exp( ik[n+1]h)を計算 行列temp(1,2)に格納
%			gsl_matrix_complex_set(temp[n],1,0,gsl_complex_mul(gsl_complex_add_real(pm,1),zm)); //(1-p)exp(-ik[n+1]h)を計算 行列temp(2,1)に格納
%			gsl_matrix_complex_set(temp[n],1,1,gsl_complex_mul(gsl_complex_add_real(pp,1),zm)); //(1+p)exp(-ik[n+1]h)を計算 行列temp(2,2)に格納

%			gsl_matrix_complex_scale(temp[n], gsl_complex_rect(0.5,0));                         //全ての要素に0.5をかける
%		}
%		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],dummy,gsl_complex_rect(0,0),trans); //行列の掛け算

%		gsl_matrix_complex_memcpy(dummy, trans);
%	}
%	if(v[0]<E)
%	{
%		if(v[N[np]]<E)
%			i=0;
%		else
%			i=0;
%	}
%	else if(v[0]>E)
%		i=1;
%	switch(i)
%	{
%	case 0:	A[0]=gsl_complex_rect(1,0);
%			B[0]=gsl_complex_negative( gsl_complex_div( gsl_matrix_complex_get(trans,1,0), gsl_matrix_complex_get(trans,1,1) ) );
%	case 1: A[0]=gsl_complex_rect(0,0);
%			B[0]=gsl_complex_rect(1,0);
%	case 2: A[0]=gsl_complex_rect(1,0);
%			B[0]=gsl_complex_rect(1,0);
%	case 3: A[0]=gsl_complex_rect(1,0);
%			B[0]=gsl_complex_rect(0,0);
%	}
%	gsl_matrix_complex *FN  = gsl_matrix_complex_calloc(2,1);
%	gsl_matrix_complex *FNN = gsl_matrix_complex_calloc(2,1);
	
%	for(n=0; n<N[np]-1; n++)
%	{
%		gsl_matrix_complex_set(FN,0,0,A[n]);
%		gsl_matrix_complex_set(FN,1,0,B[n]);
%		gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,gsl_complex_rect(1,0),temp[n],FN,gsl_complex_rect(0,0),FNN); //行列の掛け算 FNN=temp*FN
%		A[n+1]=gsl_matrix_complex_get(FNN,0,0);
%		B[n+1]=gsl_matrix_complex_get(FNN,1,0);
%	}
%	for(n=0; n<N[np]; n++){
%		gsl_matrix_complex_free(temp[n]);
%	}
%	gsl_matrix_complex_free(dummy);
%	gsl_matrix_complex_free(trans);
%}








%double func(double E, double v, int n)
%{
%	if(fabs(v-E)<1e-15)
%		return 0;
%	double t;
%	t=1.0+exp(ELEC*(v+Ef[n]-E)/KB/temp);
%	return log(t);
%}



%double calconfinedstate(double v[], int np, double E1, double E2)
%{
%	/*	E1からE2のエネルギー範囲において、最も低い準位を探す	*/
%	/*	ない場合には、EMAXを返す。								*/
%	/*	変数	v	ポテンシャル								*/
%	/*			np	0なら空間全体、1ならRTD部分のみ				*/
%	double E,T,S,dE;
%	E=E1;											/* 探索範囲の下限													*/
%	while(E<E2)
%	{											
%		dE=1e-3;									/* 刻み幅(粗い)の設定                                              */
%		T=0;
%		S=cal(E,v,np,0);
%		while(E<E2)									/* 透過率の減少区間の検索                                          */
%		{											/* 以下の準位を求めるプログラムは                                  */
%			E+=dE;									/* エネルギーが上昇していくとともに、透過率が上昇していくと仮定し、*/
%			T=cal(E,v,np,0);						/* その最大値をとる箇所が準位であるとしている                      */
%			if(S>T)									/* よって、透過率が減少していく区間は不要                          */
%				S=T;
%			else
%				break;								/* 透過率の減少区間の終了 */
%		}
%		S=T=0;										/* S、Tの初期化                          */
%		while(E<E2)									/* 基底準位の探索 刻み幅：粗い           */
%		{
%			E+=dE;									/* エネルギーの設定                      */
%			T=cal(E,v,np,0);						/* 透過率の計算                          */
%			if(S<T)									/* 増減の判断　増加関数ならばwhile文続行 */
%				S=T;
%			else
%			{
%				E-=2*dE;							/* エネルギーを少し戻して、刻み幅を細かくしていく */
%				dE/=10;
%				S=0;
%				if(dE<1e-10)						/*　*/
%					break;
%			}
%		}
%		dE/=10;										/* 刻み幅(細かい)の設定                           */
%		S=0;										/* Sの初期化                                      */
%		while(E<E2)									/* 基底準位の探索 刻み幅：細かい                  */
%		{
%			E+=dE;									/* エネルギーの設定                               */
%			T=cal(E,v,np,0);						/* 透過率の計算                                   */
%			if(S<T)
%				S=T;								/* 増減の判断　増加関数ならばwhile文続行          */
%			else
%			{
%				if(fabs(E-dE-v[0])<1e-9||fabs(E-dE-v[N[np]-1])<1e-9)
%				{									/*cal関数の影響で、E=v[0]のところが準位に見えてしまうので*/
%					E+=10*dE;						/*それを避けるために導入*/
%					break;							/*検知できない可能性あり*/
%				}
%				E-=dE;
%				return E;
%			}
%		}
%	}
%	E2=EMAX;
%	return E2;
%}

%double getconfinedstate(int n, double v[], int np)
%{											/* 第n準位を出力												*/
%	int i;
%	double E,E1,E2;
%	E1=E2=EMAX;
%	for(i=0;i<N[np]+1;i++)										
%	{
%		if(E1>v[i])							/* エネルギーの探索範囲の決定									*/
%			E1=v[i];						/* 分割されたポテンシャルの最小値を求める						*/
%	}
%	for(i=0;i<n;i++)
%	{
%		E=calconfinedstate(v,np,E1,E2);
%		E1=E+1e-9;
%	}
%	return E;
%}

%void getconfinedstates(int n, double v[], int np, double En[])
%{											/* 複数の準位を取得し、配列に出力								*/
%	int i;
%	double E1,E2;
%	E1=E2=EMAX;
%	for(i=0;i<N[np]+1;i++)										
%	{
%		if(E1>v[i])							/* エネルギーの探索範囲の決定									*/
%			E1=v[i];						/* 分割されたポテンシャルの最小値を求める						*/
%	}
%	for(i=0;i<n;i++)
%	{
%		En[i]=calconfinedstate(v,np,E1,E2);
%		E1=En[i]+1e-9;
%	}
%}

function setmaterial() % m=-100だと、出力されない。実行はされる。引き数mは、物性値を変更したい場合に使用。
    global const;
    global layer;
    global RTD_Designs;
	%base = 0;		%/* 伝導帯エネルギーの基準。大抵はSiの電子親和力	*/
	%temp=TEMP;			%/* 温度						*/
	%/* 各物質における物性値の設定。name:名前、die：誘電率、bar：バンド不連続、mass:有効質量、valley:谷の数。valley:基本1。cond:導電性。0が金属、1がn-Si、2がp-Si、3が絶縁物。*/
	
	%/* 物質を追加したい場合は、この行の上をコピーして貼り付け。iを変更するのとmtconstの配列に収まるように注意。*/

	%if(m!=flag)
	%	printf("層番号\t材料番号\t物質名\tML数\t層厚[nm]\t比誘電率\t有効質量\t障壁の高さ\t谷\t分割数\tNX\tEF\t仕事関数\t電荷量\n");
	design_data = readtable('set.csv');
	layer = length(design_data.materialName);							%層数
	RTD_Designs = Materials.empty(0, length(design_data.materialName));
	NX = cumsum(design_data.ML*const.DX);
	for i = 1:layer
		RTD_Designs(i) = Materials(design_data.materialName(i), design_data.ML(i), design_data.Q(i), NX(i));
    end
	% if(m!=flag)
	% 	printf("%d\t%d\t%s\t%d\t%.4g\t%.3g\t%.2g\t%.4g\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\n", i, smt[i], name[i], j, j*ML*1e9, die[i]/DIEELECSTAR, mass[i], bar[i], valley[i], divnum[i], NX[i],Ef[i],Work[i],Q[i]);//出力。
	% }

	% if(m!=flag){
	% 	printf("\n1ML[nm]\t分割数\t積分範囲\tΔE\tピーク時ΔE\tN\tNRTD\tTEMP\tII*temp\tbar[%d]\tbar[%d]\n%8.3g\t%6d\t%8.2g\t%8.2g\t%8.2g\t%5d\t%4d\t%4.3g\t%8.2g",layer,LBAR-2,ML*1e9,DX,VI,DELTA,DELTAE,N[0],N[1],temp,II*temp);
	% 	printf("なし\tなし\n");
	% }
end