SETFILE = "set.dat";

count = 0 %csv保存用のカウンタ　手動で変える


filename = datestr(now, 'yyyy_mm_dd_') + num2str(count);

j = 0;
n = 0;
loop = 0; %齋藤が変数loop(ループを回す為の変数)を追加 (2016.10.13)
global RTD_Designs;
global layer;
const = Constant();

setmaterial(const);

				%/* ML=1e-9にすれば、z[nm]。ML=0.31e-9ならz[ML]	*/
global N;
N = [0, 0];
N(const.ALL) = sum([RTD_Designs.divnum]);
N(const.RTD) = RTD_Designs(const.LBAR).divnum + RTD_Designs(const.WELL).divnum + RTD_Designs(const.RBAR).divnum;	%/* RTD構造のみの層数*/

v = zeros(N(const.ALL)+1, 1);
vRTD = zeros(N(const.RTD)+1, 1);				%/* ポテンシャルの格納用			*/
QW=0*const.ELEC*1e4*5e12;

VRTD=0.2;

[v, vRTD] = potential(v, vRTD, VRTD, QW);
% FILE *fp;
% fp = fopen(filename_potential, "w");
% for(j=0; j<N[0]+1; j++){
% 	fprintf(fp, "%d\t%e\t%e\n", j, j*dx*1e9, v[j]);
% }
% fclose(fp);


% //double V = v[0]+Ef[0]-(v[N[0]-1]+Ef[layer]);
% n=60;								/* 計算したい準位の数	*/
% double En[n];						/* 準位の格納			*/
% double wavestore[n][DIV*N[0]];  //波動関数格納
% getconfinedstates(n,v,0,En);
% for(j=0;j<n;j++){
% 	wavefunction(v,En[j],0, wavestore[j]); //齋藤が引数wavestore[j]を追加 (2016.10.13)
% }

% fp = fopen(filename_wave, "w");
% for(j=0; j<DIV*N[0]; j++){
% 	fprintf(fp, "%g\t", (double)j*dx*1e9/DIV);
% 	for(loop=0; loop<n; loop++)
% 		fprintf(fp, "%e\t", wavestore[loop][j]);
% 	fprintf(fp, "\n");
% }
% fclose(fp);
% //Beep( 440, 1200 );
% return 0;


function x = mt(n)
    global RTD_Designs layer;
	for x = 1 : layer
		if  n <= RTD_Designs(x).NX
            print(x)
			return 
		end
    end
end

function n = nRTD(np)
	if np == 0
		n = 0;
	elseif np == 1
		n = RTD_Designs(1).NX;
	else
		printf("想定していない値が引き数npに使われています。\n");
	end
end


function [v, vRTD] = potential(v, vRTD, VRTD, QW)

	vRTD = selfpotential(QW, VRTD, vRTD);
	
	for n = 0 : N(const.RTD)+1
		v(n + RTD_Designs(1).NX) = vRTD(n);
	end

	DL = RTD_Designs(const.LBAR).die * (vRTD(1) - vRTD(0))/const.dx;
	DR = RTD_Designs(const.RBAR).die * (vRTD(N(const.RTD)) - vRTD(N(const.RTD) - 1))/const.dx;

	for n = const.LBAR-1 : -1 : 0
		DL = calVRL(-1, n, DL, v);
	end

	for n = const.RBAR+1 : layer
		DR = calVRL(1, n, DR, v);
	end

	setpotential(v, 0);
end

function vRTD = selfpotential(Q, VRTD, vRTD)
    const = Constant();
	np=const.RTD;										%/*npはRTD構造か、全体かを表している。np=1はRTD構造。=0は素子全体 */
	aa=-1;										%/*計算終了の判断のために使用*/
	vRTD = potential0(Q, VRTD, vRTD);			%/*近似式 初期値として使用*/
	vRTD = setpotential(vRTD, np);				%/*階段近似適用*/
	A = ones(N(np), 1) * 0+0i;
	B = ones(N(np), 1) * 1+0i;
	k = ones(N(np), 1);							%/*Aは前進波、Bは後進波の係数。kは波数 */

	while aa == -1								% /*自己無撞着計算* 計算が収束しない可能性あり。*/
		E = confinedstate(1, vRTD, np);
		[A, B] = cal2(E, vRTD, A, B, np);			%/*エネルギーEにおける波動関数の係数AとBを計算*/
		k = setk(E, vRTD, k, np);					%/*波数計算*/
		[aa, vRTD] = makewave2(Q, E, VRTD, vRTD, k, A, B, np);	%/*電荷Qを考慮したポテンシャル計算 出力は、前のポテンシャルと計算後のポテンシャルの差が1e-6以下なら-1　それ以外が1*/
	end
end

function vRTD = potential0 (Q, VRTD, vRTD)
    global RTD_Designs;
    const = Constant();
	nrtd = RTD_Designs(const.LBAR-1).NX;
	
	DL = (VRTD - (RTD_Designs(const.WELL).d/2/RTD_Designs(const.WELL).die + RTD_Designs(const.RBAR).d/RTD_Designs(const.RBAR).die)*Q) / (RTD_Designs(const.LBAR).d/RTD_Designs(const.LBAR).die+RTD_Designs(const.WELL).d/RTD_Designs(const.WELL).die + RTD_Designs(const.RBAR).d/RTD_Designs(const.RBAR).die);
	DR = DL + Q;

	for n = 0 : nrtd+1	
		x = n*const.dx;
		switch mt(n + nrtd) 
			case 0
				vRTD(n) = VRTD;
			case 1
				vRTD(n) = VRTD;
			case 2
				vRTD(n) = -DL*x/RTD_Designs(const.LBAR).die + VRTD;
			case 3
				vRTD(n) = -Q * ( pow(x-RTD_Designs(const.LBAR).d, 2) + RTD_Designs(const.WELL).d*RTD_Designs(const.WELL).d*(cos(2*pi*(x-RTD_Designs(const.LBAR).d)/RTD_Designs(const.WELL).d) - 1)/2/(pi*pi) ) / 2/RTD_Designs(const.WELL).die/RTD_Designs(const.WELL).d - DL*(x - RTD_Designs(const.LBAR).d)/RTD_Designs(const.WELL).die - DL*RTD_Designs(const.LBAR).d/RTD_Designs(const.LBAR).die + VRTD;
			case 4
				vRTD(n) = -DR*( x - RTD_Designs(const.LBAR).d - RTD_Designs(const.WELL).d - RTD_Designs(const.RBAR).d ) / RTD_Designs(const.RBAR).die;
			case 5
				vRTD(n) = 0;
			case 6
				vRTD(n) = 0;
			case 7
				vRTD(n) = 0;
		end
	end
end

function v = setpotential (v, np)							%/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
	nrtd = nRTD(np);
	vr = zeros(N(np)+1);
	vl = zeros(N(np)+1);									%/* ポテンシャルの左側からの極限vrと右側からの極限			*/
	for n = 0 : N(np)+1
		vl(n) = v(n) + RTD_Designs(mt(n+nrtd)).bar;
		vr(n) = v(n) + RTD_Designs(mt(n+nrtd+1)).bar;
	end

	for n = 0 : N(np)
		v(n) = ( vr(n) + vl(n+1) )/2;						%/* 階段近似適用	*/
	end
end

function E = confinedstate(num, v, np)						%/* 引数numは、量子数nのこと。nは別のところで使われているのでnumにした。 */
															%/*バグあり。恐らく、左右の一番端のエネルギーが準位として検出されてしまう。*/
	emax = 5;												%/* whileの無限ループを避けるため */
	min = min(v);												%/* 分割されたポテンシャルの最小値を求める */
	E=min;

	for i = 1 : num+1
		dE = 1e-3;											%/* 刻み幅(粗い)の設定                                              */
		T = 0;
		j = 0;
		S = cal(E, v, np, 0);
		while E < emax                                		%/* 透過率の減少区間の検索                                          */
															%/* 以下の準位を求めるプログラムは                                  */
			E = E + dE;										%/* エネルギーが上昇していくとともに、透過率が上昇していくと仮定し、*/
			T = cal(E, v, np, 0);							%/* その最大値をとる箇所が準位であるとしている                      */
			if S > T										%/* よって、透過率が減少していく区間は不要                          */
				S = T;
			else
				break;										%/* 透過率の減少区間の終了 */
			end
		end

		S = 0;
		T = 0;												%/* S、Tの初期化                          */
		while E < emax										%/* 基底準位の探索 刻み幅：粗い           */
			E = E + dE;										%/* エネルギーの設定                      */
			T = cal(E, v, np, 0);							%/* 透過率の計算                          */
			if S < T										%/* 増減の判断　増加関数ならばwhile文続行 */
				S = T;
			else
				E = E - 2*dE;								%/* エネルギーを少し戻して、刻み幅の細かい方へ移行 */
				dE = dE/10;
				%j++;
				S = 0;
				if dE < 1e-10
					break;
				end
			end
		end

		dE = dE/10;											%/* 刻み幅(細かい)の設定                           */
		S=0;												%/* Sの初期化                                      */
		while E < emax										%/* 基底準位の探索 刻み幅：細かい                  */
			E = E + dE;										%/* エネルギーの設定                               */
			T=cal(E, v, np, 0);								%/* 透過率の計算                                   */
			if S < T
				S = T;										%/* 増減の判断　増加関数ならばwhile文続行          */
			else
				if abs(E-dE-v(0)) < 1e-9 || abs(E-dE-v(N(np)-1)) < 1e-9
															%/*cal関数の影響で、E=v[0]のところが準位に見えてしまうので*/
					E = E + 10*dE;							%/*それを避けるために導入*/
					i = i - 1;								%/*従って、準位がv[0]±1e-10[eV]の近辺に存在した場合、*/
					break;									%/*検知できない可能性あり*/
				end
				if i == num
					E = E - dE;
					return;
				end
				E = E + 10*dE;
				break;										%/* 透過率の増加区間の終了 その一個前の計算結果を  */
			end
		end													%/* 簡略化したい */
	end
	E = E - dE;												%/* num番目の準位を返す                           */
end

function t = cal(E, v, np, f) 
    MSTAR = const.MSTAR;
    ELEC = const.ELEC;
    HBAR = const.HBAR;
    dx = const.dx;
	mass = [RTD_Designs.mass]

    nrtd = nRTD(np);

    if(f == 1 && (E - v(N(np)) <= 0 || E - v(1) <= 0))
        t = 0;
        return;
    end

    temp = eye(2);
    dummy = eye(2);
    trans = eye(2);

    for n = 1:N(np)-1
        if ((E-v(n))*(E-v(n+1)) ~= 0)
            m = sqrt(mass(mt(n+nrtd+2))/mass(mt(n+nrtd+1)));
            kk = sqrt((E-v(n))/(E-v(n+1)));

            pp = kk * m;               
            pm = -kk * m;           

            kn1 = sqrt(2*mass(mt(n+nrtd+2))*(E-v(n+1))*MSTAR*ELEC);  

            zp = exp(1i * kn1 * dx/HBAR);                        
            zm = exp(-1i * kn1 * dx/HBAR);                      

            temp(1,1) = 0.5 * ((1 + pp) * zp); 
            temp(1,2) = 0.5 * ((1 + pm) * zp); 
            temp(2,1) = 0.5 * ((1 + pm) * zm); 
            temp(2,2) = 0.5 * ((1 + pp) * zm); 

        elseif (E - v(n+1) == 0 && E - v(n) ~= 0)
            kn = sqrt(2*mass(mt(n+nrtd))*(E-v(n))*MSTAR*ELEC/(HBAR*HBAR));

            m = mass(mt(n+nrtd+1))/mass(mt(n+nrtd));

            zp = 1i * kn * m;
            zm = -1i * kn * m;

            temp(1,1) = zp;
            temp(1,2) = zm;

            zp = 1i * kn * m * dx;
            zm = -1i * kn * m * dx;

            temp(2,1) = 1 + zp;
            temp(2,2) = 1 + zm;

        elseif (E - v(n) == 0 && E - v(n+1) ~= 0)
            kn1 = sqrt(2*mass(mt(n+nrtd+1))*(E-v(n+1))*MSTAR*ELEC);
            zp = exp(1i * kn1 * dx/HBAR);
            zm = exp(-1i * kn1 * dx/HBAR);

            m = mass(mt(n+nrtd+1))/mass(mt(n+nrtd));
            kn1 = sqrt(HBAR*HBAR/(2*mass(mt(n+nrtd+1))*(E-v(n+1))*MSTAR*ELEC));
            pp = 1i * kn1 * m;
            pm = -1i * kn1 * m;

            temp(1,1) = 0.5 * (pm * zp); 
            temp(1,2) = 0.5 * zp;                     
            temp(2,1) = 0.5 * (pp * zm); 
            temp(2,2) = 0.5 * zm;

        elseif (E - v(n) == 0 && E - v(n+1) == 0)
            temp(1,1) = mass(mt(n+nrtd+1))/mass(mt(n+nrtd));  
            temp(1,2) = 0;                                    
            temp(2,1) = mass(mt(n+nrtd+1))*dx/mass(mt(n+nrtd));
            temp(2,2) = 1;
        end
        trans = temp * dummy;
        dummy = trans;
    end

    if (f == 1)
        if (E - v(N(np)) > 0 && E - v(1) > 0)
            t = sqrt((mass(layer)/mass(1)) * ((E - v(1))/(E - v(N(np)))));
        else
            t = 0;
        end
    else
        t = 1;
    end
    t = t * 1/(abs(trans(2,2))^2);
end

function [A, B] = cal2(E, v, A, B, np)

	nrtd=nRTD(np);

	temp = eye(2,2);
	FN = zeros(2,1);
	FNN = zeros(2,1);

	for n=1:(N(np)-1)
		FN(1,1)=A(n);
		FN(2,1)=B(n);

		if((E-v(n))*(E-v(n+1))~=0)
			m=sqrt(RTD_Designs(n+nrtd+2).mass / RTD_Designs(n+nrtd+1).mass);
			kk=sqrt(E-v(n))/sqrt(E-v(n+1));

			pp=kk*m;
			pm=kk*(-m);

			kn1=sqrt(2*RTD_Designs(n+nrtd+2).mass*(E-v(n+1))*const.MSTAR*const.ELEC);

			zp=exp(1i*kn1*const.dx/const.HBAR);
			zm=exp(-1i*kn1*const.dx/const.HBAR);

			temp(1,1)=(pp+1)*zp;
			temp(1,2)=(pm+1)*zp;
			temp(2,1)=(pm+1)*zm;
			temp(2,2)=(pp+1)*zm;

			temp=temp*0.5;
		elseif(E-v(n+1)==0 && E-v(n)~=0)
			kn=sqrt(2*RTD_Designs(n+nrtd).mass*(E-v(n))*const.MSTAR*const.ELEC/(const.HBAR^2));

			m=RTD_Designs(n+nrtd+1).mass/RTD_Designs(n+nrtd).mass;

			zp=1i*kn*m;
			zm=-1i*kn*m;

			temp(1,1)=zp;
			temp(1,2)=zm;

			zp=1i*kn*m*const.dx;
			zm=-1i*kn*m*const.dx;

			temp(2,1)=zp+1;
			temp(2,2)=zm+1;
		elseif(E-v(n)==0 && E-v(n+1)~=0)
			kn1=sqrt(2*RTD_Designs(n+nrtd+1).mass*(E-v(n+1))*const.MSTAR*const.ELEC);
			zp=exp(1i*kn1*const.dx/const.HBAR);
			zm=exp(-1i*kn1*const.dx/const.HBAR);

			m=RTD_Designs(n+nrtd+1).mass/RTD_Designs(n+nrtd).mass;
			kn1=sqrt(const.HBAR^2/(2*RTD_Designs(n+nrtd+1).mass*(E-v(n+1))*const.MSTAR*const.ELEC));
			pp=1i*kn1*m;
			pm=-1i*kn1*m;

			temp(1,1)=pm*zp;
			temp(1,2)=zp;
			temp(2,1)=pp*zm;
			temp(2,2)=zm;

			temp=temp*0.5;
		elseif(E-v(n)==0 && E-v(n+1)==0)
			temp(1,1)=RTD_Designs(n+nrtd+1).mass/RTD_Designs(n+nrtd);
			temp(1,2)=0;
			temp(2,1)=(RTD_Designs(n+nrtd+1).mass*const.dx)/RTD_Designs(n+nrtd);
			temp(2,2)=1;
		end

		FNN=temp*FN;

		A(n+1)=FNN(1,1);
		B(n+1)=FNN(2,1);

		FN=FNN;
	end
end

function k = setk(E, v, k, np)		%/* ポテンシャルから波数kを計算して格納*/
	nrtd=nRTD(np);
	for n = 0 : N(np)
		k(n) = sqrt( 2*RTD_Designs(mt(n+nrtd+1)).mass*(E-v(n))*const.MSTAR*const.ELEC) / const.HBAR;
	end
end

function [aa, v] = makewave2(Q, E, VRTD, v, k, A, B, np)
	max=0;													%// 規格化用
	qmax=0;
	p = zeros(const.DIV*N(np), 1);
	q = zeros(const.DIV*N(np), 1);								%// pがx座標、qがy座標 qは波動関数の絶対値の2乗を出力予定
	qtemp = zeros(N(np)+1, 1)
	vnew = zeros(N(np), 1);
	%gsl_complex tempA,tempB;
	j=0;	%//	p[0]=-ML;
	for n = 0 : N(np)
		xn = (n+1)*const.dx;
		for a = 0 : const.DIV
			p(j) = j*const.dx/const.DIV;
			tempA = A(n) * exp(k(n) * ( (p(j)-xn) * 1i));	%/* tempA= A[n] × exp( ik[n](x-x[n])) */
			tempB = B(n) * exp(k(n) * ( (xn-p(j)) * 1i));	%/* tempB= B[n] × exp(-ik[n](x-x[n])) */
			q(j) = abs(tempA + tempB)^2;					%/* 波動関数の絶対値の2乗を計算        */
			qmax = qmax + q(j)/const.DIV;
			j = j+1;
		end
	end

	j=0;
	for n = 1 : N(np)
		qtemp(n) = qtemp(n) - Q*q(const.DIV*n)/qmax;
	end

	qtemp(0) = 0;
	qtemp(N(np)) = 0;
	vnew = potential0(0, VRTD, vnew);
	vnew = calpotential(vnew, qtemp, np);
	opv = zeros(N(np)+1, 1);
	for n = 0 : N(np)+1
		opv(n) = vnew(n);
	end

	vnew = setpotential(vnew,np);

	for n = 1 : N(np)-1
		if abs(v(n) - vnew(n)) > max
			max = abs(v(n) - vnew(n));
		end
	end
	for n = 0 : N(np)
		v(n) = vnew(n);
	end

	if max < 1e-6			%/* 重要　自己無撞着計算の収束判定*/
		for n = 0 : N(np) 
			v(n) = opv(n);
		end
		v(N(np)) = 0;
		aa = 1;
	else
		aa = -1;
	end
end

function vnew = calpotential(vnew, q, np)
    nrtd = nRTD(np);
    div = N(np) - 1;

    V = zeros(div, 1);
    X = zeros(div, 1);
    S = zeros(div, div);
    
    for n = 1:div
        temp = q(n + 1) * dx;
        if n + nrtd + 1 == NX(3)
            dietemp = (die(LBAR) + die(WELL)) / 2;
        elseif n + nrtd + 1 == NX(4)
            dietemp = (die(WELL) + die(RBAR)) / 2;
        else
            dietemp = die(mt(n + nrtd + 1));
        end

        if n == 1
            S(1, 1) = -2 * dietemp;
            S(2, 1) = 1 * dietemp;
            V(1) = temp;
        elseif n == div
            S(n, n - 1) = 1 * dietemp;
            S(n, n) = -2 * dietemp;
            V(n) = temp;
        else
            S(n, n - 1) = 1 * dietemp;
            S(n, n) = -2 * dietemp;
            S(n + 1, n) = 1 * dietemp;
            V(n) = temp;
        end
    end

    [L,U,P] = lu(S);  % LU decomposition
    X = P * (U \ (L \ V));  % Solve linear system

    for n = 1:div
        vnew(n + 1) = vnew(n + 1) + X(n);
    end
end

function D = calVRL(direction, n, D, v)
	
	la = RTD_Designs(n).NX;
	
	if n-1 < 0
		sm = 0;
	else
		sm = RTD_Designs(n-1).NX;
	end

	if RTD_Designs(n).smt == 7 || RTD_Designs(n).smt == 5 || RTD_Designs(n).smt == 23 || RTD_Designs(n).smt == 24 || RTD_Designs(n).smt == 22 %//i-Si smt[n]==24を追加 07/28大野 8/19ていがn=22を追加
		E = D/RTD_Designs(n).die;
		if direction == -1
			for a = la : -1 : sm
				v(a-1) = v(a) - E*const.dx;
			end
		else
			for a = sm : la
				v(a+1) = v(a) + E*const.dx;
			end
		end
		return;
	end

	if RTD_Designs(n).cond == 0 %//金属
		if direction == -1
			for a = la : -1 : sm
				v(a-1) = v(n);
			end
		else
			for a = sm : la
				v(a+1) = v(RTD_Designs(n-1).NX);
			end
		end
		D = 0;
		return;
	end
	
	if RTD_Designs(n).cond == 1		%//n-SiSub.
		if direction < 0
			j=la;
		else
			j=sm;
		end
		if  (direction < 0 && D >= 0) || (direction > 0 && D <= 0)	%//空乏層
			A = const.ELEC*RTD_Designs(n).Q*1e6/2/(RTD_Designs(n).die);
			w = abs(D/const.ELEC/(RTD_Designs(n).Q*1e6));
			wmax = sqrt(4*RTD_Designs(n).die*const.KB*const.TEMP*log(RTD_Designs(n).Q*1e6/const.NI)/(const.ELEC*const.ELEC*RTD_Designs(n).Q*1e6));
			x=0;
			%//w<d[n] || wmax<d[n])
			if abs(w) > wmax
				w=wmax;
			end

			if direction < 0
				for a = la : -1 : sm
					if x < w
						v(a) = A*x*(x-2*w) + v(j);
					else
						v(a) = -A*w*w + v(j);
					end
					x = x + const.dx;
				end
			else
				for a = sm : la
					if x < w
						v(a) = A*x*(x-2*w) + v(j);
					else
						v(a) = -A*w*w + v(j);
					end
					x = x + const.dx;
				end
			end
			% else		%//エラー処理
			% 	printf("\n%d層目:空乏層幅が膜厚以上です。w=%e\n膜厚を増やして下さい。\n空乏層幅Max=%e",n,w,wmax);
			% 	int wML=wmax/ML;
			% 	printf("=%6.4gML<%dML \n現在の膜厚d[%d]=%e\n",wmax/ML,wML+1,n,d[n]);
			% 	exit(1);
			% }
			D = 0;
			return;
		else %//蓄積層
			A = sqrt(2*const.KB*const.TEMP*RTD_Designs(n).Q*1e6/RTD_Designs(n).die);
			q = const.ELEC/(const.KB*const.TEMP);
			fai = -vacc(D, RTD_Designs(n).Q);
			E = A*sqrt(exp(-q*fai)+q*fai-1);
			tempV = 0;
			if direction > 0
				for a = sm : la
					v(a) = tempV + v(j);
					fai = fai + E*const.dx;
					tempV = tempV + direction*E*const.dx;
					E = A*sqrt(exp(-q*fai)+q*fai-1);
				end
			else
				for a = la : -1 : sm
					v(a) = tempV + v(j);
					fai = fai + E*dx;
					tempV = tempV + direction*E*const.dx;
					E = A*sqrt(exp(-q*fai)+q*fai-1);
				end
			end
		end
	end

	if RTD_Designs(n).cond == 2  %//p-Si Sub.
		if direction < 0
			j=la;
		else
			j=sm;
		end
		if  (direction < 0 && D <= 0) || (direction > 0 && D >= 0)  %//空乏層
			A = -const.ELEC*RTD_Designs(n).Q*1e6/2/RTD_Designs(n).die;
			w = abs(D/const.ELEC/(RTD_Designs(n).Q*1e6));
			wmax = sqrt(4*RTD_Designs(n).die*const.KB*const.TEMP*log(RTD_Designs(n).Q*1e6/const.NI)/(const.ELEC*const.ELEC*RTD_Designs(n).Q*1e6));
			x=0;
			if w < RTD_Designs(n).d || wmax < RTD_Designs(n).d
				if abs(w) > wmax
					w=wmax;
				end
				if direction < 0
					for a = la : -1 : sm
						if x < w
							v(a) = A*x*(x-2*w) + v(j);
						else
							v(a) = -A*w*w + v(j);
						end
						x= x + dx;
					end
				else
					for a = sm : la
						if x < w
							v(a) = A*x*(x-2*w) + v(j);
						else
							v(a) = -A*w*w + v(j);
						end
						x= x + dx;
					end
				end
			else %//エラー処理
				% printf("\n%d層目:空乏層幅が膜厚以上です。w=%e\n膜厚を増やして下さい。\n空乏層幅Max=%e",n,w,wmax);
				% int wML=wmax/ML;
				% printf("=%6.4gML<%dML\n",wmax/ML,wML+1);
				exit(1);
			end
			D = 0;
			return
		else %//蓄積層
			A = sqrt(2*const.KB*cosnt.TENP*RTD_Designs(n).Q*1e6/RTD_Designs(n).die);
			q = const.ELEC/(const.KB*const.TEMP);
			fai = -vacc(D, RTD_Designs(n).Q);
			E = A*sqrt(exp(q*fai)-q*fai-1);
			tempV = 0;
			if direction > 0
				for a = sm : la
					v(a) = tempV + v(j);
					fai = fai + E*const.dx;
					tempV = tempV - direction*E*const.dx;
					E = A*sqrt(exp(q*fai)-q*fai-1);
				end
			else
				for a = la : -1 : sm
					v(a) = tempV + v(j);
					fai = fai + E*const.dx;
					tempV = tempV + direction*E*const.dx;
					E = A*sqrt(exp(-q*fai)+q*fai-1);
				end
			end
		end
	end
	D = 0;
end

function Vacc =  vacc(D, nd)
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


%void makewave(double E, double v[], gsl_complex k[], gsl_complex A[], gsl_complex B[], int np, double wavestore[]) //齋藤が引数wavestoreを追加 (2016.10.13)
%{
%	int i,j,n;
%	double max=0;                     // 規格化用
%	double sum=0;
%	double xn,px,py, pyy[N[np]*DIV]; //齋藤が配列追加(pyy, 2016.10.13)
%	gsl_complex tempA,tempB;
%	j=max=0;//p[0]=-ML;q[0]=E;
%	for(n=0;n<N[np];n++)
%	{
%		xn=(n+1)*dx;
%		for(i=0; i<DIV; i++)
%		{
%			px=(double)j*dx/DIV;
%			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
%			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
%			py=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
%			if(mt(n)==LBAR ||mt(n)==WELL||mt(n)==RBAR)
%				sum+=py*dx*1e9/DIV;
%			j++;
%		}
%	}j=0;//2015/5/20追加　思考停止 とりあえず波動関数の面積を求める(↑のmaxに格納)。波動関数の2乗の全区間で積分すると1になるように規格化。
%	 //量子井戸の場合、膜厚が厚すぎるとexp(±ikz)の掛け算でオーバーフロー(？)を起こし、減衰項ではなく増幅項が支配的になるので注意。
%	for(n=0;n<N[np];n++)
%	{

%	//ここから先，出力の仕方を大幅に変更(旧QCL.cとほぼ同じにした) by齋藤 (2016.10.13)

%		xn=(n+1)*dx;
%		for(i=0; i<DIV; i++)
%		{
%			px=(double)j*dx/DIV;
%			tempA=gsl_complex_mul(A[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],px-xn)));
%			tempB=gsl_complex_mul(B[n],gsl_complex_exp(gsl_complex_mul_imag(k[n],xn-px)));
%			pyy[j]=gsl_complex_abs2(gsl_complex_add(tempA,tempB));
%			if(pyy[j]>max)
%				max=pyy[j];
%			j++;
%		}
%	}
%	for(i=0; i<N[np]*DIV; i++)
%		wavestore[i] = E+pyy[i]/max;
%	return;
%}




%void wavefunction(double v[], double E,int np, double wavestore[]) //齋藤が引数wavestoreを追加 (2016.10.13)
%{
%	gsl_complex A[N[np]],B[N[np]],k[N[np]];
%	cal3(E,v,A,B,np);                               /* 波動関数の計算                  */
%	setk(E,v,k,np);                                 /* 波数の計算                      */		
%	makewave(E,v,k,A,B,np, wavestore);                         /* 波動関数の計算と表示            */
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

function setmaterial(const) % m=-100だと、出力されない。実行はされる。引き数mは、物性値を変更したい場合に使用。
	%base = 0;		%/* 伝導帯エネルギーの基準。大抵はSiの電子親和力	*/
	%temp=TEMP;			%/* 温度						*/
	%/* 各物質における物性値の設定。name:名前、die：誘電率、bar：バンド不連続、mass:有効質量、valley:谷の数。valley:基本1。cond:導電性。0が金属、1がn-Si、2がp-Si、3が絶縁物。*/
	
	%/* 物質を追加したい場合は、この行の上をコピーして貼り付け。iを変更するのとmtconstの配列に収まるように注意。*/

	design_data = readtable('set.csv');
	%if(m!=flag)
	%	printf("層番号\t材料番号\t物質名\tML数\t層厚[nm]\t比誘電率\t有効質量\t障壁の高さ\t谷\t分割数\tNX\tEF\t仕事関数\t電荷量\n");
	RTD_Designs = Materials.empty(0, length(design_data.materialName));

	layer = length(design_data.materialName);							%層数
    
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