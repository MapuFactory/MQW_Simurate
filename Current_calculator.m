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

mass = zeros(1,N);
mass(1:RTD_Designs(1).NX) = RTD_Designs(1).mass*const.MSTAR;
for i = 2:layer
    mass(RTD_Designs(i-1).NX+1:RTD_Designs(i).NX) = RTD_Designs(i).mass*const.MSTAR; 
end

dv = 0.0025;
v = 0:dv:2;
for index = 1:length(v)
    pot = potential(v(index));
    zn = (1 : N)*const.dx*1e9;
    
    E1 = min(pot);
    E2 = max(pot);
    dE = 5e-4;
    E = E1+1e-9 :dE: E2-1e-9;
    i_element = zeros(1, length(E));
    T = zeros(1, length(E));
    for n = 1 : length(E)
        [i_element(n), T(n)] = cal_Current(pot, E(n), const, mass, RTD_Designs, RTD_Designs(const.InjectLayer).NX, N);
    end
    J(index) = const.II * const.TEMP * trapz(E, i_element)/1e4;
end
plot(v, J);

function v = potential(V_all)

	v = potential0(V_all);			%/*近似式 初期値として使用*/
	v = setpotential(v);				%/*階段近似適用*/
    v = v-min(v);
	
end

function v = potential0 (V_all)
    global RTD_Designs;
    global layer;
    global const;
    global N;

    D = V_all / sum([RTD_Designs.d] ./ [RTD_Designs.die]);%全体の電束密度を計算
    % 先に全ての要素を0で初期化
    die = zeros(1,N);
    d = zeros(1,N);
    % 繰り返し処理でaの要素を更新
    v = zeros(1,N);
    die(1 : RTD_Designs(1).NX) = RTD_Designs(1).die;
    d(1 : RTD_Designs(1).NX) = RTD_Designs(1).d;
    for i = 2:layer
        die(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).die;
        d(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).d;
    end
    E_point = D./die*const.dx;
    v(1) = V_all;
    for n = 2:N
        v(n) = v(n-1) - E_point(n);
    end
    if strcmp(RTD_Designs(1).name, 'n-Si.Sub.')
        die = die(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX);
        v(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX) = V_all - cumsum(D./die .* const.dx); %1つ前のvに対して、-D/ε_n*dxを計算することで蓄電状態にない系のポテンシャルを計算

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
    end
    if strcmp(RTD_Designs(layer).name, 'n-Si.Sub.')
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
end
function v = setpotential (v)							%/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
    global N;
    global RTD_Designs;
    global layer
	vr = zeros(1, N+1);
	vl = zeros(1, N+1);									%/* ポテンシャルの左側からの極限vrと右側からの極限			*/
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
    global RTD_Designs;
    E1 = min(v);
    E2 = max(v);
    %E2 = const.EMAX;
    dE = 5e-4;
    E = E1+1e-9 :dE: E2-1e-9;
    t22s = zeros(1,length(E));
    T = zeros(1,length(E));
	for i = 1 : length(E)
		[T(i), t11, t12, t21, t22] = TransMatrix(v, E(i), const, mass, N);
        t22s(i) = t22(N);
        %phi = waveFunction_zenkashiki(v, E(i), const, mass, N);
        %phi_end(i) = phi(N);
    end    
    [p, Es] = findpeaks(1./log(abs(t22s).^2), E);
    if n < length(Es)
        En = Es(1:n);
    else
        En = Es;
    end
end

function En = eig_confinedStates(n, v, mass)
    global N;
    global const;
    global RTD_Designs;
    global layer;

    L = RTD_Designs(layer).NX/const.DX*const.ML;
    H2 = const.HBAR^2;
    C_se = (-H2/2./const.MSTAR)*(1/L^2/const.ELEC);
    dx2 = (0.31/const.DX)^2;
    SD = (diag(-2*ones(1,N-2)) + diag(ones(1,N-3), -1) + diag(ones(1,N-3), 1))/dx2;
    K = C_se.*SD;
    h = K + diag(v(2:N-1));
    [u, E] = eig(h);
    E = eig(E);
    zn = (1:N).*0.31/const.DX;
    for i = 1:40
        psi(:,i) = [0;u(:,i);0];
        %area = trapz(psi(:,i).*conj(psi(:,i)));
        %psi(:,i) = psi(:,i).*conj(psi(:,i));
        plot(zn, psi(:,i) + E(i))
    end
end

function plot_eig(st, E, u)
    global N;
    global const;
    plot((0:N-1) * const.dx, E(st,st) + u(:,st));
end


function [T, t11, t12, t21, t22] = TransMatrix(v, E, const, mass, N_all)
    kn = sqrt(2.*mass.*(E-v).*const.ELEC) / const.HBAR;
    P = zeros(2);
    ex = zeros(2);
    dx = const.dx;

    trans = diag([1 1]);%通常計算
    for n = 2:N_all
            P(1,1) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(1,2) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(2,1) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(2,2) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            ex(1,1) = exp(1i * kn(n-1) * dx);
            ex(2,2) = exp(-1i * kn(n-1) * dx);
            trans = 0.5*P*ex*trans;
            t11(n) = trans(1, 1);
            t12(n) = trans(1, 2);
            t21(n) = trans(2, 1);
            t22(n) = trans(2, 2);
    end
    T = mass(n)*kn(1)/mass(1)/kn(n)./abs(trans(1,1)).^2;
end

function wavestore = wavefunction(v, E, const, mass, N) %//齋藤が引数wavestoreを追加 (2016.10.13)
	[T, t11, t12, t21, t22] = TransMatrix(v, E, const, mass, N);                               %/* 波動関数の計算                  */
    A = zeros(1,N);
    B = zeros(1,N);
    if v(1) <= E
        A(1) = 1;
        B(1) = -t21(N)/t22(N)*A(1);
    else
        A(1) = 0;
        B(1) = 1;
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
    tra = trapz(pyy);
	wavestore = pyy./max(pyy);
end

function wavestore = waveFunction_zenkashiki(v, E, const, mass, N)
    k = sqrt(2.*mass(1).*(E-v(1)).*const.ELEC) / const.HBAR;
    dx = const.dx;
    phi = zeros(1, N);
    f = zeros(1, N);
    
    phi(1) = 0;
    f(1) = 1i * k * exp(1i * k * dx);

    for n = 2:N
        phi(n) = f(n-1)*dx + phi(n-1);
        f(n) = -const.HBAR^2/2/mass(n) * (E - v(n-1))*phi(n-1)*dx + f(n-1);
    end
    wavestore = abs(phi);
end

function [i_element, T] = cal_Current(v, E, const, mass, RTD_Designs, N_all, N_Inject)
    global layer
    kn = sqrt(2.*mass.*(E-v).*const.ELEC) / const.HBAR;
    P = zeros(2);
    ex = zeros(2);
    dx = const.dx;


	sigma = 0.1;%ガウス関数の幅
    mu = 1.3;
    amp = 0.2;
    gaus = amp*exp(-1/2/sigma/sigma*(E-mu)^2);
    T_phonon = (gaus  + E*E/270 + exp(15*E-31.5)/400+ 0.004) / (exp((E-2.3)/0.02)+1);

    k0 = sqrt(2*mass(1)*(E-v(1))*const.ELEC) / const.HBAR;
    k_Inject = sqrt(2*mass(N_Inject)*(E-v(N_Inject))*const.ELEC) / const.HBAR;

    t11_0 = k0*mass(N_Inject)/k_Inject/mass(1) / T_phonon;
    t12_0 = t11_0 - k0/k_Inject;

    t11 = ones(1, N_all) * t11_0;
    t12 = ones(1, N_all) * t12_0;
    t21 = ones(1, N_all) * conj(t12_0);
    t22 = ones(1, N_all) * conj(t11_0);
    trans = [t11_0, t12_0; conj(t12_0), conj(t11_0)];%gaussプロット
    for n = N_Inject:N_all
            P(1,1) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(1,2) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(2,1) = 1 - (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            P(2,2) = 1 + (kn(n-1) * mass(n))/(kn(n) * mass(n-1));
            ex(1,1) = exp(1i * kn(n-1) * dx);
            ex(2,2) = exp(-1i * kn(n-1) * dx);
            trans = 0.5*P*ex*trans;
            t11(n) = trans(1, 1);
            t12(n) = trans(1, 2);
            t21(n) = trans(2, 1);
            t22(n) = trans(2, 2);
    end
    T = mass(N_all)*kn(1)/mass(1)/kn(N_all)./abs(trans(1,1)).^2;
    i_plus  = RTD_Designs(1).massxy     * RTD_Designs(1).valley     * abs(T)^2 * log(1 + exp( const.ELEC*(v(1) + RTD_Designs(1).Ef - E)/const.KT ));
    i_minus = RTD_Designs(layer).massxy * RTD_Designs(layer).valley * abs(T)^2 * log(1 + exp( const.ELEC*(v(N_all) + RTD_Designs(layer).Ef - E)/const.KT ));
    i_element = i_plus + i_minus;
end

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