SETFILE = "set.dat";

count = 0 %csv保存用のカウンタ　手動で変える

filename = datestr(now, 'yyyy_mm_dd_') + num2str(count);

loop = 0; %齋藤が変数loop(ループを回す為の変数)を追加 (2016.10.13)
global RTD_Designs;
global layer;
global const;
const = Constant();

setmaterial();

				%/* ML=1e-9にすれば、z[nm]。ML=0.31e-9ならz[ML]	*/
global N;
N = RTD_Designs(layer).NX;



V_all=1.0;%4.722;%2.21906;

v = potential(V_all);
zn = (1 : N)*const.dx*1e9;
plot(zn, v);
hold on

mass = zeros(1,N);
mass(1:RTD_Designs(1).NX) = RTD_Designs(1).mass*const.MSTAR;
for i = 2:layer
    mass(RTD_Designs(i-1).NX+1:RTD_Designs(i).NX) = RTD_Designs(i).mass*const.MSTAR; 
end
n = 50;								%/* 計算したい準位の数	*/
En = eig_confinedStates(n, v, mass);
hold off

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

function En = eig_confinedStates(n, v, mass)
    global N;
    global const;
    global RTD_Designs;
    global layer;

    L = 1e-9;
    H2 = const.HBAR^2;
    mass = mass'./const.MSTAR;
    C_se = (-H2/2./const.MSTAR)*(1/L^2/const.ELEC);
    dx = 0.31/const.DX;
    MD = 0.5.* ([mass(1:N) ; mass(1)] + [mass(N) ; mass(1:N)]);
    SD = (( diag(ones(1,N-1), 1) - diag(ones(1,N)) ) ./ MD(2:N+1) - ( diag(ones(1,N)) - diag(ones(1,N-1), -1) ) ./ MD(1:N)) ./dx ./dx;
    SD(N,1) = 1/MD(1);
    SD(1,N) = 1/MD(N+1);
    K = C_se .* SD;
    h = K + diag(v(1:N));
    [u, E] = eig(h);
    E = eig(E);
    En = E;
    zn = (1:N).*0.31/const.DX;
    
%     i = 12;
%     plot_wavefunc(zn, E(i), u(:,i).^2);
%     i = 14;
%     plot_wavefunc(zn, E(i), u(:,i).^2);
%     i = 17;
%     plot_wavefunc(zn, E(i), u(:,i).^2);
%     
%     for i = 20:22
%         plot_wavefunc(zn, E(i), u(:,i).^2);
%     end
    for i = 1:30
        plot_wavefunc(zn, E(i), u(:,i).^2);
    end
end

function plot_wavefunc(zn, E, psi)
    area = trapz(zn, psi);
    psi = psi./sqrt(area);
    plot(zn, psi + E)

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