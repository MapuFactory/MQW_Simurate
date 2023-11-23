SETFILE = "set.dat";

%count = 0 %csv保存用のカウンタ　手動で変える


%filename = datestr(now, 'yyyy_mm_dd_') + num2str(count);

const = Constant();

[layer, RTD_Designs] = setmaterial(const);

				%/* ML=1e-9にすれば、z[nm]。ML=0.31e-9ならz[ML]	*/
N = RTD_Designs(layer).NX;

zn = (1 : N)*const.dx*1e9;

mass = zeros(1,N);
die = zeros(1,N);
d = zeros(1,N);
bar = zeros(1,N);
mass(1:RTD_Designs(1).NX) = RTD_Designs(1).mass*const.MSTAR;
die(1 : RTD_Designs(1).NX) = RTD_Designs(1).die;
d(1 : RTD_Designs(1).NX) = RTD_Designs(1).d;
bar(1 : RTD_Designs(1).NX) = RTD_Designs(1).bar;
for i = 2:layer
    mass(RTD_Designs(i-1).NX+1:RTD_Designs(i).NX) = RTD_Designs(i).mass*const.MSTAR; 
    die(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).die;
    d(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).d;
    bar(RTD_Designs(i-1).NX+1 : RTD_Designs(i).NX) = RTD_Designs(i).bar;
end


n = 20;								%/* 計算したい準位の数	*/
V_max = 4;
V_width = 100;
V_shift = 0;
dv = V_max/V_width;
I = zeros(1, V_width+1);
V = zeros(1, V_width+1);
M = int16(const.VI/const.DELTA)+1;
I_elem = zeros(1, M);

check = [strcmp(RTD_Designs(1).name, 'n-Si.Sub.'), strcmp(RTD_Designs(layer).name, 'n-Si.Sub.')];

for i = 0 : V_width
    V_input = (i+V_shift) * dv;
    potential = calc_potential(V_input, die, bar, RTD_Designs, layer, const, N, check);
    V(i+1)=potential(1)+RTD_Designs(1).Ef-(potential(N)+RTD_Designs(layer).Ef);
    %if V(i+1) > 5
      %  break;
    %end
   I(i+1) =  current(potential, mass, layer, const, RTD_Designs, N, M);
    if round(i/V_width*100) > round((i-1)/V_width*100)
        fprintf(string(fix(i/V_width*10)))
    end
end
fprintf('\n')
J = const.II * RTD_Designs(1).massxy * const.TEMP / 1e4 .* I;

function v = calc_potential(V_all, die, bar, RTD_Designs, layer, const, N, check)

	v = potential0(V_all, die, RTD_Designs, layer, const, N, check);			%/*近似式 初期値として使用*/
	v = setpotential(v, bar, N);				%/*階段近似適用*/
    v = v-min(v);
	
end
function v = potential0 (V_all, die, RTD_Designs, layer, const, N, check)
    D = V_all / sum([RTD_Designs.d] ./ [RTD_Designs.die]);%全体の電束密度を計算
    % 先に全ての要素を0で初期化
    % 繰り返し処理でaの要素を更新
    E_point = D./die*const.dx;
    v = V_all - cumsum([0, E_point(2:N)]);
    if check(1)
        die = die(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX);
        v(RTD_Designs(1).NX+1:RTD_Designs(layer-1).NX) = V_all - cumsum(D./die .* const.dx); %1つ前のvに対して、-D/ε_n*dxを計算することで蓄電状態にない系のポテンシャルを計算

        A = sqrt(2*const.KB*const.TEMP*RTD_Designs(1).Q*1e6/RTD_Designs(1).die);
        q = const.ELEC / (const.KB * const.TEMP);
        fai = - vacc(D, RTD_Designs(1).Q, const, RTD_Designs);
        E = A * sqrt(exp(-q*fai) + q*fai - 1);
        tempV = 0;
        for i = RTD_Designs(1).NX : -1 :1
            v(i) = tempV + v(RTD_Designs(1).NX+1);
            fai = fai + E*const.dx;
            tempV= tempV + E*const.dx;
            E = A*sqrt(exp(-q*fai) + q*fai - 1);
        end
    end
    if check(2)
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
function v = setpotential (v, bar, N)							%/* 電位分布vに伝導バンド不連続を導入し、階段近似を適用する。*/
    v = v(1:N) + bar(1:N);
end
function Vacc =  vacc(D, nd, const, RTD_Designs)
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

function En = eig_confinedStates(n, v, mass, N, const)
    L = 1e-9;
    H2 = const.HBAR^2;
    C_se = (-H2/2./mass(1,N-2))*(1/L^2/const.ELEC);
    dx2 = (0.31/const.DX)^2;
    SD = (diag(-2*ones(1,N-2)) + diag(ones(1,N-3), -1) + diag(ones(1,N-3), 1))/dx2;
    K = C_se.*SD;
    h = K + diag(v(2:N-1));
    [u, E] = eig(h);
    E = eig(E);
    En = E(1:n);
end

function I = current(potential, mass, layer, const, RTD_Designs, N, M)
    T = zeros(1, M);
    Emin = potential(1);
    Emax = potential(1) + const.VI;
    Es = Emin:const.DELTA:Emax;
    ev = reshape(Es, [], 1) - potential;
    mev = mass.*ev;
    kn = sqrt(2.*mev.*const.ELEC) / const.HBAR;
    T = TransMatrix(kn, const, mass, N, M);
    I_elem = T .*  RTD_Designs(1).valley .* log( (1 + exp( const.ELEC/(const.KB*const.TEMP) .* (potential(1) - RTD_Designs(1).Ef - Es) ))./(1 + exp( const.ELEC/(const.KB*const.TEMP) .* (potential(layer) - RTD_Designs(layer).Ef - Es) ) ) );
    I = trapz(Es, I_elem);
end

function T = TransMatrix(kn, const, mass,  N, M)
    dx = const.dx;
    
    P = zeros(2, 2, M, N-1);
    ex = zeros(2, 2, M, N-1);
    T = zeros(1, M);
    P(1,1, :, :) = 1 + (sin(kn(:, 1:N-1)*dx) .* mass(:, 2:N))./(sin(kn(:, 2:N)*dx) .* mass(1:N-1));
    P(1,2, :, :) = 1 - (sin(kn(:, 1:N-1)*dx) .* mass(:, 2:N))./(sin(kn(:, 2:N)*dx) .* mass(1:N-1));
    P(2,1, :, :) = 1 - (sin(kn(:, 1:N-1)*dx) .* mass(:, 2:N))./(sin(kn(:, 2:N)*dx) .* mass(1:N-1));
    P(2,2, :, :) = 1 + (sin(kn(:, 1:N-1)*dx) .* mass(:, 2:N))./(sin(kn(:, 2:N)*dx) .* mass(1:N-1));    
    ex(1,1, :, :) = exp(1i * kn(:, 1:N-1) * dx);
    ex(2,2, :, :) = exp(-1i * kn(:, 1:N-1) * dx);
    trans = repmat(diag([1 1]), 1, 1, M);
    p_ex = pagemtimes( P, ex );
    for n = 1:N-1
        trans = 0.5* pagemtimes( p_ex(:, :, :, n), trans );
    end
    T = (mass(N).*kn(:, 1)./mass(1)./kn(:, N))./reshape(abs(trans(1,1, :)).^2, [], 1);
    T = reshape(T, 1, []);
end

function T = TransMatrix2(kn, const, mass,  N)
    P = zeros(2);
    ex = zeros(2);
    dx = const.dx;

    trans = diag([1 1]);%通常計算
    for n = 2:N
            P(1,1) = 1 + (sin(kn(n-1)*dx) * mass(n))/(sin(kn(n)*dx) * mass(n-1));
            P(1,2) = 1 - (sin(kn(n-1)*dx) * mass(n))/(sin(kn(n)*dx) * mass(n-1));
            P(2,1) = 1 - (sin(kn(n-1)*dx) * mass(n))/(sin(kn(n)*dx) * mass(n-1));
            P(2,2) = 1 + (sin(kn(n-1)*dx) * mass(n))/(sin(kn(n)*dx) * mass(n-1));
            ex(1,1) = exp(1i * kn(n-1) * dx);
            ex(2,2) = exp(-1i * kn(n-1) * dx);
            trans = 0.5*P*ex*trans;
    end
    T = (mass(n)*kn(1)/mass(1)/kn(n))./(abs(trans(1,1)).^2);
end

function [layer, RTD_Designs] = setmaterial(const) % m=-100だと、出力されない。実行はされる。引き数mは、物性値を変更したい場合に使用。
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