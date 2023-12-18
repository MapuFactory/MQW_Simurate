% パラメータの設定
% ドーピング濃度の設定
ND = 1e16; % n型シリコンのドナー濃度 (atoms/cm^3)

% 物理定数
q = 1.602e-19; % 電気素量 (C)

% 電荷密度の計算
% n型シリコンでは、ドーピング濃度が電子濃度に近似される
N = 100; % グリッド点の数
L = 1e-2; % 領域の長さ (m)
epsilon = 11.7 * 8.854e-12; % シリコンの誘電率 (F/m)
rho = q * ND; % 電荷密度 (C/m^3)
dx = L / (N-1); % グリッド間隔

% 境界条件の電位
phi_left = 5; % 左端の電位 (V)
phi_right = 0; % 右端の電位 (V)

% 電位の初期化
phi = zeros(N, 1);

% 電荷密度の設定（一定と仮定）
rho_vec = rho * ones(N, 1);

% 有限差分法による行列の構築
A = -2 * eye(N);
for i = 1:N-1
    A(i, i+1) = 1;
    A(i+1, i) = 1;
end
A = A / dx^2;

% 境界条件の適用
phi(1) = phi_left; % 左端での電位
phi(end) = phi_right; % 右端での電位
A(1, :) = 0; A(1, 1) = 1;
A(end, :) = 0; A(end, end) = 1;

% ポアソン方程式の解
B = -rho_vec / epsilon;
B(1) = phi_left;
B(end) = phi_right;
phi = A \ B;

% 結果のプロット
x = linspace(0, L, N);
plot(x, phi);
xlabel('Position (m)');
ylabel('Potential (V)');