clear;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infinity = 500000;
% input
WFP = [1, 3, 1];
robot_mode = 'dual';

w = 3;
v = 2;
rho_raw = [70; 280; 280; 280; 280; 70];
% tauI = rho_raw - [7; -7; 71; 149; 227; 65];
% tauE = rho_raw - [-16; 0; 204; 12; 108; 70];

rho = rho_raw([1, 3:6]);
% varthetaI = [7; 71; 149; 227; 65];
varthetaI = [50; 114; 192; 270; 65];
% varthetaE = [-16; 204; 12; 108; 70];
varthetaE = [40; 68; 164; 260; 70];
MT0_PM = [1; 1; 1; 1; 1];
MTn_PM = [1; 1; 1; 1; 1];
varpi_0 = [0; 0; 0; 1];
varpi_1 = [1; 0; 0; 0];
delta = ones(5, 1) * 30;

A_PM_idx = [2, 4:6, 8];
A_R_idx = [1, 3, 7, 9];
Zl_idx = [2, 4:6, 10];
Zu_idx = [3, 7:9, 11];
step = 33;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPM = sum(WFP) + 2;
dim_s = sum(WFP) + 2 + size(WFP, 2) - 1;

% M S K
MT0 = zeros(dim_s, 1);
MTn = zeros(dim_s, 1);
S = zeros(dim_s + 1, 1);
S(1) = 1;
S(end) = 1;
K = ones(dim_s + 1, 1);
curid = 1;
curid_PM = 1;
id_robot = [1];
for i = 1:size(WFP, 2)
    curid = curid + 1;
    MT0(curid:curid + WFP(i) - 1) = MT0_PM(curid_PM:curid_PM + WFP(i) - 1);
    MTn(curid:curid + WFP(i) - 1) = MTn_PM(curid_PM:curid_PM + WFP(i) - 1);
    curid_PM = curid_PM + WFP(i);
    curid = curid + WFP(i);
    S(curid) = 1;
    id_robot = [id_robot, curid];
end

if strcmpi(robot_mode, 'single')
    R = 1;
elseif strcmpi(robot_mode, 'dual')
    MT0(1) = 0;
    MTn(1) = 1;
    R = 2;
else
    disp('No matched robot mode')
end

K(id_robot) = R;
K(end) = R;
% A
dim_action = 2 * sum(WFP) + 2;
A = zeros(1, dim_s);
A(1, 1) = 1;
A_plus = -1;
curid = 1;
for i = 1:size(WFP, 2)
    curid = curid + 1;
    atmp = zeros(dim_s, 2 * WFP(i));
    atmp_plus = zeros(1, 2 * WFP(i));
    for j = 1:WFP(i)
        atmp(curid - 1, j) = -1;
        atmp(curid + j - 1, j) = 1;
        atmp(curid + j - 1, WFP(i) + j) = -1;
        atmp(curid + WFP(i), WFP(i) + j) = 1;
        
        atmp_plus(j) = 1;
        atmp_plus(WFP(i) + j) = -1;
    end
    curid = curid + WFP(i);
    A = [A; atmp'];
    A_plus = [A_plus; atmp_plus'];
end
atmp = zeros(1, dim_s);
atmp(end) = -1;
A = [A; atmp];
A = A';
A_plus = [A_plus; 1];
A_plus = A_plus';
A_plus = [A; A_plus];

A_PM = A(A_PM_idx, :);
A_R = A(A_R_idx, :);

B = 500000;

% variable
gamma = binvar(sum(WFP), step, 'full');
varpi = intvar(size(WFP, 2) + 1, step, 'full');
vartheta = sdpvar(sum(WFP), step, 'full');

Z = binvar(dim_action, step, 'full');
omega = sdpvar(1, step, 'full');

x = sdpvar(sum(WFP), step, 'full');
y = sdpvar(sum(WFP), step, 'full');
h = sdpvar(sum(WFP), step, 'full');
phi = sdpvar(sum(WFP), step, 'full');

idle_time = sdpvar(1);

lambda = sdpvar(1);

% Constraint
Constraints = [
    gamma(:, 1) == MT0_PM + A_PM * Z(:, 1);
    varpi(:, 1) == varpi_0 + A_R * Z(:, 1);
%     vartheta(:, 1) == varthetaI + Z([2, 4, 5], 1) .* rho(2:4) + Z([3, 6, 7], 1) .* (-varthetaI + w + v + omega(1)) - MT0_PM * (w + v + omega(1));
    vartheta(:, 1) == varthetaI + Z(Zl_idx, 1) .* rho + Z(Zu_idx, 1) .* (-varthetaI + w + v) + y(:, 1) - MT0_PM * (w + v + omega(1));
%     x(:, 1) <= varthetaI;
%     x(:, 1) >= varthetaI - rho(2:4) .* (1 - Z([3, 6, 7], 1));
%     -delta .* Z([3, 6, 7], 1) <= x(:, 1) <= rho(2:4) .* Z([3, 6, 7], 1);
    y(:, 1) <= omega(1);
    y(:, 1) >= omega(1) - B * (1 - Z(Zu_idx, 1));
    0 * Z(Zu_idx, j) <= y(:, 1) <= B * Z(Zu_idx, 1);
    sum(varpi, 1) <= ones(1, step) * R;
    varpi >= 0;
    varpi <= R;
    vartheta >= -repmat(delta, 1, step);
    vartheta <= repmat(rho, 1, step);
    omega >= 0;
    sum(Z, 1) == ones(1, step);
    varthetaI - omega(1) - v <= B * (1 - Z(Zu_idx, 1));
    ];

for j = 2:step
    Constraints = [Constraints;
        gamma(:, j) == gamma(:, j - 1) + A_PM * Z(:, j);
        varpi(:, j) == varpi(:, j - 1) + A_R * Z(:, j);
%         vartheta(:, j) == vartheta(:, j - 1) + Z([2, 4, 5], j) .* rho(2:4) + Z([3, 6, 7], j) .* (-vartheta(:, j - 1) + w + v + omega(j)) - gamma(:, j - 1) * (w + v + omega(j));
        vartheta(:, j) == vartheta(:, j - 1) + Z(Zl_idx, j) .* rho - x(:, j) + Z(Zu_idx, j) .* (w + v + delta) + y(:, j) - gamma(:, j - 1) * (w + v) - h(:, j);
        x(:, j) <= vartheta(:, j - 1) + delta;
        x(:, j) >= vartheta(:, j - 1) + delta - (rho + delta) .*(1 - Z(Zu_idx, j));
        0 * Z(Zu_idx, j) <= x(:, j) <= (rho + delta) .* Z(Zu_idx, j);
        y(:, j) <= omega(j);
        y(:, j) >= omega(j) - B * (1 - Z(Zu_idx, j));
        0 * Z(Zu_idx, j) <= y(:, j) <= B * Z(Zu_idx, j);
        h(:, j) <= omega(j);
        h(:, j) >= omega(j) - B * (1 - gamma(:, j - 1));
        0 * Z(Zu_idx, j) <= h(:, j) <= B * gamma(:, j - 1);
        vartheta(:, j - 1) - omega(j) - v <= B * (1 - Z(Zu_idx, j));
        ];
end

Constraints = [Constraints;
%     lambda >= 0;
    lambda == (w + v) * step + sum(omega) + idle_time;
    gamma(:, end) == MTn_PM;
    varthetaE == vartheta(:, end) - idle_time * MTn_PM;
%     -50 <= varthetaE - vartheta(:, end) <= 50;
    idle_time >= 0;
    varpi(:, end) == varpi_1;
    ];

% ops = sdpsettings('verbose', 0, 'solver', 'bnb','usex0', 1);
% ops = sdpsettings(ops,'bmibnb.uppersolver','fmincon');
ops = sdpsettings('verbose', 0, 'solver', 'gurobi', 'usex0', 1);
% ops = sdpsettings('verbose',1,'debug',1,'solver','cplex','savesolveroutput',1,'savesolverinput',1);
% ops.cplex.exportmodel='abc.lp';

Objective = lambda;
tic
result = optimize(Constraints, Objective, ops);
toc

if result.problem == 0
    disp(value(Objective));
    Z_value = value(Z);
    v_value = value(vartheta);
    o_value = value(omega);
    m_value_pm = value(gamma);
    output = value(sum(Z(end, :)))
else
    disp('No solution');
    result.info
    yalmiperror(result.problem)
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Additional cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 1

% infinity = 500000;
% % input
% WFP = [1, 3, 1];
% robot_mode = 'dual';
% 
% w = 3;
% v = 2;
% rho_raw = [70; 280; 280; 280; 280; 70];
% % tauI = rho_raw - [7; -7; 71; 149; 227; 65];
% % tauE = rho_raw - [-16; 0; 204; 12; 108; 70];
% 
% rho = rho_raw([1, 3:6]);
% % varthetaI = [7; 71; 149; 227; 65];
% varthetaI = [50; 114; 192; 270; 65];
% % varthetaE = [-16; 204; 12; 108; 70];
% varthetaE = [40; 68; 164; 260; 70];
% MT0_PM = [1; 1; 1; 1; 1];
% MTn_PM = [1; 1; 1; 1; 1];
% varpi_0 = [0; 0; 0; 1];
% varpi_1 = [1; 0; 0; 0];
% delta = ones(5, 1) * 30;
% 
% A_PM_idx = [2, 4:6, 8];
% A_R_idx = [1, 3, 7, 9];
% Zl_idx = [2, 4:6, 10];
% Zu_idx = [3, 7:9, 11];
% step = 33;



% Case 2

% infinity = 500000;
% % input
% WFP = [1, 1, 1];
% robot_mode = 'single';
% 
% w = 3;
% v = 3;
% rho_raw = [80; 189; 189; 80];
% % tauI = rho_raw - [23; 15; 120; -1];
% % tauE = rho_raw - [76; 0; 94; -27];
% 
% rho = rho_raw([1, 2, 4]);
% varthetaI = [23; 120; -1];
% varthetaE = [76; 94; -27];
% MT0_PM = [1; 1; 1];
% MTn_PM = [1; 1; 1];
% varpi_0 = [0; 0; 0; 0];
% varpi_1 = [0; 0; 0; 0];
% delta = ones(3, 1) * 120;
% 
% A_PM_idx = [2, 4, 6];
% A_R_idx = [1, 3, 5, 7];
% Zl_idx = [2, 4, 6];
% Zu_idx = [3, 5, 7];
% step = 16;


% Case 3

% infinity = 500000;
% % input
% WFP = [1, 2, 1, 1];
% robot_mode = 'dual';
% 
% w = 3;
% v = 2;
% rho_raw = [48; 130; 130; 130; 40; 40];
% % tauI = rho_raw - [2; -18; 38; 94; 14; 24];
% % tauE = rho_raw - [-11; 0; 81; 12; 2; 13];
% 
% rho = rho_raw([1, 3:6]);
% varthetaI = [2; 38; 94; 14; 24];
% varthetaE = [-11; 81; 12; 2; 13];
% MT0_PM = [1; 1; 1; 1; 1];
% MTn_PM = [1; 1; 1; 1; 1];
% varpi_0 = [0; 0; 0; 0; 0];
% varpi_1 = [1; 0; 0; 0; 0];
% delta = ones(5, 1) * 3100;
% 
% A_PM_idx = [2, 4:5, 7, 9];
% A_R_idx = [1, 3, 6, 8, 10];
% Zl_idx = [2, 4:5, 8, 10];
% Zu_idx = [3, 6:7, 9, 11];
% step = 31;



% Case 4

% infinity = 500000;
% % input
% WFP = [1, 2, 1, 1];
% robot_mode = 'single';
% 
% w = 4;
% v = 2;
% rho_raw = [80; 230; 230; 230; 73; 80];
% % tauI = rho_raw - [2; -18; 38; 94; 14; 24];
% % tauE = rho_raw - [-11; 0; 81; 12; 2; 13];
% 
% rho = rho_raw([1, 3:6]);
% varthetaI = [62; 200; -4; 14; 2];
% varthetaE = [80; 214; 88; 45; 40];
% MT0_PM = [1; 1; 1; 1; 1];
% MTn_PM = [1; 1; 1; 1; 1];
% varpi_0 = [0; 0; 0; 0; 0];
% varpi_1 = [0; 0; 0; 0; 0];
% delta = [25; 30; 30; 25; 25];
% 
% A_PM_idx = [2, 4:5, 7, 9];
% A_R_idx = [1, 3, 6, 8, 10];
% Zl_idx = [2, 4:5, 8, 10];
% Zu_idx = [3, 6:7, 9, 11];
% step = 40;
