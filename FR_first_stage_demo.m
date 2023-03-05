clear;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infinity = 500000;
% input
WFP = [1, 2];
robot_mode = 'dual';

w = 1;
v = 1;
MT0_PM = [1; 1; 1];
MTn_PM = [1; 1; 1];

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

% varaible
Steps = sdpvar(1);
Zk = intvar(dim_action, 1);
% constraints
Constraints = [
    Steps == sum(Zk);
    MTn - MT0 == A * Zk;
    Steps >= 0;
    Zk >= 0;
    Steps >= 11;
    ];
ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
Objective = Steps;
result = optimize(Constraints, Objective, ops);
if result.problem == 0
    disp(value(Objective));
    disp(value(Zk));
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
% MT0_PM = [1; 1; 1; 1; 1];
% MTn_PM = [1; 1; 1; 1; 1];
% % fPMid = 1;

% Case 2
% 
% infinity = 500000;
% % input
% WFP = [1, 1, 1];
% robot_mode = 'single';
% 
% w = 3;
% v = 3;
% MT0_PM = [1; 1; 1];
% MTn_PM = [1; 1; 1];
% 

% Case 3

% infinity = 500000;
% % input
% WFP = [1, 2, 1, 1];
% robot_mode = 'single';
% 
% w = 3;
% v = 2;
% MT0_PM = [1; 1; 1; 1; 1];
% MTn_PM = [1; 1; 1; 1; 1];
