clear;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infinity = 500000;
% input
robot_mode = 'dual';
n = 3;
w = 3;
v = 2;
r = [70, 280, 70];
residency = [30, 30, 30];
K_PM = [1, 3, 1];
assert (size(K_PM, 2) == n, 'Dimension mismatch');
assert (size(residency, 2) == n, 'Dimension mismatch');
assert(size(r, 2) == n, 'Dimension mismatch');

% constant
T = 2 * n + 2;
K = zeros(T, 1);
K(2:2:end - 1) = K_PM;
if strcmpi(robot_mode, 'single')
    K(1:2:end - 1) = 1;
    K(end) = 1;
    R = 1;
elseif strcmpi(robot_mode, 'dual')
    K(1:2:end - 1) = 2;
    K(end) = 2;
    R = 2;
else
    disp('No matched robot mode.')
end

A = zeros(T);
A(1:(T + 1):end - T) = 1;
A(2:(T + 1):end - T) = -1;
A = A';
A(end, 1:2:end) = -1;
A(end, 2:2:end) = 1;

S = zeros(T,1);
S(1:2:end - 1) = 1;
S(end) = 1;

B = 500000;

H = zeros(T-1, 1);
H(1:2:end) = w + v;
H(2:2:end) = w + r;

J = zeros(T-1, 1);
J(1:2:end) = infinity;
J(2:2:end) = residency + r + w;


% varaible
lambda = sdpvar(1);
M = intvar(T, T, 'full'); 

M_0tmp = zeros(T, 1);
M_0tmp(2:2:end - 1) = K_PM;

if strcmpi(robot_mode, 'single')
    M_0tmp(end) = 1;
elseif strcmpi(robot_mode, 'dual')
    M_0tmp(1) = 1;
    M_0tmp(end) = 1;
else
    disp('No matched robot mode.')
end

M_0 = M_0tmp;
% M_0 = intvar(T, 1);
% assign(M_0, M_0tmp);
Z = binvar(T, T, 'full');
x = sdpvar(T, 1);
y = sdpvar(T, 1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constraints

Constraints = [
    lambda >= 0;
    M_0 >= 0;
    M_0 <= K;
    M >= 0;
    M <= repmat(K, 1, T);
    M(:, 1) == M_0 + A * Z(:, 1);
    ];
for j = 2:T
    Constraints = [Constraints;
    M(:, j) == M(:, j - 1) + A * Z(:, j);    
    ] ;
end
Constraints = [Constraints;
    sum(Z, 1) == ones(1, T);
    sum(Z, 2) == ones(T, 1);
    ];

Constraints = [Constraints;
    S' * M_0 == R;
    ];

Constraints = [Constraints;
    x(2:T, 1) + lambda * M_0(1:T-1, 1) - x(1:T-1, 1) >= H;
    x(2:T, 1) + lambda * M_0(1:T-1, 1) - x(1:T-1, 1) <= J;
    ];
Constraints = [Constraints;
    repmat(y', T, 1) - repmat(x, 1, T) <= (1 - Z) * B;
    repmat(y', T, 1) - repmat(x, 1, T) >= (Z - 1) * B;
    y >= 0;
    y <= lambda;
    y(2:T) - y(1:T-1) >= zeros(T-1, 1);
    y(2:T) - y(1:T-1) >= (w+v) * ones(T-1,1);
    y(1) + lambda - y(T) >= w + v;
    ];

ops = sdpsettings('verbose', 0, 'solver', 'gurobi', 'usex0', 1);
Objective = lambda;
result = optimize(Constraints, Objective, ops);
if result.problem == 0
    disp(value(Objective));
    
    %disp(value(M_0));
    %disp(value(M));
    %disp(value(Z));
    m0 = value(M_0);
    m = value(M);
    z = value(Z);
    Y = value(y);
    X = value(x);
%     disp(r' - (value(lambda) - X(2:2:end - 1)));
else
    disp('No solution');
    result.info
    yalmiperror(result.problem)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Gantte Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cycle = max(K_PM);
% cycle = sum(K_PM);
xmax = max(value(lambda)) * cycle;
% xmax = 142;
ymax = sum(K_PM) + 1.5;
axis([0, int32(xmax), 0, ymax]);
set(gca, 'xtick', 0:1:int32(xmax));
set(gca, 'ytick', 0:1:ymax);
set(gca, 'YTickLabel', {''; num2str((1:ymax - 1.5)','PM%d');'Robot';''});
xlabel('Time');
ylabel('Robot and PMs');


start_time = value(x);
duration_robot = zeros(T, 1);
duration_robot(1:2:end - 1) = w + v;
duration_robot(2:2:end) = w + v;
duration_PM = zeros(T, 1);
duration_PM(2:2:end - 1) = r;

robot_color = [[0.3, 0.3, 0.3];
                [0.7, 0.7, 0.7]];

period_PM = zeros(T, 1);
period_PM(2:2:end - 1) = cycle;
% period_PM(1:2:end - 1) = 1;
% period_PM(2) = 0;

id_PM_temp = zeros(T, 1);
id_PM_temp(2:2:end - 1) = K_PM;

residency_PM = zeros(T, 1);
X = value(x);
residency_PM(2:2:end - 1) = X(3:2:end - 1) - X(2:2:end - 1)  + M_0(2:2:end - 1) * value(lambda) - r' - w - v;
% residency_PM(2:2:end - 1) = 20;

S_T0 = [];

rect = zeros(1,4);
for i = 1:size(value(x), 1)
    for j = 1:period_PM(i)
        % PM
        rect(1) = start_time(i) + (j - 1) * value(lambda) + w + v;
        id_PM = sum(id_PM_temp(1:i)) - id_PM_temp(i) + min(j, id_PM_temp(i));
        rect(2) = id_PM - 1 + 0.7;
%         rect(3) = min(duration_PM(i), cycle * value(lambda) - rect(1));
        rect(3) = duration_PM(i);
        rect(4) = 0.6;
        rectangle('Position', rect, 'LineWidth', 0.5, 'LineStyle','-', 'FaceColor', 'w');
        
        remain_process = max(rect(1) + duration_PM(i) - cycle * value(lambda), 0);
        S_T0 = [S_T0;remain_process];
        rect(1) = rect(1) + rect(3);
        rect(3) = max(residency_PM(i), 0);
        rectangle('Position', rect, 'LineWidth', 0.5, 'LineStyle','none', 'FaceColor', 'red');
        WRTC = cycle * value(lambda) - rect(1)
        remain_WRTC = max(rect(1) + rect(3) - cycle * value(lambda), 0);
        if (remain_WRTC == 0)
            continue;
        end
        rect(1) = 0;
        rect(3) = remain_WRTC;
        rectangle('Position', rect, 'LineWidth', 0.5, 'LineStyle','none', 'FaceColor', 'red');
        
        if remain_process == 0
            continue;
        end
        rect(1) = 0;
        rect(3) = remain_process;
        rectangle('Position', rect, 'LineWidth', 0.5, 'LineStyle','-', 'FaceColor', 'w');
    end
    for j = 1:cycle
        %robot
        rect(1) = start_time(i) + (j - 1) * value(lambda);
        rect(2) = sum(K_PM) + 1 - 1 + 0.7;
        rect(3) = duration_robot(i);
        rect(4) = 0.6;
        rectangle('Position', rect, 'LineWidth', 0.5, 'LineStyle','none', 'FaceColor', robot_color(mod(i, 2) + 1, :));
        hold on;
        if(mod(i, 2) == 0)
            plot([rect(1) + rect(3), rect(1) + rect(3)], [0, sum(K_PM) + 1.3], '--k', 'LineWidth', 0.7);
        end
        
        if(i == size(value(x), 1))
                continue;
        end
        if(mod(i, 2) == 1 && rect(1) ~= 0)
            plot([rect(1), rect(1)], [0, sum(K_PM) + 1.3], '--k', 'LineWidth', 0.7);        
        end
    end
end
S_T0

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Additional cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case 1

% infinity = 500000;
% % input
% robot_mode = 'dual';
% n = 4;
% w = 3;
% v = 2;
% r = [48, 130, 40, 40];
% residency = [20, 30, 20, 20];
% K_PM = [1, 2, 1, 1];
% assert (size(K_PM, 2) == n, 'Dimension mismatch');
% assert (size(residency, 2) == n, 'Dimension mismatch');
% assert(size(r, 2) == n, 'Dimension mismatch');

% Case 2

% infinity = 500000;
% % input
% robot_mode = 'single';
% n = 4;
% w = 4;
% v = 2;
% r = [80, 230, 73, 80];
% residency = [20, 40, 20, 20];
% K_PM = [1, 2, 1, 1];
% assert (size(K_PM, 2) == n, 'Dimension mismatch');
% assert (size(residency, 2) == n, 'Dimension mismatch');
% assert(size(r, 2) == n, 'Dimension mismatch');

% Case 3

% infinity = 500000;
% % input
% robot_mode = 'single';
% n = 3;
% w = 5;
% v = 2;
% r = [50, 280, 50];
% residency = [30, 30, 30];
% K_PM = [1, 3, 1];
% assert (size(K_PM, 2) == n, 'Dimension mismatch');
% assert (size(residency, 2) == n, 'Dimension mismatch');
% assert(size(r, 2) == n, 'Dimension mismatch');