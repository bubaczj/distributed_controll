%% script start

close all;
clear;

addpath('functions');
addpath('controls');

data.A1 = 505;
data.a1 = 22;
data.a2 = 16;
data.C2 = 0.65;
data.Fd = 13;

h20 = 33.7852;

[A, B, X0, U0] = linAB(h20, data);

step = [1; 0];
step_val = 0.1 * U0(1);

step_dis = [0; 1];
step_val_dis = 0.1 * U0(2);

T = [0, 2000];


%% step response generation

options = odeset('MaxStep', 1, 'Refine', 1);

%control
controlhandle = step * step_val;

control_delay([], [], 0);
step_sol = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle), X0, A, B), T, X0, options);

t = [0 : 1 : 2000];
X = (deval(step_sol, t))';

%disruption
controlhandle = step_dis * step_val_dis;

control_delay([], [], 0);
step_sol_dis = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle), X0, A, B), T, X0, options);

t_dis = [0 : 1 : 2000];
X_dis = (deval(step_sol_dis, t))';

%% plots

% figure(1); %h2 (X(2))
% plot(t, X(:,2));
% grid on;
% 
% figure(2); %h2 (X(2))
% plot(t_dis, X_dis(:,2));
% grid on;


%% matrix

%control

s = X(:,2);
s = (s - s(1)) / (s(end) - s(1));

D = 1600; %horyzont dynamiki
N = 1600; %horyzont predykcji
Nu = 120; %horyzont sterowania

Mp = zeros(N, D-1);
M = zeros(N, Nu);

lambda = 0.01;
I = lambda * eye(Nu);

for i = 1 : N
    for j = 1 : D - 1
        Mp(i, j) = s(min(i + j, D)) - s(j);
    end
    
    for j = 1 : Nu
        if i - j + 1 > 0
            M(i, j) = s(i - j + 1);
        end    
    end
end

K = (M' * M + I) \ M';
K1 = K(1,:);

mpc.ku = K1 * Mp;
mpc.ke = sum(K1);
mpc.D = D;

%disruption

s_dis = X_dis(:,2);
s_dis = (s_dis - s_dis(1)) / (s_dis(end) - s_dis(1));

D_dis = 1100; %horyzont dynamiki
N_dis = 1100; %horyzont predykcji
Nu_dis = 180; %horyzont sterowania

Mp_dis = zeros(N_dis, D_dis-1);
M_dis = zeros(N_dis, Nu_dis);

lambda_dis = 1;
I_dis = lambda_dis * eye(Nu_dis);

for i = 1 : N_dis
    for j = 1 : D_dis - 1
        Mp_dis(i, j) = s_dis(min(i + j, D_dis)) - s_dis(j);
    end
    
    for j = 1 : Nu_dis
        if i - j + 1 > 0
            M_dis(i, j) = s_dis(i - j + 1);
        end    
    end
end

K_dis = (M_dis' * M_dis + I_dis) \ M_dis';
K1_dis = K_dis(1,:);

mpc.ku_dis = K1_dis * Mp_dis;
%mpc.ke = sum(K1);
mpc.D_dis = D_dis;

save('controls\mpc.mat', 'mpc');



