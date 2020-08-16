%% script start

close all;
clear;

addpath('functions');
addpath('controls');
addpath('plants');

data.A1 = 505;
data.a1 = 22;
data.a2 = 16;
data.C2 = 0.65;
data.Fd = 13;
data.F10 = 80;
data.h20 = 33.7852;

n = 5;
rng = 0.5;
[A_TS_U, B_TS_U, X0_TS_U, U0_TS_U, range_TS_U] = linAB_TS_U(data, n, rng);
[A_TS_X, B_TS_X, X0_TS_X, U0_TS_X, range_TS_X] = linAB_TS_X(data, n, rng);

[A, B, X0, U0] = linAB(data.h20, data);

mpc_TS.U0_TS_U = U0_TS_U;
mpc_TS.U0_TS_X = U0_TS_X;
mpc_TS.X0_TS_U = X0_TS_U;
mpc_TS.X0_TS_X = X0_TS_X;
mpc_TS.U0 = U0;
mpc_TS.X0 = X0;
mpc_TS.range_U = range_TS_U;
mpc_TS.range_X = range_TS_X;

step = [1; 0];
step_val_U = 0.1 * U0(1);
step_val_X = 0.1 * X0(2);

% step_dis = [0; 1];
% step_val_dis = 0.1 * U0(2);

T = [0, 1600];


%% step response generation

options = odeset('MaxStep', 1, 'Refine', 1);

% %control
% %TS_U
% for i = 1 : n
%     controlhandle = step * step_val_U;
% 
%     control_delay([], [], 0);
%     step_sol = ode45(@(t, x)plant_lin(t, x, ...
%         control_delay(t, controlhandle), X0_TS_U(:,i), A_TS_U(:,:,i), B_TS_U(:,:,i)), T, X0_TS_U(:,i), options);
% 
%     t = [0 : 1 : 1600];
%     X_U(:,:,i) = (deval(step_sol, t))';
% end


%TS_X
for i = 1 : n
    controlhandle = step * step_val_X;

    control_delay([], [], 0);
    step_sol = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_TS_X(:,i), A_TS_X(:,:,i), B_TS_X(:,:,i)), T, X0_TS_X(:,i), options);

    t = [0 : 1 : 1600];
    X_X(:,:,i) = (deval(step_sol, t))';
end


% %disruption
% controlhandle = step_dis * step_val_dis;
% 
% control_delay([], [], 0);
% step_sol_dis = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle), X0, A, B), T, X0, options);
% 
% t_dis = [0 : 1 : 1600];
% X_dis = (deval(step_sol_dis, t))';

%% plots
% figure(1);
% for i = 1 : n
% %figure(1); %h2 (X(2))
% plot(t, X(:,2,i));
% grid on;
% hold on;
% end
% 
% figure(2); %h2 (X(2))
% plot(t_dis, X_dis(:,2));
% grid on;


%% matrix

%control


D = 1600;
N = 1600;
Nu = 40;
lambda = [10 1 0.5 0.1 0.01];

% %TS_U
% for k = 1 : n
% 
%     s = X_U(:,2,k);
%     s = (s - s(1)) / (s(end) - s(1));
% 
%     D = D; %horyzont dynamiki
%     N = N; %horyzont predykcji
%     Nu = Nu; %horyzont sterowania
% 
%     Mp = zeros(N, D-1);
%     M = zeros(N, Nu);
% 
%     lambda = lambda;
%     I = lambda * eye(Nu);
% 
%     for i = 1 : N
%         for j = 1 : D - 1
%             Mp(i, j) = s(min(i + j, D)) - s(j);
%         end
% 
%         for j = 1 : Nu
%             if i - j + 1 > 0
%                 M(i, j) = s(i - j + 1);
%             end    
%         end
%     end
% 
%     K = (M' * M + I) \ M';
%     K1 = K(1,:);
% 
%     mpc_TS.ku_U(:,k) = K1 * Mp;
%     mpc_TS.ke_U(k) = sum(K1);
%     mpc_TS.D_U(k) = D;
% 
% end
% mpc_TS.n_U = n;

%TS_X
    Mp = zeros(N, D-1, n);
    M = zeros(N, Nu, n);
for k = 1 : n

    s = X_X(:,2,k);
    s = (s - s(1)) / (s(end) - s(1));

    mpl_TS.D_X(i) = D; %horyzont dynamiki
    mpl_TS.N_X(i) = N; %horyzont predykcji
    mpl_TS.Nu_X(i) = Nu; %horyzont sterowania

    lambdaTS = lambda(k);
    I = lambdaTS * eye(Nu);

    for i = 1 : N
        for j = 1 : D - 1
            Mp(i, j, k) = s(min(i + j, D)) - s(j);
        end

        for j = 1 : Nu
            if i - j + 1 > 0
                M(i, j, k) = s(i - j + 1);
            end    
        end
    end
    
    

    K = (M(:,:,k)' * M(:,:,k) + I) \ M(:,:,k)';
    K1 = K(1,:);

    mpc_TS.ku_X(:,k) = K1 * Mp(:,:,k);
    mpc_TS.ke_X(k) = sum(K1);
    mpc_TS.D_X(k) = D;
    mpc_TS.Nu_X(k) = Nu;
    mpc_TS.N_X(k) = N;
    mpc_TS.lambda_X(k) = lambdaTS;
    

end
mpc_TS.Mp_X(:,:, :) = Mp(:,:,:);
mpc_TS.M_X(:,:, :) = M(:,:,:);
mpc_TS.n_X = n;

% %disruption
% for k = 1 : n
% s_dis = X_dis(:,2);
% s_dis = (s_dis - s_dis(1)) / (s_dis(end) - s_dis(1));
% 
% D_dis = 1100; %horyzont dynamiki
% N_dis = 1100; %horyzont predykcji
% Nu_dis = 180; %horyzont sterowania
% 
% Mp_dis = zeros(N_dis, D_dis-1);
% M_dis = zeros(N_dis, Nu_dis);
% 
% lambda_dis = 3e6;
% I_dis = lambda_dis * eye(Nu_dis);
% 
% for i = 1 : N_dis
%     for j = 1 : D_dis - 1
%         Mp_dis(i, j) = s_dis(min(i + j, D_dis)) - s_dis(j);
%     end
%     
%     for j = 1 : Nu_dis
%         if i - j + 1 > 0
%             M_dis(i, j) = s_dis(i - j + 1);
%         end    
%     end
% end
% 
% K_dis = (M_dis' * M_dis + I_dis) \ M_dis';
% K1_dis = K_dis(1,:);
% 
% mpc_TS.ku_dis = K1_dis * Mp_dis;
% %mpc_TS.ke = sum(K1);
% mpc_TS.D_dis = D_dis;


save('controls\mpc_TS.mat', 'mpc_TS');