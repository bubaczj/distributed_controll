%% script start

close all;
clear;

addpath('functions');
addpath('controls');
addpath('plants');

load('controls\mpc.mat');
load('controls\mpc_TS.mat');

data.A1 = 505;
data.a1 = 22;
data.a2 = 16;
data.C2 = 0.65;
data.Fd = 13;
data.F10 = 80;
data.h20 = 33.7852;

[A, B, X0, U0] = linAB(data.h20, data);

n = 5;
rng = 0.5;
[A_TS_U, B_TS_U, X0_TS_U, U0_TS_U, range_TS_U] = linAB_TS_U(data, n, rng);
[A_TS_X, B_TS_X, X0_TS_X, U0_TS_X, range_TS_X] = linAB_TS_X(data, n, rng);

step = [1; 0];
step_val = 20;
T = [0, 900];

des_h2 = 20;

%% simulation

%controlhandle = step * step_val;
%controlhandle = @P;
controlhandle = @DMC_TS_X;
%controlhandle = @DMC_dis;
%controlhandle = @DMC;


options = odeset('MaxStep', 1, 'Refine', 1);

%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[sim_1, t_U, U] = control_delay([], [], 1);
t_U = t_U';
U = U';

controlhandle = @DMC_TS_U;%---------------------------------------------------------------------------
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[tU, XU] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[sim_1, t_UU, UU] = control_delay([], [], 1);
t_UU = t_UU';
UU = UU';

controlhandle = @DMC;
controlhandle([], [], [], mpc, 0);
control_delay([], [], 0);
[tUU, XUU] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[sim_1, t_UUU, UUU] = control_delay([], [], 1);
t_UU = t_UU';
UUU = UUU';

controlhandle = @DMC_TS_X;%-----------------------------------------------------------------------------------


% figure(6);
% plot(t_U);

%lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t_lin, X_lin] = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), X0, A, B), T, X0, options);
%[t_lin, X_lin] = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle), X0, A, B), T, X0, options);
[sim_2, t_U_lin, U_lin] = control_delay([], [], 1);
t_U_lin = t_U_lin';
U_lin = U_lin';

%TS_U
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t_TS_U, X_TS_U] = ode45(@(t, x)plant_TS_U(t, x, ...
    control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), X0_TS_U, U0_TS_U, A_TS_U, B_TS_U, n, U0, range_TS_U), ...
    T, X0, options);
[sim_3, t_U_TS_U, U_TS_U] = control_delay([], [], 1);
t_U_TS_U = t_U_TS_U';
U_TS_U = U_TS_U';

%TS_X
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t_TS_X, X_TS_X] = ode45(@(t, x)plant_TS_X(t, x, ...
    control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), X0_TS_X, U0_TS_X, A_TS_X, B_TS_X, n, U0, range_TS_X), ...
    T, X0, options);
[sim_4, t_U_TS_X, U_TS_X] = control_delay([], [], 1);
t_U_TS_X = t_U_TS_X';
U_TS_X = U_TS_X';

%% post-processing
[t_imax, ~] = size(t);
[t_lin_imax, ~] = size(t_lin);

% for i = 1 : t_imax
%     U(i, :) = controlhandle(X(i, :), des_h2, mpc);
% end
% for i = 1 : t_lin_imax
%     U_lin(i, :) = controlhandle(X_lin(i, :), des_h2, mpc);
% end

figure(1); %h1 (X(1))
plot(t, X(:, 1), t_lin, X_lin(:,1), t_TS_U, X_TS_U(:, 1), t_TS_X, X_TS_X(:, 1));
grid on;
legend('non lin', 'lin','TS','Location','southeast');
title('h1(t)');

figure(2); %h2 (X(2))
plot(t, X(:, 2), t_lin, X_lin(:,2), t_TS_U, X_TS_U(:, 2), t_TS_X, X_TS_X(:, 2));
legend('non lin', 'lin','TS','Location','southeast');
grid on;
title('h2(t)');

figure(3); %F1 (U(1))
plot(t_U, U(:, 1) + U0(1), t_U_lin,  U_lin(:, 1) + U0(1), t_U_TS_U, U_TS_U(:, 1) + U0(1), t_U_TS_X, U_TS_X(:, 1) + U0(1));
legend('non lin', 'lin', 'TS','Location','southeast');
grid on;
title('F1(t)');

figure(4); %Fd (U(1))
plot(t_U, U0(2) + U(:, 2), t_U_lin, U0(2) + U_lin(:, 2), t_U_TS_U, U_TS_U(:, 2) + U0(2), t_U_TS_X, U_TS_X(:, 2) + U0(2));
legend('non lin', 'lin');
grid on;
title('Fd(t)');

figure(5); %h2 (X(2))
plot(t, X(:, 2), tU, XU(:, 2),tUU, XUU(:, 2));
legend('TS_X', 'TS_U', 'DMC');
grid on;
title('h2(t)');


figure(6); %h2 (X(2))
plot(t_U, U(:, 1) + U0(1), t_UU, UU(:, 1) + U0(1), t_UUU, UUU(:, 1) + U0(1));
legend('TS_X', 'TS_U', 'DMC');
grid on;
title('F1(t)');
