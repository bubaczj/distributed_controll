%% script start

%1. run
%2. pray to God

close all;
clear;

addpath('../functions');
addpath('../controls');
addpath('../plants');

load('controls/mpc.mat');
load('controls/mpc_TS.mat');

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
T = [0, 600]; %500

des_h2 = 45;

%% simulation

controlhandle = @SL_TS_X;

disp('step 1');
options = odeset('MaxStep', 1, 'Refine', 1);


des_h2 = data.h20 *1.5;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
[~, t_U, U] = control_delay([], [], 1);
t_U = t_U';
U = U';
save('SL_results');

disp('step 2');
des_h2 = data.h20 *1.25;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t2, X2] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[~, t_U2, U2] = control_delay([], [], 1);
t_U2 = t_U2';
U2 = U2';
save('SL_results');

disp('step 3');
des_h2 = data.h20 *1.1;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t3, X3] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[~, t_U3, U3] = control_delay([], [], 1);
t_U3 = t_U3';
U3 = U3';
save('SL_results');

disp('step 4');
des_h2 = data.h20 *0.9;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t4, X4] = ode15s(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[~, t_U4, U4] = control_delay([], [], 1);
t_U4 = t_U4';
U4 = U4';
save('SL_results');

disp('step 5');
des_h2 = data.h20 *0.75;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t5, X5] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[~, t_U5, U5] = control_delay([], [], 1);
t_U5 = t_U5';
U5 = U5';
save('SL_results');


disp('step 6');
des_h2 = data.h20 *0.5;
%non-lin
controlhandle([], [], [], mpc_TS, 0);
control_delay([], [], 0);
[t6, X6] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, mpc_TS)), U0, data), T, X0, options);
%[t, X] = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle), U0, data), T, X0, options);
[~, t_U6, U6] = control_delay([], [], 1);
t_U6 = t_U6';
U6 = U6';


%% MY LEGACY
 figure(21); %h2 (X(2))
 hold on;
 plot(t_U, U(:, 1) + U0(1));
 plot(t_U2, U2(:, 1) + U0(1));
 plot(t_U3, U3(:, 1) + U0(1));
 plot(t_U4, U4(:, 1) + U0(1));
 plot(t_U5, U5(:, 1) + U0(1));
 plot(t_U6, U6(:, 1) + U0(1));
 grid on
 
 figure(22); %h2 (X(2))
 hold on;
 plot(t, X(:, 2));
 plot(t2, X2(:, 2));
 plot(t3, X3(:, 2));
 plot(t4, X4(:, 2));
 plot(t5, X5(:, 2));
 plot(t6, X6(:, 2));
 grid on
 
 save('SL_results');