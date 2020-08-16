%%%%% Porównanie modelu nieliniowego, liniowego, rozmytego równomiernie i
%%%%% nierównomiernie wzglêdem h

close all;
clear;
warning ('off','all');   % LaTeX interpreter warnings

addpath('../functions');
addpath('../controls');
addpath('../plant');


%% Parametry modelu 
data.A1 = 505;
data.a1 = 22;
data.a2 = 16;
data.C2 = 0.65;
data.Fd = 13;
data.F10 = 80;
data.h20 = 33.7852;


% modele s¹ wyznaczane w n równomiernie rozmieszczonych punktach tak, ¿e
% skrajne modele s¹ wyznaczone wokó³ punktów p_pracy_min/max = p0 -/+ rng * p0
n = 5;                                                      % liczba punktów pracy
rng = 0.5;                                                  % zakres rozmywania
[A_X, B_X, X0_X, U0_X, range_X] = linAB_TS_X(data, n, rng); % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów po h
q = 1.3;
[A_X_q, B_X_q, X0_X_q, U0_X_q, range_X_q, activ_X_q] = linAB_TS_X_q(data, n, rng, q); % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów po h
[A_, B_, X0_, U0_] = linAB(data.h20, data);                     % wektor sterowania i stanu w zadanym punkcie pracy

%% Parametry symulacji
T = [0, 9000];                                  % Czas symulacji
controlhandle = @(t)F_steps(t);                  % Obs³uga sterowania
options = odeset('MaxStep', 1, 'Refine', 1);    % Ustawienia solvera


%% Symulacja
% Obiekt nieliniowy
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
sol_nonlin = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t)), U0_, data), T, X0_, options);
[~, t_nonlin, U_nonlin] = control_delay([], [], 1);% Przebieg sterowania

% Obiekt liniowy
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
sol_lin = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(t)), X0_, A_, B_), T, X0_, options);
[~, t_lin, U_lin] = control_delay([], [], 1);   % Przebieg sterowania

% Obiekt rozmyty równomiernie
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
sol_TS = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(t)), X0_X, U0_X, A_X, B_X, n, U0_, range_X), T, X0_, options);
[~, t_TS, U_TS] = control_delay([], [], 1);         % Przebieg sterowania

% Obiekt rozmyty nierównomiernie
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
sol_TS_q = ode45(@(t, x)plant_TS_X_q(t, x, control_delay(t, controlhandle(t)), X0_X_q, U0_X_q, A_X_q, B_X_q, n, U0_, activ_X_q), T, X0_, options);
[~, t_TS_q, U_TS_q] = control_delay([], [], 1);         % Przebieg sterowania

t = T(1) : T(2);                                % Wektor czasu symulacji
X_nonlin = (deval(sol_nonlin, t))';
X_lin = (deval(sol_lin, t))';
X_TS = (deval(sol_TS, t))';
X_TS_q = (deval(sol_TS_q, t))';


%% Wykresy
figure('NumberTitle', 'off', 'Name', 'Porównanie odpowiedzi na skok ró¿nych obiektów');
subplot(1,2,1);
plot(t, X_nonlin(:,2),'b',t, X_lin(:,2),'r',t, X_TS(:,2),'c',t, X_TS_q(:,2),'m');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2[cm]$', 'interpreter', 'latex');
title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
legend('obiekt nieliniowy','obiekt liniowy','obiekt rozmyty równomiernie','obiekt rozmyty nierównomiernie','location', 'northwest');
%legend('off');
grid on;

subplot(1,2,2);
plot(t_nonlin, U_nonlin(1,:) + 80);
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$F[\frac{cm^3}{s}]$', 'interpreter', 'latex');
title('Przebieg sterowania', 'interpreter', 'latex');
grid on;



