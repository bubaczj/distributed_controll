%%%%% Odpowied� skokowa modelu rozmytego wzgl�dem warto�ci sterowania F1

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

% modele s� wyznaczane w n r�wnomiernie rozmieszczonych punktach tak, �e
% skrajne modele s� wyznaczone wok� punkt�w F_pracy_min/max = F_0 -/ +rng * F0
n = 5;                                                  % liczba punkt�w pracy
rng = 0.5;                                              % zakres rozmywania
[A, B, X0, U0, range] = linAB_TS_U(data, n, rng);       % A, B, X0, U0 - parametry n zlinearyzowanych uk�ad�w
[~,~,X0_, U0_] = linAB(data.h20, data);                 % wektor sterowania i stanu w zadanym punkcie pracy

%% Parametry symulacji
step = [1; 0];                                  % Skok jedynie warto�ci steruj�cej
step_val = 40*[-1,-0.5, 0.5, 1];                % Przyrosty warto�ci steruj�cej F = F1 + FD + step_val 
T = [0, 1400];                                  % Czas symulacji
controlhandle = step * step_val;                % Obs�uga sterowania
options = odeset('MaxStep', 1, 'Refine', 1);    % Ustawienia solvera

%% Symulacja
control_delay([], [], 0);                       % Reset funkcji realizuj�cej op�nienie
solution1 = ode45(@(t, x)plant_TS_U(t, x, control_delay(t, controlhandle(:,1)), X0, U0, A, B, n, U0_, range), T, X0_, options);
[~, t_U1, U1] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj�cej op�nienie
solution2 = ode45(@(t, x)plant_TS_U(t, x, control_delay(t, controlhandle(:,2)), X0, U0, A, B, n, U0_, range), T, X0_, options);
[~, t_U2, U2] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj�cej op�nienie
solution3 = ode45(@(t, x)plant_TS_U(t, x, control_delay(t, controlhandle(:,3)), X0, U0, A, B, n, U0_, range), T, X0_, options);
[~, t_U3, U3] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj�cej op�nienie
solution4 = ode45(@(t, x)plant_TS_U(t, x, control_delay(t, controlhandle(:,4)), X0, U0, A, B, n, U0_, range), T, X0_, options);
[~, t_U4, U4] = control_delay([], [], 1);         % Przebieg sterowania

t = T(1) : T(2);                                % Wektor czasu symulacji
X1 = (deval(solution1, t))';                    % Wektor rozwi�za�
X2 = (deval(solution2, t))';
X3 = (deval(solution3, t))';
X4 = (deval(solution4, t))';


%% Wykresy
figure('NumberTitle', 'off', 'Name', 'Odpowied� modelu rozmytego wzgl�dem h2 na skok');
subplot(2,2,1);
plot(t, X1(:,1),t, X2(:,1),t, X3(:,1),t, X4(:,1));
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_1[cm]$', 'interpreter', 'latex');
title('Przebieg wysoko�ci p�ynu w pierwszym zbiorniku', 'interpreter', 'latex');
legend('skok 1','skok 2','skok 3','skok 4', 'location', 'east');
grid on;

subplot(2,2,2);
plot(t, X1(:,2),t, X2(:,2),t, X3(:,2),t, X4(:,2));
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2[cm]$', 'interpreter', 'latex');
title('Przebieg wysoko�ci p�ynu w drugim zbiorniku', 'interpreter', 'latex');
legend('skok 1','skok 2','skok 3','skok 4', 'location', 'east');
grid on;

subplot(2,2,3);
U_s = U0(1) * ones(1, T(2) - T(1) + 1);
U_s1 = [U_s(T(1)+1:T(1) + 79),  (U_s(T(1) + 80 : T(2) + 1) + step_val(1))];
U_s2 = [U_s(T(1)+1:T(1) + 79),  (U_s(T(1) + 80 : T(2) + 1) + step_val(2))];
U_s3 = [U_s(T(1)+1:T(1) + 79),  (U_s(T(1) + 80 : T(2) + 1) + step_val(3))];
U_s4 = [U_s(T(1)+1:T(1) + 79),  (U_s(T(1) + 80 : T(2) + 1) + step_val(4))];
plot(t, U_s1, t, U_s2, t, U_s3, t, U_s4 );                                                          % op�nione 
%%plot(t_U1, U0(1) + U1(1,:), t_U2, U0(1) + U2(1,:), t_U3, U0(1) + U3(1,:), t_U4, U0(1) + U4(1,:))  
                                                                                                    % nieop�nione
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$U[\frac{cm^3}{s}]$', 'interpreter', 'latex');
title('Przebieg warto�ci sterowania i dop�ywu', 'interpreter', 'latex');
legend('skok 1','skok 2','skok 3','skok 4', 'location', 'east');
grid on;

subplot(2,2,4);
plot(t_U1, U0(2) + U1(2,:));
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$F_D[\frac{cm^3}{s}]$', 'interpreter', 'latex');
title('Przebieg warto�ci zak��cenia', 'interpreter', 'latex');
grid on;


