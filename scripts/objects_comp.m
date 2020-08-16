%%%%% Porównanie odpowiedzi modelu rozmytego wzglêdem h, F i obiektu
%%%%% nieliniowego

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
[A_X2, B_X2, X0_X2, U0_X2, range_X2] = linAB_TS_X(data, 2, rng); % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów po h
[A_X3, B_X3, X0_X3, U0_X3, range_X3] = linAB_TS_X(data, 3, rng);
[A_X4, B_X4, X0_X4, U0_X4, range_X4] = linAB_TS_X(data, 4, rng);
[A_X5, B_X5, X0_X5, U0_X5, range_X5] = linAB_TS_X(data, 5, rng);

[A,B,X0, U0] = linAB(data.h20, data);                     % wektor sterowania i stanu w zadanym punkcie pracy

%% Parametry symulacji
step = [1; 0];                                  % Skok jedynie wartoœci steruj¹cej
step_val = 40*[-1,-0.5, 0.5, 1];                % Przyrost wartoœci steruj¹cej F = F1 + FD + step_val 
T = [0, 600];                                  % Czas symulacji
controlhandle = step * step_val;                % Obs³uga sterowania
options = odeset('MaxStep', 1, 'Refine', 1);    % Ustawienia solvera


%% Symulacja
% Obiekt nieliniowy
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution1_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,1)), U0, data), T, X0, options);
[~, t_U1_nl, U1_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution2_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,2)), U0, data), T, X0, options);
[~, t_U2_nl, U2_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution3_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,3)), U0, data), T, X0, options);
[~, t_U3_nl, U3_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution4_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,4)), U0, data), T, X0, options);
[~, t_U4_nl, U4_nl] = control_delay([], [], 1);         % Przebieg sterowania

% Obiekt rozmyty wzglêdem h 2
n = 2;
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_X2 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,1)), X0_X2, U0_X2, A_X2, B_X2, n, U0, range_X2), T, X0, options);
[~, t_U1_X2, U1_X2] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_X2 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,2)), X0_X2, U0_X2, A_X2, B_X2, n, U0, range_X2), T, X0, options);
[~, t_U2_X2, U2_X2] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_X2 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,3)), X0_X2, U0_X2, A_X2, B_X2, n, U0, range_X2), T, X0, options);
[~, t_U3_X2, U3_X2] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_X2 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,4)), X0_X2, U0_X2, A_X2, B_X2, n, U0, range_X2), T, X0, options);
[~, t_U4_X2, U4_X2] = control_delay([], [], 1);         % Przebieg sterowania

% Obiekt rozmyty wzglêdem h 3
n = 3;
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_X3 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,1)), X0_X3, U0_X3, A_X3, B_X3, n, U0, range_X3), T, X0, options);
[~, t_U1_X3, U1_X3] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_X3 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,2)), X0_X3, U0_X3, A_X3, B_X3, n, U0, range_X3), T, X0, options);
[~, t_U2_X3, U2_X3] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_X3 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,3)), X0_X3, U0_X3, A_X3, B_X3, n, U0, range_X3), T, X0, options);
[~, t_U3_X3, U3_X3] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_X3 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,4)), X0_X3, U0_X3, A_X3, B_X3, n, U0, range_X3), T, X0, options);
[~, t_U4_X3, U4_X3] = control_delay([], [], 1);         % Przebieg sterowania


% Obiekt rozmyty wzglêdem h 4
n = 4;
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_X4 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,1)), X0_X4, U0_X4, A_X4, B_X4, n, U0, range_X4), T, X0, options);
[~, t_U1_X4, U1_X4] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_X4 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,2)), X0_X4, U0_X4, A_X4, B_X4, n, U0, range_X4), T, X0, options);
[~, t_U2_X4, U2_X4] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_X4 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,3)), X0_X4, U0_X4, A_X4, B_X4, n, U0, range_X4), T, X0, options);
[~, t_U3_X4, U3_X4] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_X4 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,4)), X0_X4, U0_X4, A_X4, B_X4, n, U0, range_X4), T, X0, options);
[~, t_U4_X4, U4_X4] = control_delay([], [], 1);         % Przebieg sterowania


% Obiekt rozmyty wzglêdem h 5
n = 5;
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_X5 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,1)), X0_X5, U0_X5, A_X5, B_X5, n, U0, range_X5), T, X0, options);
[~, t_U1_X5, U1_X5] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_X5 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,2)), X0_X5, U0_X5, A_X5, B_X5, n, U0, range_X5), T, X0, options);
[~, t_U2_X5, U2_X5] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_X5 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,3)), X0_X5, U0_X5, A_X5, B_X5, n, U0, range_X5), T, X0, options);
[~, t_U3_X5, U3_X5] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_X5 = ode45(@(t, x)plant_TS_X(t, x, control_delay(t, controlhandle(:,4)), X0_X5, U0_X5, A_X5, B_X5, n, U0, range_X5), T, X0, options);
[~, t_U4_X5, U4_X5] = control_delay([], [], 1);         % Przebieg sterowania


% Obiekt liniowy
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_lin = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(:,1)), X0,  A, B), T, X0, options);
[~, t_U1_lin, U1_lin] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_lin = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(:,2)), X0,  A, B), T, X0, options);
[~, t_U2_lin, U2_lin] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_lin = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(:,3)), X0,  A, B), T, X0, options);
[~, t_U3_lin, U3_lin] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_lin = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle(:,4)), X0,  A, B), T, X0, options);
[~, t_U4_lin, U4_lin] = control_delay([], [], 1);         % Przebieg sterowania




t = T(1) : T(2);                                % Wektor czasu symulacji
X1_nl = (deval(solution1_nl, t))';
X2_nl = (deval(solution2_nl, t))';
X3_nl = (deval(solution3_nl, t))';
X4_nl = (deval(solution4_nl, t))';

X1_X2 = (deval(solution1_X2, t))';
X2_X2 = (deval(solution2_X2, t))';
X3_X2 = (deval(solution3_X2, t))';
X4_X2 = (deval(solution4_X2, t))';

X1_X3 = (deval(solution1_X3, t))';
X2_X3 = (deval(solution2_X3, t))';
X3_X3 = (deval(solution3_X3, t))';
X4_X3 = (deval(solution4_X3, t))';

X1_X4 = (deval(solution1_X4, t))';
X2_X4 = (deval(solution2_X4, t))';
X3_X4 = (deval(solution3_X4, t))';
X4_X4 = (deval(solution4_X4, t))';

X1_X5 = (deval(solution1_X5, t))';
X2_X5 = (deval(solution2_X5, t))';
X3_X5 = (deval(solution3_X5, t))';
X4_X5 = (deval(solution4_X5, t))';

X1_lin = (deval(solution1_lin, t))';
X2_lin = (deval(solution2_lin, t))';
X3_lin = (deval(solution3_lin, t))';
X4_lin = (deval(solution4_lin, t))';


%% Wykresy
figure('NumberTitle', 'off', 'Name', 'Ró¿ne obiekty');
hold on;
%plot(t, X1_nl(:,2),'b',t, X2_nl(:,2),'b',t, X3_nl(:,2),'b',t, X4_nl(:,2),'b');
R1 = X4_X2(:,2)-X4_nl(:,2);   %% zmieniaæ pierwsze cyfry 1-4
 plot(t, log(abs(R1) + 1),'r');
s1 = cumsum(abs(R1));
R2 = X4_X3(:,2)-X4_nl(:,2);
 plot(t, log(abs(R2) + 1),'g');
s2 = cumsum(abs(R2));
R3 = X4_X4(:,2)-X4_nl(:,2);
 plot(t, log(abs(R3) + 1),'black');
s3 = cumsum(abs(R3));
R4 = X4_X5(:,2)-X4_nl(:,2);
 plot(t, log(abs(R4) + 1),'y');
s4 = cumsum(abs(R4));
R5 = X4_lin(:,2)-X4_nl(:,2);
 plot(t, log(abs(R5) + 1), 'm');
s5 = cumsum(abs(R5));
grid on;
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2 - h_{nonlin}[cm]$', 'interpreter', 'latex');
%title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
% legend('skok 1, rozmycie h','skok 2, rozmycie h','skok 3, rozmycie h','skok 4, rozmycie h',...
%     'skok 1, rozmycie F','skok 2, rozmycie F','skok 3, rozmycie F','skok 4, rozmycie F',...
%     'skok 1, obiekt nieliniowy','skok 2, obiekt nieliniowy','skok 3, obiekt nieliniowy','skok 4, nieobiekt liniowy',...
%     'location', 'east');
%legend('off');
legend('n=2','n=3', 'n=4', 'n=5', 'n=1');
grid on;
figure();
% plot(t, s1, t, s2, t, s3, t, s4, t, s5);
% legend('n=2','n=3', 'n=4', 'n=5', 'n=1');
% grid on;
% subplot(1,2,2);
% plot(t, X1_X(:,2),':b',t, X2_X(:,2),':r',t, X3_X(:,2),':c',t, X4_X(:,2),':m',...
%     t, X1_U(:,2),'--b',t, X2_U(:,2),'--r',t, X3_U(:,2),'--c',t, X4_U(:,2),'--m',...
%     t, X1_nl(:,2),  'b',t, X2_nl(:,2), 'r',t, X3_nl(:,2), 'c',t, X4_nl(:,2), 'm');
% xlabel('$t[s]$', 'interpreter', 'latex');
% ylabel('$h_2[cm]$', 'interpreter', 'latex');
% title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
% legend('skok 1, rozmycie h','skok 2, rozmycie h','skok 3, rozmycie h','skok 4, rozmycie h',...
%     'skok 1, rozmycie F','skok 2, rozmycie F','skok 3, rozmycie F','skok 4, rozmycie F',...
%     'skok 1, obiekt nieliniowy','skok 2, obiekt nieliniowy','skok 3, obiekt nieliniowy','skok 4, nieobiekt liniowy',...
%     'location', 'east');
% %legend('off');
% grid on;
