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
q = 1.3;
[A_X, B_X, X0_X, U0_X, range_X, activ_X] = linAB_TS_X_q(data, n, rng, q); % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów po h
[A_U, B_U, X0_U, U0_U, range_U, activ_U] = linAB_TS_U_q(data, n, rng, q); % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów po F
[~,~,X0_, U0_] = linAB(data.h20, data);                     % wektor sterowania i stanu w zadanym punkcie pracy

%% Parametry symulacji
step = [1; 0];                                  % Skok jedynie wartoœci steruj¹cej
step_val = 40*[-1,-0.5, 0.5, 1];                % Przyrost wartoœci steruj¹cej F = F1 + FD + step_val 
T = [0, 1400];                                  % Czas symulacji
controlhandle = step * step_val;                % Obs³uga sterowania
options = odeset('MaxStep', 1, 'Refine', 1);    % Ustawienia solvera


%% Symulacja
% Obiekt nieliniowy
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution1_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,1)), U0_, data), T, X0_, options);
[~, t_U1_nl, U1_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution2_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,2)), U0_, data), T, X0_, options);
[~, t_U2_nl, U2_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution3_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,3)), U0_, data), T, X0_, options);
[~, t_U3_nl, U3_nl] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                         % Reset funkcji realizuj¹cej opóŸnienie
solution4_nl = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(:,4)), U0_, data), T, X0_, options);
[~, t_U4_nl, U4_nl] = control_delay([], [], 1);         % Przebieg sterowania

% Obiekt rozmyty wzglêdem h
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_X = ode45(@(t, x)plant_TS_X_q(t, x, control_delay(t, controlhandle(:,1)), X0_X, U0_X, A_X, B_X, n, U0_, activ_X), T, X0_, options);
[~, t_U1_X, U1_X] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_X = ode45(@(t, x)plant_TS_X_q(t, x, control_delay(t, controlhandle(:,2)), X0_X, U0_X, A_X, B_X, n, U0_, activ_X), T, X0_, options);
[~, t_U2_X, U2_X] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_X = ode45(@(t, x)plant_TS_X_q(t, x, control_delay(t, controlhandle(:,3)), X0_X, U0_X, A_X, B_X, n, U0_, activ_X), T, X0_, options);
[~, t_U3_X, U3_X] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_X = ode45(@(t, x)plant_TS_X_q(t, x, control_delay(t, controlhandle(:,4)), X0_X, U0_X, A_X, B_X, n, U0_, activ_X), T, X0_, options);
[~, t_U4_X, U4_X] = control_delay([], [], 1);         % Przebieg sterowania


% Obiekt rozmyty wzglêdem F
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution1_U = ode45(@(t, x)plant_TS_U_q(t, x, control_delay(t, controlhandle(:,1)), X0_U, U0_U, A_U, B_U, n, U0_, activ_U), T, X0_, options);
[~, t_U1_U, U1_U] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution2_U = ode45(@(t, x)plant_TS_U_q(t, x, control_delay(t, controlhandle(:,2)), X0_U, U0_U, A_U, B_U, n, U0_, activ_U), T, X0_, options);
[~, t_U2_U, U2_U] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution3_U = ode45(@(t, x)plant_TS_U_q(t, x, control_delay(t, controlhandle(:,3)), X0_U, U0_U, A_U, B_U, n, U0_, activ_U), T, X0_, options);
[~, t_U3_U, U3_U] = control_delay([], [], 1);         % Przebieg sterowania
control_delay([], [], 0);                       % Reset funkcji realizuj¹cej opóŸnienie
solution4_U = ode45(@(t, x)plant_TS_U_q(t, x, control_delay(t, controlhandle(:,4)), X0_U, U0_U, A_U, B_U, n, U0_, activ_U), T, X0_, options);
[~, t_U4_U, U4_U] = control_delay([], [], 1);         % Przebieg sterowania





t = T(1) : T(2);                                % Wektor czasu symulacji
X1_nl = (deval(solution1_nl, t))';              % Wektor rozwi¹zañ obiektu nieliniowego
X2_nl = (deval(solution2_nl, t))';
X3_nl = (deval(solution3_nl, t))';
X4_nl = (deval(solution4_nl, t))';

X1_X = (deval(solution1_X, t))';              % Wektor rozwi¹zañ obiektu rozmytego po h
X2_X = (deval(solution2_X, t))';
X3_X = (deval(solution3_X, t))';
X4_X = (deval(solution4_X, t))';

X1_U = (deval(solution1_U, t))';              % Wektor rozwi¹zañ obiektu rozmytego po F
X2_U = (deval(solution2_U, t))';
X3_U = (deval(solution3_U, t))';
X4_U = (deval(solution4_U, t))';

%% Wykresy
figure('NumberTitle', 'off', 'Name', 'Porównanie odpowiedzi na skok obiektu rozmytego z nieliniowym');
subplot(1,2,1);
plot(t, X1_X(:,1),':b',t, X2_X(:,1),':r',t, X3_X(:,1),':c',t, X4_X(:,1),':m',...
    t, X1_U(:,1),'--b',t, X2_U(:,1),'--r',t, X3_U(:,1),'--c',t, X4_U(:,1),'--m',...
    t, X1_nl(:,1),  'b',t, X2_nl(:,1), 'r',t, X3_nl(:,1), 'c',t, X4_nl(:,1), 'm');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_1[cm]$', 'interpreter', 'latex');
title('Przebieg wysokoœci p³ynu w pierwszym zbiorniku', 'interpreter', 'latex');
legend('skok 1, rozmycie h','skok 2, rozmycie h','skok 3, rozmycie h','skok 4, rozmycie h',...
    'skok 1, rozmycie F','skok 2, rozmycie F','skok 3, rozmycie F','skok 4, rozmycie F',...
    'skok 1, obiekt nieliniowy','skok 2, obiekt nieliniowy','skok 3, obiekt nieliniowy','skok 4, obiekt nieliniowy',...
    'location', 'east');
%legend('off');
grid on;

subplot(1,2,2);
plot(t, X1_X(:,2),':b',t, X2_X(:,2),':r',t, X3_X(:,2),':c',t, X4_X(:,2),':m',...
    t, X1_U(:,2),'--b',t, X2_U(:,2),'--r',t, X3_U(:,2),'--c',t, X4_U(:,2),'--m',...
    t, X1_nl(:,2),  'b',t, X2_nl(:,2), 'r',t, X3_nl(:,2), 'c',t, X4_nl(:,2), 'm');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2[cm]$', 'interpreter', 'latex');
title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
legend('skok 1, rozmycie h','skok 2, rozmycie h','skok 3, rozmycie h','skok 4, rozmycie h',...
    'skok 1, rozmycie F','skok 2, rozmycie F','skok 3, rozmycie F','skok 4, rozmycie F',...
    'skok 1, obiekt nieliniowy','skok 2, obiekt nieliniowy','skok 3, obiekt nieliniowy','skok 4, obiekt nieliniowy',...
    'location', 'east');
%legend('off');
grid on;




figure('NumberTitle', 'off', 'Name', 'Funkcje aktywacji');
h = 10:0.1:70;
[~,s] = size(h);
v = zeros(n,s);
for j = 1:s
    total_act = 0;
    for i = 1 : n
        if i == 1
            total_act = total_act + activation(activ_X(1, i), activ_X(2,i), h(j), -inf);
        elseif i == n
            total_act = total_act + activation(activ_X(1, i), activ_X(2,i), h(j), inf);
        else
            total_act = total_act + activation(activ_X(1, i), activ_X(2,i), h(j));
        end
    end
     for i = 1 : n
        if i == 1
            v(i,j) = activation(activ_X(1, i), activ_X(2,i), h(j), -inf)/total_act;
        elseif i == n
            v(i,j) = activation(activ_X(1, i), activ_X(2,i), h(j), inf)/total_act;
        else
            v(i,j) = activation(activ_X(1, i), activ_X(2,i), h(j))/total_act;
        end
    end   
end
for i = 1 : n
   plot(h, v(i,:));
   hold on;
end
grid on;