%%%%% Regulacja rozmytym r�wnomiernie DMC ze wzgl�du na parametry lambda i horyzont sterowania
%%%%% oraz warto�� skoku

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
% skrajne modele s� wyznaczone wok� punkt�w p_pracy_min/max = p0 -/+ rng * p0
n = 5;                                                      % liczba punkt�w pracy
rng = 0.5;                                                  % zakres rozmywania
[A_X, B_X, X0_X, U0_X, range_X] = linAB_TS_X(data, n, rng); % A, B, X0, U0 - parametry n zlinearyzowanych uk�ad�w po h
[A_U, B_U, X0_U, U0_U, range_U] = linAB_TS_U(data, n, rng); % A, B, X0, U0 - parametry n zlinearyzowanych uk�ad�w po F
[~,~,X0, U0] = linAB(data.h20, data);                     % wektor sterowania i stanu w zadanym punkcie pracy

%% Punkty r�wnowagi dla kolejnych regulator�w sk�adowych
dmc_TS.U0_U = U0_U;
dmc_TS.U0_X = U0_X;
dmc_TS.X0_U = X0_U;
dmc_TS.X0_X = X0_X;
dmc_TS.U0 = U0;
dmc_TS.X0 = X0;
dmc_TS.range_U = range_U;
dmc_TS.range_X = range_X;


%% Parametry symulacji uk�adu z regulatorem
%h2_des = 
T_sym = [0, 1200];                                          % Czas symulacji regulacji
dis_ampl = 15;                                             % amplituda zak�uce�

%% Synteza regulator�w DMC na podstawie modelu liniowego
% Parametry regulatora
N       = 1600;                                             % Horyzont predykcji[s]   
Nu      = 40;                                              % Horyzont sterowania[s] 
D       = 1600;                                             % Horyzont dynamiki[s]
lambda  = 2e-1;                                             % Parametr lambda 

% Wyznaczenie odpowiedzi skokowej z rz�dnymi co 1s
step = [1; 0];                                              % Skok warto�ci sterowania                        
step_val = 0.1 * U0(1);                                   % Warto�� skoku jako 10%sterowania w stanie ustalonym (nie ma znaczenia, bo uk�ad jest zlinearyzowany)
T = [0, N];                                                 % Czas symulacji do horyzontu predykcji
options = odeset('MaxStep', 1);                             % Ustawienia solvera
controlhandle = step * step_val;                            % Sterowanie skokiem
control_delay([], [], 0);                                   % Reset funkcji generuj�cej op�nienie

X_U = zeros(N + 1, 2, n);
X_X = zeros(N + 1, 2, n);
t = 0 : 1 : N;                                              % Chwile rozwi�za�
for i = 1 : n
    control_delay([], [], 0);
    step_sol_U = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_U(:,i), A_U(:,:,i), B_U(:,:,i)), T, X0_U(:,i), options);
    control_delay([], [], 0);
    step_sol_X = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_X(:,i), A_X(:,:,i), B_X(:,:,i)), T, X0_X(:,i), options);
    X_U(:,:,i) = (deval(step_sol_U, t))';                   % Warto�ci rozwi�za�
    X_X(:,:,i) = (deval(step_sol_X, t))';                   % Warto�ci rozwi�za�
end             

Mp_X = zeros(N, D-1);                                         % Macierz Mp jest wymiaru NxD-1
Mp_U = zeros(N, D-1);
M_X = zeros(N, Nu);                                           % Macierz M jest wymiaru NxNu
M_U = zeros(N, Nu);
L = lambda * eye(Nu);                                       % Macierz lambda * I

dmc_TS.n = n;
dmc_TS.D = D;
for reg = 1 : n
    s_U = X_U(:,2, reg);                                        % Odpowied� skokowa regulatora reg
    s_U = (s_U - s_U(1)) / (s_U(end) - s_U(1));                           % Znormalizowane rz�dne odpowiedzi skokowej
    s_X = X_X(:,2, reg);
    s_X = (s_X - s_X(1)) / (s_X(end) - s_X(1));
    % Wyznaczanie macierzy M i Mp
    for i = 1 : N
        for j = 1 : D - 1
            Mp_X(i, j) = s_X(min(i + j, D)) - s_X(j);
            Mp_U(i, j) = s_U(min(i + j, D)) - s_U(j);
        end
        for j = 1 : Nu
            if i - j + 1 > 0
                M_X(i, j) = s_X(i - j + 1);
                M_U(i, j) = s_U(i - j + 1);
            end    
        end
    end

    % Wyznaczanie wektora K1
    K_X = (M_X' * M_X + L) \ M_X';
    K_U = (M_U' * M_U + L) \ M_U';
    K1_X = K_X(1,:);
    K1_U = K_U(1,:);

    % Struktura zawieraj�ca wszystkie dane potrzebne do realizacji regulacji
    dmc_TS.ku_X(:,reg) = K1_X * Mp_X;
    dmc_TS.ku_U(:,reg) = K1_U * Mp_U;
    dmc_TS.ke_X(reg) = sum(K1_X);
    dmc_TS.ke_U(reg) = sum(K1_U);
end

%% Symulacja przebiegu regulacji
controlhandle_U = @DMC_TS_U;                                       % Regulacja regulatorem DMC
controlhandle_X = @DMC_TS_X;

rang = 0.5;                                                  % Zakres zadanych warto�ci +/- 50% h0

des_h2 = data.h20 * (1 - rang);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol1_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U1_U, U1_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol1_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U1_X, U1_X] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 - rang*2/3);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol2_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U2_U, U2_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol2_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U2_X, U2_X] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 - rang/6);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol3_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U3_U, U3_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol3_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U3_X, U3_X] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rang/6);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol4_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U4_U, U4_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol4_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U4_X, U4_X] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rang*2/3);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol5_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U5_U, U5_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol5_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U5_X, U5_X] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rang);                              % ��dana warto�� wysoko�ci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_U([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol6_U = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_U(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U6_U, U6_U] = control_delay([], [], 1);                    % Przebieg sterowania
control_delay([], [], 0);                                   % Reset funkcji realizuj�cej op�nienie
controlhandle_X([],[], [], dmc_TS, 0);                               % Reset historii regulatora
sol6_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS) + szszsz(t, dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U6_X, U6_X] = control_delay([], [], 1);                    % Przebieg sterowania


t_sym = 0:1:T_sym(2);                                       % Wektor czasu symulacji
X1_X = (deval(sol1_X, t_sym))';                                     % Wektory rozwi�za�
X2_X = (deval(sol2_X, t_sym))';
X3_X = (deval(sol3_X, t_sym))';
X4_X = (deval(sol4_X, t_sym))';
X5_X = (deval(sol5_X, t_sym))';
X6_X = (deval(sol6_X, t_sym))';

X1_U = (deval(sol1_U, t_sym))';                                     % Wektory rozwi�za�
X2_U = (deval(sol2_U, t_sym))';
X3_U = (deval(sol3_U, t_sym))';
X4_U = (deval(sol4_U, t_sym))';
X5_U = (deval(sol5_U, t_sym))';
X6_U = (deval(sol6_U, t_sym))';


%% Wykres
figure('NumberTitle', 'off', 'Name', 'Odpowiedzi regulator�w DMC rozmytych na skokowe zmiany warto�ci zadanej');
hz = ones(1, (T_sym(2) - T_sym(1)) + 1);
hz1 = hz * ( data.h20 * (1 - rang) );
hz2 = hz * ( data.h20 * (1 - rang*2/3) );
hz3 = hz * ( data.h20 * (1 - rang/6) );
hz4 = hz * ( data.h20 * (1 + rang/6) );
hz5 = hz * ( data.h20 * (1 + rang*2/3) );
hz6 = hz * ( data.h20 * (1 + rang) );
plot(t_sym, X1_X(:,2),':b',t_sym, X2_X(:,2),':r',t_sym, X3_X(:,2),':c',...
    t_sym, X4_X(:,2),':m',t_sym, X5_X(:,2),':g',t_sym, X6_X(:,2),':black',...
    t_sym, X1_U(:,2),'--b',t_sym, X2_U(:,2),'--r',t_sym, X3_U(:,2),'--c',...
    t_sym, X4_U(:,2),'--m',t_sym, X5_U(:,2),'--g',t_sym, X6_U(:,2),'--black',...
    t_sym, hz1,  'b',t_sym, hz2, 'r',t_sym, hz3, 'c',...
    t_sym, hz4, 'm',t_sym, hz5, 'g',t_sym, hz6, 'black');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2[cm]$', 'interpreter', 'latex');
title('Przebieg wysoko�ci p�ynu w drugim zbiorniku', 'interpreter', 'latex');
hold on;
p1 = plot(t_sym, X1_X(:,2),':b');
p2 = plot(t_sym, X1_U(:,2),'--b');
p3 = plot(t_sym, hz1,  'b');
legend([p1, p2, p3], {'rozmycie po stanie','rozmycie po sterowaniu', 'warto�� zadana'});
% legend('DMC_TS_X 1','DMC_TS_X 2','DMC_TS_X 3','DMC_TS_X 4','DMC_TS_X 5','DMC_TS_X 6',...
%     'DMC_TS_U 1','DMC_TS_U 2','DMC_TS_U 3','DMC_TS_U 4','DMC_TS_U 5','DMC_TS_U 6',...
%     'warto�� zadana 1','warto�� zadana 2','warto�� zadana 3','warto�� zadana 4','warto�� zadana 5','warto�� zadana 6',...
%     'location', 'northwest');
grid on;

figure('NumberTitle', 'off', 'Name', 'Przebiegi sterowa� odpowiedzi na wymuszenie skokowe');
plot(t_U1_X,U0(1) + U1_X(1,:),':b',t_U2_X,U0(1) +  U2_X(1,:),':r',t_U3_X,U0(1) + U3_X(1,:),':c',...
    t_U4_X,U0(1) + U4_X(1,:),':m',t_U5_X,U0(1) + U5_X(1,:),':g',t_U6_X,U0(1) + U6_X(1,:),':black',...
    t_U1_U,U0(1) + U1_U(1,:),'--b',t_U2_U,U0(1) +  U2_U(1,:),'--r',t_U3_U,U0(1) + U3_U(1,:),'--c',...
    t_U4_U,U0(1) + U4_U(1,:),'--m',t_U5_U,U0(1) + U5_U(1,:),'--g',t_U6_U,U0(1) + U6_U(1,:),'--black');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$U[\frac{cm^3}{s}]$', 'interpreter', 'latex');
title('Przebieg warto�ci sterowania', 'interpreter', 'latex');
% legend('DMC_TS_X 1','DMC_TS_X 2','DMC_TS_X 3','DMC_TS_X 4','DMC_TS_X 5','DMC_TS_X 6',...
%     'DMC_TS_U 1','DMC_TS_U 2','DMC_TS_U 3','DMC_TS_U 4','DMC_TS_U 5','DMC_TS_U 6','location', 'northwest');
grid on;

