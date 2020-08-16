%%%%% Porównanie ró¿nych regulatorów analitycznych

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




% modele s¹ wyznaczane w n punktach
n = 5;                                                  % liczba punktów pracy
rng = 0.5;                                              % zakres rozmywania
q = 1.3;                                                % krok rozmywania

% parametry do regulacji rozmytej, nierównomiernej
[A_X_q, B_X_q, X0_X_q, U0_X_q, range_X_q, activ_X_q] = linAB_TS_X_q(data, n, rng, q);       % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów
% parametry uk³adu liniowego
[A,B,X0,U0] = linAB(data.h20, data);                 % wektor sterowania i stanu w zadanym punkcie pracy
% parametry do regulacji rozmytej równomiernie
[A_X, B_X, X0_X, U0_X, range_X] = linAB_TS_X(data, n, rng);


%% Punkty równowagi dla kolejnych regulatorów sk³adowych DMC równomiernego
dmc_TS.U0_X = U0_X;
dmc_TS.X0_X = X0_X;
dmc_TS.U0 = U0;
dmc_TS.X0 = X0;
dmc_TS.range_X = range_X;

%% Punkty równowagi dla kolejnych regulatorów sk³adowych DMC nierównomiernego
dmc_TS_q.U0_X = U0_X_q;
dmc_TS_q.X0_X = X0_X_q;
dmc_TS_q.X0 = X0;
dmc_TS_q.U0 = U0;
dmc_TS_q.range_X = range_X_q;


%% Parametry symulacji
T_sym = [0, 3000];                                          % Czas symulacji regulacji
disruption_type = 0;                                        % 0 - brak dodatkowych zak³óceñ

%% Synteza regulatorów

% Parametry regulatorów
N       = 1600;                                             % Horyzont predykcji[s]   
Nu      = 20;                                              % Horyzont sterowania[s] 
D       = 1600;                                             % Horyzont dynamiki[s]
lambda  = 1e-1;                                               % Parametr lambda 

%% klasyczny regulator DMC 

% Wyznaczenie odpowiedzi skokowej z rzêdnymi co 1s
step = [1; 0];                                              % Skok wartoœci sterowania
step_val = 0.1 * U0(1);                                     % Wartoœæ skoku jako 10%sterowania w stanie ustalonym (nie ma znaczenia, bo uk³ad jest zlinearyzowany)
T = [0, N];                                                 % Czas symulacji do horyzontu predykcji
options = odeset('MaxStep', 1, 'Refine',1);                 % Ustawienia solvera
controlhandle = step * step_val;                            % Sterowanie skokiem
control_delay([], [], 0);                                   % Reset funkcji generuj¹cej opóŸnienie

step_sol = ode45(@(t, x)plant_lin(t, x, control_delay(t, controlhandle), X0, A, B), T, X0, options);
t = 0 : 1 : N;                                              % Chwile rozwi¹zañ
X = (deval(step_sol, t))';                                  % Wartoœci rozwi¹zañ

s = X(:,2);                                                 % OdpowiedŸ skokowa
s = (s - s(1)) / (s(end) - s(1));                           % Znormalizowane rzêdne odpowiedzi skokowej

Mp = zeros(N, D-1);                                         % Macierz Mp jest wymiaru NxD-1
M = zeros(N, Nu);                                           % Macierz M jest wymiaru NxNu
L = lambda * eye(Nu);                                       % Macierz lambda * I

% Wyznaczanie macierzy M i Mp
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

% Wyznaczanie wektora K1
K = (M' * M + L) \ M';
K1 = K(1,:);

% Struktura zawieraj¹ca wszystkie dane potrzebne do realizacji regulacji
dmc.ku = K1 * Mp;
dmc.ke = sum(K1);
dmc.D = D;

%% rozmyty regulator DMC (równomiernie)

X_X = zeros(N + 1, 2, n);
t = 0 : 1 : N;                                              % Chwile rozwi¹zañ
for i = 1 : n
    control_delay([], [], 0);
    step_sol_X = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_X(:,i), A_X(:,:,i), B_X(:,:,i)), T, X0_X(:,i), options);
    X_X(:,:,i) = (deval(step_sol_X, t))';                   % Wartoœci rozwi¹zañ
end             


lambda = [5 1 0.5 0.05 0.001];
dmc_TS.n = n;
dmc_TS.D = D;
for reg = 1 : n
    s_X = X_X(:,2, reg);
    s_X = (s_X - s_X(1)) / (s_X(end) - s_X(1));
    Mp_X = zeros(N, D-1);   
    M_X = zeros(N, Nu);
    % Wyznaczanie macierzy M i Mp
    for i = 1 : N
        for j = 1 : D - 1
            Mp_X(i, j) = s_X(min(i + j, D)) - s_X(j);
        end
        for j = 1 : Nu
            if i - j + 1 > 0
                M_X(i, j) = s_X(i - j + 1);
            end    
        end
    end

    L = lambda(reg)*eye(Nu);
    % Wyznaczanie wektora K1
    K_X = (M_X' * M_X + L) \ M_X';
    K1_X = K_X(1,:);

    % Struktura zawieraj¹ca wszystkie dane potrzebne do realizacji regulacji
    dmc_TS.ku_X(:,reg) = K1_X * Mp_X;
    dmc_TS.ke_X(reg) = sum(K1_X);
end

%% nierównomiernie rozmyty regulator DMC

for i = 1 : n
    control_delay([], [], 0);
    step_sol_X = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_X_q(:,i), A_X_q(:,:,i), B_X_q(:,:,i)), T, X0_X_q(:,i), options);
    X_X(:,:,i) = (deval(step_sol_X, t))';                   % Wartoœci rozwi¹zañ
end            
for reg = 1 : n 
    s_X = X_X(:,2, reg);
    s_X = (s_X - s_X(1)) / (s_X(end) - s_X(1));
    Mp_X = zeros(N, D-1);   
    M_X = zeros(N, Nu);
    % Wyznaczanie macierzy M i Mp
    for i = 1 : N
        for j = 1 : D - 1
            Mp_X(i, j) = s_X(min(i + j, D)) - s_X(j);
        end
        for j = 1 : Nu
            if i - j + 1 > 0
                M_X(i, j) = s_X(i - j + 1);
            end    
        end
    end

    L = lambda(reg)*eye(Nu);
    % Wyznaczanie wektora K1
    K_X = (M_X' * M_X + L) \ M_X';
    K1_X = K_X(1,:);

    % Struktura zawieraj¹ca wszystkie dane potrzebne do realizacji regulacji
    dmc_TS_q.n = n;
    dmc_TS_q.D = D;
    dmc_TS_q.ku_X(:,reg) = K1_X * Mp_X;
    dmc_TS_q.ke_X(reg) = sum(K1_X);
end

%% symulacja 
controlhandle_X_q = @DMC_TS_X_q;
controlhandle_X = @DMC_TS_X;
controlhandle = @DMC;

des_h2 = @(t)step_steering(t);

control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([],[], [], dmc, [] );                               % Reset historii regulatora
sol_DMC = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2(t), dmc)), U0, data), T_sym, X0, options);
[~, t_DMC, U_DMC] = control_delay([], [], 1);                    % Przebieg sterowania


control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle_X([],[], [], dmc_TS, [] );                               % Reset historii regulatora
sol_DMC_TS = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2(t), dmc_TS)), U0, data), T_sym, X0, options);
[~, t_DMC_TS, U_DMC_TS] = control_delay([], [], 1);                    % Przebieg sterowania

control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle_X_q([],[], [], dmc_TS_q,[], [] );                               % Reset historii regulatora
sol_DMC_TS_q = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X_q(t, x, des_h2(t), dmc_TS_q, activ_X_q)), U0, data), T_sym, X0, options);
[~, t_DMC_TS_q, U_DMC_TS_q] = control_delay([], [], 1);                    % Przebieg sterowania


%solutions
t = T_sym(1) : T_sym(2);
X_DMC = (deval(sol_DMC, t))';
X_DMC_TS = (deval(sol_DMC_TS, t))';
X_DMC_TS_q = (deval(sol_DMC_TS_q, t))';


%% wykres

figure('NumberTitle', 'off', 'Name', 'ró¿ne regulatory');
[~,s] = size(t);
steer = zeros(1,s);
for i = 1 : s
    steer(i) = step_steering(i);
end
plot(t, X_DMC(:,2),t, X_DMC_TS(:,2), t, X_DMC_TS_q(:,2),t, steer );
legend('DMC', 'DMC rozmyty regularnie', 'DMC rozmyty nieregularnie', 'wartoœæ zadana');
grid on;





