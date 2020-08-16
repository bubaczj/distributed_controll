%%%%% Regulacja DMC ze wzglêdu na parametry lambda i horyzont sterowania
%%%%% oraz wartoœæ skoku

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

[A, B, X0, U0] = linAB(data.h20, data);                     % Model liniowy

%% Parametry symulacji uk³adu z regulatorem
%h2_des = 
T_sym = [0, 1200];                                          % Czas symulacji regulacji
disruption_type = 0;                                        % 0 - brak dodatkowych zak³óceñ

%% Synteza regulatora DMC na podstawie modelu liniowego
% Parametry regulatora
N       = 1600;                                             % Horyzont predykcji[s]   
Nu      = 40;                                              % Horyzont sterowania[s] 
D       = 1600;                                             % Horyzont dynamiki[s]
lambda  = 2e-1;                                               % Parametr lambda 

% Wyznaczenie odpowiedzi skokowej z rzêdnymi co 1s
step = [1; 0];                                              % Skok wartoœci sterowania
step_val = 0.1 * U0(1);                                     % Wartoœæ skoku jako 10%sterowania w stanie ustalonym (nie ma znaczenia, bo uk³ad jest zlinearyzowany)
T = [0, N];                                                 % Czas symulacji do horyzontu predykcji
options = odeset('MaxStep', 1);                             % Ustawienia solvera
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

%% Symulacja przebiegu regulacji
controlhandle = @DMC;                                       % Regulacja regulatorem DMC
dis_ampl = 15;
rng = 0.5;                                                  % Zakres zadanych wartoœci +/- 50% h0

des_h2 = data.h20 * (1 - rng);                              % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([],[], [], dmc, 0);                               % Reset historii regulatora
sol1 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U1, U1] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 - rng*2/3);                            % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([],[], [], dmc, 0);                               % Reset historii regulatora
sol2 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U2, U2] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 - rng/6);                            % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([], [], [], dmc, 0);                               % Reset historii regulatora
sol3 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U3, U3] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rng/6);                              % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([], [], [], dmc, 0);                               % Reset historii regulatora
sol4 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U4, U4] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rng*2/3);                              % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([], [], [], dmc, 0);                               % Reset historii regulatora
sol5 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U5, U5] = control_delay([], [], 1);                    % Przebieg sterowania

des_h2 = data.h20 * (1 + rng);                              % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle([], [], [], dmc, 0);                               % Reset historii regulatora
sol6 = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle(t, x, des_h2, dmc) + szszsz(t,dis_ampl)), U0, data), T_sym, X0, options);
[~, t_U6, U6] = control_delay([], [], 1);                    % Przebieg sterowania


t_sym = 0:1:T_sym(2);                                       % Wektor czasu symulacji
X1 = (deval(sol1, t_sym))';                                     % Wektory rozwi¹zañ
X2 = (deval(sol2, t_sym))';
X3 = (deval(sol3, t_sym))';
X4 = (deval(sol4, t_sym))';
X5 = (deval(sol5, t_sym))';
X6 = (deval(sol6, t_sym))';


%% Wykres
figure('NumberTitle', 'off', 'Name', 'Odpowiedzi regulatora DMC na skokowe zmiany wartoœci zadanej');
hz = ones(1, (T_sym(2) - T_sym(1)) + 1);
hz1 = hz * ( data.h20 * (1 - rng) );
hz2 = hz * ( data.h20 * (1 - rng*2/3) );
hz3 = hz * ( data.h20 * (1 - rng/6) );
hz4 = hz * ( data.h20 * (1 + rng/6) );
hz5 = hz * ( data.h20 * (1 + rng*2/3) );
hz6 = hz * ( data.h20 * (1 + rng) );
plot(t_sym, X1(:,2),'b',t_sym, X2(:,2),'r',t_sym, X3(:,2),'c',...
    t_sym, X4(:,2),'m',t_sym, X5(:,2),'g',t_sym, X6(:,2),'y',...
    t_sym, hz1,  '--b',t_sym, hz2, '--r',t_sym, hz3, '--c',...
    t_sym, hz4, '--m',t_sym, hz5, '--g',t_sym, hz6, '--y');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$h_2[cm]$', 'interpreter', 'latex');
title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
hold on;
p1 = plot(t_sym, X1(:,2),'b');
p2 = plot(t_sym, hz1,  '--b');
legend([p1, p2], {'poziom p³ynu w zbiorniku drugim', 'wartoœæ zadana'});
% legend('regulacja 1','regulacja 2','regulacja 3','regulacja 4','regulacja 5','regulacja 6',...
%     'wartoœæ zadana 1','wartoœæ zadana 2','wartoœæ zadana 3','wartoœæ zadana 4','wartoœæ zadana 5','wartoœæ zadana 6',...
%     'location', 'northwest');
grid on;

figure('NumberTitle', 'off', 'Name', 'Przebiegi sterowañ odpowiedzi na wymuszenie skokowe');
plot(t_U1,U0(1) + U1(1,:),'b',t_U2,U0(1) +  U2(1,:),'r',t_U3,U0(1) + U3(1,:),'c',...
    t_U4,U0(1) + U4(1,:),'m',t_U5,U0(1) + U5(1,:),'g',t_U6,U0(1) + U6(1,:),'y');
xlabel('$t[s]$', 'interpreter', 'latex');
ylabel('$U[\frac{cm^3}{s}]$', 'interpreter', 'latex');
title('Przebieg wartoœci sterowania', 'interpreter', 'latex');
legend('sterowanie 1','sterowanie 2','sterowanie 3','sterowanie 4','sterowanie 5','sterowanie 6','location', 'northwest');
grid on;

figure();
plot(t_U1,U1(2,:),'b',t_U2,U2(2,:),'r',t_U3, U3(2,:),'c',...
    t_U4,U4(2,:),'m',t_U5,U5(2,:),'g',t_U6,U6(2,:),'y');
grid on;
