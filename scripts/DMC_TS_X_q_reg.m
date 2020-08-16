%%%%% Regulacja rozmytym nierównomiernie DMC ze wzglêdu na parametry lambda i horyzont sterowania
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

% modele s¹ wyznaczane w n punktach
n = 5;                                                  % liczba punktów pracy
rng = 0.5;                                              % zakres rozmywania
q = 1.15;                                                % krok rozmywania
[A_X, B_X, X0_X, U0_X, range_X, activ_X] = linAB_TS_X_q(data, n, rng, q);       % A, B, X0, U0 - parametry n zlinearyzowanych uk³adów
[A_X_nq, B_X_nq, X0_X_nq, U0_X_nq, range_X_nq] = linAB_TS_X(data, n, rng);
[~,~,X0_, U0_] = linAB(data.h20, data);                 % wektor sterowania i stanu w zadanym punkcie pracy

%% Punkty równowagi dla kolejnych regulatorów sk³adowych
dmc_TS.U0_X = U0_X;
dmc_TS.X0_X = X0_X;
dmc_TS.X0 = X0_;
dmc_TS.range_X = range_X;

dmc_TS_nq.U0_X = U0_X_nq;
dmc_TS_nq.X0_X = X0_X_nq;
dmc_TS_nq.X0 = X0_;
dmc_TS_nq.range_X = range_X_nq;


%% Parametry symulacji uk³adu z regulatorem
%h2_des = 
T_sym = [0, 1200];                                          % Czas symulacji regulacji
dis_ampl = 0;

%% Synteza regulatorów DMC na podstawie modelu liniowego
% Parametry regulatora
N       = 1600;                                             % Horyzont predykcji[s]   
Nu      = 15;                                              % Horyzont sterowania[s] 
D       = 1600;                                             % Horyzont dynamiki[s]
lambda  = 1e-1;                                             % Parametr lambda 

% Wyznaczenie odpowiedzi skokowej z rzêdnymi co 1s
step = [1; 0];                                              % Skok wartoœci sterowania                        
step_val = 0.1 * U0_(1);                                    % Wartoœæ skoku jako 10%sterowania w stanie ustalonym (nie ma znaczenia, bo uk³ad jest zlinearyzowany)
T = [0, N];                                                 % Czas symulacji do horyzontu predykcji
options = odeset('MaxStep', 1, 'Refine',1);                             % Ustawienia solvera
controlhandle = step * step_val;                            % Sterowanie skokiem

X_X = zeros(N + 1, 2, n);
t = 0 : 1 : N;                                              % Chwile rozwi¹zañ
for i = 1 : n
    control_delay([], [], 0);
    step_sol_X = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_X(:,i), A_X(:,:,i), B_X(:,:,i)), T, X0_X(:,i), options);
    X_X(:,:,i) = (deval(step_sol_X, t))';                   % Wartoœci rozwi¹zañ
end             

Mp_X = zeros(N, D-1);                                         % Macierz Mp jest wymiaru NxD-1
M_X = zeros(N, Nu);                                           % Macierz M jest wymiaru NxNu
%L = lambda * eye(Nu);                                       % Macierz lambda * I

lambda = [0.8 0.004 0.0032 0.00049 0.0001];
dmc_TS.n = n;
dmc_TS.D = D;
for reg = 1 : n 
    s_X = X_X(:,2, reg);
    s_X = (s_X - s_X(1)) / (s_X(end) - s_X(1));
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

    L = lambda(reg) * eye(Nu);
    % Wyznaczanie wektora K1
    K_X = (M_X' * M_X + L) \ M_X';
    K1_X = K_X(1,:);

    % Struktura zawieraj¹ca wszystkie dane potrzebne do realizacji regulacji
    dmc_TS.ku_X(:,reg) = K1_X * Mp_X;
    dmc_TS.ke_X(reg) = sum(K1_X);
end


lambda = [10 1 0.1 0.01 0.001];
X_X = zeros(N + 1, 2, n);
t = 0 : 1 : N;                                              % Chwile rozwi¹zañ
for i = 1 : n
    control_delay([], [], 0);
    step_sol_X = ode45(@(t, x)plant_lin(t, x, ...
        control_delay(t, controlhandle), X0_X_nq(:,i), A_X_nq(:,:,i), B_X_nq(:,:,i)), T, X0_X_nq(:,i), options);
    X_X(:,:,i) = (deval(step_sol_X, t))';                   % Wartoœci rozwi¹zañ
end             

Mp_X = zeros(N, D-1);                                         % Macierz Mp jest wymiaru NxD-1
M_X = zeros(N, Nu);                                           % Macierz M jest wymiaru NxNu
L = lambda(reg) * eye(Nu);                                       % Macierz lambda * I

dmc_TS_nq.n = n;
dmc_TS_nq.D = D;
for reg = 1 : n 
    s_X = X_X(:,2, reg);
    s_X = (s_X - s_X(1)) / (s_X(end) - s_X(1));
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

    L = lambda(reg) * eye(Nu);
    % Wyznaczanie wektora K1
    K_X = (M_X' * M_X + L) \ M_X';
    K1_X = K_X(1,:);

    % Struktura zawieraj¹ca wszystkie dane potrzebne do realizacji regulacji
    dmc_TS_nq.ku_X(:,reg) = K1_X * Mp_X;
    dmc_TS_nq.ke_X(reg) = sum(K1_X);
end



%% Symulacja przebiegu regulacji
controlhandle_X_q = @DMC_TS_X_q;                                       % Regulacja regulatorem DMC
controlhandle_X = @DMC_TS_X;


des_h2 = data.h20 * (1.25);                                  % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle_X_q([],[], [], dmc_TS, [], 0);                               % Reset historii regulatora
sol_X = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X_q(t, x, des_h2, dmc_TS, activ_X) + szszsz(t, dis_ampl)), U0_, data), T_sym, X0_, options);
[~, t_X, U_X] = control_delay([], [], 1);                    % Przebieg sterowania
                                  % ¯¹dana wartoœæ wysokoœci h2
control_delay([], [], 0);                                   % Reset funkcji realizuj¹cej opóŸnienie
controlhandle_X([],[], [], dmc_TS_nq, []);                               % Reset historii regulatora
sol_X_nq = ode45(@(t, x)plant(t, x, control_delay(t, controlhandle_X(t, x, des_h2, dmc_TS_nq) + szszsz(t, dis_ampl)), U0_, data), T_sym, X0_, options);
[~, t_X_nq, U_X_nq] = control_delay([], [], 1);                    % Przebieg sterowania







t_sym = 0:1:T_sym(2);                                       % Wektor czasu symulacji
X = (deval(sol_X, t_sym))';  
X_nq = (deval(sol_X_nq, t_sym))'; 



%% Wykres

figure();
hold on;
plot(t_sym, X(:,2));
plot(t_sym, X_nq(:,2));
plot(t_sym, des_h2 * t_sym./t_sym);
plot(t_sym, (1.05*(des_h2 - data.h20)+data.h20) * t_sym./t_sym, ':black');
plot(t_sym, (0.95*(des_h2 - data.h20)+data.h20) * t_sym./t_sym, ':black');
grid on;
legend('zmodyfikowane funkcje aktywacji', 'równomierne funkcje aktywacji', 'wartoœæ zadana');
xlabel('czas [s]');
ylabel('poziom p³ynu w drugim zbiorniku[cm]');




% figure('NumberTitle', 'off', 'Name', 'Odpowiedzi regulatorów DMC rozmytych nierównomiernie na skokowe zmiany wartoœci zadanej');
% hz = ones(1, (T_sym(2) - T_sym(1)) + 1);
% hz1 = hz * ( data.h20 * (1 - rang) );
% hz2 = hz * ( data.h20 * (1 - rang*2/3) );
% hz3 = hz * ( data.h20 * (1 - rang/6) );
% hz4 = hz * ( data.h20 * (1 + rang/6) );
% hz5 = hz * ( data.h20 * (1 + rang*2/3) );
% hz6 = hz * ( data.h20 * (1 + rang) );
% plot(t_sym, X1_X(:,2),'--b',t_sym, X2_X(:,2),'--r',t_sym, X3_X(:,2),'--c',...
%     t_sym, X4_X(:,2),'--m',t_sym, X5_X(:,2),'--g',t_sym, X6_X(:,2),'--y',...
%     t_sym, hz1,  'b',t_sym, hz2, 'r',t_sym, hz3, 'c',...
%     t_sym, hz4, 'm',t_sym, hz5, 'g',t_sym, hz6, 'y');
% xlabel('$t[s]$', 'interpreter', 'latex');
% ylabel('$h_2[cm]$', 'interpreter', 'latex');
% title('Przebieg wysokoœci p³ynu w drugim zbiorniku', 'interpreter', 'latex');
% legend('DMC_TS_X 1','DMC_TS_X 2','DMC_TS_X 3','DMC_TS_X 4','DMC_TS_X 5','DMC_TS_X 6',...
%     'wartoœæ zadana 1','wartoœæ zadana 2','wartoœæ zadana 3','wartoœæ zadana 4','wartoœæ zadana 5','wartoœæ zadana 6',...
%     'location', 'northwest');
% grid on;
% 
% figure('NumberTitle', 'off', 'Name', 'Przebiegi sterowañ odpowiedzi na wymuszenie skokowe');
% plot(t_U1_X,U0_(1) + U1_X(1,:),'b',t_U2_X,U0_(1) +  U2_X(1,:),'r',t_U3_X,U0_(1) + U3_X(1,:),'c',...
%     t_U4_X,U0_(1) + U4_X(1,:),'m',t_U5_X,U0_(1) + U5_X(1,:),'g',t_U6_X,U0_(1) + U6_X(1,:),'y');
% xlabel('$t[s]$', 'interpreter', 'latex');
% ylabel('$U[\frac{cm^3}{s}]$', 'interpreter', 'latex');
% title('Przebieg wartoœci sterowania', 'interpreter', 'latex');
% legend('DMC_TS_X 1','DMC_TS_X 2','DMC_TS_X 3','DMC_TS_X 4','DMC_TS_X 5','DMC_TS_X 6','location', 'northwest');
% grid on;

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
