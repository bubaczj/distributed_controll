%%plik informacyjny
% Na wszystko, co œwiête i zbo¿ne: nie usuwaæ plikó mpc i mpc_TS z folderu
% controls

%% Odpowiedzi ukladow

%lin_nonlin_stepresp.m 
%porownanie odpowiedzi ukladu liniowego i nieliniowego 
% linia 25 step_val - wartoœæ skoków sterowania
% linia 26 T - czas symulacji

%obj_compare_stepresp.m 
%porownanie obiektow rozmytych rownomiernie i nierownomiernie z liniowym i nieliniowym
% linia 25 n - liczba punktów rozmycia
% linia 28 q - podstawa ci¹gu geometrycznego
% linia 33 T - czas symulacji

%obj_lin_stepresp
%odpowiedz ukladu liniowego na skok
% linia 25 step_val - wartoœci skoków sterowania
% linia 26 T - czas symulacji


%obj_TS_X_stepresp
%odpowiedz rozmytego rownomiernie na skok
% linia 23 n - liczba punktów rozmycia ( n > 1 )
% linia 30 step_val - wartoœæ skoków sterowania
% linia 31 T - czas symulacji


%obj_TX_X_q_stepresp
%odpowiedz rozmytego nierownomiernie na skok
% linia 24 n - liczba punktów rozmycia ( n > 1 )
% linia 26 q - podstawa ci¹gu geometrycznego
% linia 32 step_val - wartoœci skoków sterowania
% linia 33 T - czas symulacji

%TS_nonlin_stepresp
%porownanie rozmycia po stanach i sterowaniu z nieliniowym
% linia 24 n - liczba punktów rozmycia (n > 1 )
% linia 32 step_val - wartoœci skoków sterowania
% linia 33 T - czas symulacji

%TS_q_nonlin_stepresp
%porownanie rozmycia nierownomiernego po stanach i sterowaniu z nieliniowym
%(+ przebiegi funkcji akywacji)
% linia 24 n - liczba punktów rozmycia (n > 1 )
% linia 26 q - podstawa ci¹gu geometrycznego
% linia 34 T - czas symulacji


%% Regulacja

%DMC_reg 
%symulacje regulacji DMC
% linia 35 T_sym - ustawienie casu sumylacji
% linia 78 dis_ampl - zakres zak³óceñ 
% linia 79 rng - ustawienie zakresu skoków; 
% wartoœæ zadana = h20 +/- rng * h20

%DMC_TS_reg
%symulacje regulacji DMC rozmytym po stanach i sterowaniu
% linia 23 n - liczba punktów rozmycia (n > 1)
% linia 42 T_sym - czas symulacji
% linia 43 dis_ampl - zakres zak³óceñ
% linia 118 rang - ustawienie zakresu skoków, analogicznie j/w

%DMC_TS_X_q_reg
%symulacje regulacji nierownomiernie rozmytym DMC
% linia 22 n - liczba punktów rozmycia (n > 1)
% linia 24 q - podstawa ci¹gu geometrycznego
% linia 43 T_sym - czas symulacji
% linia 44 dis_ampl - zakres zak³óceñ
% linia 150 des_h2 - wartoœæ zadana regulacji

% SL_script
% symulacja uk³adu z regulatorem SL
% linia 26 n - liczba punktów rozmycia (n > 1)
% linia 33 T - czas symulacji
% linie 45, 56, 68, 80, 92, 105 des_h2 - wartoœci zadane dla kolejnych
% symulacji


