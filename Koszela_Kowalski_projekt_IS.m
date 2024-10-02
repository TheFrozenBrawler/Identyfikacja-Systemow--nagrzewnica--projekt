close all; clear variables;
data = load('dane.mat');
Tp = 0.08;
u = data.in;
y = data.out;
D = [y u];
M = 150;

%% ----------------------------------------------------
%% METODA ANALIZY KORELACYJNEJ
%% ----------------------------------------------------

%% 1. Podział danych na estymacyjne i weryfikacyjne 
smpl_size = size(u);        %liczba dostępnych probek sygnalu wejsciowego
trshld = 0.5*smpl_size(1);  %Podział zbioru danych na dwie czesci
u_est = u(1:trshld);
u_wer = u(trshld+1:end);
y_est = y(1:trshld);
y_wer = y(trshld+1:end);

% wykreślenie wektorów czasu dla próbek estymaty i weryfikacji
t = 0:Tp:999*Tp;
t_wer = 0:Tp:(length(u_wer)-1)*Tp;
t_est = 0:Tp:(trshld-1)*Tp;

%% 2. Obliczanie autokorelacji sygnału wejściowego
r_uulong = [];
l = 1;

%Dla ujemnych wartości przesuniecia tau:
for tau = M-1:-1:1 
    r_uulong(l) = Covar([u u], tau);
    l = l+1;
end

%Dla tau = 0
r_uulong(l) = Covar([u u], 0);
l = l+1;

%Dla tau>0
for tau = 1:1:M-1 
    r_uulong(l) = Covar([u u], tau);
    l = l+1;
end

%Stworzenie macierzy R_uu
R_uu = [];
k = M;
for i = 1:M  % iteracja po rzędach
    R_uu(i, :) = r_uulong(k:k+M-1);
    k = k-1;
end

%% 3. Stworzenie wektora korelacji wzajemnej
r_yu = [];
for i = 1:1:M
    r_yu(i, 1) = Covar([y u], i-1);
end

%% 4. Obliczenie estymowanej odpowiedzi impulsowej
g_M = 1/Tp*pinv(R_uu)*r_yu;
t_gM = 0:Tp:M*Tp-Tp;

%% 5. Obliczenie estymowanej odpowiedzi skokowej 
for i = 1:M
    h_M(i) = Tp* sum(g_M(1:i));
end

%% 6. Plot 
figure;
plot(t_gM, h_M, 'r')
title("Estymowana odpowiedź skokowa");
legend("odp. skokowa");
xlabel("Czas [s]");
ylabel("odpowiedź [V]");

figure;
plot(t_gM, g_M);
title("Estymowana odpowiedź impulsowa");
legend("odp. impulsowa");
xlabel("Czas [s]");
ylabel("odpowiedź [V]");

%% 7. Wyznaczenie stałej czasowej oraz czasu opóźnienia
Au = 1;
y_ust = 0.974;
K = y_ust/Au;
%max(estymowana odpowiedz impulsowa) = 1.75786
ag = 1.75786;
tg = 0.4;
% sg = s(tg) w h(t)
sg = 0.5;
Td = (ag*tg-sg)/ag;
T = Au*K/ag;

%% -----------------------------------------------------------------------------
%% IDENTYFIKACJA POŚREDNIA SYSTEMU DYNAMICZNEGO CZASU CIĄGŁEGO METODAMI LS I IV
%% -----------------------------------------------------------------------------


%% IDENTYFIKACJA PARAMETRÓW METODĄ LS
%% 1. Definicja transmitancji: 
G_s = tf(K, [T 1], 'InputDelay', Td)
G_z = c2d(G_s, Tp, 'zoh')

%% 2. Podział danych na estymacyjne oraz weryfikacyjne
smpl_size = size(u);
trshld = 0.5 * smpl_size(1);

u_est = u(1:trshld);
u_wer = u(trshld:smpl_size);

y_est = y(1:trshld);
y_wer = y(trshld:smpl_size);

%% 2. Stworzenie macierzy regresji
Phi = zeros(trshld, 3);
Phi(1, :) = [0 0 0];             %zerowe warunki poczatkowe
Phi(2, :) = [y_est(1) 0 0];       %zerowe warunki poczatkowe v2
Phi(3, :) = [y_est(2) u_est(1) 0];%zerowe warunki poczatkowe v3 
for i = 4:trshld
    Phi(i,:) = [y_est(i-1) u_est(i-2) u_est(i-3)]; %wektor regresji
end


p_hatLS = pinv(Phi) * y_est;    %Estymator parametrów LS 

%Przypisanie estymaty do parametrow
p3 = p_hatLS(1);
p1 = p_hatLS(2);
p2 = p_hatLS(3);

Gm = tf([p1 p2], [1 -p3], Tp, 'InputDelay', 2);
Gm_ss = ss(Gm);     %konwersja

%% 3. PLOT LS
tm = 0:Tp:(length(u_wer)-1)*Tp;
ym = lsim(Gm_ss, u_wer, tm, 13.21);
figure;
hold on;
plot(tm, y_wer, 'b');
plot(tm, ym, 'r--');
hold off;
title("Odpowiedź modelu z parametrami LS");
legend("odp. obiektu", "odp. modelu");
xlabel("Czas [s]");
ylabel("odpowiedź [V]");

%% 4. JFIT dla LS
ymean = mean(y_wer);
Jfitls = (1 -(( sum ( abs ( y_wer - ym ) ) ) /( sum ( abs ( y_wer - ymean ) ) ) ) ) *100;
disp("Jfit dla metody LS:")
disp(Jfitls);

%% Przedziały ufności dla LS
% Phi_p = zeros(trshld, 3);
% Phi_p(1, :) = [0 0 0];             %zerowe warunki poczatkowe
% Phi_p(2, :) = [y(1) 0 0];       %zerowe warunki poczatkowe v2
% Phi_p(3, :) = [y(2) u_est(1) 0];%zerowe warunki poczatkowe v3 
% for i = 4:length(y)
%     Phi_p(i,:) = [y(i-1) u(i-2) u(i-3)]; 
% end
% p_hatLS_przedz = pinv(Phi_p) * y;
% eps = y-Phi_p * p_hatLS_przedz;
% sum = 0;
% for i = 1:length(eps)
%     sum = sum + (eps(i)*eps(i));
% end
% 
% N = length(y);
% d = 3;
% sigm2 = 1 / (N-d) * sum;

% macierz kowariancji
% cov = sigm2* inv( (transpose(Phi_p)*Phi_p) );
% Pinf = N*sigm2* inv( (transpose(Phi_p)*Phi_p) );
% Macierz kowariancji estymat - estymuje rzeczywistą macierz kowariancji która opisuje zmienność i zależność między zmiennymi
%% przedziały ufności
% p_ufn_min = p_hatLS_przedz - 1.96*sqrt(diag(Pinf)/N);
% p_ufn_max = p_hatLS_przedz + 1.96*sqrt(diag(Pinf)/N);

% wyniki przedziałów ufności:
%  p_ufn_min =
%     0.8812
%    -0.0087
%     0.0965
% 
% p_ufn_max =
%     0.8996
%     0.0109
%     0.1173
%


%% 5. Predyktor LS
y_pred = Phi * p_hatLS;
figure;
hold on;
plot(t_est, y_est);
plot(t_est, y_pred);
hold off;
title("Wykres predyktora dla metody LS")
legend("odpowiedź obiektu rzeczywistego", "odpowiedź predyktora LS");
xlabel("Czas [s]");
ylabel("odpowiedź [V]")

ymean = mean(y_est);
Jfitls = (1 -(( sum ( abs ( y_est - y_pred ) ) ) /( sum ( abs ( y_est - ymean ) ) ) ) ) *100;
disp("Jfit predyktora LS:");
disp(Jfitls);

%% IDENTYFIKACJA PARAMETROW METODA IV

%% 1. WEKTOR X
G_LS = tf([p1 p2], [1 -p3], Tp, 'InputDelay', 2);
x = lsim(G_LS, u);

%% 2. MACIERZ ZMIENNYCH INSTRUMENTALNYCH
Z = zeros(trshld, 3);
Z(1, :) = [0 0 0];
Z(2, :) = [x(1) 0 0];
Z(3, :) = [x(2) u(1) 0];

for n=4 : trshld
    Z(n, :) = [x(n-1) u_est(n-2) u_est(n-3)];
end

p_hatIV = pinv(Z'*Phi) * Z' * y_est;  %Estymator parametrów IV

%Przypisanie wartości estymat do parametrów
p3 = p_hatIV(1);
p1 = p_hatIV(2);
p2 = p_hatIV(3);

%% 3. STWORZENIE ESTYMOWANEJ TRANSMITANCJI Z PARAMETRAMI IV
Gmiv = tf([p1 p2], [1 -p3], Tp, 'InputDelay', 2);
Gmiv_ss = ss(Gmiv); 

%% 4. PLOT IV

tmiv = 0:Tp:(length(u_wer)-1)*Tp;
ymiv = lsim(Gmiv_ss, u_wer, tm, 13.31);
figure;
hold on;
plot(tmiv, y_wer, 'b');
plot(tmiv, ymiv, 'r--');
hold off;
title("Odpowiedź modelu z parametrami IV");
legend("odp. obiektu", "odp. modelu");
xlabel("Czas [s]");
ylabel("odpowiedź [V]");

%% 5. Jfit dla IV
ymean = mean(y_wer);
Jfitls = (1 -(( sum ( abs ( y_wer - ymiv ) ) ) /( sum ( abs ( y_wer - ymean ) ) ) ) ) *100;
disp("Jfit dla metody IV");
disp(Jfitls);

%% 6. PREDYKTOR IV
y_prediv = Phi * p_hatIV;
figure;
hold on;
plot(t_est, y_est);
plot(t_est, y_prediv);
title("Wykres predyktora dla metody IV")
legend("odpowiedź obiektu rzeczywistego", "odpowiedź predyktora IV");
xlabel("Czas [s]");
ylabel("odpowiedź [V]")
hold off;

% Jfit predyktora IV
ymean = mean(y_est);
Jfitiv = (1 -(( sum ( abs ( y_est - y_prediv ) ) ) /( sum ( abs ( y_est - ymean ) ) ) ) ) *100;
disp("Jfit predykatora IV:");
disp(Jfitiv);