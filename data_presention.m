close all; clear variables;

%% prezentacja danych
data = load('dane.mat');
Tp = 0.08;
u = data.in;
y = data.out;
t = 0:Tp:999*Tp;

figure;
hold on;
plot(t, u, 'r');
plot(t, y, 'b');
hold off;
% 
% figure;
% plot(u_o, y_o, 'g');
% xlabel("u_o");
% ylabel("y_o");

%% określenie transmitancji dyskretnej
p3 = 0.8823;
p1 = 0.0051;
p2 = 0.1097;


k = 0.9456;
T = 0.;
Td = 2;

Gm_z = tf([p1 p2], [1 p3], Tp, 'InputDelay', Td);


G = tf([k], [T 1]);
Gz = c2d(G, Tp, 'zoh');

ym = lsim(G, u, t, [0 4]);

figure;
hold on;
ym_wer = ym(800:1000,1);
y_wer = y(800:1000,1);
t_wer = t(1,800:1000);
plot(t_wer, ym_wer, 'r--');
plot(t_wer, y_wer, 'b');

ym_wer = ym(800:1000,1);
y_wer = y(800:1000,1);
ymean = mean(y_wer);
Jfit = (1 -(( sum ( abs ( y_wer - ym_wer ) ) ) /( sum ( abs ( y_wer - ymean ) ) ) ) ) *100;
disp(Jfit);
