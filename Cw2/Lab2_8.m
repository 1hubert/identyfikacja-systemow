Tp = 0.5;
TF = 2.5*Tp/(2*pi);
sigma2v = 0.001;
N = 1000;     
t = (0:N-1) * Tp;

u = Zdata(1:N, 1);
y = Zdata(1:N, 2);

figure;

subplot(2,1,1);
plot(t, u, 'b');
title('Sygnał pobudzający u(t) po filtracji');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(t, y, 'r');
title('Odpowiedź systemu y(t) po filtracji');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;