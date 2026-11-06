Tp = 1;
sigma2v = 0.1;
M = 50;
D = [Zdata(:,2), Zdata(:,1)]; 

g12_uproszczony = zeros(M, 1);
ruu0 = Covar([D(:,2), D(:,2)], 0);

for i = 0:M-1
    g12_uproszczony(i+1) = (1/Tp) * (Covar(D, i) / ruu0);
end

ryu = zeros(M, 1);
Ruu = zeros(M, M);

for i = 0:M-1
    ryu(i+1) = Covar(D, i);
    for j = 0:M-1
        Ruu(i+1, j+1) = Covar([D(:,2), D(:,2)], abs(i-j));
    end
end

g12_pelny = (1/Tp) * pinv(Ruu) * ryu; 

t = (0:M-1) * Tp;
figure;
stem(t, g12_pelny, 'b', 'LineWidth', 1.5); hold on;
stem(t, g12_uproszczony, 'r--', 'LineWidth', 1.5);
legend('Wzór (12) - pełny', 'Wzór (12) - uproszczony');
title('Estymaty odpowiedzi impulsowej systemu G_o(s)');
xlabel('Czas [s]');
ylabel('Amplituda');

licznik_s = 0.5;
mianownik_s = [5, 11, 7, 1];
opoznienie = 3; 

G_s = tf(licznik_s, mianownik_s, 'InputDelay', opoznienie);

t = (0:M-1) * Tp; 

g_o_t = impulse(G_s, t);

figure;
stem(t, g_o_t, 'g', 'LineWidth', 1.5);
title('Rzeczywista odpowiedź impulsowa g_o(t) w chwilach próbkowania');
xlabel('Czas [s]');
ylabel('Amplituda')



h_est = Tp * cumsum(g12_pelny);


figure;
stairs(t, h_est, 'b', 'LineWidth', 1.5);
title('Estymata odpowiedzi skokowej h_M(nT_p)');
xlabel('Czas [s]');
ylabel('Amplituda');