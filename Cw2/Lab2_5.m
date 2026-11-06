M = 50;
sigma2v = 0.0;
D = [Zdata(:,2), Zdata(:,1)];

g11 = zeros(M, 1);
ruu0 = Covar([D(:,2), D(:,2)], 0); 
for i = 0:M-1
    g11(i+1) = Covar(D, i) / ruu0;
end

ryu = zeros(M, 1);
Ruu = zeros(M, M);

for i = 0:M-1
    ryu(i+1) = Covar(D, i);
    for j = 0:M-1
        Ruu(i+1, j+1) = Covar([D(:,2), D(:,2)], abs(i-j)); 
    end
end

g10 = pinv(Ruu) * ryu; 

figure;
stem(0:M-1, g10, 'b', 'LineWidth', 1.5); hold on;
stem(0:M-1, g11, 'r--', 'LineWidth', 1.5);
legend('Wzór (10)', 'Wzór (11)');
title('Estymaty odpowiedzi impulsowej');
xlabel('Próbka (i)');
ylabel('Amplituda');


licznik = 2;
mianownik = [1, -0.7];

g_o = impz(licznik, mianownik, M);

figure;
stem(0:M-1, g_o, 'g', 'LineWidth', 1.5);
title('Rzeczywista odpowiedź impulsowa g_o(n) systemu G_o(z)');
xlabel('Próbka (i)');
ylabel('Amplituda');