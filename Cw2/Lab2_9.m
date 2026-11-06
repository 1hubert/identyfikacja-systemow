N = 1000;
Tp = 0.5;
Mw = 200;
sigma2v = 0.001;

u = Zdata(1:N, 1);
y = Zdata(1:N, 2);


U_jw = fft(u);
Y_jw = fft(y);
G15 = Y_jw ./ U_jw;

w = 0.5 * (1 + cos((0:Mw) * pi / Mw)); 

ruuP = zeros(1, Mw+1);
ryuP_dod = zeros(1, Mw+1);
ryuP_uje = zeros(1, Mw+1);

for tau = 0:Mw
    ruuP(tau+1) = Covar([u, u], tau) * w(tau+1);
    ryuP_dod(tau+1) = Covar([y, u], tau) * w(tau+1); 
    ryuP_uje(tau+1) = Covar([u, y], tau) * w(tau+1); 
end

Ruu = [ruuP(1:Mw+1), zeros(1, 2*N - 2*Mw - 2), ruuP(Mw+1:-1:2)];
Ryu = [ryuP_dod(1:Mw+1), zeros(1, 2*N - 2*Mw - 2), ryuP_uje(Mw+1:-1:2)];

Phi_uu = Tp * fft(Ruu);
Phi_yu = Tp * fft(Ryu);

G16 = Phi_yu ./ Phi_uu; 

D1 = N;                     
D1_star = D1 / 2 - 1;      
k1 = 0:D1_star;
omega1 = (2*pi / (D1*Tp)) * k1; 
Lm15 = 20 * log10(abs(G15(k1+1)));
Arg15 = unwrap(angle(G15(k1+1))) * 180/pi;

D2 = length(G16);           
D2_star = (D2 - 1) / 2;    
k2 = 0:D2_star;
omega2 = (2*pi / (D2*Tp)) * k2;
Lm16 = 20 * log10(abs(G16(k2+1)));
Arg16 = unwrap(angle(G16(k2+1))) * 180/pi;

figure;

subplot(2,1,1);
semilogx(omega1, Lm15, 'b', omega2, Lm16, 'r', 'LineWidth', 1.5);
grid on;
title('Charakterystyka modułu Lm(\omega)');
xlabel('\omega [rad/s]');
ylabel('Lm [dB]');
legend('Wzór (15) - estymator transmitancji', 'Wzór (16) - iloraz estymatorów');

subplot(2,1,2);
semilogx(omega1, Arg15, 'b', omega2, Arg16, 'r', 'LineWidth', 1.5);
grid on;
title('Charakterystyka fazowa Arg(G(j\omega))');
xlabel('\omega [rad/s]');
ylabel('Arg [deg]');
legend('Wzór (15) - estymator transmitancji', 'Wzór (16) - iloraz estymatorów');