v_row = v';

N = 2000;
Tp = 0.001;
sigma = 0.8;

t = 0:N-1;
H = tf([0.1], [1 -0.9],Tp);

x = sin(2 * pi * 5 * t * Tp) + 0.5 * sin(2 * pi * 10 * t * Tp) + 0.25 * sin(2 * pi * 30 * t * Tp);
e = sigma * randn(1, N);
v = lsim(H, e, t*Tp);

% plotting
figure;
plot(t, x);
title('x');
figure;
plot(t, e);
title('e');
figure;
plot(t, v);
title('v');

% fft

X = fft(x) * Tp;
spectrum_x = abs(X) * (2 / (N * Tp));

f = t / (N * Tp);

figure;
plot(f, spectrum_x);
title('spectrum x (wycinek 0 do 50)');
xlim([0 50]);

figure;
plot(f, spectrum_x);
title('spectrum x (pełne)');

% w okresie sygnału x(nTp) mieści się 200 próbek
% dlatego N<200 sprawi że doświadczymy wycieku widma

E_czas = Tp * sum(x.^2);

E_czest = sum((1 / (N * Tp)) * (abs(X).^2));

% Porównanie energii, czyli sprawdzenie twierdzenia Parsevala
fprintf('Energia w dziedzinie czasu: %f\n', E_czas);
fprintf('Energia w dziedzinie częstotliwości: %f\n', E_czest);



% estymata gęstości widmowej bezpośrednio i pośrednio (sygnał e)

E = fft(e) * Tp; 

Phi_ee_1 = (1 / (N * Tp)) * abs(E).^2; 

Mw = 100; % Przykładowa szerokość okna

[r_ee, lags] = xcorr(e, 'biased'); 


idx_zero = find(lags == 0);

r_ee_cut = r_ee(idx_zero - Mw : idx_zero + Mw); 


w = hamming(2 * Mw + 1)'; 


r_ee_w = r_ee_cut .* w;


r_ee_padded = zeros(1, N);
r_ee_padded(1:Mw+1) = r_ee_w(Mw+1:end); % Dodatnie opóźnienia
r_ee_padded(N-Mw+1:N) = r_ee_w(1:Mw);   % Ujemne opóźnienia


Phi_ee_2 = Tp * abs(fft(r_ee_padded));

figure;
plot(f, Phi_ee_1, 'b');
hold on;
plot(f, Phi_ee_2, 'r', 'LineWidth', 1.5);
title('Estymaty gęstości widmowej mocy sygnału e');
xlabel('Częstotliwość [Hz]');
ylabel('\Phi_{ee}');
legend('Estymator ze wzoru (14) - Periodogram', 'Estymator ze wzoru (16) - Korelogram');
xlim([0 (1/(2*Tp))]); % Wyświetlanie do częstotliwości Nyquista

% estymata gęstości widmowej bezpośrednio i pośrednio (sygnał v)
% --- Estymator 1: Periodogram dla v ---
V = fft(v_row) * Tp; 
Phi_vv_1 = (1 / (N * Tp)) * abs(V).^2; 

% --- Estymator 2: Korelogram dla v ---
Mw_v = 100;
[r_vv, lags_v] = xcorr(v_row, 'biased'); 
idx_zero_v = find(lags_v == 0);
r_vv_cut = r_vv(idx_zero_v - Mw_v : idx_zero_v + Mw_v); 
w_v = hamming(2 * Mw_v + 1)'; 
r_vv_w = r_vv_cut .* w_v;

r_vv_padded = zeros(1, N);
r_vv_padded(1:Mw_v+1) = r_vv_w(Mw_v+1:end);
r_vv_padded(N-Mw_v+1:N) = r_vv_w(1:Mw_v);
Phi_vv_2 = Tp * abs(fft(r_vv_padded));

% Rysowanie wyników dla sygnału v
figure;
plot(f, Phi_vv_1, 'b');
hold on;
plot(f, Phi_vv_2, 'r', 'LineWidth', 1.5);
title('Estymaty gęstości widmowej mocy sygnału v (szum kolorowy)');
xlabel('Częstotliwość [Hz]');
ylabel('\Phi_{vv}');
legend('Estymator (14) - Periodogram', 'Estymator (16) - Korelogram');
xlim([0 (1/(2*Tp))]);
