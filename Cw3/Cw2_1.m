

%% 1. Wczytanie danych i parametry z instrukcji
load IdentWsadowaDyn.mat;   
Tp = 0.01;                  
N  = size(DaneDynW, 1); 

N_est = floor(N / 2);

u_est = DaneDynW(1:N_est, 1);
y_est = DaneDynW(1:N_est, 2);

u_wer = DaneDynW(N_est+1:end, 1);
y_wer = DaneDynW(N_est+1:end, 2);

%% 2. Wpisanie modelu w postaci regresji liniowej: Y = Phi * p
Y = y_est(2:end);

Phi = [-y_est(1:end-1), u_est(1:end-1)];

%% 3. Estymacja parametrów modelu dyskretnego metodą LS
p_LS = (Phi' * Phi) \ (Phi' * Y);   % wzór (12): (Phi^T * Phi)^(-1) * Phi^T * Y

a1 = p_LS(1);
b1 = p_LS(2);

%% 4. Rekonstrukcja parametrów systemu czasu ciągłego (k oraz T)
T_hat = -Tp / log(-a1);
k_hat = b1 / (1 + a1);

fprintf('Wyestymowane parametry: k = %.3f (prawdziwe: 2.0), T = %.3f (prawdziwe: 0.5)\n', k_hat, T_hat);

%% 5. Zapisanie zidentyfikowanego modelu jako transmitancja w MATLAB-ie
Gz = tf(b1, [1 a1], Tp);     
Gs = tf(k_hat, [T_hat 1]);   

fprintf('====================================================\n');
display(Gz)

fprintf('====================================================\n');
display(Gs)



%% Weryfikacja modelu na zbiorze Z_wer (przykład dla DaneDynW)
Nv    = length(u_wer);
t_wer = (0:Nv-1)' * Tp;

% 1. Niezakłócona odpowiedź prawdziwego obiektu Go(s) = 2.0 / (0.5s + 1)
Go     = tf(2.0, [0.5 1]);
yo_wer = lsim(Go, u_wer, t_wer, 'zoh');

% 2. Odpowiedź predyktora jednokrokowego y_hat(n|n-1)
Phi_wer = [-y_wer(1:end-1), u_wer(1:end-1)];
y_pred  = [y_wer(1); Phi_wer * p_LS];

% 3. Odpowiedź modelu symulowanego y_m(n)
ym_wer = lsim(Gz, u_wer, t_wer);

%% Wskaźniki jakości (wzór 21)
Vp = mean((y_wer - y_pred).^2);
Vm = mean((yo_wer - ym_wer).^2);

fprintf('====================================================\n');
fprintf('Wskaźniki jakości: Vp = %.6f, Vm = %.6f\n', Vp, Vm);

%% Wspólny wykres 4 przebiegów w dziedzinie czasu
figure;
plot(t_wer, y_wer,  'k', 'LineWidth', 1.1); hold on;
plot(t_wer, y_pred, 'g--', 'LineWidth', 1.1);
plot(t_wer, yo_wer, 'b',   'LineWidth', 1.8);
plot(t_wer, ym_wer, 'r--', 'LineWidth', 1.5);
grid on;
grid on;
xlabel('Czas t [s]');
ylabel('Amplituda');
legend('y(n) - zmierzone', 'y\_hat(n|n-1) - predyktor', ...
       'y_o(n) - niezakłócone', 'y_m(n) - model symulowany');
title('Weryfikacja identyfikacji metodą LS na zbiorze Z_{wer}');


%% Wyznaczenie macierzy kowariancji (wzór 14) oraz przedziałów ufności (Uwaga 1)
N_p = length(Y);          % liczba równań użytych do estymacji (N_est - 1)
dp  = length(p_LS);       % liczba estymowanych parametrów (dp = 2: a1, b1)

% 1. Błędy resztowe i estymata wariancji szumu (wzór 14)
eps_res   = Y - Phi * p_LS;
sigma2_hat = sum(eps_res.^2) / (N_p - dp);

% 2. Macierz kowariancji estymat Cov[p_LS] (wzór 14)
Cov_p = sigma2_hat * inv(Phi' * Phi);     % lub: sigma2_hat * ((Phi' * Phi) \ eye(dp))

% 3. Macierz P_inf dla regresora stochastycznego (Uwaga 1, str. 7)
P_inf = N_p * Cov_p;

% 4. 95% przedziały ufności PU_95% (wzór 25)
% Uwaga: P_inf(i,i) / N_p jest równe Cov_p(i,i)
std_p = sqrt(diag(Cov_p));                % odchylenia standardowe estymat
PU_95 = [p_LS - 1.96 * std_p, p_LS + 1.96 * std_p];

fprintf('====================================================\n');

disp('Macierz kowariancji Cov[p_LS]:');
disp(Cov_p);


fprintf('====================================================\n');
fprintf('Przedziały ufności: \n');
fprintf('a1 = %.4f, PU_95%% = [%.4f ; %.4f]\n', p_LS(1), PU_95(1,1), PU_95(1,2));
fprintf('b1 = %.4f, PU_95%% = [%.4f ; %.4f]\n', p_LS(2), PU_95(2,1), PU_95(2,2));