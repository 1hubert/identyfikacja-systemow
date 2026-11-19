

%% 1. Wczytanie danych DaneDynC (zakłócenie kolorowe) i podział 50% / 50%
load IdentWsadowaDyn.mat;
Tp = 0.01;
N  = size(DaneDynC, 1);
N_est = floor(N / 2);

u_est = DaneDynC(1:N_est, 1);
y_est = DaneDynC(1:N_est, 2);

u_wer = DaneDynC(N_est+1:end, 1);
y_wer = DaneDynC(N_est+1:end, 2);

%% 2. Etap I: Wstępna estymacja metodą LS (wzór 12)
Y   = y_est(2:end);
Phi = [-y_est(1:end-1), u_est(1:end-1)];

p_LS = (Phi' * Phi) \ (Phi' * Y);

%% 3. Generacja zmiennych instrumentalnych z(n) (wzory 17-18)
Gz_LS = tf(p_LS(2), [1 p_LS(1)], Tp);
x_est = lsim(Gz_LS, u_est);             % symulowana odpowiedź x(n)

Z = [-x_est(1:end-1), u_est(1:end-1)];  % macierz zmiennych instrumentalnych

%% 4. Etap II: Estymacja metodą zmiennych instrumentalnych IV (wzór 16)
p_IV = (Z' * Phi) \ (Z' * Y);

a1_IV = p_IV(1);
b1_IV = p_IV(2);

%% 5. Rekonstrukcja parametrów ciągłych k oraz T
T_IV = -Tp / log(-a1_IV);
k_IV = b1_IV / (1 + a1_IV);

fprintf('Metoda IV (DaneDynC): k = %.4f (prawdziwe: 2.0), T = %.4f (prawdziwe: 0.5)\n', k_IV, T_IV);

%% 6. Weryfikacja na zbiorze Z_wer i wskaźniki jakości (wzór 21)
Nv    = length(u_wer);
t_wer = (0:Nv-1)' * Tp;

% Niezakłócona odpowiedź prawdziwego obiektu yo(n)
Go     = tf(2.0, [0.5 1]);
yo_wer = lsim(Go, u_wer, t_wer, 'zoh');

% Predyktor jednokrokowy y_hat(n|n-1)
Phi_wer = [-y_wer(1:end-1), u_wer(1:end-1)];
y_pred  = [y_wer(1); Phi_wer * p_IV];

% Model symulowany ym(n)
Gz_IV  = tf(b1_IV, [1 a1_IV], Tp);
ym_wer = lsim(Gz_IV, u_wer, t_wer);

% Wskaźniki jakości Vp oraz Vm
Vp_IV = mean((y_wer - y_pred).^2);
Vm_IV = mean((yo_wer - ym_wer).^2);

fprintf('Wskaźniki IV: Vp = %.6f, Vm = %.6f\n', Vp_IV, Vm_IV);

%% 7. Wspólny wykres 4 przebiegów w dziedzinie czasu
figure;
plot(t_wer, y_wer,  'k', 'LineWidth', 1.1); hold on;
plot(t_wer, y_pred, 'g--', 'LineWidth', 1.1);
plot(t_wer, yo_wer, 'b',   'LineWidth', 1.8);
plot(t_wer, ym_wer, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Czas t [s]');
ylabel('Amplituda');
legend('y(n) - zmierzone', 'y\_hat(n|n-1) - predyktor', ...
       'y_o(n) - niezakłócone', 'y_m(n) - model symulowany (IV)');
title('Weryfikacja identyfikacji metodą IV dla DaneDynC');