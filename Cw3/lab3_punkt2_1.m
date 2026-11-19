load("IdentWsadowaDyn.mat");

Tp = 0.01;
N = 4001;
t = 0:Tp:(N-1)*Tp;

u_DaneDynW = DaneDynW(:, 1);
y_DaneDynW = DaneDynW(:, 2);

u_DaneDynC = DaneDynC(:, 1);
y_DaneDynC = DaneDynC(:, 2);

figure;
subplot(2, 1, 1);
plot(t, u_DaneDynW);
title("DaneDynW: Sygnał wejściowy - u(t)")
xlabel('Czas [s]'); ylabel('u'); grid on;

subplot(2, 1, 2);
plot(t, y_DaneDynW);
title('DaneDynW: Sygnał wyjściowy zakłócony szumem białym');
xlabel('Czas [s]'); ylabel('y'); grid on;


figure;
subplot(2, 1, 1);
plot(t, u_DaneDynC);
title("DaneDynC: Sygnał wejściowy - u(t)")
xlabel('Czas [s]'); ylabel('u'); grid on;

subplot(2, 1, 2);
plot(t, y_DaneDynC);
title('DaneDynC: Sygnał wyjściowy zakłócony szumem kolorowym');
xlabel('Czas [s]'); ylabel('y'); grid on;


% -------------------------------------------------------------------------
% Podział na dane estymacyjne (Z_est) i dane weryfikacyjne (Z_wer)
% -------------------------------------------------------------------------

t_est = t(1:2000);
t_wer = t(2001:4001);

u_DaneDynW_est = u_DaneDynW(1:2000);
y_DaneDynW_est = y_DaneDynW(1:2000);
u_DaneDynC_est = u_DaneDynC(1:2000);
y_DanDynC_est = y_DaneDynC(1:2000);

u_DaneDynW_wer = u_DaneDynW(2001:4001);
y_DaneDynW_wer = y_DaneDynW(2001:4001);
u_DaneDynC_wer = u_DaneDynC(2001:4001);
y_DaneDynC_wer = y_DaneDynC(2001:4001);

% -------------------------------------------------------------------------
% Teraz zajmujemy się tylko `_est`. Osobno dla szumu białego i kolorowego.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% 1. ESTYMACJA (na zbiorze _est)
% -------------------------------------------------------------------------
Phi_W_est = [-y_DaneDynW_est(1:end-1), u_DaneDynW_est(1:end-1)];
Y_W_est   = y_DaneDynW_est(2:end);

p_LS_W = (Phi_W_est' * Phi_W_est) \ (Phi_W_est' * Y_W_est);

a_1_W = p_LS_W(1);
b_1_W = p_LS_W(2);

T_W = -Tp / log(-a_1_W);
k_W = b_1_W / (1 + a_1_W);

% -------------------------------------------------------------------------
% 2. WERYFIKACJA (na zbiorze _wer)
% -------------------------------------------------------------------------
t_wer = t(2001:4001)'; % wektor czasu dla danych weryfikacyjnych (kolumnowy)

% a) Zmierzona odpowiedź ze zbioru wer:
y_meas_W = y_DaneDynW_wer;

% b) Niezakłócona odpowiedź prawdziwego obiektu Go(s):
s = tf('s');
G_o = 2 / (0.5 * s + 1);
y_o_W = lsim(G_o, u_DaneDynW_wer, t_wer, 'zoh');

% c) Odpowiedź predyktora jednokrokowego y_hat(n|n-1):
Phi_W_wer = [-y_DaneDynW_wer(1:end-1), u_DaneDynW_wer(1:end-1)];
y_pred_W  = Phi_W_wer * p_LS_W; % uwaga: ma długość N_wer - 1 (od 2. próbki)

% d) Odpowiedź modelu symulowanego y_m(n) (bez użycia pomiarów y_wer!):
z = tf('z', Tp);
G_m_W = b_1_W / (z + a_1_W);
y_m_W = lsim(G_m_W, u_DaneDynW_wer, t_wer);

% -------------------------------------------------------------------------
% 3. WYKRES I WSKAŹNIKI JAKOŚCI
% -------------------------------------------------------------------------
figure;
hold on; grid on;
plot(t_wer, y_meas_W, 'Color', [0.7 0.7 0.7]);      % 1. y(n) - zmierzone (szare tło)
plot(t_wer, y_o_W, 'k', 'LineWidth', 1.5);         % 2. y_o(n) - prawdziwe niezakłócone
plot(t_wer(2:end), y_pred_W, 'b--');               % 3. y(n|n-1) - predyktor 1-krokowy
plot(t_wer, y_m_W, 'r', 'LineWidth', 1.2);         % 4. y_m(n) - model symulowany
legend('y(n) zmierzone', 'y_o(n) prawdziwe', 'y(n|n-1) predyktor', 'y_m(n) model symulowany');
xlabel('Czas [s]'); ylabel('Wyjście y');
title('Weryfikacja modelu LS dla danych DaneDynW');

% Wskaźniki jakości ze wzoru (21) - dopasowane wymiarami (od 2. próbki)
Vp_W = mean((y_meas_W(2:end) - y_pred_W).^2);
Vm_W = mean((y_o_W(2:end) - y_m_W(2:end)).^2);
disp("Vp_W: ");
disp(Vp_W);
disp("Vm_W: ");
disp(Vm_W);