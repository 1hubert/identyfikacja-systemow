load("IdentWsadowaDyn.mat");

Tp = 0.01;
N = 4001;

u_DaneDynW = DaneDynW(:, 1);
y_DaneDynW = DaneDynW(:, 2);

M = 3200; % wybór liczby danych do zbioru Ze
tE = Tp*(N-M:N)'; % wektor próbek chwil czasowych dla danych estymujących
uE = DaneDynW(N-M:N,1); % wybór wektora próbek sygnału pobudzającego u do zbioru Ze
yE = DaneDynW(N-M:N,2); % wybór wektora próbek sygnału wyjściowego y do zbioru Ze
s = tf('s'); % zmienna operatorowa Laplace'a
TF = 50*Tp; % wybór wartości stałej czasowej dla filtrów SVF
nF = 2; % wybór rzędu dynamiki dla filtrów SVF
F0 = 1/(1+s*TF)^nF; % definicja filtru SFV typu F^0
F1 = s/(1+s*TF)^nF; % definicja filtru SFV typu F^1

yF = lsim(F0,yE,tE,'foh'); % filtracja SVF filtrem F^0 sekwencji yE z ekstrapolacją 'foh'
ypF = lsim(F1,yE,tE,'foh'); % filtracja SVF filtrem F^1 sekwencji yE z ekstrapolacją 'foh'
uF = lsim(F0,uE,tE,'foh'); % filtracja SVF filtrem F^0 sekwencji uE z ekstrapolacją 'foh'

% LS

Phi_W_est = [-ypF(1:end-1), uF(1:end-1)];
Y_W_est   = yF(2:end);

p_LS_W = (Phi_W_est' * Phi_W_est) \ (Phi_W_est' * Y_W_est);
disp(p_LS_W);


% Przygotowanie zbioru weryfikującego Zw (pozostałe próbki od 1 do N-M-1)
tW = Tp*(0:N-M-2)';
uW = u_DaneDynW(1:N-M-1);
yW = y_DaneDynW(1:N-M-1);

% Funkcje pomocnicze do oceny jakości (MSE oraz FIT w %)
calc_MSE = @(y, y_hat) mean((y - y_hat).^2);
calc_FIT = @(y, y_hat) (1 - norm(y - y_hat) / norm(y - mean(y))) * 100;

% -------------------------------------------------------------------------
% 1. Wpływ stanów przejściowych filtrów SVF na jakość identyfikacji
% -------------------------------------------------------------------------
% Przyjmujemy czas zaniku stanów przejściowych jako ok. 5 * TF
k_odrzucone = 2000; % liczba początkowych próbek do odrzucenia

% Model A: Z uwzględnieniem stanów przejściowych (p_LS_W policzone wyżej)
T_z_przejsc = p_LS_W(1);
k_z_przejsc = p_LS_W(2);
G_z_przejsc = tf(k_z_przejsc, [T_z_przejsc, 1]);
y_sym_z_przejsc = lsim(G_z_przejsc, uW, tW, 'foh');

% Model B: Bez stanów przejściowych (odrzucenie pierwszych k_odrzucone próbek)
Phi_bez_przejsc = [-ypF(k_odrzucone+1:end-1), uF(k_odrzucone+1:end-1)];
Y_bez_przejsc   = yF(k_odrzucone+2:end);

p_bez_przejsc = (Phi_bez_przejsc' * Phi_bez_przejsc) \ (Phi_bez_przejsc' * Y_bez_przejsc);
T_bez_przejsc = p_bez_przejsc(1);
k_bez_przejsc = p_bez_przejsc(2);
G_bez_przejsc = tf(k_bez_przejsc, [T_bez_przejsc, 1]);
y_sym_bez_przejsc = lsim(G_bez_przejsc, uW, tW, 'foh');

% -------------------------------------------------------------------------
% 2. Porównanie dwóch alternatywnych postaci regresji liniowej
%    (Wariant 1: Y = yF vs Wariant 2: Y = ypF, czyli yF^(1))
% -------------------------------------------------------------------------
% Wariant 1 (Y = yF) - korzystamy z wyników bez stanów przejściowych
T_w1 = T_bez_przejsc;
k_w1 = k_bez_przejsc;
y_sym_w1 = y_sym_bez_przejsc;

% Wariant 2 (Y = yF^(1)): ypF = -a * yF + b * uF, gdzie a = 1/T, b = k/T
Phi_w2 = [-yF(k_odrzucone+1:end-1), uF(k_odrzucone+1:end-1)];
Y_w2   = ypF(k_odrzucone+2:end);

p_w2 = (Phi_w2' * Phi_w2) \ (Phi_w2' * Y_w2);
T_w2 = 1 / p_w2(1);
k_w2 = p_w2(2) / p_w2(1);
G_w2 = tf(k_w2, [T_w2, 1]);
y_sym_w2 = lsim(G_w2, uW, tW, 'foh');

% -------------------------------------------------------------------------
% Wyświetlenie wyników i wykresów
% -------------------------------------------------------------------------
fprintf('\n--- 1. WPŁYW STANÓW PRZEJŚCIOWYCH SVF ---\n');
fprintf('Ze stanami przejściowymi:  T = %.4f, k = %.4f | MSE = %.6f | FIT = %.2f%%\n', ...
    T_z_przejsc, k_z_przejsc, calc_MSE(yW, y_sym_z_przejsc), calc_FIT(yW, y_sym_z_przejsc));
fprintf('Bez stanów przejściowych: T = %.4f, k = %.4f | MSE = %.6f | FIT = %.2f%%\n', ...
    T_bez_przejsc, k_bez_przejsc, calc_MSE(yW, y_sym_bez_przejsc), calc_FIT(yW, y_sym_bez_przejsc));

fprintf('\n--- 2. PORÓWNANIE POSTACI REGRESJI LINIOWEJ ---\n');
fprintf('Wariant 1 (Y = yF):       T = %.4f, k = %.4f | MSE = %.6f | FIT = %.2f%%\n', ...
    T_w1, k_w1, calc_MSE(yW, y_sym_w1), calc_FIT(yW, y_sym_w1));
fprintf('Wariant 2 (Y = yF^(1)):   T = %.4f, k = %.4f | MSE = %.6f | FIT = %.2f%%\n', ...
    T_w2, k_w2, calc_MSE(yW, y_sym_w2), calc_FIT(yW, y_sym_w2));

figure('Name', 'Analiza identyfikacji SVF');
subplot(2,1,1);
plot(tW, yW, 'k', 'LineWidth', 1.5); hold on;
plot(tW, y_sym_z_przejsc, 'r--', 'LineWidth', 1.2);
plot(tW, y_sym_bez_przejsc, 'b-.', 'LineWidth', 1.2);
grid on; title('1. Wpływ stanów przejściowych filtrów SVF (weryfikacja na Zw)');
xlabel('Czas [s]'); ylabel('Wyjście y');
legend('Dane rzeczywiste (Zw)', 'Ze stanami przejściowymi', 'Bez stanów przejściowych');

subplot(2,1,2);
plot(tW, yW, 'k', 'LineWidth', 1.5); hold on;
plot(tW, y_sym_w1, 'b--', 'LineWidth', 1.2);
plot(tW, y_sym_w2, 'g:', 'LineWidth', 1.5);
grid on; title('2. Porównanie postaci regresji: Y = y_F vs Y = y_F^{(1)} (weryfikacja na Zw)');
xlabel('Czas [s]'); ylabel('Wyjście y');
legend('Dane rzeczywiste (Zw)', 'Wariant 1: Y = y_F', 'Wariant 2: Y = y_F^{(1)}');



% dla nF = 1;
% dla TF=50*tp
% otrzymaliśmy k_O = 1.5166, spodziewane 2.0
% otrzymaliśmy T_O = 0.2904, spodziewane 0.5

% dla TF=50*tp
% otrzymaliśmy k_O = 1.9272, spodziewane 2.0
% otrzymaliśmy T_O = 0.4781, spodziewane 0.5

% dla TF=50*tp
% otrzymaliśmy k_O = 1.8831, spodziewane 2.0
% otrzymaliśmy T_O = 0.4660, spodziewane 0.5

%{
dla nF = 2

    0.3872
    1.7311

    0.5011      <-- najdokładniejsze
    1.9684      <-- najdokładniejsze

   0.2033
    0.9458
%}

%{
dla nF = 3
    0.4190
    1.8011

   0.5057
    1.9662

   -1.1380
   -6.4137
%}
