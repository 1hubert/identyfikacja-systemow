% Aproksymacja systemu modelem inercyjnym z opóźnieniem Gm2(s)

% --- CZĘŚĆ 1: Dane bez zakłóceń ---
load('ProcessStepResponse.mat'); % Wczytanie danych (macierz S)
t = S(:,1);                      % Wektor czasu
h = S(:,2);                      % Odpowiedź skokowa systemu

Au = 1;                          % Amplituda wymuszenia skokowego

% 1. Oszacowanie wzmocnienia statycznego (Wzór 2)
K_hat = h(end) / Au;             

% 2. Oszacowanie odpowiedzi impulsowej za pomocą różnic skończonych (Wzór 6)
dt = diff(t);
dh = diff(h);
g_hat = dh ./ dt;
g_hat = [0; g_hat];              % Wyrównanie długości wektora do t

% 3. Znalezienie punktu przegięcia (maksimum odpowiedzi impulsowej)
[a_g, idx_g] = max(g_hat);       % a_g to współczynnik kierunkowy stycznej
t_g = t(idx_g);                  % Współrzędna czasowa punktu przegięcia
s_g = h(idx_g);                  % Wartość odpowiedzi skokowej w pkt. przegięcia

% 4. Oszacowanie stałych czasowych T0 (opóźnienie) oraz T (inercja) (Wzory 4 i 5)
T0_hat = (a_g * t_g - s_g) / a_g;
T_hat = (Au * K_hat) / a_g;

% Wyświetlenie wyników dla danych czystych
fprintf('--- Dane bez zakłóceń ---\n');
fprintf('K = %.4f\n', K_hat);
fprintf('T0 = %.4f\n', T0_hat);
fprintf('T = %.4f\n', T_hat);

% Definicja modelu Gm2(s)
s = tf('s');
Gm2 = (K_hat / (T_hat * s + 1)) * exp(-T0_hat * s);

% Symulacja modelu
[y_Gm2, t_Gm2] = step(Gm2, t(end));
[g_Gm2, ~] = impulse(Gm2, t(end));

% Wykresy dla danych bez zakłóceń
figure('Name', 'Punkt 1.2 - Dane bez zakłóceń');
subplot(2,1,1);
plot(t, h, 'b', t_Gm2, y_Gm2, 'r--');
title('Odpowiedź skokowa - dane czyste');
legend('System h(t)', 'Model G_{m2}(s)');
grid on;

subplot(2,1,2);
plot(t, g_hat, 'b', t_Gm2, g_Gm2, 'r--');
title('Odpowiedź impulsowa - dane czyste');
legend('Oszacowanie \approx g(t)', 'Model G_{m2}(s)');
grid on;


% --- CZĘŚĆ 2: Dane zakłócone ---
load('NoisyProcessStepResponse.mat'); % Wczytanie danych zakłóconych (macierz nS)
t_n = nS(:,1);
h_n = nS(:,2);

% 1. Oszacowanie wzmocnienia statycznego (Wzór 2)
% UWAGA: W przypadku zaszumionych danych warto użyć średniej z ostatnich próbek
K_hat_n = mean(h_n(end-50:end)) / Au; 

% 2. Oszacowanie odpowiedzi impulsowej (Wzór 6)
dt_n = diff(t_n);
dh_n = diff(h_n);
g_hat_n = dh_n ./ dt_n;
g_hat_n = [0; g_hat_n];

% 3. Znalezienie punktu przegięcia dla zaszumionych danych
% Z powodu zakłóceń prosta różniczka bardzo zniekształca wynik
[a_g_n, idx_g_n] = max(g_hat_n); 
t_g_n = t_n(idx_g_n);
s_g_n = h_n(idx_g_n);

% 4. Oszacowanie stałych czasowych dla danych zaszumionych (Wzory 4 i 5)
T0_hat_n = (a_g_n * t_g_n - s_g_n) / a_g_n;
T_hat_n = (Au * K_hat_n) / a_g_n;

% Wyświetlenie wyników dla danych zakłóconych
fprintf('\n--- Dane zakłócone ---\n');
fprintf('K = %.4f\n', K_hat_n);
fprintf('T0 = %.4f\n', T0_hat_n);
fprintf('T = %.4f\n', T_hat_n);

% Wykresy dla danych zakłóconych
figure('Name', 'Punkt 1.2 - Dane zakłócone');
subplot(2,1,1);
plot(t_n, h_n, 'b');
title('Odpowiedź skokowa - dane zakłócone');
grid on;

subplot(2,1,2);
plot(t_n, g_hat_n, 'b');
title('Odpowiedź impulsowa oszacowana różnicami skończonymi');
grid on;