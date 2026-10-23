load ProcessStepResponse.mat
load NoisyProcessStepResponse.mat

t = S(:, 1);
y = S(:, 2);

war50 = 0.5 * y(end);
war90 = 0.9 * y(end);

idx50 = find(y >= war50, 1);
idx90 = find(y >= war90, 1);

T50 = t(idx50);
T90 = t(idx90);

T90_50 = T90 / T50;
disp('T90/T50: ');

% Przyjmujemy p=4
T = ((T50 / 3.67) + (T90/6.68)) / 2;

% Mamy już wymagane dane do modelowania transmitancji
K = y(end);
p = 4;
s = tf('s');
Gm1 = K / (T*s + 1)^p;
disp(Gm1);

% imp.

t_imp = t(1:600);
y_imp = diff(y) ./ diff(t);

figure;
plot(t, y);
hold on;
plot(t_imp, y_imp);
hold off;
title('S')

% ----- teraz dla noisy (nS) ----


t_noise = nS(:, 1);
y_noise = nS(:, 2);

war50_noise = 0.5 * y_noise(end);
war90_noise = 0.9 * y_noise(end);

idx50_noise = find(y >= war50_noise, 1);
idx90_noise = find(y >= war90_noise, 1);

T50_noise = t_noise(idx50_noise);
T90_noise = t_noise(idx90_noise);

T90_50_noise = T90_noise / T50_noise;  % 1.7177


% Przyjmujemy p=6
T_noise = ((T50_noise / 5.67) + (T90_noise/9.27)) / 2;

% Mamy już wymagane dane do modelowania transmitancji
K_noise = y_noise(end);
p_noise = 6;
Gm_noise = K_noise / (T_noise*s + 1)^p_noise;
disp(Gm_noise);

% imp.

t_imp_2 = t_noise(1:1200);
y_imp_2 = diff(y_noise) ./ diff(t_noise);

figure;
plot(t_noise, y_noise);
hold on;
plot(t_imp_2, y_imp_2);
hold off;
title('noisy S')

% punkt 2

