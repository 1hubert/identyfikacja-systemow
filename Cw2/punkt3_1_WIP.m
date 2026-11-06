% WIP - wersja Huberta

% semilogx - wymuszenie log skali pulsacji
load("./tabledata.txt")
tabledata = tabledata(2:end, :);

u = tabledata(:, 1);
y = tabledata(:, 2);

figure;
plot(u);
title("u(t)")

figure;
plot(y);
title("y(t)")

% estymator transmitancji ze wzoru (15)
N = 1000;
omega = (2 * pi) / N;
omega_array = (0:1:N-1) * omega;

frequencyResponse = fft(y) ./ fft(u);

figure;
semilogx(omega_array, abs(frequencyResponse));

% estymator ze wzoru (16)
gamma = N - 1;
M_w = 200;

omega_array = (2 * pi * k) / N;

yu_suma = 0;
for tau = -gamma:1:gamma
    yu_suma = yu_suma + xcorr(y, u) * 
end
