load("StochasticProcess.mat")

M = 1001;
N = 500;

t = StochasticProcess(1, :);
y1 = StochasticProcess(2, :);
y2 = StochasticProcess(3, :);
y3 = StochasticProcess(4, :);

figure(1);
plot(t, y1);
hold on;
plot(t, y2);
hold on;
plot(t, y3);
hold off;
title('Trzy wybrane realizacje w funkcji czasu');

% estymaty m_n oraz sigma_n^2 po realizacjach
A = StochasticProcess(2:end, :);

lista_m_n = [];
lista_sigma_n_2 = [];
for i=1:M
    m_n = mean(A(:, i));
    sigma_n_2 = 0;
    for j=1:N
        sigma_n_2  = sigma_n_2 + (A(j, i) - m_n)^2;
    end
    sigma_n_2 = sigma_n_2 / N;

    lista_m_n = [lista_m_n; m_n];
    lista_sigma_n_2 = [lista_sigma_n_2; sigma_n_2];
end

% estymaty m oraz sigma^2 po czasie
lista_m = [];
lista_sigma_2 = [];
for i=1:N
    m = mean(A(i, :));
    sigma_2 = 0;
    for j=1:M
        sigma_2  = sigma_2 + (A(i, j) - m)^2;
    end
    sigma_2 = sigma_2 / M;

    lista_m = [lista_m; m];
    lista_sigma_2 = [lista_sigma_2; sigma_2];
end


% plotting
figure(2);
plot(t, lista_m_n);
hold on;
plot(0:499, lista_m);
hold off;
title('Średnie')

figure(3);
plot(t, lista_sigma_n_2);
hold on;
plot(0:499, lista_sigma_2);
hold off;
title('Wariancje')

% porównanie średnich z list
disp('Średnia średnich po realizacjach: ');
disp(mean(lista_m_n));
disp('Średnia średnich po czasie: ');
disp(mean(lista_m));
