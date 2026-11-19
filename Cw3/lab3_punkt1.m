load("IdentWsadowaStat.mat");

u_DaneStatW = DaneStatW(:, 1);
y_DaneStatW = DaneStatW(:, 2);

u_DaneStatC = DaneStatC(:, 1);
y_DaneStatC = DaneStatC(:, 2);

figure;
plot(u_DaneStatW, y_DaneStatW);
title('DaneStatW: pomiary zakłócone zakłóceniem nieskorelowanym e');

figure;
plot(u_DaneStatC, y_DaneStatC);
title('DaneStatC: pomiary zakłócone zakłóceniem skorelowanym v');


% -------------------------------------------------------------------------
% 1) dla DaneStatW
% -------------------------------------------------------------------------

u_DaneStatW = u_DaneStatW(1:1:end);
y_DaneStatW = y_DaneStatW(1:1:end);

% N - liczba próbek, u i y to wektory kolumnowe danych
N = length(u_DaneStatW);
disp("N dla DaneStatW: ")
disp(N);
% 1. Konstrukcja macierzy regresji Phi
Phi_W = [ones(N, 1), 1./u_DaneStatW, 1./(u_DaneStatW.^2), 1./(u_DaneStatW.^3)];



% 2. Obliczenie parametrów estymatora LS wbudowanym operatorem
p_LS = Phi_W \ y_DaneStatW;

y_m_W = zeros(1, N);
for j = 1:N
    p_sum = p_LS(1);
    for i = 2:4
        p_sum = p_LS(i) / u_DaneStatW(j)^(i-1);
    end
    y_m_W(j) = p_sum;
end


figure;
hold on;
plot(u_DaneStatW, y_DaneStatW);
plot(u_DaneStatW, y_m_W);
hold off;
title("DaneStatW");


% -------------------------------------------------------------------------
% 2) dla DaneStatC
% -------------------------------------------------------------------------

u_DaneStatC = u_DaneStatC(1:1:end);
y_DaneStatC = y_DaneStatC(1:1:end);

% N - liczba próbek, u i y to wektory kolumnowe danych
N = length(u_DaneStatC);
disp("N dla DaneStatC: ")
disp(N);
% 1. Konstrukcja macierzy regresji Phi
Phi = [ones(N, 1), 1./u_DaneStatC, 1./(u_DaneStatC.^2), 1./(u_DaneStatC.^3)];



% 2. Obliczenie parametrów estymatora LS wbudowanym operatorem
p_LS = Phi \ y_DaneStatC;

y_m = zeros(1, N);
for j = 1:N
    p_sum = p_LS(1);
    for i = 2:4
        p_sum = p_LS(i) / u_DaneStatC(j)^(i-1);
    end
    y_m(j) = p_sum;
end


figure;
hold on;
plot(u_DaneStatC, y_DaneStatC);
plot(u_DaneStatC, y_m);
hold off;
title("DaneStatC");

% -------------------------------------------------------------------------
% 3) oszacować macierz kowariancji dla DaneStatW
% -------------------------------------------------------------------------
N = length(u_DaneStatW);

epsilon = y_DaneStatW - y_m_W;
sigma_2 = 1 / (N - 4) * epsilon.^2;
disp(size(inv(Phi_W'*Phi_W)));
cov = sigma_2 .* inv(Phi_W'*Phi_W);