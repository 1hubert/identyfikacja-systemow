figure;
subplot(2,1,1);
stairs(Zdata(:,1));
title('Sygnał wejściowy u(n)');

subplot(2,1,2);
stairs(Zdata(:,2));
title('Sygnał wyjściowy y(n)');