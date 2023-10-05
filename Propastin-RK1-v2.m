clear force all
close all

c = 3e8;
f = 1e10;           % Частота
lambda = c/f;       % Длина волны
K = 10;             % Количество элементов АР
k = 2*pi/lambda;    % Волновое число
d = lambda/2;       % Шаг элементов АР

n = 1:K;
x = (2/(K-1))*(n-1) - 1;
y = 1 + 5*cos(pi/2*x).^2;

phi = n./n;
otk = 15; % Отклонение луча АР

phi(1) = 0;
% Фазовое распределения
for i = 2:K
    phi(i) = phi(i - 1) + 2*pi*d*sin(0.017*otk)/lambda;
end

A = y./(5 + 1).*exp(-1.i*phi); % Амплитудное распределение
plot(1:K, abs(A)); % График амплитудного распределения
hold on;

% dz - на каком расстоянии строить поле
dz1 = lambda/2;
dz2 = lambda*5;

% Подготовка к расчётам
x = n*d;
xz = x;

% Первый рассчёт на расстоянии dz1
A1 = A;
x = n*d;
xz = x;
temp = 0;
for i_xz = 1:length(xz)
    for i_x = 1:length(x)
        R = sqrt((xz(i_xz) - x(i_x))^2 + dz1^2);
        temp = temp - A(i_x)*exp(1i*k*R)/(4*pi*R);
    end
    A1(i_xz) = temp;
    temp = 0;
end

% Первый рассчёт на расстоянии dz2
A2 = A;
temp = 0;
for i_xz = 1:length(xz)
    for i_x = 1:length(x)
        R = sqrt((xz(i_xz) - x(i_x))^2 + dz2^2);
        temp = temp - A(i_x)*exp(1i*k*R)/(4*pi*R);
    end
    A2(i_xz) = temp;
    temp = 0;
end

plot(n, abs(A1)./max(abs(A1))); % Нормированный график для dz1
hold on;
plot(n, abs(A2)./max(abs(A2))); % Нормированный график для dz2
hold on;
grid on;

% Проверка правильности расчёта через нахождение отклонения луча
izm_otk = 2.5; % Отклонение луча по графику dz2 в элементах
test = sqrt((izm_otk*d)^2 + dz2^2);
disp(asin((izm_otk*d)/test)*57);