clc; close all; clear all;

model = 'zad6_sim.slx';

%Uloha 2 - Prevodova charakteristika
y_ciel = 0:0.5:4;
u_hodnoty = 10*y_ciel + 7.5*y_ciel.^3;
y_ust = zeros(size(y_ciel));

for i = 1:length(u_hodnoty)
    U = u_hodnoty(i);
    out = sim(model, 'StopTime', '20');
    y = out.get('yout');
    y_ust(i) = y.Data(end);
end

figure;
plot(u_hodnoty, y_ust, 'o-');
xlabel('Vstup u');
ylabel('Ustaleny vystup y');
title('Prevodova charakteristika');
grid on;

% Uloha 3 - Casove odozvy
y_ciel = [0.5 1 2 3 4];

figure;
hold on;

for i = 1:5
    U = 10*y_ciel(i) + 7.5*y_ciel(i)^3;
    out = sim(model, 'StopTime', '20');
    y = out.get('yout');

    plot(y.Time, y.Data, ...
        'DisplayName', sprintf('u = %.2f', U));
end

xlabel('Cas [s]');
ylabel('Vystup y');
title('Casove odozvy na 5 skokov');
legend;
grid on;
