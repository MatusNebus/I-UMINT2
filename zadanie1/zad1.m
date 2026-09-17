clc;
clear;
close all;

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 3;                                      % pocet spusteni
popsize = 50;                                  % velkost populacie
ngen = 200;                                    % pocet generacii
eliteCount = 2;                                % elitizmus

% Vyber pripadu
variant = 'e';

% Nastavenie pripadov a-e
switch variant
    case 'a' %velky selektivny tlak, velka diverzita
        selectionName = "selsus";
        pmut_global = 0.10;
        pmut_local = 0.10;
        Amp_local = 100 * ones(1,nvars);

    case 'b' % velky selektivny tlak + mala diverzita
        selectionName = "selsus";
        pmut_global = 0.01;
        pmut_local = 0.01;
        Amp_local = 10 * ones(1,nvars);

    case 'c' % maly selektivny tlak + velka diverzita
        selectionName = "seltourn";
        pmut_global = 0.10;
        pmut_local = 0.10;
        Amp_local = 100 * ones(1,nvars);

    case 'd' % maly selektivny tlak + mala diverzita
        selectionName = "seltourn";
        pmut_global = 0.01;
        pmut_local = 0.01;
        Amp_local = 10 * ones(1,nvars);

    case 'e' % kompromis
        selectionName = "selsus";
        pmut_global = 0.03;
        pmut_local = 0.05;
        Amp_local = 50 * ones(1,nvars);
end

bestHistAll = zeros(runs, ngen); % Ulozenie priebehov

%samotny geneticky algoritmus:
for r = 1:runs
    pop = genrpop(popsize, Space);  % vytvorenie pociatocnej populacie
    f = zeros(popsize,1);           % vyhodnotenie populacie

    for i = 1:popsize
        f(i) = testfn3(pop(i,:));
    end

    for g = 1:ngen
        % vyber eliteCount najlepsich jedincov
        elites = selsort(pop, f, eliteCount);
        % elites = selbest(pop, f, eliteCount);

        nKids = popsize - eliteCount; % pocet novych jedincov

        % selekcia rodicov
        switch selectionName
            case "selsus"
                parents = selsus(pop, f, nKids);
            case "seltourn"
                parents = seltourn(pop, f, nKids);
        end

        kids = crossov(parents, 2, 0);                  % krizenie
        kids = mutx(kids, pmut_global, Space);          % globalna mutacia
        kids = muta(kids, pmut_local, Amp_local, Space); % lokalna mutacia
        pop = [elites; kids];                            % nova populacia

        % vyhodnotenie novej populacie
        for i = 1:popsize
            f(i) = testfn3(pop(i,:));
        end

        bestHistAll(r,g) = min(f); %najlepsia hodnota v generacii
    end

    fprintf("Run %d: best = %.6f\n", r, min(f));
end

averageHist = mean(bestHistAll, 1); % Priemer troch behov


% Graf
figure;
hold on;
grid on;

plot(1:ngen, bestHistAll(1,:), 'LineWidth', 1);
plot(1:ngen, bestHistAll(2,:), 'LineWidth', 1);
plot(1:ngen, bestHistAll(3,:), 'LineWidth', 1);

plot(1:ngen, averageHist, 'k', 'LineWidth', 2.5);

xlabel('Generacia');
ylabel('Najlepsia hodnota funkcie');

title(['Pripad ', upper(variant), ' - Schwefelova funkcia']);

legend('Beh 1', 'Beh 2', 'Beh 3', 'Priemer');

hold off;