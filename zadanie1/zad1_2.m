clc;clear;close all;

global evals
evals = 0; %to pre eggholder, idk

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 3;                                      % pocet spusteni
popsize = 50;                                  % velkost populacie
ngen = 200;                                    % pocet generacii

% Parametre variantu E podla GA01
bestNums = [1,1];
oldCount = 10;
work1Count = 18;
work2Count = 20;
pmut_global = 0.15;
pmut_local = 0.15;
Amp_local = 5 * ones(1,nvars);

variants = ['e','f','g','h','i'];
averageHistAll = zeros(5,ngen);                % priemer kazdeho variantu

figure;
tiledlayout(2,3);

% Spustenie variantov e-i
for v = 1:5

    variant = variants(v);

    % Zapnutie/vypnutie genetickych operacii
    crossoverOn = true;
    mutxOn = true;
    mutaOn = true;

    switch variant
        case 'e'                                % povodny GA01
        case 'f'
            crossoverOn = false;                % bez krizenia
        case 'g'
            mutxOn = false;                     % bez globalnej mutacie
        case 'h'
            mutaOn = false;                     % bez lokalnej mutacie
        case 'i'
            mutxOn = false;                     % bez oboch mutacii
            mutaOn = false;
    end

    bestHistAll = zeros(runs,ngen);

    % Samotny geneticky algoritmus
    for r = 1:runs

        pop = genrpop(popsize,Space);           % vytvorenie pociatocnej populacie

        for g = 1:ngen

            f = testfn3(pop);                   % vyhodnotenie populacie
            %f = eggholder(pop); %EGGHOLDER
            bestHistAll(r,g) = min(f);          % najlepsia hodnota v generacii

            Best = selbest(pop,f,bestNums);     % 2 najlepsi jedinci
            Old = selrand(pop,f,oldCount);      % nahodne zachovani jedinci
            Work1 = selsus(pop,f,work1Count);   % skupina pre krizenie
            Work2 = selsus(pop,f,work2Count);   % skupina pre mutacie

            if crossoverOn
                Work1 = crossov(Work1,1,0);     % krizenie
            end

            if mutxOn
                Work2 = mutx(Work2,pmut_global,Space);          % globalna mutacia
            end

            if mutaOn
                Work2 = muta(Work2,pmut_local,Amp_local,Space); % lokalna mutacia
            end

            pop = [Best;Old;Work1;Work2];       % nova populacia

        end

        fprintf("Variant %c, Run %d: best = %.6f\n",variant,r,min(f));

    end

    averageHist = mean(bestHistAll,1);          % priemer troch behov
    averageHistAll(v,:) = averageHist;

    % Graf konkretneho variantu
    titles = ["Original", ...
        "Bez krizenia", ...
        "Bez globalnej mutacie", ...
        "Bez lokalnej mutacie", ...
        "Bez oboch mutacii"];

    nexttile;
    hold on;grid on;
    plot(1:ngen,bestHistAll(1,:),'LineWidth',1);
    plot(1:ngen,bestHistAll(2,:),'LineWidth',1);
    plot(1:ngen,bestHistAll(3,:),'LineWidth',1);
    plot(1:ngen,averageHist,'k','LineWidth',2.5);
    xlabel('Generacia');
    ylabel('Najlepsia hodnota');
    title(titles(v));
    legend('Beh 1','Beh 2','Beh 3','Priemer');
    hold off;

end

% Porovnanie priemerov vsetkych variantov
nexttile;
hold on;grid on;
plot(1:ngen,averageHistAll(1,:),'LineWidth',1.5);
plot(1:ngen,averageHistAll(2,:),'LineWidth',1.5);
plot(1:ngen,averageHistAll(3,:),'LineWidth',1.5);
plot(1:ngen,averageHistAll(4,:),'LineWidth',1.5);
plot(1:ngen,averageHistAll(5,:),'LineWidth',1.5);
xlabel('Generacia');
ylabel('Priemerna najlepsia hodnota');
title('Porovnanie');
legend('Original', ...
       'Bez krizenia', ...
       'Bez globalnej mutacie', ...
       'Bez lokalnej mutacie', ...
       'Bez oboch mutacii');
hold off;