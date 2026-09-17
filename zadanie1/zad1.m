clc;clear;close all;

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 3;                                      % pocet spusteni
popsize = 50;                                  % velkost populacie
ngen = 200;                                    % pocet generacii

variants = ['a','b','c','d','e'];
averageHistAll = zeros(5,ngen);                % priemer kazdeho variantu

figure;
tiledlayout(2,3);

% Spustenie vsetkych variantov
for v = 1:5
    variant = variants(v);

    % Nastavenie pripadov a-e
    switch variant
        case 'a' % velky selektivny tlak + velka diverzita
            bestNums = [6,4];                    % viac kopii najlepsich -> vyssi selektivny tlak
            oldCount = 1;                        % ziadny nahodny vyber -> vyssi selektivny tlak
            work1Count = 17;
            work2Count = 22;
            selectionName = "selsus";            % vyssi selektivny tlak
            pmut_global = 0.30;                  % silna globalna mutacia -> vyssia diverzita
            pmut_local = 0.30;                   % viac lokalnych mutacii -> vyssia diverzita
            Amp_local = 100 * ones(1,nvars);     % vacsie zmeny pri muta -> vyssia diverzita
    
        case 'b' % velky selektivny tlak + mala diverzita
            bestNums = [6,4];                    % rovnaky vysoky selektivny tlak ako A
            oldCount = 1;
            work1Count = 17;
            work2Count = 22;
            selectionName = "selsus";
            pmut_global = 0.005;                 % skoro ziadna globalna mutacia -> mala diverzita
            pmut_local = 0.005;                  % skoro ziadna lokalna mutacia -> mala diverzita
            Amp_local = 1 * ones(1,nvars);       % velmi male lokalne zmeny -> mala diverzita
    
        case 'c' % maly selektivny tlak + velka diverzita
            bestNums = [1,1];                    % iba 2 elity -> nizsi selektivny tlak
            oldCount = 20;                       % vela nahodne vybranych -> nizsi tlak a vyssia diverzita
            work1Count = 14;
            work2Count = 14;
            selectionName = "seltourn";          % nizsi selektivny tlak ako selsus
            pmut_global = 0.30;                  % silna globalna mutacia -> vyssia diverzita
            pmut_local = 0.30;                   % silna lokalna mutacia -> vyssia diverzita
            Amp_local = 100 * ones(1,nvars);     % vacsie lokalne zmeny -> vyssia diverzita
    
        case 'd' % maly selektivny tlak + mala diverzita
            bestNums = [1,1];                    % rovnaky nizky selektivny tlak ako C
            oldCount = 20;
            work1Count = 14;
            work2Count = 14;
            selectionName = "seltourn";
            pmut_global = 0.005;                 % skoro ziadna globalna mutacia -> mala diverzita
            pmut_local = 0.005;                  % skoro ziadna lokalna mutacia -> mala diverzita
            Amp_local = 1 * ones(1,nvars);       % velmi male lokalne zmeny -> mala diverzita
    
        case 'e' % kompromis - podla GA01
            bestNums = [1,1];
            oldCount = 10;
            work1Count = 18;
            work2Count = 20;
            selectionName = "selsus";
            pmut_global = 0.15;
            pmut_local = 0.15;
            Amp_local = 5 * ones(1,nvars);
    end

    bestHistAll = zeros(runs,ngen);

    % Samotny geneticky algoritmus
    for r = 1:runs

        pop = genrpop(popsize,Space);           % vytvorenie pociatocnej populacie

        for g = 1:ngen

            f = testfn3(pop);                   % vyhodnotenie populacie
            bestHistAll(r,g) = min(f);          % najlepsia hodnota v generacii
            Best = selbest(pop,f,bestNums);     % najlepsi jedinci
            Old = selrand(pop,f,oldCount);      % nahodne zachovani jedinci

            % selekcia pracovnych skupin
            %work1 prechadza iba cez krizenie, work2 iba cez mutaciu.
            switch selectionName
                case "selsus"
                    Work1 = selsus(pop,f,work1Count);
                    Work2 = selsus(pop,f,work2Count);

                case "seltourn"
                    Work1 = seltourn(pop,f,work1Count);
                    Work2 = seltourn(pop,f,work2Count);
            end

            Work1 = crossov(Work1,1,0);                         % krizenie
            Work2 = mutx(Work2,pmut_global,Space);              % globalna mutacia
            Work2 = muta(Work2,pmut_local,Amp_local,Space);     % lokalna mutacia

            pop = [Best;Old;Work1;Work2];                       % nova populacia

        end

        fprintf("Variant %c, Run %d: best = %.6f\n",variant,r,min(f));

    end

    averageHist = mean(bestHistAll,1);           % priemer troch behov
    averageHistAll(v,:) = averageHist;

    % Graf konkretneho variantu
    nexttile;
    hold on;grid on;
    plot(1:ngen,bestHistAll(1,:),'LineWidth',1);
    plot(1:ngen,bestHistAll(2,:),'LineWidth',1);
    plot(1:ngen,bestHistAll(3,:),'LineWidth',1);
    plot(1:ngen,averageHist,'k','LineWidth',2.5);
    xlabel('Generacia');
    ylabel('Najlepsia hodnota');
    title(['Pripad ',upper(variant)]);
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
title('Porovnanie A-E');
legend('A','B','C','D','E');
hold off;