
clc;clear;close all;

% Uvodne nastavenia
nvars = 10;                     % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];   % definicny obor
runs = 5;                       % pocet behov
seeds = [1 2 3 4 5];            % rovnake seedy pre oba algoritmy
ngen = 2000;                    % pocet iteracii
popSize = 50;                   % pocet jedincov

% PSO
c0 = 0.85;                      % zotrvacnost castice
c1 = 1.5;                       % individualny faktor
c2 = 1.5;                       % socialny faktor
vMax = 200*ones(1,nvars);       % max. rychlost

% SGA
pmut_global = 0.3;
pmut_local = 0.3;
Amp_local = 5*ones(1,nvars);
sgaSize = 50;
bestNumsSGA = ones(1,2);
RandOldCountSGA = 10;
work1CountSGA = 18;
work2CountSGA = 20;

% Ulozenie vysledkov
% Prvy stlpec je inicializacia, dalsie stlpce su iteracie

histPSO = zeros(runs,ngen+1);
histSGA = zeros(runs,ngen+1);
finalPSO = zeros(runs,1);
finalSGA = zeros(runs,1);
bestVectorsPSO = zeros(runs,nvars);
bestVectorsSGA = zeros(runs,nvars);

% 5 behov

for r = 1:runs

    fprintf('\nBeh %d, seed %d\n',r,seeds(r));

    % Spolocna pociatocna populacia
    rng(seeds(r));
    initialPop = genrpop(popSize,Space);

    % Ulozenie stavu generatora po inicializacii
    initialRngState = rng;

    % PSO
    rng(initialRngState);
    x = initialPop;                             % aktualne polohy castic

    % Nahodne pociatocne rychlosti v intervale [-200,200]
    v = (2*rand(popSize,nvars)-1).*vMax;

    % Inicializacia pbest a gbest
    fPSO = testfn3(x);
    fPSO = fPSO(:);

    pbest = x;
    pbestFitness = fPSO;

    [gbestFitness,bestIndex] = min(pbestFitness);
    gbest = pbest(bestIndex,:);

    histPSO(r,1) = gbestFitness;

    for g = 1:ngen

        % Aktualizacia rychlosti podla prednasky
        % c0 - zotrvacnost, c1 - individualny, c2 - socialny faktor

        v = c0*v ...
            + c1*rand(popSize,nvars).*(pbest-x) ...
            + c2*rand(popSize,nvars).*(gbest-x);

        % Obmedzenie rychlosti na interval [-vMax,vMax]
        v = max(min(v,vMax),-vMax);

        % Aktualizacia polohy
        x = x + v;

        % Osetrenie hranic
        outside = x < Space(1,:) | x > Space(2,:);
        x = max(min(x,Space(2,:)),Space(1,:));
        v(outside) = -v(outside);

        % Vyhodnotenie novych poloh
        fPSO = testfn3(x);
        fPSO = fPSO(:);

        % Aktualizacia pbest
        % Castica si zapamata svoju novu polohu, iba ak je lepsia
        better = fPSO < pbestFitness;
        pbest(better,:) = x(better,:);
        pbestFitness(better) = fPSO(better);

        % Aktualizacia gbest
        % Najlepsia doteraz najdena poloha celeho krdla
        [currentBest,bestIndex] = min(pbestFitness);

        if currentBest < gbestFitness
            gbestFitness = currentBest;
            gbest = pbest(bestIndex,:);
        end

        % Najlepsia doteraz dosiahnuta fitness
        histPSO(r,g+1) = gbestFitness;
    end

    finalPSO(r) = gbestFitness;
    bestVectorsPSO(r,:) = gbest;
    fprintf('PSO = %.6f\n',finalPSO(r));

    % SGA

    rng(initialRngState);
    popSGA = initialPop;
    fSGA = testfn3(popSGA);
    fSGA = fSGA(:);

    [bestFitnessSGA,bestIndex] = min(fSGA);
    bestVectorSGA = popSGA(bestIndex,:);

    histSGA(r,1) = bestFitnessSGA;

    for g = 1:ngen

        % Selekcia
        Best = selbest(popSGA,fSGA,bestNumsSGA);
        RandOld = selrand(popSGA,fSGA,RandOldCountSGA);
        Work1 = selsus(popSGA,fSGA,work1CountSGA);
        Work2 = selsus(popSGA,fSGA,work2CountSGA);

        % Krizenie
        Work1 = crossov(Work1,1,0);

        % Mutacia
        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        % Nova populacia
        popSGA = [Best;RandOld;Work1;Work2];

        % Vyhodnotenie novej populacie
        fSGA = testfn3(popSGA);
        fSGA = fSGA(:);

        % Ulozenie najlepsieho doteraz dosiahnuteho riesenia
        [currentBest,bestIndex] = min(fSGA);

        if currentBest < bestFitnessSGA
            bestFitnessSGA = currentBest;
            bestVectorSGA = popSGA(bestIndex,:);
        end

        histSGA(r,g+1) = bestFitnessSGA;
    end

    finalSGA(r) = bestFitnessSGA;
    bestVectorsSGA(r,:) = bestVectorSGA;
    fprintf('SGA = %.6f\n',finalSGA(r));
end

% Priemerne priebehy 5 behov

averagePSO = mean(histPSO,1);
averageSGA = mean(histSGA,1);
iterations = 0:ngen;

% Spolocna mierka grafov
allFitness = [
    histPSO(:);
    histSGA(:)
];

commonYLim = [min(allFitness) max(allFitness)];

% Farby algoritmov
% PSO - modra, SGA - oranzova
% Jednotlive behy svetle, priemery vyrazne

colorPSO = [0 0.4470 0.7410];
colorSGA = [0.8500 0.3250 0.0980];

lightPSO = 0.65*colorPSO + 0.35*[1 1 1];
lightSGA = 0.65*colorSGA + 0.35*[1 1 1];

% Graf 1 - PSO
fig1 = figure;

hold on;grid on;

for r = 1:runs
    plot(iterations,histPSO(r,:), ...
        'Color',lightPSO,'LineWidth',1);
end

plot(iterations,averagePSO, ...
    'Color',colorPSO,'LineWidth',2.5);
xlabel('Iteracia');
ylabel('Fitness');
title('PSO: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    sprintf('Seed %d: %.2f',seeds(1),finalPSO(1)), ...
    sprintf('Seed %d: %.2f',seeds(2),finalPSO(2)), ...
    sprintf('Seed %d: %.2f',seeds(3),finalPSO(3)), ...
    sprintf('Seed %d: %.2f',seeds(4),finalPSO(4)), ...
    sprintf('Seed %d: %.2f',seeds(5),finalPSO(5)), ...
    sprintf('Priemer: %.2f',averagePSO(end)), ...
    'Location','best');
ylim(commonYLim);

hold off;

% Graf 2 - porovnanie PSO a SGA

fig2 = figure;
hold on;grid on;

for r = 1:runs
    plot(iterations,histPSO(r,:), ...
        'Color',lightPSO,'LineWidth',1);
end

plot(iterations,averagePSO, ...
    'Color',colorPSO,'LineWidth',2.5);

for r = 1:runs
    plot(iterations,histSGA(r,:), ...
        'Color',lightSGA,'LineWidth',1);
end

plot(iterations,averageSGA, ...
    'Color',colorSGA,'LineWidth',2.5);

xlabel('Iteracia');
ylabel('Fitness');
title('Porovnanie PSO a SGA na Schwefelovej funkcii');

legend( ...
    sprintf('PSO seed %d: %.2f',seeds(1),finalPSO(1)), ...
    sprintf('PSO seed %d: %.2f',seeds(2),finalPSO(2)), ...
    sprintf('PSO seed %d: %.2f',seeds(3),finalPSO(3)), ...
    sprintf('PSO seed %d: %.2f',seeds(4),finalPSO(4)), ...
    sprintf('PSO seed %d: %.2f',seeds(5),finalPSO(5)), ...
    sprintf('PSO priemer: %.2f',averagePSO(end)), ...
    sprintf('SGA seed %d: %.2f',seeds(1),finalSGA(1)), ...
    sprintf('SGA seed %d: %.2f',seeds(2),finalSGA(2)), ...
    sprintf('SGA seed %d: %.2f',seeds(3),finalSGA(3)), ...
    sprintf('SGA seed %d: %.2f',seeds(4),finalSGA(4)), ...
    sprintf('SGA seed %d: %.2f',seeds(5),finalSGA(5)), ...
    sprintf('SGA priemer: %.2f',averageSGA(end)), ...
    'Location','best');

ylim(commonYLim);
hold off;

% Tabulka vysledkov 5 behov + priemer + smerodajna odchylka
Beh = [
    compose("Seed %d",seeds(:));
    "Priemer";
    "Smerodajna odchylka"
];

PSO = [
    finalPSO;
    mean(finalPSO);
    std(finalPSO)
];

SGA = [
    finalSGA;
    mean(finalSGA);
    std(finalSGA)
];

resultsTable = table(Beh,PSO,SGA);

fprintf('\n');
fprintf('============================================\n');
fprintf('TABULKA VYSLEDKOV\n');
fprintf('============================================\n');

disp(resultsTable);

% Najlepsi vektor kazdeho algoritmu

[bestFitnessPSO,bestRunPSO] = min(finalPSO);
[bestFitnessSGA,bestRunSGA] = min(finalSGA);

bestVectors = [
    bestVectorsPSO(bestRunPSO,:);
    bestVectorsSGA(bestRunSGA,:)
];

Algoritmus = ["PSO";"SGA"];

Seed = [
    seeds(bestRunPSO);
    seeds(bestRunSGA)
];

Fitness = [
    bestFitnessPSO;
    bestFitnessSGA
];

vectorNames = cellstr(compose("x%d",1:nvars));

vectorTable = array2table(bestVectors,'VariableNames',vectorNames);

bestTable = [
    table(Algoritmus,Seed,Fitness) ...
    vectorTable
];

fprintf('\n');
fprintf('============================================\n');
fprintf('NAJLEPSIE NAJDENE VECTORY\n');
fprintf('============================================\n');

disp(bestTable);
