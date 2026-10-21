
clc;clear;close all;

% Uvodne nastavenia
nvars = 10;                     % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];   % definicny obor
runs = 5;                       % pocet behov
seeds = [1 2 3 4 5];            % rovnake seedy pre oba algoritmy
ngen = 200;                    % pocet generacii
popSize = 50;                   % pocet jedincov

% DE
F = 0.25;                        % faktor diferencialnej mutacie
CR = 0.05;                       % pravdepodobnost krizenia

%nastavenia SGA
pmut_global = 0.3;
pmut_local = 0.3;
Amp_local = 5*ones(1,nvars);
bestNumsSGA = ones(1,2);
RandOldCountSGA = 10;
work1CountSGA = 18;
work2CountSGA = 20;
popSize = 50;
ngen = 200;

% Ulozenie vysledkov
% Prvy stlpec je inicializacia, dalsie stlpce su generacie

histDE = zeros(runs,ngen+1);
histSGA = zeros(runs,ngen+1);
finalDE = zeros(runs,1);
finalSGA = zeros(runs,1);
bestVectorsDE = zeros(runs,nvars);
bestVectorsSGA = zeros(runs,nvars);

% 5 behov

for r = 1:runs

    fprintf('\nBeh %d, seed %d\n',r,seeds(r));

    % Spolocna pociatocna populacia
    rng(seeds(r));
    initialPop = genrpop(popSize,Space);

    initialRngState = rng; % Ulozenie stavu generatora po inicializacii

    % DE
    rng(initialRngState);
    popDE = initialPop;
    fDE = testfn3(popDE);
    fDE = fDE(:);

    [bestFitnessDE,bestIndex] = min(fDE);
    bestVectorDE = popDE(bestIndex,:);
    histDE(r,1) = bestFitnessDE;

    for g = 1:ngen

        % Populacia skusobnych jedincov
        trialPop = popDE;

        for i = 1:popSize

            % Vyber troch nahodnych jedincov
            candidates = [1:i-1 i+1:popSize];
            selected = candidates(randperm(popSize-1,3));

            r1 = selected(1);
            r2 = selected(2);
            r3 = selected(3);

            % vytvorim mutanta podla vzorca z prednasky:
            v = popDE(r1,:) + F*(popDE(r2,:) - popDE(r3,:));

            % Binomicke krizenie:
            mask = rand(1,nvars) < CR;
            jrand = randi(nvars);         % Aspon jeden gen musi pochadzat z mutanta podla zadania
            mask(jrand) = true;

            trialPop(i,mask) = v(mask);

            % Osetrenie hranic
            trialPop(i,:) = max(min(trialPop(i,:),Space(2,:)),Space(1,:));

        end

        % Vyhodnotenie skusobnych jedincov
        fTrial = testfn3(trialPop);
        fTrial = fTrial(:);

        % Vyber lepsieho jedinca
        better = fTrial <= fDE;
        popDE(better,:) = trialPop(better,:);
        fDE(better) = fTrial(better);

        % Ulozenie najlepsieho doteraz dosiahnuteho riesenia
        [currentBest,bestIndex] = min(fDE);

        if currentBest < bestFitnessDE
            bestFitnessDE = currentBest;
            bestVectorDE = popDE(bestIndex,:);
        end

        histDE(r,g+1) = bestFitnessDE;
    end

    finalDE(r) = bestFitnessDE;
    bestVectorsDE(r,:) = bestVectorDE;
    fprintf('DE = %.6f\n',finalDE(r));

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

averageDE = mean(histDE,1);
averageSGA = mean(histSGA,1);
iterations = 0:ngen;

% Spolocna mierka grafov
allFitness = [
    histDE(:);
    histSGA(:)
];

commonYLim = [min(allFitness) max(allFitness)];

% Farby algoritmov
% DE - modra, SGA - oranzova
% Jednotlive behy svetle, priemery vyrazne

colorDE = [0 0.4470 0.7410];
colorSGA = [0.8500 0.3250 0.0980];

lightDE = 0.65*colorDE + 0.35*[1 1 1];
lightSGA = 0.65*colorSGA + 0.35*[1 1 1];

% Graf 1 - DE

fig1 = figure;
hold on;grid on;

for r = 1:runs
    plot(iterations,histDE(r,:), ...
        'Color',lightDE,'LineWidth',1);
end

plot(iterations,averageDE, ...
    'Color',colorDE,'LineWidth',2.5);

xlabel('Generacia');
ylabel('Fitness');
title('DE: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    sprintf('Seed %d: %.2f',seeds(1),finalDE(1)), ...
    sprintf('Seed %d: %.2f',seeds(2),finalDE(2)), ...
    sprintf('Seed %d: %.2f',seeds(3),finalDE(3)), ...
    sprintf('Seed %d: %.2f',seeds(4),finalDE(4)), ...
    sprintf('Seed %d: %.2f',seeds(5),finalDE(5)), ...
    sprintf('Priemer: %.2f',averageDE(end)), ...
    'Location','best');

xlim([0 ngen]);
ylim(commonYLim);
hold off;

% Graf 2 - porovnanie DE a SGA

fig2 = figure;
hold on;grid on;

for r = 1:runs
    plot(iterations,histDE(r,:), ...
        'Color',lightDE,'LineWidth',1);
end

plot(iterations,averageDE, ...
    'Color',colorDE,'LineWidth',2.5);

for r = 1:runs
    plot(iterations,histSGA(r,:), ...
        'Color',lightSGA,'LineWidth',1);
end

plot(iterations,averageSGA, ...
    'Color',colorSGA,'LineWidth',2.5);

xlabel('Generacia');
ylabel('Fitness');
title('Porovnanie DE a SGA na Schwefelovej funkcii');

legend( ...
    sprintf('DE seed %d: %.2f',seeds(1),finalDE(1)), ...
    sprintf('DE seed %d: %.2f',seeds(2),finalDE(2)), ...
    sprintf('DE seed %d: %.2f',seeds(3),finalDE(3)), ...
    sprintf('DE seed %d: %.2f',seeds(4),finalDE(4)), ...
    sprintf('DE seed %d: %.2f',seeds(5),finalDE(5)), ...
    sprintf('DE priemer: %.2f',averageDE(end)), ...
    sprintf('SGA seed %d: %.2f',seeds(1),finalSGA(1)), ...
    sprintf('SGA seed %d: %.2f',seeds(2),finalSGA(2)), ...
    sprintf('SGA seed %d: %.2f',seeds(3),finalSGA(3)), ...
    sprintf('SGA seed %d: %.2f',seeds(4),finalSGA(4)), ...
    sprintf('SGA seed %d: %.2f',seeds(5),finalSGA(5)), ...
    sprintf('SGA priemer: %.2f',averageSGA(end)), ...
    'Location','best');

xlim([0 ngen]);
ylim(commonYLim);
hold off;

% Tabulka pouzitych parametrov DE
parametersTable = table(F,CR);
fprintf('\n');
fprintf('============================================\n');
fprintf('PARAMETRE DE\n');
fprintf('============================================\n');
disp(parametersTable);

% Tabulka vysledkov 5 behov + priemer + smerodajna odchylka

Beh = [
    compose("Seed %d",seeds(:));
    "Priemer";
    "Smerodajna odchylka"
];

DE = [
    finalDE;
    mean(finalDE);
    std(finalDE)
];

SGA = [
    finalSGA;
    mean(finalSGA);
    std(finalSGA)
];

resultsTable = table(Beh,DE,SGA);

fprintf('\n');
fprintf('============================================\n');
fprintf('TABULKA VYSLEDKOV\n');
fprintf('============================================\n');

disp(resultsTable);

% Najlepsi vektor kazdeho algoritmu

[bestFitnessDE,bestRunDE] = min(finalDE);
[bestFitnessSGA,bestRunSGA] = min(finalSGA);

bestVectors = [
    bestVectorsDE(bestRunDE,:);
    bestVectorsSGA(bestRunSGA,:)
];

Algoritmus = ["DE";"SGA"];

Seed = [
    seeds(bestRunDE);
    seeds(bestRunSGA)
];

Fitness = [
    bestFitnessDE;
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
fprintf('NAJLEPSIE NAJDENE VEKTORY\n');
fprintf('============================================\n');

disp(bestTable);
