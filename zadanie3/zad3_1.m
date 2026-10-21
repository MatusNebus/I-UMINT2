clc;clear;close all;

global evals
evals = 0;

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych Eggholder funkcie
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 5;                                      % pocet behov
seeds = [1 2 3 4 5];                           % rovnake seedy pre vsetky algoritmy
ngen = 2000;                                   % pocet generacii
pmut_global = 0.3;
pmut_local = 0.3;
Amp_local = 5 * ones(1,nvars);

% CPGA
nRows = 30;
nCols = 15;
cpgaSize = nRows*nCols;

% Mooreovo susedstvo - 8 susedov
neighbourOffsets = [
    -1 -1
    -1  0
    -1  1
     0 -1
     0  1
     1 -1
     1  0
     1  1
];

% PGA
numIslands = 9;
islandSize = 50;
migrationPeriod = 100;
bestNumsPGA = [1,1];
RandOldCountPGA = 10;
work1CountPGA = 18;
work2CountPGA = 20;

% SGA
sgaSize = 450;
bestNumsSGA = ones(1,18);
RandOldCountSGA = 90;
work1CountSGA = 162;
work2CountSGA = 180;

% Ulozenie vysledkov
histCPGA = zeros(runs,ngen);
histPGA = zeros(runs,ngen);
histSGA = zeros(runs,ngen);

finalCPGA = zeros(runs,1);
finalPGA = zeros(runs,1);
finalSGA = zeros(runs,1);

bestVectorsCPGA = zeros(runs,nvars);
bestVectorsPGA = zeros(runs,nvars);
bestVectorsSGA = zeros(runs,nvars);


% 5 behov
for r = 1:runs
    fprintf('\nBeh %d, seed %d\n',r,seeds(r));

    % CPGA
    rng(seeds(r));
    popCPGA = genrpop(cpgaSize,Space);
    fCPGA = eggholder(popCPGA);
    fCPGA = fCPGA(:);

    for g = 1:ngen

        histCPGA(r,g) = min(fCPGA);

        if g == ngen
            break;
        end

        partnerIndex = zeros(cpgaSize,1);
        % Pre kazdu bunku sa vyberie jeden nahodny sused
        for row = 1:nRows
            for col = 1:nCols
                cellIndex = sub2ind([nRows,nCols],row,col);
                neighbourNumber = randi(8);
                dr = neighbourOffsets(neighbourNumber,1);
                dc = neighbourOffsets(neighbourNumber,2);

                % Toroidna mriezka
                neighbourRow = mod(row - 1 + dr,nRows) + 1;
                neighbourCol = mod(col - 1 + dc,nCols) + 1;

                partnerIndex(cellIndex) = ...
                    sub2ind([nRows,nCols],neighbourRow,neighbourCol);
            end
        end

        % Vytvorenie dvojic rodicov
        parents = zeros(2*cpgaSize,nvars);
        parents(1:2:end,:) = popCPGA;
        parents(2:2:end,:) = popCPGA(partnerIndex,:);

        children = crossov(parents,1,0);    % Krizenie

        % Z kazdej dvojice pouzijeme jedneho potomka
        offspring = children(1:2:end,:);

        % Mutacia
        offspring = mutx(offspring,pmut_global,Space);
        offspring = muta(offspring,pmut_local,Amp_local,Space);

        fOffspring = eggholder(offspring);
        fOffspring = fOffspring(:);

        % Synchronna aktualizacia
        % Potomok nahradi bunku iba ak je lepsi
        better = fOffspring < fCPGA;

        popCPGA(better,:) = offspring(better,:);
        fCPGA(better) = fOffspring(better);

    end

    [finalCPGA(r),bestIndex] = min(fCPGA);
    bestVectorsCPGA(r,:) = popCPGA(bestIndex,:);

    fprintf('CPGA = %.6f\n',finalCPGA(r));


    % PGA
    rng(seeds(r));
    popPGA = cell(1,numIslands);

    for k = 1:numIslands
        popPGA{k} = genrpop(islandSize,Space);
    end

    fitPGA = cell(1,numIslands);
    for g = 1:ngen

        % Vyhodnotenie ostrovov
        for k = 1:numIslands
            fitPGA{k} = eggholder(popPGA{k});
            fitPGA{k} = fitPGA{k}(:);
        end

        % Periodicka migracia 2-9 -> 1
        if mod(g,migrationPeriod) == 0
            migrants = zeros(8,nvars);

            % Najlepsi jedinec z ostrovov 2 az 9
            for k = 2:numIslands

                [~,bestIndex] = min(fitPGA{k});
                migrants(k-1,:) = popPGA{k}(bestIndex,:);

            end

            % 8 najhorsich jedincov centralneho ostrova
            [~,worstCenter] = sort(fitPGA{1},'descend');

            % Migranti nahradia najhorsich
            popPGA{1}(worstCenter(1:8),:) = migrants;

            % Nove vyhodnotenie ostrova 1
            fitPGA{1} = eggholder(popPGA{1});
            fitPGA{1} = fitPGA{1}(:);

        end

        % Globalne najlepsia fitness zo vsetkych ostrovov
        globalBestPGA = inf;

        for k = 1:numIslands

            islandBest = min(fitPGA{k});

            if islandBest < globalBestPGA
                globalBestPGA = islandBest;
            end

        end

        histPGA(r,g) = globalBestPGA;

        if g == ngen
            break;
        end

        % GA na jednotlivych ostrovoch
        for k = 1:numIslands

            Best = selbest(popPGA{k},fitPGA{k},bestNumsPGA);
            RandOld = selrand(popPGA{k},fitPGA{k},RandOldCountPGA);

            Work1 = selsus(popPGA{k},fitPGA{k},work1CountPGA);
            Work2 = selsus(popPGA{k},fitPGA{k},work2CountPGA);

            Work1 = crossov(Work1,1,0);

            Work2 = mutx(Work2,pmut_global,Space);
            Work2 = muta(Work2,pmut_local,Amp_local,Space);

            popPGA{k} = [Best;RandOld;Work1;Work2];

        end

    end

    % Najlepsi jedinec zo vsetkych ostrovov
    finalPGA(r) = inf;

    for k = 1:numIslands

        [islandBest,bestIndex] = min(fitPGA{k});

        if islandBest < finalPGA(r)

            finalPGA(r) = islandBest;
            bestVectorsPGA(r,:) = popPGA{k}(bestIndex,:);

        end

    end

    fprintf('PGA  = %.6f\n',finalPGA(r));


    % SGA

    rng(seeds(r));

    popSGA = genrpop(sgaSize,Space);

    for g = 1:ngen

        fSGA = eggholder(popSGA);
        fSGA = fSGA(:);

        histSGA(r,g) = min(fSGA);

        if g == ngen
            break;
        end

        Best = selbest(popSGA,fSGA,bestNumsSGA);
        RandOld = selrand(popSGA,fSGA,RandOldCountSGA);

        Work1 = selsus(popSGA,fSGA,work1CountSGA);
        Work2 = selsus(popSGA,fSGA,work2CountSGA);

        Work1 = crossov(Work1,1,0);

        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        popSGA = [Best;RandOld;Work1;Work2];

    end

    [finalSGA(r),bestIndex] = min(fSGA);
    bestVectorsSGA(r,:) = popSGA(bestIndex,:);

    fprintf('SGA  = %.6f\n',finalSGA(r));

end


% Priemerne priebehy 5 behov
averageCPGA = mean(histCPGA,1);
averagePGA = mean(histPGA,1);
averageSGA = mean(histSGA,1);


% Spolocna mierka grafov
allFitness = [
    histCPGA(:);
    histPGA(:);
    histSGA(:)
q];

commonYLim = [min(allFitness) max(allFitness)];

% Graf 1 - CPGA

fig1 = figure;
hold on;grid on;
a
for r = 1:runs
    plot(1:ngen,histCPGA(r,:),'LineWidth',1);
end

plot(1:ngen,averageCPGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('CPGA: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    'Seed 1','Seed 2','Seed 3','Seed 4','Seed 5','Priemer', ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 2 - PGA

fig2 = figure;
hold on;grid on;

for r = 1:runs
    plot(1:ngen,histPGA(r,:),'LineWidth',1);
end

plot(1:ngen,averagePGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('PGA: priebehy globalne najlepsej fitness z 5 behov a ich priemer');

legend( ...
    'Seed 1','Seed 2','Seed 3','Seed 4','Seed 5','Priemer', ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 3 - SGA

fig3 = figure;
hold on;grid on;

for r = 1:runs
    plot(1:ngen,histSGA(r,:),'LineWidth',1);
end

plot(1:ngen,averageSGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('SGA: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    'Seed 1','Seed 2','Seed 3','Seed 4','Seed 5','Priemer', ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 1 - CPGA

fig1 = figure;
hold on;grid on;

for r = 1:runs
    plot(1:ngen,histCPGA(r,:),'LineWidth',1);
end

plot(1:ngen,averageCPGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('CPGA: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    sprintf('Seed %d: %.2f',seeds(1),histCPGA(1,end)), ...
    sprintf('Seed %d: %.2f',seeds(2),histCPGA(2,end)), ...
    sprintf('Seed %d: %.2f',seeds(3),histCPGA(3,end)), ...
    sprintf('Seed %d: %.2f',seeds(4),histCPGA(4,end)), ...
    sprintf('Seed %d: %.2f',seeds(5),histCPGA(5,end)), ...
    sprintf('Priemer: %.2f',averageCPGA(end)), ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 2 - PGA

fig2 = figure;
hold on;grid on;

for r = 1:runs
    plot(1:ngen,histPGA(r,:),'LineWidth',1);
end

plot(1:ngen,averagePGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('PGA: priebehy globalne najlepsej fitness z 5 behov a ich priemer');

legend( ...
    sprintf('Seed %d: %.2f',seeds(1),histPGA(1,end)), ...
    sprintf('Seed %d: %.2f',seeds(2),histPGA(2,end)), ...
    sprintf('Seed %d: %.2f',seeds(3),histPGA(3,end)), ...
    sprintf('Seed %d: %.2f',seeds(4),histPGA(4,end)), ...
    sprintf('Seed %d: %.2f',seeds(5),histPGA(5,end)), ...
    sprintf('Priemer: %.2f',averagePGA(end)), ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 3 - SGA

fig3 = figure;
hold on;grid on;

for r = 1:runs
    plot(1:ngen,histSGA(r,:),'LineWidth',1);
end

plot(1:ngen,averageSGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('SGA: priebehy najlepsiej fitness z 5 behov a ich priemer');

legend( ...
    sprintf('Seed %d: %.2f',seeds(1),histSGA(1,end)), ...
    sprintf('Seed %d: %.2f',seeds(2),histSGA(2,end)), ...
    sprintf('Seed %d: %.2f',seeds(3),histSGA(3,end)), ...
    sprintf('Seed %d: %.2f',seeds(4),histSGA(4,end)), ...
    sprintf('Seed %d: %.2f',seeds(5),histSGA(5,end)), ...
    sprintf('Priemer: %.2f',averageSGA(end)), ...
    'Location','best');

ylim(commonYLim);

hold off;


% Graf 4 - porovnanie priemerov

fig4 = figure;
hold on;grid on;

plot(1:ngen,averageCPGA,'LineWidth',2);
plot(1:ngen,averagePGA,'LineWidth',2);
plot(1:ngen,averageSGA,'LineWidth',2);

xlabel('Generacia');
ylabel('Fitness');
title('Porovnanie priemernych priebehov CPGA, PGA a SGA');

legend( ...
    sprintf('CPGA: %.2f',averageCPGA(end)), ...
    sprintf('PGA: %.2f',averagePGA(end)), ...
    sprintf('SGA: %.2f',averageSGA(end)), ...
    'Location','best');

ylim(commonYLim);

hold off;


% Export grafov
exportgraphics(fig1,'CPGA.png','Resolution',300);
exportgraphics(fig2,'PGA.png','Resolution',300);
exportgraphics(fig3,'SGA.png','Resolution',300);
exportgraphics(fig4,'Porovnanie.png','Resolution',300);


% Tabulka vysledkov 5 behov + priemer + smerodajna odchylka

Beh = [
    compose("Seed %d",seeds(:));
    "Priemer";
    "Smerodajna odchylka"
];

CPGA = [
    finalCPGA;
    mean(finalCPGA);
    std(finalCPGA)
];

PGA = [
    finalPGA;
    mean(finalPGA);
    std(finalPGA)
];

SGA = [
    finalSGA;
    mean(finalSGA);
    std(finalSGA)
];

resultsTable = table(Beh,CPGA,PGA,SGA);

fprintf('\n');
fprintf('============================================\n');
fprintf('TABULKA VYSLEDKOV\n');
fprintf('============================================\n');

disp(resultsTable);


% Najlepsi vektor kazdeho algoritmu

[bestFitnessCPGA,bestRunCPGA] = min(finalCPGA);
[bestFitnessPGA,bestRunPGA] = min(finalPGA);
[bestFitnessSGA,bestRunSGA] = min(finalSGA);

bestVectors = [
    bestVectorsCPGA(bestRunCPGA,:);
    bestVectorsPGA(bestRunPGA,:);
    bestVectorsSGA(bestRunSGA,:)
];

Algoritmus = ["CPGA";"PGA";"SGA"];

Seed = [
    seeds(bestRunCPGA);
    seeds(bestRunPGA);
    seeds(bestRunSGA)
];

Fitness = [
    bestFitnessCPGA;
    bestFitnessPGA;
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


% Tabulka nastaveni

Nastavenie = [
    "Pocet premennych";
    "Rozsah premennych";
    "Pocet generacii";
    "Pocet behov";
    "Seedy";
    "Pocet jedincov";
    "CPGA mriezka";
    "CPGA susedstvo";
    "CPGA hranice";
    "CPGA vyber suseda";
    "CPGA aktualizacia";
    "PGA pocet ostrovov";
    "PGA velkost ostrova";
    "PGA perioda migracie";
    "PGA migracia";
    "SGA velkost populacie";
    "Globalna mutacia";
    "Lokalna mutacia";
    "Amplituda lokalnej mutacie"
];

Hodnota = [
    string(nvars);
    "[-500, 500]";
    string(ngen);
    string(runs);
    mat2str(seeds);
    string(cpgaSize);
    string(nRows) + " x " + string(nCols);
    "Mooreovo - 8 susedov";
    "Toroidne";
    "Nahodny z 8 susedov";
    "Synchronna, lepsi potomok nahradi bunku";
    string(numIslands);
    string(islandSize);
    string(migrationPeriod);
    "Najlepsi z ostrovov 2-9 do ostrova 1";
    string(sgaSize);
    string(pmut_global);
    string(pmut_local);
    string(Amp_local(1))
];

settingsTable = table(Nastavenie,Hodnota);

fprintf('\n');
fprintf('============================================\n');
fprintf('NASTAVENIA\n');
fprintf('============================================\n');

disp(settingsTable);