clc;clear;close all;
global evals
evals = 0;

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych eggholder funkcie
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 5;                                      % pocet behov
ngen = 3000;                                   % pocet generacii
islandSize = 100;                              % 3 ostrovy po 100 = 300
sgaSize = 300;                                 % rovnaky pocet jedincov ako cele PGA
stagnationLimit = 50;                          % po 50 generaciach bez zlepsenia zasah
pmut_global = 0.3;
pmut_local = 0.3;
Amp_local = 5 * ones(1,nvars);
resetCount = 50;                               % reset polovice ostrova A
warmupAmp = 0.30 * (Space(2,:) - Space(1,:));  % warm-up +-30 % rozsahu

% Rozdelenie populacie na jednom ostrove:
bestNumsPGA = [1,1];
RandOldCountPGA = 21;
work1CountPGA = 36; %tito idu na krizenie
work2CountPGA = 41; %tito idu na mutaciu

% Rozdelenie populacie SGA:
bestNumsSGA = ones(1,8);
RandOldCountSGA = 60;
work1CountSGA = 108;
work2CountSGA = 124;

histA = zeros(runs,ngen);
histB = zeros(runs,ngen);
histC = zeros(runs,ngen);
histSGA = zeros(runs,ngen);

migrationHist = cell(runs,1);

figure;
tiledlayout(3,2);

% 5 behov
for r = 1:runs

    migrationGenerations = [];

    popA = genrpop(islandSize,Space);
    popB = genrpop(islandSize,Space);
    popC = genrpop(islandSize,Space);
    popSGA = genrpop(sgaSize,Space);

    stagnationA = 0;
    stagnationB = 0;
    stagnationC = 0;

    prevBestA = inf;
    prevBestB = inf;
    prevBestC = inf;

    for g = 1:ngen

        % Vyhodnotenie
        fA = eggholder(popA);
        fB = eggholder(popB);
        fC = eggholder(popC);
        fSGA = eggholder(popSGA);

        bestAvalue = min(fA);
        bestBvalue = min(fB);
        bestCvalue = min(fC);

        histA(r,g) = bestAvalue;
        histB(r,g) = bestBvalue;
        histC(r,g) = bestCvalue;
        histSGA(r,g) = min(fSGA);

        % Kontrola stagnacie
        if bestAvalue < prevBestA - 1e-6
            stagnationA = 0;
            prevBestA = bestAvalue;
        else
            stagnationA = stagnationA + 1;
        end

        if bestBvalue < prevBestB - 1e-6
            stagnationB = 0;
            prevBestB = bestBvalue;
        else
            stagnationB = stagnationB + 1;
        end

        if bestCvalue < prevBestC - 1e-6
            stagnationC = 0;
            prevBestC = bestCvalue;
        else
            stagnationC = stagnationC + 1;
        end

        % Najlepsi jedinci pred pripadnym zasahom
        [~,bestAindex] = min(fA);
        [~,bestBindex] = min(fB);

        bestAindividual = popA(bestAindex,:);
        bestBindividual = popB(bestBindex,:);

        migrationNow = false;

        % A - RESET + migrant z B
        if stagnationA >= stagnationLimit

            [~,worstA] = sort(fA,'descend');

            popA(worstA(1:resetCount),:) = ...
                genrpop(resetCount,Space);      % reset najhorsej polovice

            popA(worstA(1),:) = bestBindividual; % B -> A

            fA = eggholder(popA);
            prevBestA = min(fA);
            stagnationA = 0;
            migrationNow = true;
        end

        % B - WARM-UP + migrant z A
        if stagnationB >= stagnationLimit

            noise = (2*rand(size(popB))-1) .* warmupAmp;
            popB = popB + noise;                % nahodny sum

            popB = max(popB,Space(1,:));        % ohranicenie na <-500,500>
            popB = min(popB,Space(2,:));

            popB(1,:) = bestBindividual;         % zachovanie najlepsieho B
            popB(2,:) = bestAindividual;         % A -> B

            fB = eggholder(popB);
            prevBestB = min(fB);
            stagnationB = 0;
            migrationNow = true;
        end

        % C - migranti z A aj B
        if stagnationC >= stagnationLimit

            [~,worstC] = sort(fC,'descend');

            popC(worstC(1),:) = bestAindividual; % A -> C
            popC(worstC(2),:) = bestBindividual; % B -> C

            fC = eggholder(popC);
            prevBestC = min(fC);
            stagnationC = 0;
            migrationNow = true;
        end

        % Zapamatanie generacie zasahu
        if migrationNow
            migrationGenerations(end+1) = g;
        end

        % OSTROV A
        Best = selbest(popA,fA,bestNumsPGA);
        RandOld = selrand(popA,fA,RandOldCountPGA);
        Work1 = selsus(popA,fA,work1CountPGA);
        Work2 = selsus(popA,fA,work2CountPGA);

        Work1 = crossov(Work1,1,0);
        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        popA = [Best;RandOld;Work1;Work2];

        % OSTROV B
        Best = selbest(popB,fB,bestNumsPGA);
        RandOld = selrand(popB,fB,RandOldCountPGA);
        Work1 = selsus(popB,fB,work1CountPGA);
        Work2 = selsus(popB,fB,work2CountPGA);

        Work1 = crossov(Work1,1,0);
        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        popB = [Best;RandOld;Work1;Work2];

        % OSTROV C
        Best = selbest(popC,fC,bestNumsPGA);
        RandOld = selrand(popC,fC,RandOldCountPGA);
        Work1 = selsus(popC,fC,work1CountPGA);
        Work2 = selsus(popC,fC,work2CountPGA);

        Work1 = crossov(Work1,1,0);
        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        popC = [Best;RandOld;Work1;Work2];

        % SGA
        Best = selbest(popSGA,fSGA,bestNumsSGA);
        RandOld = selrand(popSGA,fSGA,RandOldCountSGA);
        Work1 = selsus(popSGA,fSGA,work1CountSGA);
        Work2 = selsus(popSGA,fSGA,work2CountSGA);

        Work1 = crossov(Work1,1,0);
        Work2 = mutx(Work2,pmut_global,Space);
        Work2 = muta(Work2,pmut_local,Amp_local,Space);

        popSGA = [Best;RandOld;Work1;Work2];

    end

    migrationHist{r} = migrationGenerations;

    fprintf("Run %d: A = %.2f, B = %.2f, C = %.2f, SGA = %.2f\n", ...
        r,histA(r,end),histB(r,end),histC(r,end),histSGA(r,end));

    % Graf konkretneho behu
    nexttile;
    hold on;grid on;

    plot(1:ngen,histA(r,:),'LineWidth',1);
    plot(1:ngen,histB(r,:),'LineWidth',1);
    plot(1:ngen,histC(r,:),'LineWidth',1);
    plot(1:ngen,histSGA(r,:),'k','LineWidth',1.5);

    % Oznacenie migracie/resetu/warm-upu
    for m = 1:length(migrationHist{r})
        xline(migrationHist{r}(m),'--','HandleVisibility','off');
    end

    xlabel('Generacia');
    ylabel('Najlepsia hodnota');
    title(['Beh ',num2str(r)]);
    legend('Ostrov A','Ostrov B','Ostrov C','SGA');
    hold off;

end

% Priemery 5 behov
averageA = mean(histA,1);
averageB = mean(histB,1);
averageC = mean(histC,1);
averageSGA = mean(histSGA,1);

nexttile;
hold on;grid on;

plot(1:ngen,averageA,'LineWidth',1.5);
plot(1:ngen,averageB,'LineWidth',1.5);
plot(1:ngen,averageC,'LineWidth',1.5);
plot(1:ngen,averageSGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Priemerna najlepsia hodnota');
title('Priemer 5 behov');
legend('Ostrov A','Ostrov B','Ostrov C','SGA');
hold off;