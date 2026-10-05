clc;clear;close all;
global evals
evals = 0;

% Uvodne nastavenia
nvars = 10;                                    % pocet premennych eggholder funkcie
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];                  % definicny obor
runs = 5;                                      % pocet behov
ngen = 3000;                                   % pocet generacii
numIslands = 9;                                % 9 ostrovov
islandSize = 50;                               % 50 jedincov na ostrove
sgaSize = numIslands * islandSize;             % SGA = 450 jedincov
migrationPeriod = 100;                         % migracia kazdych 100 generacii
pmut_global = 0.3;
pmut_local = 0.3;
Amp_local = 5 * ones(1,nvars);

% Rozdelenie populacie na ostrove
bestNumsPGA = [1,1];
RandOldCountPGA = 10;
work1CountPGA = 18;
work2CountPGA = 20;

% Rozdelenie populacie SGA
bestNumsSGA = ones(1,18);
RandOldCountSGA = 90;
work1CountSGA = 162;
work2CountSGA = 180;

histPGA = zeros(runs,numIslands,ngen);
histSGA = zeros(runs,ngen);

migrationHist = cell(runs,1);

figure;
tiledlayout(3,2);

% 5 behov
for r = 1:runs

    migrationGenerations = [];

    % Pociatocne populacie 9 ostrovov
    pop = cell(1,numIslands);

    for k = 1:numIslands
        pop{k} = genrpop(islandSize,Space);
    end

    popSGA = genrpop(sgaSize,Space);

    for g = 1:ngen

        % Vyhodnotenie ostrovov
        fit = cell(1,numIslands);

        for k = 1:numIslands
            fit{k} = eggholder(pop{k});
            histPGA(r,k,g) = min(fit{k});
        end

        % Vyhodnotenie SGA
        fSGA = eggholder(popSGA);
        histSGA(r,g) = min(fSGA);

        % Periodicka migracia 2-9 -> 1
        if mod(g,migrationPeriod) == 0

            migrants = zeros(8,nvars);

            % Najlepsi jedinci z ostrovov 2 az 9
            for k = 2:9
                [~,bestIndex] = min(fit{k});
                migrants(k-1,:) = pop{k}(bestIndex,:);
            end

            % 8 najhorsich jedincov centralneho ostrova
            [~,worstCenter] = sort(fit{1},'descend');

            % Migranti nahradia najhorsich v ostrove 1
            pop{1}(worstCenter(1:8),:) = migrants;

            % Nove vyhodnotenie centra po migracii
            fit{1} = eggholder(pop{1});

            migrationGenerations(end+1) = g;
        end

        % Geneticky algoritmus na vsetkych 9 ostrovoch
        for k = 1:numIslands

            Best = selbest(pop{k},fit{k},bestNumsPGA);
            RandOld = selrand(pop{k},fit{k},RandOldCountPGA);
            Work1 = selsus(pop{k},fit{k},work1CountPGA);
            Work2 = selsus(pop{k},fit{k},work2CountPGA);

            Work1 = crossov(Work1,1,0);
            Work2 = mutx(Work2,pmut_global,Space);
            Work2 = muta(Work2,pmut_local,Amp_local,Space);

            pop{k} = [Best;RandOld;Work1;Work2];

        end

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

    fprintf("Run %d: Ostrov 1 = %.2f, SGA = %.2f\n", ...
        r,histPGA(r,1,end),histSGA(r,end));

    % Graf konkretneho behu
    nexttile;
    hold on;grid on;

    plot(1:ngen,squeeze(histPGA(r,1,:)),'LineWidth',1.5);
    plot(1:ngen,squeeze(histPGA(r,2,:)),'LineWidth',1);
    plot(1:ngen,squeeze(histPGA(r,5,:)),'LineWidth',1);
    plot(1:ngen,squeeze(histPGA(r,9,:)),'LineWidth',1);
    plot(1:ngen,histSGA(r,:),'k','LineWidth',1.5);

    % Oznacenie periodickych migracii
    for m = 1:length(migrationHist{r})
        xline(migrationHist{r}(m),'--','HandleVisibility','off');
    end

    xlabel('Generacia');
    ylabel('Najlepsia hodnota');
    title(['Beh ',num2str(r)]);
    legend('Ostrov 1','Ostrov 2','Ostrov 5','Ostrov 9','SGA');
    hold off;

end

% Priemery 5 behov
averagePGA = squeeze(mean(histPGA,1));
averageSGA = mean(histSGA,1);

nexttile;
hold on;grid on;

plot(1:ngen,averagePGA(1,:),'LineWidth',1.5);
plot(1:ngen,averagePGA(2,:),'LineWidth',1);
plot(1:ngen,averagePGA(5,:),'LineWidth',1);
plot(1:ngen,averagePGA(9,:),'LineWidth',1);
plot(1:ngen,averageSGA,'k','LineWidth',2);

xlabel('Generacia');
ylabel('Priemerna najlepsia hodnota');
title('Priemer 5 behov');
legend('Ostrov 1','Ostrov 2','Ostrov 5','Ostrov 9','SGA');
hold off;