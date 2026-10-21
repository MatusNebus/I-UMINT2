
clc;clear;close all;

% Uvodne nastavenia
nvars = 10;                     % pocet premennych
Space = [-500*ones(1,nvars);
          500*ones(1,nvars)];   % definicny obor
ngen = 2000;                    % pocet generacii
popSize = 50;                   % pocet jedincov
seed = 1;

% Testovane parametre DE
Fvalues = [0.25 0.5 0.75 1 1.25];
CRvalues = [0.01 0.05 0.1 0.25 0.5 0.75];

% Spolocna pociatocna populacia
rng(seed);
initialPop = genrpop(popSize,Space);
initialRngState = rng;

% Ulozenie vysledkov
% Prvy prvok historie je inicializacia, dalsie su generacie

histDE = zeros(5,6,ngen+1);
finalFitness = zeros(5,6);

% Testovanie 25 kombinacii

for a = 1:5
    for b = 1:6

        F = Fvalues(a);
        CR = CRvalues(b);

        rng(initialRngState);
        popDE = initialPop;
        fDE = testfn3(popDE);
        fDE = fDE(:);

        histDE(a,b,1) = min(fDE);

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

                % Diferencialna mutacia
                v = popDE(r1,:) + F*(popDE(r2,:) - popDE(r3,:));

                % Binomicke krizenie
                mask = rand(1,nvars) < CR;

                % Aspoň jeden gen musi pochadzat z mutanta
                jrand = randi(nvars);
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

            % Najlepsia doteraz dosiahnuta fitness
            histDE(a,b,g+1) = min(fDE);

        end

        finalFitness(a,b) = min(fDE);

        fprintf('F = %.2f, CR = %.2f, Fitness = %.6f\n', ...
            F,CR,finalFitness(a,b));

    end
end

% Graf 1 - priebehy konvergencie pre rozne F a CR

fig1 = figure;
tiledlayout(2,3);

iterations = 0:ngen;
commonYLim = [min(histDE(:)) max(histDE(:))];

colors = lines(5);

for b = 1:6

    nexttile;
    hold on;grid on;

    for a = 1:5

        plot(iterations,squeeze(histDE(a,b,:)), ...
            'Color',colors(a,:),'LineWidth',1.5, ...
            'DisplayName',sprintf('F = %.2f',Fvalues(a)));

    end

    xlabel('Generacia');
    ylabel('Fitness');
    title(sprintf('CR = %.2f',CRvalues(b)));
    xlim([0 ngen]);
    ylim(commonYLim);
    legend('Location','best');

    hold off;
end

sgtitle('DE - porovnanie parametrov F a CR');

% Tabulka vysledkov

[Fgrid,CRgrid] = ndgrid(Fvalues,CRvalues);

resultsTable = table( ...
    Fgrid(:),CRgrid(:),finalFitness(:), ...
    'VariableNames',{'F','CR','Fitness'});

resultsTable = sortrows(resultsTable,'Fitness');

fprintf('\n');
fprintf('============================================\n');
fprintf('VYSLEDKY 30 KOMBINACII\n');
fprintf('============================================\n');

disp(resultsTable);

% Najlepsia kombinacia parametrov

fprintf('\n');
fprintf('============================================\n');
fprintf('NAJLEPSIE PARAMETRE DE\n');
fprintf('============================================\n');

fprintf('F = %.2f\n',resultsTable.F(1));
fprintf('CR = %.2f\n',resultsTable.CR(1));
fprintf('Fitness = %.6f\n',resultsTable.Fitness(1));
