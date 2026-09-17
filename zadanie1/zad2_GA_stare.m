clc; clear; close all;

%nastavenia, parametre:
nvars = 100;                         % pocet premennych
Space = [-500*ones(1,nvars);500*ones(1,nvars)]; % definicny obor pre eggholder
%Space = [-1000*ones(1,nvars);1000*ones(1,nvars)]; % definicny obor pre schwefelovu
%target_success = -7920;             % hranica uspechu pre schwefelovu 10D
%target_success = -80000; %-79200;           % hranica uspechu pre schwefelovu 100D
target_success = -100000;           % hranica uspechu pre eggholder

runs    = 5;                        % počet behov
popsize = 1000;                      % veľkosť populácie
ngen    = 1000;                       % počet generácií

eliteCount = 100;                     % b = elitizmus (koľko najlepších ide priamo ďalej)

selectionName = "selsus";           % metoda vyberu: selsus, selbest, selrand, seltourn
crossoverName = "crossov";          % metoda krizenia: crossov, among (among asi nemame ale)
pcross        = 0.85;               % pravdepodobnosť kríženia v generácii
nCrossPoints  = 2;                  % len pre crossov: počet bodov (1..4)

pmut_global = 0.15;                 % pravdepodobnost mutacie pre mutx
pmut_local  = 0.2;                  % pravdepodobnost mutacie pre muta
Amp_local = 500 * ones(1, nvars);    % velkost zmeny (+-) pre muta (teraz 50 pre kazdy parameter)
globalPhaseRatio = 0.6;             % prvých 60% generácií mutx, zvyšok muta

%Úložiská na výsledky (kvoli grafu a porovnaniu behov)
bestHistAll = nan(runs, ngen);      % každý riadok = jeden beh, stĺpce = generácie
bestXAll    = nan(runs, nvars);     % najlepší jedinec z každého behu
bestFAll    = nan(runs, 1);         % najlepšia fitness z každého behu


%HLAVNY CYKLUS:
for r = 1:runs
    pop = genrpop(popsize, Space);  % Inicializacia populácie
    f = eval_pop(pop);              % vyhodnotenie fitness pre celú populáciu

    for g = 1:ngen
        elites = selbest(pop, f, eliteCount); % eliteCount najlepších ide priamo ďalej

        nKids = popsize - eliteCount;                           % pocet deti n-b
        
        %selekcia
        parents = do_selection(selectionName, pop, f, nKids);

        %krizenie
        kids = parents;             % ak by krizenie nenastalo
        if rand < pcross            % nastane s pravdepodobnostou pcross
            kids = do_crossover(crossoverName, parents, nCrossPoints);
        end

        %oba typy mutacie naraz:
        kidsPoMutx = mutx(kids, pmut_global, Space); 
        kids = muta(kidsPoMutx, pmut_local, Amp_local, Space);

        %iba globalna:
        %kids = mutx(kids, pmut_global, Space);

        %iba lokalna s apmplitudou 50:
        %kids = muta(kids, pmut_local, Amp_local, Space);


        %mutacia najskor globalna, potom lokalna:
        % if g <= round(globalPhaseRatio * ngen)                  %najskor globalna, potom lokalna
        %     kids = mutx(kids, pmut_global, Space);              % globalna
        % else
        %     kids = muta(kids, pmut_local, Amp_local, Space);    % lokalna
        % end

        pop = [elites; kids];       %nova populacia
        f = eval_pop(pop);          %vyhodnotenie novej populacie

        bestF = min(f);             %ulozenie najlepsich hodnot
        bestHistAll(r, g) = bestF;

        if bestF <= target_success             % ak dosiahneme target, ukoncime
            bestHistAll(r, g:end) = bestF;     % nech graf pekne dobehne rovno
            break;
        end
    end

    [bestF, bestIdx] = min(f);      %najlepsi jedinec z konca behu, z poslednej populacie
    bestFAll(r) = bestF;
    bestXAll(r,:) = pop(bestIdx,:);
    fprintf("Run %02d: bestF = %.6f\n", r, bestF);
end

% vyhodnotenie, graf:
[globalBestF, ridx] = min(bestFAll);
globalBestX = bestXAll(ridx,:);

fprintf("\n====================\n");
fprintf("GLOBAL BEST (z %d behov):\n", runs);
fprintf("bestF = %.6f\n", globalBestF);
fprintf("bestX = ["); fprintf(" %.4f", globalBestX); fprintf(" ]\n");
fprintf("====================\n");

figure; hold on; grid on;
plot(1:ngen, bestHistAll', 'LineWidth', 1);
xlabel("Generácia");
ylabel("Najlepšia fitness v populácii");
title("Genetický algoritmus");
yline(target_success, "--", sprintf("hranica úspechu %d", target_success));
hold off;

% ==========pomocne funkcie=============
% vyhodnotenie populácie:
function f = eval_pop(pop)
    n = size(pop, 1);
    f = zeros(n, 1);
    for i = 1:n
        %f(i) = testfn3c(pop(i, :)); %schwefelova funkcia
        f(i) = eggholder(pop(i, :)); %eggholder
    end
end

% vyber rodičov podľa zvoleného typu:
function parents = do_selection(name, pop, f, nSel)
    switch lower(name)
        case "selsus"
            parents = selsus(pop, f, nSel);
        case "selbest"
            parents = selbest(pop, f, nSel);
        case "selrand"
            parents = selrand(pop, f, nSel);
        case "seltourn"
            parents = seltourn(pop, f, nSel);
        otherwise
            error("Neznáma selekcia: %s", name);
    end
end

% kríženie podľa zvoleného
% typu 
function kids = do_crossover(name, parents, nCrossPoints)
    switch lower(name)
        case "crossov"
            kids = crossov(parents, nCrossPoints, 0);
        case "among"
            kids = among(parents, 0);
        otherwise
            error("Neznáme kríženie: %s", name);
    end
end