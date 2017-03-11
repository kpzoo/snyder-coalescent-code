% Script to test heterochronous performance
clearvars
clc
close all

% Tree no. tips fixed for all tests
n = 200;
% No. repetitions for averaging
M = 100;
% Function of interest and its grid size
fnid = 5;
numRV = 4;
mi = 15*ones(1, numRV);

% Sample times to run against
nSampTimes = [1 7 15 30 50];
nRuns = length(nSampTimes);

% Creat a label for legend
hetLab = cell(1, nRuns);
for i = 1:nRuns
    hetLab{i} = num2str(nSampTimes(i));
end

% Variables to store results
x = cell(1, nRuns);
xhat = cell(1, nRuns);
Nhatmean = cell(1, nRuns);
Nhatlb = cell(1, nRuns);
Nhatub = cell(1, nRuns);
nLin = cell(1, nRuns);
tLin = cell(1, nRuns);
sampSny = cell(1, nRuns);
tN = cell(1, nRuns);
Nt = cell(1, nRuns);

% Run in a loop
for i = 1:nRuns
    % Main code runs M sims of n tip trees over grid mi for nSampTimes
    [x{i}, xhat{i}, Nhatmean{i}, Nhatlb{i}, Nhatub{i}, nLin{i}, tLin{i},...
        sampSny{i}, tN{i}, Nt{i}] = runSnyCoalHetFn(M, n, fnid, mi, nSampTimes(i));
    disp(['Finished run ' num2str(i) ' of ' num2str(nRuns)]);
    disp('***************************************************************');
end

% Save data and clear other variables
save(['testHet_' num2str(fnid) '.mat']);

% Ensure all true values are the same as use x{1} to compare in plots
dx = zeros(1, nRuns-1);
for i = 2:nRuns
    dx(i) = any(x{i} - x{1} ~= 0);
end
if any(dx)
    error('The true values are not the same across runs');
end

% Boxplot of the relative errors
figure;
% For each parameter
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    % Error storage variable
    ehat = zeros(M, nRuns);
    % For each sampling scheme
    hold on
    for j = 2:nRuns
        % Error for a given parameter for a given run
        ehat(:, j) = (xhat{j}(i, :)' - x{j}(i))/x{j}(i);
        [~, ~, bw] = ksdensity(ehat(:, j));
        newBw = 3*bw;
        [fsny, xsny] = ksdensity(ehat(:, j), 'width', newBw);
        plot(xsny, fsny, 'linewidth', 2);
    end
    %boxplot(ehat, 'labels', nSampTimes);
    legend(hetLab{2:end}, 'location', 'best');
    box on
    xlabel('no. sample times');
    ylabel('relative error');
    title(['Parameter: x_' num2str(i)]);
    grid;
end

% Smoothed histograms of each 
figure;
postLab = [hetLab {'true'}];
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    hold on
    for j = 1:nRuns
        % Main density of each parameter for both methods
        [~, ~, bw] = ksdensity(sampSny{j}(:, i));
        newBw = 4*bw;
        [fsny, xsny] = ksdensity(sampSny{j}(:, i), 'width', newBw);
        plot(xsny, fsny, 'linewidth', 2);
    end
    % Get current axes and plot true value if not empirical data
    h = gca;
    plot(x{1}(i)*ones(1, 2), h.YLim, 'k', 'linewidth', 2);
    hold off
    grid;
    xlabel('parameter space');
    ylabel('smoothed pdf');
    hold off
    title(['Posterior: x_' num2str(i) ' = ' num2str(x{1}(i))]);
    legend(postLab, 'location', 'best');
    box on
end

% Recalculate sample times
if fnid == 5
    maxSampTimes = 2*x{1}(3) + x{1}(4);
end
svec = cell(1, nRuns);
for i = 1:nRuns
    if nSampTimes(i) > 1
        svec{i} = linspace(0, maxSampTimes, nSampTimes(i));
    else
        % Because linspace picks maxSampTime if nSampTimes = 1
        svec{i} = 0;
    end
end

% Combined trajectories without homochronous and with linked axes
figure;
for j = 1:nRuns-1
    i = j+1;
    subplot(ceil((nRuns-1)/2), 2, j);
    plot(tN{i}, Nt{i}, 'k', 'linewidth', 2);
    hold on
    plot(tN{i}, Nhatmean{i}, 'r--', 'linewidth', 2);
    plot(tN{i}, Nhatlb{i},'g:', tN{i}, Nhatub{i}, 'g:', 'linewidth', 2);
    hAx = gca;
    % Add sample times
    %h = stem(svec{i}, hAx.YLim(2)*ones(size(svec{i})), 'g:', 'filled');
    h = stem(svec{i}, (4.8*10^4)*ones(size(svec{i})), 'g:', 'filled');
    set(h, 'color', 0.5*ones(1, 3));
    hold off
    xlim([tN{i}(1), tN{i}(end)]);
    legend('true trajectory', 'mean', '2.5% percentile', '97.5% percentile', 'sample times', 'location', 'best');
    xlabel('time from present');
    ylabel('N(t) estimate');
    if fnid == 5
        title(['Con-exp-con: [n, M, K] = [' [num2str(n) ' ' num2str(M) ' ' num2str(nSampTimes(i))] ']']);
    end
end
% Link the axes (in case of 4)
if nRuns-1 == 4
    h1 = subplot(2, 2, 1);
    h2 = subplot(2, 2, 2);
    h3 = subplot(2, 2, 3);
    h4 = subplot(2, 2, 4);
    linkaxes([h1 h2 h3 h4], 'y');
end
        
        
        
        
