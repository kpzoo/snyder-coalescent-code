%% Perform Snyder filtering for multi-variable N(t) - main script

% Usage:
% - write parametric function as a new fnid
% - simulates a coalescent and then infers its parameters
% - multiple runs uses distribution of estimates from runs
% - single run uses the posterior from that run

% Assumptions and Modifications
% - x(3) is start of exp phase and x(3) + x(4) the end
% - only coded for con-exp-con and exponential at moment
% - allows for heterochronous sampling
% - requires a parametric form for the demographic function
% - checked numerical effects for multimodal exponential
% - used rejection sampling to simulate NHPP coalescent
% - the data is k-coalescent times and the process is self-exciting

% Clean the workspace, start timer
clc
close all
clearvars
tic;

%% Set simulation variables

% Perform a batch of M runs with trees of n tips
M = 1;
disp(['Simulating ' num2str(M) ' trees']);
n = 200;
% Get function and grid limits
fnid = 5;
[xmin, xmax, fnname, fn1, fn2, numRV] = setDemogFn2(fnid, n);

% Function identifier structure
fn.id = fnid;
fn.name = fnname;
fn.fn1 = fn1;
fn.fn2 = fn2;

% Variable to store execution times
tIter = zeros(1, M);
% Decide if true value in median of param space
fn.midpt = 1;
% Choose if will do linear or non-linear Snyder ODES
fn.linBool = 1;

% Set grid dimensions (mi) and filter complexity m
mi = 15*ones(1, numRV);
m = prod(mi);

% Check for pre-specified parameters not to be estimated - mi values of 1
specParam = find(mi == 1);
nEstRV = numRV - length(specParam);
% Modify spaces if specify parameters so that they are set to min values
if nEstRV > 0
    xmin(specParam) = xmax(specParam);
end

% Add more fields to function struc
fn.mi = mi;
fn.m = m;
fn.numRV = numRV;
fn.nEstRV = nEstRV;

% Create a matrix of identifiers to tell which xset{i} values are used for
% each entry in N(t) and lam(t) calculations
[xset, m, xsetMx, IDMx] = getxsetMxCoal(numRV, xmin, xmax, mi);
fn.xset = xset;
fn.IDMx = IDMx;
fn.xsetMx = xsetMx;

% Decompose xsetMx into column and row cell arrays for getCoalRate
fn.gridSz = size(xsetMx);
xsetRow = mat2cell(xsetMx, ones(1, fn.gridSz(1)), fn.gridSz(2));

% Pre-calculate the grid with time as the only input
fn.fnTime = @(tx) fn.fn2(xsetRow, tx);

% Choose true value from space    
if ~fn.midpt
    % Random sample from grid
    x = cellfun(@datasample, xset, num2cell(ones(1, numRV)));
else
    % Median of grid
    x = cellfun(@median, xset);
end
fn.param = x;

% Set uniform prior q0
q0 = ones(1, m)/m;


% Set number of sample times, assume uniformly divided in time
nSampTimes = 2;
disp(['Sampling at ' num2str(nSampTimes) ' different times']);
if fnid == 5
    % Max sample time to be symmetric on hetero
    maxSamp = 2*x(3) + x(4);
    pyb = 0; % boolean to control simulation of coalescent
else
    maxSamp = 100;
    pyb = 1;
end
if nSampTimes > 1
    svec = linspace(0, maxSamp, nSampTimes);
else
    % Because linspace picks maxSampTime if nSampTimes = 1
    svec = 0;
end

% Samples introduced at each sample time and n re-calculation
nvec = round(n/nSampTimes)*ones(1, nSampTimes);
n = sum(nvec);

% Set no. data points as distinct coalescent events and sample times -
% becomes n-1 for isochronous when nSampTimes = 1
nData = n + nSampTimes - 2;

% Get coalescent binomial factors for lineages and make no. lineages an
% index so facSort(1) = 0, facSort(k) = k choose 2
nset = n:-1:2;
fac = nset.*(nset-1)/2;
facSort = sort([0 fac]);
fn.facSort = facSort;
fn.nData = nData;
fn.svec = svec;

% Check as only fnid = 1, 5 allowed at the
if ~ismember(fnid, [1 5])
    error('Demographic function not specified for input fnid')
end

%% Main loop doing simulation and filtering on same tree

% Store parameter estimates
xhat = zeros(numRV, M);

% Times for evaluating population of each tree
nPts = 5000;
tset = zeros(nPts, M);
Nhat = zeros(nPts, M);
NhatSq = zeros(nPts, M);

% Function to calculate Nhat for a given time and posterior
NhatFn = @(q,tx)q*fn.fnTime(tx)'; 
NhatFnSq = @(q,tx)q*((fn.fnTime(tx)').^2); 

% Coalescent times, lineage times (includes coalescents and sample times)
tcoal = cell(1, M);
tLin = cell(1, M);
nLin = cell(1, M);

% Variables for batch results
qnlast = cell(1, M);
qmarg = cell(1, M);
xhatEv = cell(1, M);

for ii = 1:M
    % Generate the coalescent data using rejection sampling
    [nLin{ii}, tLin{ii}, tcoal{ii}] = simExpPyb(svec, nvec, n, nData, x, pyb, facSort, nSampTimes);
    
    % Filter the simulated event times - main filter function
    % Input is lineage times not just coalescent times
    [qev, ~] = snyderCoalHet(fn, q0, tLin{ii}, nLin{ii});
    
    % Last posterior and marginals and times
    qnlast{ii} = qev(end, :);
    [qmarg{ii}, ~] = marginalise(numRV, IDMx, qnlast{ii}, mi);
    
    % Conditional estimates, xhatEv gives estimates at event times
    xhatEv{ii} = qev*xsetMx';
    
    % Time points at which population size will be estimated
    tset(:, ii) = linspace(tcoal{ii}(1), tcoal{ii}(end), nPts);

    % Store batch results with xhat having rows for numRVs
    xhat(:, ii) = xhatEv{ii}(end, :);
    disp('***************************************************************');
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('***************************************************************');
    tIter(ii) = toc; % set to 0 only for parfor 
end

% True population, min(max) as smallest tree length wanted
tN = linspace(min(min(tset)), min(max(tset)), nPts);
Nt = fn.fn1(x, tN);

% Get all population estimates on same time points
for ii = 1:M
    % Collapse Nhat function with correct posterior
    Nq = @(tx) NhatFn(qnlast{ii}, tx);
    % Calculate N(t) posterior estimate across time set
    Nhat(:, ii) = arrayfun(Nq, tN);
end

% Single run uses std deviations
if M == 1
    Nq2 = @(tx) NhatFnSq(qnlast{1}, tx);
    NhatSq(:, 1) = arrayfun(Nq2, tN);
    Nstd = sqrt(NhatSq - Nhat.^2);
    % Upper and lower bounds with lower not < 0
    Nlb = Nhat - 2*Nstd;
    Nlb(Nlb < 0) = 0;
    Nub = Nhat + 2*Nstd;
end

% Statistics of demographics
Nhatmean = mean(Nhat, 2);
Nhatlb = quantile(Nhat', 0.025);
Nhatub = quantile(Nhat', 1-0.025);


% Clock time set to minutes
tIter = tIter/60;
disp(['Total time for ' num2str(M) ' iterations is ' num2str(tIter(end))]);
tavgIter = mean(diff([0, tIter]));
disp(['Average time with [n m M] = ' [num2str(n) ' ' num2str(m) ' ' num2str(M)] ' is ' num2str(tavgIter)]);

%% Analysing and plotting of data

% Variables for samples from posteriors
if M == 1
    % If want an estimate based on a single run
    nSamps = 10000;
else
    % Less samples per run as more data
    nSamps = 1000;
end
sampSny = zeros(M*nSamps, numRV);

% Sample from the discrete posteriors
for ii = 1:M
    % Get start and end indices in sample array
    id1 = (ii-1)*nSamps + 1;
    id2 = id1 + nSamps - 1;
    for i = 1:numRV
        % Frequentist approximation to probabilities
        sampSny(id1:id2, i) = datasample(xset{i}, nSamps, 'Weights', qmarg{ii}{i});
    end
end

% Get estimate relative errors
xmat = repmat(x, [M 1])';
ehat = 1 - xhat./xmat;
    
% Plot distribution of conditional mean estimates
if M > 1
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        % Main density of each parameter for both methods
        [fsny, xsny] = ksdensity(xhat(i, :));
        plot(xsny, fsny, 'linewidth', 2);
        hold on
        % Get current axes and plot true value if not empirical data
        h = gca;
        plot(x(i)*ones(1, 2), h.YLim, 'k', 'linewidth', 2);
        hold off
        grid;
        xlabel('parameter estimates');
        ylabel('pdf of conditional means');
        hold off
        title(['Estimates: x_' num2str(i) ' = ' num2str(x(i))]);
        legend('snyder', 'true', 'location', 'best');
    end
    
    % Distribution of error
    figure;
    for i = 1:numRV
        subplot(ceil(numRV/2), 2, i);
        % Main density of each parameter for both methods
        [fsny, xsny] = ksdensity(ehat(i, :));
        plot(xsny, fsny, 'linewidth', 2);
        grid;
        xlabel('parameter estimate error');
        ylabel('pdf of errors');
        hold off
        title(['Errors around: x_' num2str(i) ' = ' num2str(x(i))]);
        legend('snyder error', 'location', 'best');
    end
    
    % Boxplot of error
    figure;
    h = boxplot(ehat');
    set(h, 'linewidth', 2)
    grid;
    xlabel('dempographic parameters');
    ylabel('relative error');
    title(['Boxplot of errors: ' fnname]);
    
end

if M > 1
    % Plot all estimated demographic functions against the true one
    figure;
    plot(tN, Nt, 'k', 'linewidth', 2);
    hold on
    ha = plot(tN, Nhat, 'r');
    if M > 1
        set(ha, 'Color', [0.8 0.8 0.8]);
    else
        set(ha, 'linewidth', 2);
    end
    hold off
    xlim([tN(1), tN(end)]);
    grid;
    legend('true trajectory', 'estimates', 'location', 'best');
    xlabel('time');
    ylabel('N(t)');
    title([fnname ': [n, m, M] = [' [num2str(n) ' ' num2str(m) ' ' num2str(M)] ']']);
    
    % Plot demographic function summary
    figure;
    plot(tN, Nt, 'k', 'linewidth', 2);
    hold on
    plot(tN, Nhatmean, 'r--', 'linewidth', 2);
    plot(tN, Nhatlb,'g:', tN, Nhatub, 'g:', 'linewidth', 2);
    hAx = gca;
    % Add sample times
    h = stem(svec, hAx.YLim(2)*ones(size(svec)), 'g:');
    set(h, 'color', 0.8*ones(1, 3));
    hold off
    xlim([tN(1), tN(end)]);
    legend('true trajectory', 'mean', '2.5% percentile', '97.5% percentile', 'sample times', 'location', 'best');
    xlabel('time from present');
    ylabel('N(t) reconstructions');
    title([fnname ': [n, m, M] = [' [num2str(n) ' ' num2str(m) ' ' num2str(M)] ']']);
else
    % Single run cannot get quantiles so use std deviation

    % Plot demographic function summary
    figure;
    plot(tN, Nt, 'k', 'linewidth', 2);
    hold on
    plot(tN, Nhat, 'r', 'linewidth', 2);
    plot(tN, Nlb,'g:', tN, Nub, 'g:', 'linewidth', 2);
    hAx = gca;
    % Add sample times
    h = stem(svec, hAx.YLim(2)*ones(size(svec)), 'g:');
    set(h, 'color', 0.8*ones(1, 3));
    hold off
    xlim([tN(1), tN(end)]);
    legend('true trajectory', 'mean', 'mean - 2*std', 'mean + 2*std', 'sample times', 'location', 'best');
    xlabel('time');
    ylabel('N(t)');
    title([fnname ': [n, m, M] = [' [num2str(n) ' ' num2str(m) ' ' num2str(M)] ']']);
    
    % Also plot lineages with time
    figure;
    hold on
    stairs(tLin{1}, nLin{1}, 'linewidth', 2);
    hAx = gca;
    % Add sample times
    stem(svec, hAx.YLim(2)*ones(size(svec)), 'k');
    % Add coalescent times
    h = stem(tcoal{1}, hAx.YLim(2)*ones(size(tcoal{1})), 'r:');
    set(h, 'color', 0.8*ones(1, 3));
    hold off
    xlabel('time from present');
    ylabel('lineage count');
    title('Lineage through time plot');
    legend('lineages', 'sample times', 'coalescent times', 'location', 'best');
end


% Plot combined smoothed posteriors
figure;
for i = 1:numRV
    subplot(ceil(numRV/2), 2, i);
    % Main density of each parameter for both methods
    [~, ~, bw] = ksdensity(sampSny(:, i));
    newBw = 4*bw;
    [fsny, xsny] = ksdensity(sampSny(:, i), 'width', newBw);
    plot(xsny, fsny, 'linewidth', 2);
    hold on
    % Get current axes and plot true value if not empirical data
    h = gca;
    plot(x(i)*ones(1, 2), h.YLim, 'k', 'linewidth', 2);
    hold off
    grid;
    xlabel('parameter space');
    ylabel('smoothed pdf');
    hold off
    title(['Posterior: x_' num2str(i) ' = ' num2str(x(i))]);
    legend('snyder', 'true', 'location', 'best');
end


% Joint final posterior if bivariate
if numRV == 2 && nEstRV == 2 && M == 1
    % Reshape last posterior to 2D form for bivariate plotting
    qtemp = reshape(qnlast{1}, mi);
    if ~isnan(qtemp)
        % Classic joint posterior
        figure;
        meshz(xset{1}, xset{2}, qtemp');
        xlabel('x_1');
        ylabel('x_2');
        zlabel('P(x_1, x_2 | data)');
        title(['Joint posterior with true [x1 x2] = ' num2str(x(1)) ' ' num2str(x(2))]);
        
        % Contour plot
        figure;
        contour(xset{1}, xset{2}, qtemp');
        h = gca;
        xLim = h.XLim;
        yLim = h.YLim;
        hold on
        plot(x(1)*[1 1], yLim, 'k');
        plot(xLim, x(2)*[1 1], 'k');
        xlabel('x_1');
        ylabel('x_2');
        legend('contours', 'true values', 'location', 'best');
        title('Joint posterior contour plot');
        grid;
        
    else
        warning('Mat:qNaN', 'qtemp has NaN values');
    end
end

%% Final storage of data and saving of dependencies

% Get all m files called in creating this simulated data
currName = dbstack; % gives a struct with current m file name
currName = currName.file;
[fList,pList] = matlab.codetools.requiredFilesAndProducts(currName);

% Remove unnecessary variables
vlist = {'qev', 'tset', 'Nhat', 'NhatSq', 'tree', 'sampSny', 'qnlast', 'xhatEv', 'xsetRow', 'IDMx'};

% Store all data
clear(vlist{:});
save([fnname '_' num2str(n) '_' num2str(M) '_' nSampTimes]);
