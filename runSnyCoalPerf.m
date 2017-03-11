%% Perform Snyder filtering for multi-variable N(t) - performance function
function [tIter, demog, est] = runSnyCoalPerf(n, M, mi, m, fnid)
% Usage:
% - write parametric function as a new fnid
% - simulates a coalescent and then infers its parameters

% Assumptions and Modifications
% - requires a parametric form for the demographic function
% - functional form of runSnyCoal, no plots but variables for calculation
% - used rejection sampling to simulate NHPP coalescent
% - the data is k-coalescent times and the process is self-exciting


%% Set simulation variables

% Start timing
tic;

% Perform a batch of M runs with trees of n tips
nData = n - 1; % no. coalescent events

%[xmin, xmax, fnname, fnForm, numRV] = setDemogFn(fnid, n);
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
fn.linBool = 0;

% Check for pre-specified parameters not to be estimated - mi values of 1
specParam = find(mi == 1);
nEstRV = numRV - length(specParam);
% Modify spaces if specify parameters so that they are set to min values
if nEstRV > 0
    xmin(specParam) = xmax(specParam);
end

% Get coalescent binomial factors for lineages and change notation
nset = n:-1:2;
fac = nset.*(nset-1)/2;

% Add more fields to function struc
fn.nData = nData;
fn.mi = mi;
fn.m = m;
fn.fac = fac;
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
%xsetCol  = mat2cell(xsetMx', fn.gridSz(2), ones(1, fn.gridSz(1)));
xsetRow = mat2cell(xsetMx, ones(1, fn.gridSz(1)), fn.gridSz(2));

% Pre-calculate the grid with time as the only input
%fn.fnTime = @(tx) fn.fn2(xsetCol, tx);
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


%% Main loop doing simulation and filtering on same tree

% Store parameter estimates
xhat = zeros(numRV, M);

% Times for evaluating population of each tree
nPts = 1000;
tset = zeros(nPts, M);
Nhat = zeros(nPts, M);

% Function to calculate Nhat for a given time and posterior
NhatFn = @(q,tx)q*fn.fnTime(tx)'; 

% Coalescent times, tree and corresponding waiting times
twait = cell(1, M);
tcoal = cell(1, M);
tree = cell(1, M);

% Variables for batch results
qnlast = cell(1, M);
qmarg = cell(1, M);
xhatEv = cell(1, M);

for ii = 1:M
    % Generate the coalescent data using rejection sampling
    [twait{ii}, tcoal{ii}] = getCoalData(fn);
    C = cumsum(twait{ii});
    % Coalescent time based tree
    [~, ~, tree{ii}] = getKingmanTree(n, C');
    
    % Filter the simulated event times - main filter function
    [qev, ~] = snyderCoal(fn, q0, tcoal{ii});
    
    % Last posterior and marginals and times
    qnlast{ii} = qev(end, :);
    [qmarg{ii}, ~] = marginalise(numRV, IDMx, qnlast{ii}, mi);
    % Conditional estimates, xhatEv gives estimates at event times
    xhatEv{ii} = qev*xsetMx';
    
    % Time points at which population size will be estimated
    tset(:, ii) = linspace(tcoal{ii}(1), tcoal{ii}(end), nPts);

    % Store batch results with xhat having rows for numRVs
    xhat(:, ii) = xhatEv{ii}(end, :);
    
    disp(['Finished event ' num2str(ii) ' of ' num2str(M)]);
    disp('...............................................................');
    tIter(ii) = toc; 
end

% True population, min(max) as smallest tree length wanted
demog.tN = linspace(min(min(tset)), min(max(tset)), nPts);
demog.Nt = fn.fn1(x, demog.tN);

% Get all population estimates on same time points
for ii = 1:M
    % Collapse Nhat function with correct posterior
    Nq = @(tx) NhatFn(qnlast{ii}, tx);
    % Calculate N(t) posterior estimate across time set
    Nhat(:, ii) = arrayfun(Nq, demog.tN);
end

% Statistics of demographics
demog.Nhatmean = mean(Nhat, 2);
demog.Nhatlb = quantile(Nhat', 0.025);
demog.Nhatub = quantile(Nhat', 1-0.025);

% Variables for samples from posteriors
nSamps = 1000;
sampSny = zeros(M*nSamps, numRV);

% Samople from the discrete posteriors
for ii = 1:M
    % Get start and end indices in sample array
    id1 = (ii-1)*nSamps + 1;
    id2 = id1 + nSamps - 1;
    for i = 1:numRV
        % Frequentist approximation to probabilities
        sampSny(id1:id2, i) = datasample(xset{i}, nSamps, 'Weights', qmarg{ii}{i});
    end
end

% Save estimation information in structure
%est.qnlast = qnlast;
%est.qmarg = qmarg;
est.xset = xset;
est.x = x;
est.xhat = xhat;
est.sampSny = sampSny;

% Clock time set to minutes
tIter = tIter/60;
disp(['Total time for ' num2str(M) ' iterations is ' num2str(tIter(end))]);
tavgIter = mean(diff([0, tIter]));
disp(['Average time with [n m M] = ' [num2str(n) ' ' num2str(m) ' ' num2str(M)] ' is ' num2str(tavgIter)]);


