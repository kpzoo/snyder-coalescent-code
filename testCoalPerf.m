% Script to test timing performance of Snyder filter on coalescent
clearvars
clc
close all

% Tree no. tips fixed for all tests
n = 100;
% Always test for 100 runs and average
M = 20;
% Variable list to clear between tests
vlist = {'tIter', 'demog', 'est', 'mi', 'm', 'numRV', 'fnid'};


%% Test 1: exponential - increase grid resolution and get estimates

% Exponential id
fnid = 2;
numRV = 4;

% Increase grid size components
nRuns = 5;
mcomp = linspace(5, 20, nRuns);
m = zeros(1, nRuns);

% Variables to store results
tIter = cell(1, nRuns);
demog = cell(1, nRuns);
est = cell(1, nRuns);
tavg = zeros(1, nRuns);

% Run in a loop
for i = 1:nRuns
    % Calculate mi and m for each 
    mi = mcomp(i)*ones(1, numRV);
    m(i) = prod(mi);
    % Main code runs M sims of n tip trees over grid size m
    [tIter{i}, demog{i}, est{i}] = runSnyCoalPerf(n, M, mi, m(i), fnid);
    % Average time for run
    tavg(i) = mean(diff([0, tIter{i}]));
    disp(['Finished run ' num2str(i) ' of ' num2str(nRuns)]);
    disp('***************************************************************');
end

% Save data and clear other variables
save(['test1_' num2str(fnid) '.mat']);
clear(vlist{:}); 

%% Plotting and Analysis for Test 1

% Timing statistics
tdiff = cell(1, nRuns);
tmed = zeros(1, nRuns);
tub = zeros(1, nRuns);
tlb = zeros(1, nRuns);
for i = 1:nRuns
    % Time of each iteration
    tdiff{i} = diff([0, tIter{i}]);
    tmed(i) = mean(tdiff{i});
    tub(i) = quantile(tdiff{i}, 0.975);
    tlb(i) = quantile(tdiff{i}, 1-0.975);
end

% Parametric estimate statistics
J1 = zeros(M, nRuns);
J2 = zeros(M, nRuns);
for i = 1:nRuns
    % Parametric info for M trees in each i
    x = est{i}.x;
    xhat = est{i}.xhat;
    sampSny = est{i}.sampSny;
    % Error metric of relative MMSE
    xx = repmat(x, M, 1);
    J = 100*(1 - xhat'./xx).^2;
    J1(1:M, i) = J(:, 1);
    J2(1:M, i) = J(:, 2);
end

% Demographic statistics
%Nval = zeros(1000, nRuns);
Nmean = zeros(1, nRuns);
Nub = zeros(1, nRuns);
Nlb = zeros(1, nRuns);

for i = 1:nRuns
    % Demographic info for M trees in each i
    N = demog{i}.Nt';
    Nhat = demog{i}.Nhatmean;
    
    % Get trajectory square error
    %Nval(1:1000, i) = 100*((1-Nhat./N).^2);
    Nsq = 100*((1-Nhat./N)).^2;
    Nmean(i) = mean(Nsq);
    Nlb(i) = quantile(Nsq, 0.025);
    Nub(i) = quantile(Nsq, 1-0.025);
end

% Figure with all information about test
figure;
subplot(2, 2, 1);
% Timing information with grid size
plot(mcomp, tmed, 'r', mcomp, tlb, 'g:', mcomp, tub, 'g:', 'linewidth', 2);
grid;
xlabel('grid size per dimension');
ylabel('execution time (mins)');
title('Exponential function grid complexity');
legend('median', '2.5% quantile', '97.5% quantile', 'location', 'best');
subplot(2, 2, 2);
% Demographic accuracy
plot(mcomp, Nmean, 'r', mcomp, Nlb, 'g:', mcomp, Nub, 'g:', 'linewidth', 2);
grid;
xlabel('grid size per dimension');
ylabel('% relative square error N(t)');
title('Exponential function demographic accuracy');
legend('median', '2.5% quantile', '97.5% quantile', 'location', 'best');
subplot(2, 2, 3);
% Parameter 1 accuracy
boxplot(J1, mcomp);
h = gca;
h.XTickLabelRotation = 45;
xlabel('grid size per dimension');
ylabel('% relative square error');
title('Exponential function x_1 accuracy');
subplot(2, 2, 4);
% Parameter 2 accuracy
boxplot(J2, mcomp);
h = gca;
h.XTickLabelRotation = 45;
xlabel('grid size per dimension');
ylabel('% relative square error');
title('Exponential function x_2 accuracy');
    

% %% Test 2: logistic with parameters fixed on a space of m points
% 
% % Logistic id
% fnid = 4;
% % Fixed grid size
% numRV = 4;
% m = 20^numRV;
% % Different no. params
% mi{1} = [m ones(1, numRV-1)];
% mi{2} = [sqrt(m)*ones(1, numRV-2) ones(1, numRV-2)];
% mi{3} = [floor(m^(1/3))*ones(1, numRV-2) ceil(m^(1/3)) 1];
% mi{4} = m^(1/4)*ones(1, numRV);
% 
% % Variables to store results
% tIter = cell(1, numRV);
% demog = cell(1, numRV);
% est = cell(1, numRV);
% tavg = zeros(1, numRV);
% 
% % Run in a loop
% for i = 1:numRV
%     % Main code runs M sims of n tip trees over grid size m
%     [tIter{i}, demog{i}, est{i}] = runSnyCoalPerf(n, M, mi{i}, m, fnid);
%     % Average time for run
%     tavg(i) = mean(diff([0, tIter{i}]));
%     disp(['Finished run ' num2str(i) ' of ' num2str(numRV)]);
%     disp('***************************************************************');
% end
% 
% % Save data and clear other variables
% save('test2.mat');
% clear(vlist{:}); 
