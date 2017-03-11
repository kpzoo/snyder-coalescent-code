% Function to perform Snyder filtering and estimation for coalescent
function [qev, tev] = snyderCoalHet(fn, q0, tLin, nLin)

% Assumptions and modifications
% - only for fnid = 1 and 5 (exponential and con-ecp-con)
% - for heterochronous sampling so input tLin not tcoal
% - output is the posterior at each event time
% - allows linear vs non-linear
% - outputs qn, joint posterior at each tn
% - works for multivariable functions defined in fn
% - includes new functions for calculating lamdiag and N
% - uses matrices xsetMx etc defined in fn

% Extract inputs related to grid and no. observations
facSort = fn.facSort;
nData = fn.nData;
m = fn.m;

% Vectorised demographic function that is a scalar of time
Nt = fn.fnTime;

% Vector of sampled times
svec = fn.svec;

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nData+1, m);
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);

% As functions are on probabilities constrain to be non-negative
options = odeset('NonNegative', 1:m);
%options= odeset('RelTol', 1e-9, 'NonNegative', 1:m);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain appropriate binomial factor
    binfac = facSort(nLin(i));
    
    % Solve posterior if nLin >= 2 else do nothing
    if binfac == 0
        % Maintain posterior at last value if lineages fall to 1
        qev(i+1, :) = qev(i, :);
        disp(['Lineages prematurely fell to 1 at ' num2str(i)]);
    else
        % Solve Snyder ODEs continuously with setting of options
        [tsol, qsol] = ode113(@(ts, y) odeSnyderCoal(ts, y, binfac, fn),...
            [tLin(i) tLin(i+1)], qev(i, :)', options);
        
        % Normalise the posterior probabilities - only if linear ODE form
        if fn.linBool
            % No. of solution elements between coalescent events is variable
            % hence recalculate size(qsol, 1)
            for j = 1:size(qsol, 1)
                qsol(j, :) = qsol(j, :)/(sum(qsol(j, :)));
            end
        end
        
        % Assign the output values of time and posteriors
        qset{i} = qsol;
        tset{i} = tsol;
        qev(i+1, :) = qsol(end, :);
        tev(i+1) = tsol(end);
        
        % Only perturb if there is a new event as opposed to sample
        if ~any(tLin(i+1) == svec)
            % Perturb the q posterior for the new event
            pert = binfac./Nt(tsol(end));
            pert(pert < 0) = 0;
            qev(i+1, :) = qev(i+1, :).*pert./(qev(i+1, :)*pert');
        else
            % Maintain posterior since new samples only
           qev(i+1, :) = qsol(end, :);
        end
    end
end

        