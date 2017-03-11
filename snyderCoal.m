% Function to perform Snyder filtering and estimation for coalescent
function [qev, tev] = snyderCoal(fn, q0, tcoal)

% Assumptions and modifications
% - output is the posterior at each event time
% - allows linear vs non-linear
% - outputs qn, joint posterior at each tn
% - works for multivariable functions defined in fn
% - includes new functions for calculating lamdiag and N
% - uses matrices xsetMx etc defined in fn

% Extract inputs related to grid and no. observations
fac = fn.fac;
nData = fn.nData;
m = fn.m;

% Vectorised demographic function that is a scalar of time
Nt = fn.fnTime;

% Posterior vectors on events, qnoev is for non-updated event q's
qev = zeros(nData+1, m);
qev(1, :) = q0;
tev = zeros(nData, 1);

% Cell to save output of ODE solver and set options
qset = cell(1, 1);
tset = cell(1, 1);

% As functions are on probabilities constrain to be non-negative
options = odeset('NonNegative', 1:m);

% Run Snyder filter across the time series
for i = 1:nData
    % Obtain appropriate binomial factor
    binfac = fac(i);
    
    % Solve Snyder ODEs continuously with setting of options
    if fn.id ~= 2
        % This solver generally works well
        [tsol, qsol] = ode113(@(ts, y) odeSnyderCoal(ts, y, binfac, fn),...
            [tcoal(i) tcoal(i+1)], qev(i, :)', options);
    else
        % Sinusoidal function works better with this solver
        [tsol, qsol] = ode45(@(ts, y) odeSnyderCoal(ts, y, binfac, fn),...
            [tcoal(i) tcoal(i+1)], qev(i, :)', options);
    end
    
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
    
    % Perturb the q posterior for the new event
    pert = binfac./Nt(tsol(end));
    %[~, pert] = getTimeVaryingN(fn, tsol(end), binfac, fn.xsetMx);
    pert(pert < 0) = 0;
    qev(i+1, :) = qsol(end, :);
    tev(i+1) = tsol(end);
    qev(i+1, :) = qev(i+1, :).*pert./(qev(i+1, :)*pert');
end

        