% Function to set the parametric form of N(t) for coalescent analysis
function [xmin, xmax, fnname, fn1, fn2, numRV] = setDemogFn2(fnid, n)

% Usage
% - write parametric form of functions of interest here and define the grid
% space limits on each parameter as well as assign a name
% - function indexed by fnid
% - fnForm defined for vector of params x and possible time vector tx
% - fnForm converts x to a scalar for each tx entry


% Assumptions and modifications
% - the x(i) version is allows for vector tx
% - the y{i} version allows for vector parameter sets, scalar tx
% - spaces defined as multiples of n to help N(t) >> n assumption


%% Main definitions 

% Set of function strings
fnstr = {'exponential', 'sinusoidal', 'kingman',...
    'logistic growth', 'piecewise exponential'};
fnname = fnstr{fnid};
disp(['Simulation of function: ' fnname]);

% Choose function and define settings
switch(fnid)    
    case 1
        % Exponential
        xmin = [100*n 0.1];
        xmax = [1000*n 10];
        fn1 = @(x, tx) x(1)*exp(-x(2)*tx);
        fn2 = @(y, tx) y{1}.*exp(-y{2}*tx);
    case 2
        % Sinusoidal
        xmin = [100*n 0.001 0 1100*n];
        xmax = [1000*n 0.01 pi/200 1200*n];
        fn1 = @(x, tx) x(1)*(sin(x(2)*tx + x(3))) + x(4);
        fn2 = @(y, tx) y{1}.*sin(y{2}*tx + y{3}) + y{4};
    case 3
        % Constant
        xmin = 100*n;
        xmax = 1000*n;
        % Output always of length equal to tx
        fn1 = @(x, tx) x(1)*ones(size(tx));
        fn2 = @(y, tx) y{1};
    case 4
        % Logistic
        xmin = [100*n 0.1 1000 50*n];
        xmax = [1000*n 10 5000 100*n];
        fn1 = @(x, tx) x(1)*(1 + exp(-x(2)*x(3)))./(1 + exp(-x(2)*(x(3) - tx))) + x(4);
        fn2 = @(y, tx) y{1}.*(1 + exp(-y{2}.*y{3}))./(1 + exp(-y{2}.*(y{3} - tx))) + y{4}; 
    case 5
        % Pybus piecewise exponential
        %xmin = [1000*n 0.1 5 40];   % for hcv
        %xmax = [2000*n 0.75 30 70];
        xmin = [10000 0 50 20];
        xmax = [50000 0.75 75 50];
        fn1 = @con_exp_con1;
        fn2 = @con_exp_con2;

    otherwise
        error('Incorrect function index specified');
end

% No. parameters for each
numRV = length(xmin);
numCheck = length(xmax);
if numRV ~= numCheck
    error('Incorrect range specification');
end



%% Sub functions

% The con-exp-con parametrised with x(4) as length of exp phase
function fnval = con_exp_con1(x, tx)
% Con-exp-con function definition for vector tx
I1 = tx <= x(3);
I2 = tx > x(3) & tx < x(3) + x(4);
I3 = tx >= x(3) + x(4);
fnval = x(1)*I1 + x(1)*exp(-x(2)*(tx - x(3))).*I2 + ...
    x(1)*exp(-x(2)*x(4)).*I3;

function fnval = con_exp_con2(y, tx)
% Con-exp-con function definition for vector parameters, scalar tx
I1 = tx <= y{3};
I2 = tx > y{3} & tx < y{3} + y{4};
I3 = tx >= y{3} + y{4};
fnval = y{1}.*I1 + y{1}.*exp(-y{2}.*(tx - y{3})).*I2 + ...
    y{1}.*exp(-y{2}.*y{4}).*I3;

