% Function to setup linear Snyder ODE set for use with ODE solvers for RVs
function dy = odeSnyderCoal(ts, y, binfac, fn)

% Assumptions
% - includes non-linear option - linBool = 0
% - to improve speed replaced y'*diag(lamtdiag) with y'.*lamtdiag
% - generalised for non-homogeneous simulations of N(t)
% - y is a column vector, ts a vector just included for ode113
% - binfac is the appropriate binomial factor based on no. events


% Ensure input dimension correct
sz = size(y);
if sz(2) ~= 1
    error('y is not a column vector');
end

% Diagonal of time dependent rate matrix for various functional forms
Nt = fn.fnTime;
rateDiag = binfac./Nt(ts);
%[~, rateDiag] = getTimeVaryingN(fn, ts, binfac, fn.xsetMx);
rateDiag(rateDiag < 0) = 0;

% Solve linear differential equation set - RV filtering
if fn.linBool
    dy = y'.*(-rateDiag);
else
    nonLinDiag = rateDiag*y;
    dy = y'.*(-rateDiag + nonLinDiag);
end

% Ensure output is column vector assuming input was
dy = dy';