function [CurrentLoop,OutputLineLength] = Deflated_Newton(obj_fun,x0,DeflatedPts,params,RecordIterates)
% Deflated_Newton_Method
% Performs deflated Newton's method on the equation grad_f(x) = 0, with
% deflation operator mu(x), where
%    - obj_fun(x) returns f(x)=(1/2)||r(x)||^2, r(x), J_r(x) and Hess_r(x)
%    - x0 is the initial guess
%    - deflation returns mu(x) and grad_mu(x), which deflates out the
%      points stored in DeflatedPts
%    - params contains tuneables such as linesearch merit function and
%      hyperparameters
%    - RecordIterates is a Boolean

% See Algorithm ##### in Deflation Techniques for Finding Multiple Local
% Minima of a Nonlinear Least Squares Problem, Riley, Webb, Baker, 2024.

% Initial Values
NIter = 0;  FuncCount = 0;   x = x0;   stop = false;  OutputLineLength =[];
if RecordIterates, CurrentLoop.Iterates = x0; end
Solver = str2func(params.method.LinearSolver);
lambda = params.method.Regularisation;
if isempty(lambda); regularise = false; else
    regularise = true;
    if isnumeric(lambda); lambda = @(varargin)lambda; end
end
% Main loop
while stop == false
    % Calculate residual, Jacobian and Hessian of R
    [X.F, X.R, X.J, X.H] = obj_fun(x, params.method.constants);   FuncCount = FuncCount +1;
    % Calculate the second derivative terms in the Hessian of F:
    X.S = zeros(length(x));
    for i = 1:length(X.R)
        X.S = X.S + X.H(:,:,i)*X.R(i);
    end
    % Calculate deflation operators
    [X.Mu,X.gradMu] = deflation(DeflatedPts, x, params.deflation);

    % Calculate undeflated Newton Step: p = - Hess_f \ grad_f
    Hess = X.J'*X.J+X.S;
    if ~regularise
        p = - Solver(Hess, X.J'*X.R);
    else
        p = - Solver(lambda(NIter)*diag(diag(Hess))+Hess, X.J'*X.R);
    end

    % Calculate effect on log(Mu) in p direction
    pTgradNu = X.gradMu*p/X.Mu;

    if (pTgradNu < params.deflation.epsilon) || (abs(pTgradNu - 1)<1e-8) % Zone 1 (or deflation blowup)
        % Perform a line search on the undeflated problem
        [alpha, FCount] = Linesearch(x,p,DeflatedPts, params, params.linesearch.merit, obj_fun, X);
        FuncCount = FuncCount+FCount;
        if alpha <= params.linesearch.merit.minalpha
            if rank(X.J) < length(x)
                CurrentLoop.ConvergenceFlag = "Merit line search terminated with rank deficient Jacobian";
            else
                CurrentLoop.ConvergenceFlag = "Merit line search terminated";
            end
            break
        end
    elseif pTgradNu > 1 % Zone 3
        % Could use a linesearch on mu here. Not necessary in practice.
        alpha = params.linesearch.merit.alpha0/(1-pTgradNu);
    else % Zone 2
        % Perform a deflated step (no line search, but possibly alpha0 damping)
        alpha = params.linesearch.merit.alpha0/(1-pTgradNu);
    end
    % Take the deflated Newton step:
    x = x + alpha*p;
    NIter = NIter + 1;

    % Apply stopping criteria
    [stop,CurrentLoop.ConvergenceFlag] = ismin(X.F, p, X.J'*X.R, NIter, params.convergence);
    % Print out convergence info
    if params.method.Verbose
        if NIter>1
            fprintf(repmat('\b',1,OutputLineLength))
        end
        OutputLineLength = fprintf('k = %d; f(x) = %d; |gradf(x)| = %d; alpha = %d; \n', NIter,X.F,norm(X.J'*X.R),alpha);
    end
    % Save iterates for plotting
    if RecordIterates, CurrentLoop.Iterates = [CurrentLoop.Iterates, x]; end
end
CurrentLoop.DeflatedPoint = x;
CurrentLoop.ErrorAtDeflatedPoint = obj_fun(x, params.method.constants);
CurrentLoop.NIter = NIter;
CurrentLoop.FuncCount= FuncCount;

end
