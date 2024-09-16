function [CurrentLoop,OutputLineLength] = Bad_Deflated_GaussNewton(obj_fun,x0,DeflatedPts,params,RecordIterates)
% Deflated_GaussNewton_Type2_Method
% Performs deflated Gauss-Newton method of type 2 on nonlinear least squares
% problem for r(x), with deflation operator mu(x), where
%    - obj_fun(x) returns F(x)=(1/2)||r(x)||^2, r(x), J_r(x)
%    - x0 is the initial guess
%    - deflation returns mu(x) and grad_mu(x), which deflates out the
%      points stored in DeflatedPts
%    - params contains tuneables such as linesearch merit function and
%      hyperparameters
%    - RecordIterates is a Boolean

% See Algorithm ##### in Deflation Techniques for Finding Multiple Local
% Minima of a Nonlinear Least Squares Problem, Riley, Webb, Baker, 2024.

% Initial Values
NIter = 0;  FuncCount = 0; x = x0;   stop = false;   OutputLineLength =[];
if RecordIterates, CurrentLoop.Iterates = x0; end
lambda = params.method.Regularisation;
if isempty(lambda); regularise = false; else
    regularise = true;
    if isnumeric(lambda); lambda = @(varargin)lambda; end
end
% Main Loop
while stop == false
    % Calculate residual, Jacobian of R
    [X.F,X.R,X.J] = obj_fun(x, params.method.constants);  FuncCount = FuncCount +1;
    % Calculate deflation operators
    [X.Mu,X.gradMu] = deflation(DeflatedPts, x, params.deflation);

    % compute p and pTgradEta, while saving the QR factorisation
    % if regularise
    %     [Q,R] = qr([X.J;lambda(NIter)*diag(diag(X.J'*X.J))],0);
    % else
        [Q,R] = qr(X.J,0);
    % end
    X.gradEta = X.gradMu.'/X.Mu;
    p = -R \ (Q'* X.R );
    pTgradEta = dot(X.gradEta,p);


    if (pTgradEta <= params.deflation.epsilon) || (abs(pTgradEta - 1)<1e-8) % Zone 1 (or deflation blowup)
        % Perform a line search on the undeflated problem
        [alpha, FCount] = Linesearch(x,p,DeflatedPts, params, params.linesearch.merit, obj_fun, X);
        FuncCount = FuncCount+FCount;
        if (alpha<=params.linesearch.merit.minalpha)
            if rank(X.J)<length(x)
                CurrentLoop.ConvergenceFlag = "Merit line search terminated with rank deficient Jacobian";
            else
                CurrentLoop.ConvergenceFlag = "Merit line search terminated";
            end
            break
        end
    else
        uPu = dot(X.R,X.R - Q*(Q'*X.R));
        beta = 1 - pTgradEta;
        w = R'*R \ X.gradEta;
        omega = uPu * dot(X.gradEta,w) + beta^2;
        deflated_p = (beta * p - uPu * w) / omega;

        if pTgradEta > 1 % Zone 3
            p = deflated_p;
            % By default this linesearch just returns alpha0.
            [alpha,FCount] = Linesearch(x, p, DeflatedPts, params, params.linesearch.Mu, obj_fun, X);
            FuncCount = FuncCount+FCount;
            if (alpha<=params.linesearch.merit.minalpha)
                if rank(X.J)<length(x)
                    CurrentLoop.ConvergenceFlag = "Mu line search terminated with rank deficient Jacobian";
                else
                    CurrentLoop.ConvergenceFlag = "Mu line search terminated";
                end
                break
            end
        else % Zone 2
            p = deflated_p;
            % By default this linesearch just returns alpha0.
            [alpha,FCount] = Linesearch(x, p, DeflatedPts, params, params.linesearch.deflatedmerit, obj_fun, X);
            FuncCount = FuncCount+FCount;
            if (alpha<=params.linesearch.merit.minalpha)
                if rank(X.J)<length(x)
                    CurrentLoop.ConvergenceFlag = "Deflated Merit line search terminated with rank deficient Jacobian";
                else
                    CurrentLoop.ConvergenceFlag = "Deflated Merit line search terminated";
                end
                break
            end
        end
    end

    % Take the deflated Gauss-Newton Type 2 step:
    x = x + alpha*p;
    NIter = NIter + 1;

    % Apply stopping criteria
    [stop, CurrentLoop.ConvergenceFlag] = ismin(X.F, p, X.J'*X.R, NIter, params.convergence);
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
CurrentLoop.FuncCount = FuncCount;

end