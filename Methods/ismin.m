function [test,flag] = ismin(f,p,gradf,NIter,ConvergenceParams)
%ISMIN Tests stopping criteria at current iteration
%
% stop = ismin(f,p,x,gradf,NIter,ConvergenceParams) outputs true if
% stopping critera reached, and false if else. f is the current value of
% the function evaluation, p is the current update step, x is the current
% iteration, gradf is the gradient of f at the current iteration, NIter is
% the iteration number, and ConvergenceParams is a structure containing the
% tolerences for each stopping criterion.
%
% [stop,flag] = ismin(f,p,gradf,NIter,ConvergenceParams) also outputs
% which stopping criteria was reached.

test = false;
flag = "";
if norm(p)<ConvergenceParams.StepTolerance
    flag = 'Step Size below tolerance';
    test=true;
elseif isnan(f)
    test = true ;
    flag = 'NaN';
elseif f <ConvergenceParams.ObjectiveTolerance
    flag = 'Objective less than tolerance';
    test = true ;
elseif norm(gradf)<ConvergenceParams.GradientTolerance
    flag = 'Gradient less than tolerance';
    test=true;
elseif NIter>=ConvergenceParams.MaximumIterations
    flag ='Max Iterations reached';
    test = true  ;
end
end
