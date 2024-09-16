function [Iterations,res] = DMin(obj_fun,x0,varargin)
%DMIN Finds multiple minima of the problem using deflation.
%
% Y= DMIN(obj_fun,x0)   Outputs a structure that contains the points that
% minimise the problem 'obj_fun' initialised as x0and the number of Iterations till
% convergence; the type of 'convergence' achieved; all the iterations
% made until convergence; the finall error at convergence. Where
% obj_fun(x,constants) is a function handle that returns the current
%    error of the least squares formulation, given the current guess, x,
%    and optional parameters, constants. It must also  output the
%    Residual and the Jacbian and (optionally) Hessian (if using the Newton
% method) of said residual of the least squares problem: [F,R,J,H] =
% obj_fun(x,constants)
%
%
%DMIN PARAMETERS
%NDeflations - The number of local minima you wish to find, the output -
% [ {1} | positive integer ]
%Method - The optimisation method desired - [ {Newton} | GaussNewtonT1 |
%           GaussNewtonT2 ]
%MaxIter - The integer value  of Maximum iterations per deflation -
%           [ {1000} | positive integer ]
%Linesearch - line search method - [ No | Basic | {Armijo} | Quadratic ]
%theta - The deflation exponent - [ {2} | positive integer | exp ]
%c1 - the armijo line search parameter - [ {1e-4} | positive scalar ]
%alpha0 - Initial value of the line search parameter each iteration - [ 1
%           | positive scalar ]
%Tau - The value of the decrease in the line search parameter - [ {0.5} |
%         scalar in [0,1] ]
%Minalpha - The minimum value of the line search parameter - [ 1e-10 |
%       positive scalar ]
%Verbose - Output value of objective function and gradient at each
%            iteration - [ true | {false} ]
%Scaled - Scale the variables by the value of the initial guess - [ true |
%       {false} ]
%Constants - Any constants the objective funtion provided requires - [ {}
%        | struct ]
%Regularisation - The value regularing parameter - [ {0} ]
%LinearSolver - The linear solver used - [ {mldivide} | lsqminnorm ]
%PreviouslyFoundIterations - A struct of prevously found points to be deflated
%        - [ {} | l \times x matrix ], where x is the number of previous
%        deflations
%ConvergenceFlags - The flags that are considered to mean convergence - [
%        {Objective less than tolerance} | {Gradient less than tolerance} |
%       Step Size too small | Line search failed | Max Iterations reached |
%        Divergence Detected | NaN ]
%DeflateLinesearch - Line searches should be applied to the deflated
%        system - [ {true} | false ]
%Sigma - The value of the shift of the deflation - [ {1} | scalar ]
%SingleShift - Deflation shift should only be applied once [ true |
%        {false} ]
%Epsilon - Tolerence on the application of deflation or line search - [
%        {0.01} | scalar ]
%NormWeighting - Weighting applied to the norms in calculation of
%           deflation operators - [ {} | l \times l matrix ]
%GradientTolerance - Stopping criterion based on gradient of function -
%        [ {0} | positive scalar ]
%StepTolerance - Stopping criterion based on difference of consecutive
%        steps - [ {1e-8} | positive scalar ]
%ObjectiveTolerance - Stopping criterion based on value of objective
%        function - [ {0} | positive scalar ]

%
% Setup Default Parameters
%
defaultMethod ='Good_GN';
defaultNDeflations = 1;
defaultVerbose = false;
defaultScaled = false;
defaultConstants=[];
defaultRegularisation = [];
defaultLinearSolver = "mldivide";

defaultPreviouslyFoundIterations = [];
defaultConvergenceFlags = ["Objective less than tolerance","Gradient less than tolerance","Step Size below tolerance","Merit line search terminated"];
defaultRecordIterates = true;
defaultMaxNonMinima = 10;
defaultRecordTimes = false;

defaultLinesearch = 'Armijo';
defaultC1 = 1e-4;
defaultTau = 0.5;
defaultAlpha0 = 1;
defaultMinalpha = 1e-18;
defaultMuLinesearch = "No";
defaultDeflatedLinesearch = "No";

defaultTheta = 2;
defaultSigma = 1;
defaultSingleShift = false;
defaultEpsilon = 0.01;
defaultNormWeighting = speye(length(x0));

defaultStepTolerance = 1e-6;
defaultGradientTolerance = 0;
defaultObjectiveTolerance = 0;
defaultMaxIter = 1000;

%
% Parse Input Parameters
%
isnumericscalar = @(x) isscalar(x) && isnumeric(x);
isstringorchar = @(x) isstring(x) || ischar(x);
isvalidstruct =@(x) (isstruct(x)&&4==sum(contains(fieldnames(x),[{'DeflatedPoint'}, ...
    {'ErrorAtDeflatedPoint'},{'NIter'},{'ConvergenceFlag'}])))||isempty(x);
IP = inputParser;
addRequired(IP,'obj_fun')
addRequired(IP,'x0')

addParameter(IP,'Method',defaultMethod,isstringorchar)
addParameter(IP,'NDeflations',defaultNDeflations,isnumericscalar)
addParameter(IP,'Verbose',defaultVerbose,@islogical)
addParameter(IP,'Scaled',defaultScaled,@islogical)
addParameter(IP,'constants',defaultConstants)
addParameter(IP,'Regularisation',defaultRegularisation)
addParameter(IP,'LinearSolver',defaultLinearSolver,@(x)contains(x,["mldivide","lsqminnorm"]))

addParameter(IP,'PreviouslyFoundIterations',defaultPreviouslyFoundIterations,isvalidstruct)
addParameter(IP,'ConvergenceFlags',defaultConvergenceFlags,isstringorchar)
addParameter(IP,'RecordIterates',defaultRecordIterates,@islogical)
addParameter(IP,'MaxNonMinima',defaultMaxNonMinima,isnumericscalar)
addParameter(IP,'RecordTimes',defaultRecordTimes,@islogical)

addParameter(IP,'Linesearch',defaultLinesearch,isstringorchar)
addParameter(IP,'c1',defaultC1,isnumericscalar)
addParameter(IP,'tau',defaultTau,isnumericscalar)
addParameter(IP,'alpha0',defaultAlpha0,isnumericscalar)
addParameter(IP,'minalpha',defaultMinalpha,isnumericscalar)
addParameter(IP,'MuLinesearch',defaultMuLinesearch,isstringorchar)
addParameter(IP,'Mualpha0',defaultAlpha0,isnumericscalar)

addParameter(IP,'DeflatedLinesearch',defaultDeflatedLinesearch,isstringorchar)
addParameter(IP,'Deflatedalpha0',[],isnumericscalar)
addParameter(IP,'Deflatedc1',[],isnumericscalar)
addParameter(IP,'Deflatedminalpha',[],isnumericscalar)
addParameter(IP,'Deflatedtau',[],isnumericscalar)


addParameter(IP,'theta',defaultTheta,@(x)isnumericscalar(x)||isstringorchar(x))
addParameter(IP,'sigma',defaultSigma,isnumericscalar)
addParameter(IP,'SingleShift',defaultSingleShift,@islogical)
addParameter(IP,'epsilon',defaultEpsilon,isnumericscalar)
addParameter(IP,'NormWeighting',defaultNormWeighting)

addParameter(IP,'ObjectiveTolerance',defaultObjectiveTolerance,isnumericscalar)
addParameter(IP,'GradientTolerance',defaultGradientTolerance,isnumericscalar)
addParameter(IP,'StepTolerance',defaultStepTolerance,isnumericscalar)
addParameter(IP,'MaxIter',defaultMaxIter,isnumericscalar)

IP.parse(obj_fun,x0,varargin{:})
res = IP.Results;

%
% Setup Method Parameters
%
params.method.StepMethod = res.Method;
params.method.Verbose = res.Verbose;
params.method.Scaled = res.Scaled;
params.method.constants = res.constants;
params.method.Regularisation = res.Regularisation;
params.method.LinearSolver = res.LinearSolver;

if length(res.PreviouslyFoundIterations) >= res.NDeflations
    error('Please request at least one more minima than input')
end

%Line search Parameters for objective function
params.linesearch.merit.method = res.Linesearch;
params.linesearch.merit.c1 = res.c1;
params.linesearch.merit.tau = res.tau;
params.linesearch.merit.alpha0 = res.alpha0;
params.linesearch.merit.minalpha = res.minalpha;


if params.linesearch.merit.method =="Quadratic"
    if isfield(varargin{1},"tau")
        warning("The optional input Tau has no effect when using the Quadratic linesearch")
    end
end


%Line search Parameters for deflated objective function
if res.Method=="Bad_GN"
    if isfield(varargin{1},"LinearSolver")
        warning("The optional input LinearSolver has no effect when using the ""Bad"" Deflated Gauss Newton method");
    end
    if isempty(res.DeflatedLinesearch)
        params.linesearch.deflatedmerit.method = res.Linesearch;
    else
        params.linesearch.deflatedmerit.method = res.DeflatedLinesearch;
        if params.linesearch.deflatedmerit.method =="Quadratic"
            if ~isempty(res.Deflatedtau)
                warning("The optional input DeflatedTau has no effect when using the deflated Quadratic linesearch")
            end
        end
    end
    for i = ["c1","tau","alpha0","minalpha"]
        if isempty(res.(strcat("Deflated",i)))
            params.linesearch.deflatedmerit.(i) = res.(i);
        else
            params.linesearch.deflatedmerit.(i) = res.(strcat("Deflated",i));
        end
    end
elseif (res.DeflatedLinesearch)~="No"
    warning("Deflated linesearch is only applied to the ""Bad"" deflated Gauss-Newton method")
end

%parameters for line search on Mu, not all adaptable currently
params.linesearch.Mu.method = res.MuLinesearch;
params.linesearch.Mu.c1 = res.c1;
params.linesearch.Mu.tau = res.tau;
params.linesearch.Mu.alpha0 = res.Mualpha0;
params.linesearch.Mu.minalpha = res.minalpha;



%Deflation Parameters
params.deflation.theta = res.theta;
params.deflation.sigma = res.sigma;
params.deflation.singleshift = res.SingleShift;
params.deflation.epsilon = res.epsilon;
params.deflation.NormWeighting = res.NormWeighting;

%Convergence Parameters
params.convergence.ObjectiveTolerance = res.ObjectiveTolerance;
params.convergence.GradientTolerance = res.GradientTolerance;
params.convergence.StepTolerance = res.StepTolerance;
params.convergence.MaximumIterations = res.MaxIter;


% Set line search objective/merit function and derivatives - phi/gradphi
%First set deflation linesearch
params.linesearch.Mu.phi = @(~,x,constants,DeflatedPts,DeflationParameters) deflation(DeflatedPts,x,DeflationParameters);
params.linesearch.Mu.gradphi = @(X) X.gradMu;

%
% Set objective merit functions and optimisation methods
%
switch res.Method
    case "Newton"
        params.linesearch.merit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters);
        params.linesearch.merit.gradphi = @(X) (X.J'*X.J+X.S)*(X.J'*X.R);

        if isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"

            params.linesearch.deflatedmerit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) Deflated_Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters);
            params.linesearch.deflatedmerit.gradphi = @(X) 2*X.Mu*X.gradMu'*(X.J'*X.R)'*(X.J'*X.R)+ 2*X.Mu^2*(X.J'*X.J+X.S)*(X.J'*X.R);

        end
        Fun  = str2func('Deflated_Newton');

    case "Good_GN"
        if isfield(params.linesearch,"deflatedmerit")&&params.linesearch.deflatedmerit.method ~= "No"
            warning("Note there is no deflated line search for any method other than for the '''Bad''' deflated Gauss-Newton method.")
        end
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        Fun = str2func('Good_Deflated_GaussNewton');

    case "Bad_GN"
        params.linesearch.merit.phi = @(objective_function,x,constants,~,~) objective_function(x,constants);
        params.linesearch.merit.gradphi = @(X) 2*(X.J'*X.R);
        if params.linesearch.deflatedmerit.method ~= "No"
            params.linesearch.deflatedmerit.phi = @(objective_function,x,constants,DeflatedPts,DeflationParameters) deflation(DeflatedPts,x,DeflationParameters)^2*objective_function(x,constants);
            params.linesearch.deflatedmerit.gradphi = @(X) 2*X.Mu^2*(X.J'*X.R)+2*X.Mu*X.R'*X.R;
        end
        Fun = str2func('Bad_Deflated_GaussNewton');
    otherwise
        error('Please select a valid method - Newton/Good_GN/Bad_GN')
end


%
% Process previously found iterations
%
if isempty(res.PreviouslyFoundIterations)
    notfirstiteration = 0;

    Iterations = struct("DeflatedPoint",[],"ErrorAtDeflatedPoint",[],"NIter",[],"ConvergenceFlag",[],"FuncCount",[]);
    if res.RecordIterates == true
        Iterations.Iterates=[];
    end


else
    if res.RecordTimes == true
        % Iterations.times = [];
        Times = res.PreviouslyFoundIterations.Times;
        res.PreviouslyFoundIterations = rmfield(res.PreviouslyFoundIterations,'Times');
    end
    notfirstiteration = 1;
    if (res.RecordIterates == false && ~isfield(res.PreviouslyFoundIterations,'Iterates')) || ...
            res.RecordIterates == true && isfield(res.PreviouslyFoundIterations,'Iterates')%||...

        Iterations = res.PreviouslyFoundIterations;
    else
        error("Please input previously found Iterations with or without Iterates and Times to match the current options")
    end

end


% Main Loop
j = 0;
if res.RecordTimes
    tic
end
for  i=length(Iterations)+notfirstiteration:res.NDeflations
    if ~isempty([Iterations.DeflatedPoint])&&any(all(abs([Iterations.DeflatedPoint] - x0)<1e-16))
        warning('The intial vector is a deflated point, stopping method...')
        break
    end

    % Main optimisation algorithm
    [Iterations(i)] = Fun(obj_fun,x0,[Iterations.DeflatedPoint],params,res.RecordIterates);

    %Record Time to compute if requested
    if res.RecordTimes
        Times{:,:,i} = toc;
        tic
    end

    % Check for NaNs
    if Iterations(i).ConvergenceFlag == "NaN"
        warning("Method diverging to NaN values, stopped deflations early.")
        return
    end
    % Print desired info
    if any(contains(res.ConvergenceFlags, Iterations(i).ConvergenceFlag))
        j = j+1;
        if params.method.Verbose

        end
        OutputNumMinimaFound(j,i-length(res.PreviouslyFoundIterations))
    else
        if params.method.Verbose
        end
        OutputNumMinimaFound(j,i-length(res.PreviouslyFoundIterations))
    end
    if (i-length(res.PreviouslyFoundIterations)-1-j)>res.MaxNonMinima
        warning("Stopping deflations as the maximum number of non minima have been deflated")
        break
    end
end
if res.RecordTimes
    Iterations = cell2struct([struct2cell(Iterations);Times],[fieldnames(Iterations);'Times']);
end
end

% Auxiliary functions

function f_out = Deflated_Gradient(objective_function,x,constants,DeflatedPts,DeflationParameters)
[~,R,J] = objective_function(x,constants);
f_out = deflation(DeflatedPts,x,DeflationParameters)^2*(J'*R)'*(J'*R);
end
function f_out = Gradient(objective_function,x,constants,~,~)
[~,R,J] = objective_function(x,constants);
f_out = (J'*R)'*(J'*R);
end

function OutputNumMinimaFound(NumMinFound,OuterIterations)
Deflations = OuterIterations-1;
if Deflations ==1
    if NumMinFound ==1
        disp('Found 1 local minimum in 1 deflation.')
    else
        disp(['Found ',num2str(NumMinFound),' local minima after 1 deflation.'])
    end
else
    if NumMinFound ==1
        disp(['Found 1 local minimum after ',num2str(Deflations),' deflations.'])
    else
        disp(['Found ',num2str(NumMinFound),' local minima after ',num2str(Deflations),' deflations.'])
    end
end
end

