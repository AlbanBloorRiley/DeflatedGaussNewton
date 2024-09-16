%Figures for  paper
warning('off','MATLAB:rankDeficientMatrix');
%% Section 2 figure
clear all
obj_fun = @Himmelblau;
x0=[0;-1];
method = 'Good_GN';
Opt = struct('NDeflations',4,'Method',method,'epsilon',0);
[Iterations,options] = DMin(obj_fun,x0,Opt);
%
f = figure(1);
clf
options.xylim = 6;   options.NPoints = 100; 
options.constants = [];  options.ShowLegend = false; %options.FontSize = 8;
options.plotlines = false; options.ShowDeflations = 1:length(Iterations)-1;

subplot(2,2,1)
options.ShowDeflations = 1:1; 
options.plotlines = 2;
% PlotFContours(Iterations,options,obj_fun)
PlotBetaContours(Iterations,options,obj_fun)

title('Deflation 1')
xlabel('x_1')
ylabel('x_2')
subplot(2,2,2)
options.ShowDeflations = 1:2; 
options.plotlines = 3;
% PlotFContours(Iterations,options,obj_fun)
PlotBetaContours(Iterations,options,obj_fun)

title('Deflation 2')
xlabel('x_1')
ylabel('x_2')
subplot(2,2,3)
options.ShowDeflations = 1:3; 
options.plotlines = 4;
% PlotFContours(Iterations,options,obj_fun)
PlotBetaContours(Iterations,options,obj_fun)
title('Deflation 3')
xlabel('x_1')
ylabel('x_2')


subplot(2,2,4)
options.ShowDeflations=[];options.xylim = 1e-20;
options.ShowLegend = true;options.FontSize = 8;
options.ShowDeflations = [];
Iterationss = Iterations(1);
Iterationss.ConvergenceFlag = "";
PlotBetaContours(Iterations,options,obj_fun)
axis off
xlim([100,101])
ylim([100,101])
legend('location','northeastoutside')
ah1 = axes('position',get(gca,'position'),'visible','off');
options.ShowDeflations = 1:4;
PlotxConvergence(Iterations,options)
axis on
xlim([0,20])
ylim([1e-14,1e3])
grid on
xlabel('k')
ylabel('error')
legend('location','eastoutside')

f.Units = 'centimeters';
% f.Position = [-50 10 20 14];
linestyleorder('default')
%
% print(f, 'sec2fig.eps', '-depsc')


%% Section 3.3 Figure
clear all
clf
obj_fun = @Himmelblau;
x0=[-1;1];
epsilon = 0;  method = 'Good_GN';
Opt = struct('NDeflations',4,'Method',method,'epsilon',epsilon);
[Iterations,options] = DMin(obj_fun,x0,Opt);
options.xylim =6;   options.NPoints = 100; 
options.constants = [];  options.ShowLegend = false; options.FontSize = 8;
options.plotlines = false; options.ShowDeflations = 1:length(Iterations)-1; 
f = figure(1);
subplot(1,2,1)
options.ShowDeflations = 1:length(Iterations)-1; PlotBetaContours(Iterations,options,obj_fun)
title("\epsilon = 0")
xlabel('x_1')
ylabel('x_2')
%
subplot(1,2,2)
options.epsilon = 0.4; options.ShowLegend = true;
options.ShowDeflations = 1:3; PlotBetaContours(Iterations,options,obj_fun)
title("\epsilon = 0.4")
xlabel('x_1')
ylabel('x_2')
legend('location',"southeast")
% f.Units = 'centimeters';
% f.Position = [-50 10 20 8];
% print(f, 'sec3fig.eps', '-depsc')

%% Section 4.2 figure
clear all
obj_fun = @FTrig;
x0=[1;3];

Opt = struct('NDeflations',42,'Method','Good_GN','epsilon',0.01);
[GoodIterations,options] = DMin(obj_fun,x0,Opt);
% Opt = struct('NDeflations',42,'Method','Bad_GN');
% [BadIterations] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',143,'Method','Newton','Regularisation',1e-4,'epsilon',0.01);
[NewtonIterations] = DMin(obj_fun,x0,Opt);
clf
f = figure(1);
options.xylim =10;   options.NPoints = 50; 
options.constants = [];  options.ShowLegend = true; options.FontSize = 7;
options.plotlines = false; options.ShowDeflations = 1:length(GoodIterations)-1;
subplot(1,2,1) 
PlotBetaContours(GoodIterations,options,obj_fun)
legend('location','SE')
title("Gauss-Newton")
xlabel('x_1')
ylabel('x_2')
subplot(1,2,2)
options.Method = 'Newton';
options.ShowLegend = false; 
options.ShowDeflations = 1:length(NewtonIterations)-1; PlotBetaContours(NewtonIterations,options,obj_fun)
title("Newton")
xlabel('x_1')
ylabel('x_2')
f.Units = 'centimeters';
% f.Position = [-50 10 20 8];
% print(f, 'sec4Ftrigfig.eps', '-depsc')
%%
 clear all
 x0=[1;3];
obj_fun1 = @FTrig;
Opt = struct('NDeflations',42,'Method','Good_GN','Linesearch','No','RecordTimes',true);
tic
[GoodIterations,~] = DMin(obj_fun1,x0,Opt);

Opt = struct('NDeflations',42,'Method','Bad_GN','Linesearch','No','RecordTimes',true);
[BadIterations] = DMin(obj_fun1,x0,Opt);

global times
constants=[];xy = 10;
obj_fun = @(x)fun(obj_fun1,x,constants);
opts = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true);
problem = createOptimProblem('lsqnonlin','x0',x0,'objective',obj_fun , ...
    'lb',[-xy ;-xy],'ub',[xy;xy],'options',opts);
ms = MultiStart('StartPointsToRun','bounds','Display','final', 'OutputFcn',{@tictoc,@savelocalminima});
    [~,~,~,~,~] = run(ms,problem,300);
    localSolTableJ = localSolTable;
[~,IDJ] = uniquetol(localSolTableJ.X,1e-3,'ByRows',true);
timesJ = times;
times = [];

opts = optimoptions('lsqnonlin','SpecifyObjectiveGradient',false);
problem = createOptimProblem('lsqnonlin','x0',x0,'objective',obj_fun , ...
    'lb',[-xy ;-xy],'ub',[xy;xy],'options',opts);
ms = MultiStart('StartPointsToRun','bounds','Display','final', 'OutputFcn',{@tictoc,@savelocalminima});
    [~,~,~,~,~] = run(ms,problem,300);
    localSolTableNoJ = localSolTable;
[~,IDnoJ] = uniquetol(localSolTable.X,1e-3,'ByRows',true);
timesnoJ = times;
times = [];
xy = 20; 
problem = createOptimProblem('lsqnonlin','x0',x0,'objective',obj_fun , ...
    'lb',[-xy ;-xy],'ub',[xy;xy],'options',opts);
ms = MultiStart('StartPointsToRun','bounds','Display','final', 'OutputFcn',{@tictoc,@savelocalminima});
    [~,~,~,~,~] = run(ms,problem,300);
    localSolTableJ20 = localSolTable;
[~,IDJ20] = uniquetol(localSolTableJ20.X,1e-3,'ByRows',true);
timesJ20 = times;
times = [];
opts = optimoptions('lsqnonlin','SpecifyObjectiveGradient',false);
problem = createOptimProblem('lsqnonlin','x0',x0,'objective',obj_fun , ...
    'lb',[-xy ;-xy],'ub',[xy;xy],'options',opts);
ms = MultiStart('StartPointsToRun','bounds','Display','final', 'OutputFcn',{@tictoc,@savelocalminima});
    [~,~,~,~,~] = run(ms,problem,300);
    localSolTableNoJ20 = localSolTable;
[~,IDnoJ20] = uniquetol(localSolTableJ20.X,1e-3,'ByRows',true);
timesnoJ20 = times;



%

f = figure(1);
clf
subplot(1,2,1)
plot(localSolTableNoJ.FuncCount(sort(IDnoJ)),'linewidth',1)
hold on
plot(localSolTableJ.FuncCount(sort(IDJ)),'linewidth',1)
plot(localSolTableNoJ20.FuncCount(sort(IDnoJ20)),'linewidth',1)
plot(localSolTableJ20.FuncCount(sort(IDJ20)),'linewidth',1)
plot(cumsum([GoodIterations.FuncCount]),'linewidth',1)
plot(cumsum([BadIterations.FuncCount]),'linewidth',1)
hold off
legend('MultiStart 10x10','MultiStart with Jacobian 10x10','MultiStart 20x20','MultiStart with Jacobian 20x20','Good Gauss-Newton','Bad Gauss-Newton','location','NW')
xlabel('Minima')
xlim([0,42])
ylabel('Function Evaluations')
%
linesnotmarkers = false;      %Change to switch between line and marker styles
if linesnotmarkers
    linestyleorder('Mixedstyles');
else
    linestyleorder('MixedMarkers');
    setMarkerNumber(f.Children(2),10)
end

subplot(1,2,2)
temp = cumsum(timesnoJ);
plot(temp(sort(IDnoJ)),'linewidth',1)
hold on
temp = cumsum(timesJ);
plot(temp(sort(IDJ)),'linewidth',1)
temp = cumsum(timesnoJ20);
plot(temp(sort(IDnoJ20)),'linewidth',1)
temp = cumsum(timesJ20);
plot(temp(sort(IDJ20)),'linewidth',1)

plot(cumsum([GoodIterations.Times]),'linewidth',1)
plot(cumsum([BadIterations.Times]),'linewidth',1)
hold off
setMarkerNumber(f.Children(1),10)
xlabel('Minima')
xlim([0,42])
ylabel('Time in seconds')


f.Units = 'centimeters';
% f.Position = [-50 20 30 9];
% print(f, 'sec4MultiStartfig.eps', '-depsc')


%% Section 4.3.1 figure comparison
clear all
obj_fun = @FTrig;
x0=[1;3];
Opt = struct('NDeflations',42,'Method','Good_GN');
[GoodIterations,options] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',42,'Method','Bad_GN');
[BadIterations] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',143,'Method','Newton','Regularisation',1e-4,'epsilon',0.1);
[NewtonIterations] = DMin(obj_fun,x0,Opt);
%
clf
f = figure(1);

subplot(3,1,1)
options.ShowLegend = false; options.ShowDeflations = 1:length(GoodIterations);
PlotxConvergence(GoodIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])
rearangelegend(GoodIterations,options)

subplot(3,1,2)
options.ShowDeflations = 1:length(BadIterations);
PlotxConvergence(BadIterations,options)
rearangelegend(BadIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])
subplot(3,1,3)

options.ShowDeflations = 1:length(NewtonIterations);
PlotxConvergence(NewtonIterations,options)
rearangelegend(NewtonIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])

linestyleorder('default')
f.Units = 'centimeters';
% f.Position = [-50 10 20 13];
% print(f, 'sec4FtrigfigCompare.eps', '-depsc')

%% Section 4.3.1 figure Bratu
clear all

n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = exp(1i*constants.xi'*(-n:n)*2*pi/T);
constants.D2 = spdiags(-((-n:n)'*2*pi/T).^2,0,2*n+1,2*n+1);
constants.e0 = ones(1,2*n+1);
constants.et = exp(1i*t*(-n:n)*2*pi/T);
obj_fun=@(x,~)EvaluateBratu(x,constants);
evalfun = @(c,x) exp(1i*x*(-n:n)*2*pi/T)*c;

x0 = zeros(2*n+1,1);

Method = 'Good_GN';
NDeflations = 2;

Opt = struct('NDeflations',NDeflations,'Method',Method,...
    'MaxIter',400,'NormWeighting',constants.A,'ObjectiveTolerance',1e-12,...
     'LinearSolver', 'lsqminnorm','Linesearch','Quadratic');
Iterations = DMin(obj_fun,x0,Opt);
%
clf
f = figure(1);
subplot(2,1,1)
options.ShowNonMinima = false;  options.ShowLegend = true;
PlotFE(xi,Iterations,obj_fun,evalfun,constants,options)
xlabel('x')
ylabel('u(x)')
ylim([0,3])
subplot(2,1,2)
options.ShowDeflations = 1:length(Iterations); options.ShowLegend = true;
PlotFConvergence(Iterations,options,obj_fun)
grid on
xlabel('k')
ylabel('residual error')
ylim([1e-15,1e10])
subplot(2,1,1)
rearangelegend(Iterations,options)

subplot(2,1,2)
legend("off")

linesnotmarkers = true;       %Change to switch between line and marker styles
if linesnotmarkers
    linestyleorder('Mixedstyles');
else
    linestyleorder('MixedMarkers');
    setMarkerNumber(f.Children(1),100)
    setMarkerNumber(f.Children(3),21)
end
f.Units = 'centimeters';
% f.Position = [-50 10 20 10];
% print(f, 'sec4Bratufig.eps', '-depsc')

%% Section 4.3.2 Figure
clear all
n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = exp(1i*constants.xi'*(-n:n)*2*pi/T);
constants.D2 = spdiags(-((-n:n)'*2*pi/T).^2,0,2*n+1,2*n+1);
constants.e0 = ones(1,2*n+1);
constants.et = exp(1i*t*(-n:n)*2*pi/T);
obj_fun=@(x,~)EvaluateCarrier(x,constants);
evalfun = @(c,x) exp(1i*x*(-n:n)*2*pi/T)*c;

x0 = zeros(2*n+1,1); 

Method = 'Good_GN';
NDeflations = 5;

Opt = struct('NDeflations',NDeflations,'Method',Method,'theta',2,...
    'MaxIter',200,'NormWeighting',constants.A,'ObjectiveTolerance',1e-10,...
    'LinearSolver','lsqminnorm','Linesearch','Quadratic');
[Iterations,options] = DMin(obj_fun,x0,Opt);
%
clf
f = figure(1);
subplot(2,1,1)
options.ShowNonMinima = true;  options.ShowLegend = false;
PlotFE(xi,Iterations,obj_fun,evalfun,constants,options)
xlabel('x')
ylabel('u(x)')
ylim([-3,5])
subplot(2,1,2)
options.ShowDeflations = 1:length(Iterations); options.ShowLegend = true;
PlotFConvergence(Iterations,options,obj_fun)

grid on
xlabel('k')
ylabel('residual error')
% xlim([0,70])
% ylim([1e-15,1e10])

subplot(2,1,1)
rearangelegend(Iterations,options)
% legend(lgnd)
% legend('location','westoutside')
subplot(2,1,2)
legend("off")
linesnotmarkers = false;       %Change to switch between line and marker styles
if linesnotmarkers
    linestyleorder('Mixedstyles');
else
    linestyleorder('MixedMarkers');
    setMarkerNumber(f.Children(1),100)
    setMarkerNumber(f.Children(3),21)
end
f.Units = 'centimeters';
% f.Position = [-50 10 20 10];
% print(f, 'sec4Carrierfig.eps', '-depsc')


%% Section 4.3.2 comparison figure part 1

clear all
n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = exp(1i*constants.xi'*(-n:n)*2*pi/T);
constants.D2 = spdiags(-((-n:n)'*2*pi/T).^2,0,2*n+1,2*n+1);
constants.e0 = ones(1,2*n+1);
constants.et = exp(1i*t*(-n:n)*2*pi/T);
obj_fun=@(x,~)EvaluateCarrier(x,constants);
evalfun = @(c,x) exp(1i*x*(-n:n)*2*pi/T)*c;

x0 = zeros(2*n+1,1);


NDeflations = 5;

Opt = struct('NDeflations',NDeflations,...
    'MaxIter',400,'NormWeighting',constants.A,'ObjectiveTolerance',1e-10,...
    'LinearSolver','lsqminnorm','linesearch','Quadratic');

Opt.Method = "Good_GN";
GoodIterations = DMin(obj_fun,x0,Opt);

vals = xi'.*(xi'-1); 
betterx0 = constants.A\vals;

Opt.Method = "Good_GN";
BetterIterations = DMin(obj_fun,betterx0,Opt);

Opt.Method = "Bad_GN";
BadIterations = DMin(obj_fun,x0,Opt);

Opt.Method = "Bad_GN";
BadBetterIterations = DMin(obj_fun,betterx0,Opt);

%% Section 4.3.2 comparison figure part 2

clf

options.ShowLegend = true;
linesnotmarkers = true;
if options.ShowLegend
    if linesnotmarkers
        linestyleorder('Mixedstyles');
    else
        linestyleorder('MixedMarkers');
    end
    loc = "eastoutside";
else
    linestyleorder('default');
end


f = figure(1);
subplot(4,1,1)

options.ShowDeflations = 1:length(GoodIterations);
options.ShowNonMinima = false;
options.fontsize = 7;
PlotFConvergence(GoodIterations,options,obj_fun)
grid on
ylabel('residual error')
xlim([0,110])
ylim([1e-15,1e10])
setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(GoodIterations,options,loc);end

subplot(4,1,2)

options.ShowDeflations = 1:length(BadIterations); 
PlotFConvergence(BadIterations,options,obj_fun)
grid on
ylabel('residual error')
xlim([0,110])
ylim([1e-15,1e10])
setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BadIterations,options,loc);end

subplot(4,1,3)

options.ShowDeflations = 1:length(BetterIterations); 
PlotFConvergence(BetterIterations,options,obj_fun)
grid on
xlabel('k')
ylabel('residual error')
xlim([0,110])
ylim([1e-15,1e10])
setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BetterIterations,options,loc);end

subplot(4,1,4)

options.ShowDeflations = 1:length(BadBetterIterations);
PlotFConvergence(BadBetterIterations,options,obj_fun)
grid on
xlabel('k')
ylabel('residual error')
xlim([0,110])
ylim([1e-15,1e10])
setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BadBetterIterations,options,loc);end

f.Units = 'centimeters';
% f.Position = [0 0 20 18];
% print(f, 'sec4CarrierfigCompareLines.eps', '-depsc')


%%
clear Sys Exp
[Sys1,Exp] = Mn12_Spin_Sys_3(1,1);
Sys = Sys1;
Vary=Sys;

Opt = struct('NDeflations',4,'Method','Good_GN','Linesearch','Quadratic',...
    'IEPType','Classic','Verbose',false,'scaled',false,'c1',1e-16);
[SysOutGood]= INS_IEP(Sys,Vary,Exp,Opt);


Opt = struct('NDeflations',4,'Method','Bad_GN','Linesearch','No',...
    'IEPType','Classic','Verbose',false,'scaled',false,'c1',1e-9,'theta',2);
[SysOutBad]= INS_IEP(Sys,Vary,Exp,Opt);


Opt = struct('NDeflations',4,'Method','Newton','Linesearch','Armijo',...
    'IEPType','Classic','Verbose',false,'scaled',false,'c1',1e-4);
[SysOutNewton]= INS_IEP(Sys,Vary,Exp,Opt);
%%
clf
f=figure(1);
subplot(3,1,1)
options.ShowLegend = true; options.ShowDeflations = 1:length(SysOutGood);
PlotxConvergence(SysOutGood,options)
grid on
xlabel('k')
ylabel('error')
legend('location','eastoutside')
ylim([1e-10,1e5])

subplot(3,1,2)
options.ShowDeflations = 1:length(SysOutBad);
PlotxConvergence(SysOutBad,options)
legend('location','eastoutside')
grid on
xlabel('k')
ylabel('error')
ylim([1e-10,1e5])

subplot(3,1,3)
options.ShowDeflations = 1:length(SysOutNewton);
PlotxConvergence(SysOutNewton,options)
legend('location','eastoutside')


linesnotmarkers = true;       %Change to switch between line and marker styles
if linesnotmarkers
    linestyleorder("mixedstyles")
else
    linestyleorder("mixedmarkers")
    n = 10;
    setMarkerNumber(f.Children(2),n)
    setMarkerNumber(f.Children(4),n)
    setMarkerNumber(f.Children(6),n)
end

grid on
xlabel('k')
ylabel('error')
ylim([1e-10,1e5])
f.Units = 'centimeters';
% f.Position = [10 10 20 12];
% print(f, 'sec4Mn12.eps', '-depsc')

%%


function [f,r,Jr] = EvaluateBratu(c,constants)
m = constants.m;
A = constants.A;
D2 = constants.D2;
e0 = constants.e0;
et = constants.et;

% Residual for the BVP: u'' + 3*exp(u) = 0, u(0) = 0, u(1) = 0:
r =  [A*(D2*c) + 3*exp(A*c); ...
     e0*c;
     et*c];

f =  dot(r,r)/2;

if nargout>2
    Jr =  [A*D2 + 3*spdiags(exp(A*c),0,m,m)*A; ...
        e0;
        et];
end
end

function [f,r,Jr] = EvaluateCarrier(c,constants)
m = constants.m;
xi = constants.xi;
A = constants.A;
D2 = constants.D2;
e0 = constants.e0;
et = constants.et;

% Residual for the BVP: 0.05u'' + u^2 + 8x(1-x)u = 1, u(0) = 0, u(1) = 0:
r =  [0.05*A*(D2*c) + (A*c).^2 + 8*(xi.*(1-xi))'.*(A*c) - 1; ...
     e0*c;
     et*c];
f =  dot(r,r)/2;

if nargout>2
    Jr =  [0.05*A*D2 + 2*spdiags(A*c,0,m,m)*A+ 8*(xi.*(1-xi))'.*A; ...
        e0;
        et];
end
end


function ferrs = format_errs(errs)
    ferrs = string([]);
    for k = 1:length(errs)
        ferrs(k) = sprintf('%0.1e', errs(k));
    end
end

function R = Rfun(fun1,x,constants)
[~,R] = fun1(x,constants);
end
function [R,J] = fun(fun1,x,constants)
[~,R,J] = fun1(x,constants);
end
function [R,J] = Jfun(fun1,x,constants)
[~,R,J] = fun1(x,constants);
end

function rearangelegend(Iterations,options,varargin)
if isempty(varargin)
    loc = 'westoutside';
else
    loc = varargin{1};
end
hLegend = findobj(gcf, 'Type', 'Legend');
if options.ShowLegend&&options.ShowNonMinima
    lgnd = strcat(hLegend(1).String, string(newline), format_errs([Iterations.ErrorAtDeflatedPoint]));
legend(lgnd,'location',loc)
elseif options.ShowLegend
str = convertCharsToStrings({Iterations(:).ConvergenceFlag});
lgnd = strcat(hLegend(1).String, string(newline), format_errs([Iterations(str~="Max Iterations reached").ErrorAtDeflatedPoint]));
legend(lgnd,'location',loc)
end
end

function setMarkerNumber(f,n)
    for i = 1:length(f.Children)
        f.Children(i).MarkerIndices = 1:(floor(length(f.Children(i).MarkerIndices)/n)):f.Children(i).MarkerIndices(end);
    end
end

function stop = tictoc(~, state)
global times
stop = false;
switch state
    case 'init'
        tic
    case 'iter'
        times(end+1) = toc;
        tic
end
end



function stop = savelocalminima(optimValues,state)
% Adapted from savelocalsolutions (Copyright 2023 The MathWorks, Inc.)
stop = false;
switch state
  case 'init'
    assignin('base','localSolTable',[])
  case 'iter'
    if ~isempty(optimValues.localsolution.X)
        t = table(...
            optimValues.localsolution.X(:)', ...
            optimValues.funccount,...
            optimValues.localsolution.Fval, ...
            optimValues.localsolution.Exitflag, ...
            'VariableNames', ["X", "FuncCount","fval", "exitflag"]);
        if isfield(optimValues.localsolution, "Constraintviolation") && ...
                ~isempty(optimValues.localsolution.Constraintviolation)
            tCon = table(...
                optimValues.localsolution.Constraintviolation, ...
                'VariableNames', "constrviolation");
            t = [t, tCon];
        end
      sols = evalin('base','localSolTable');
      sols = [sols; t];
      assignin('base','localSolTable',sols)      
    end
  case 'done'   
end
end



function [F,R,J,H] = FTrig(X,varargin)
x = X(1); y = X(2);
    R = [-10*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1);
        -10*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1);
        x^2/100 + y^2/100 + 10];
    F = (sum(R.^2));
    if nargout>2
        J = [- 10*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 10*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 10*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 10*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015), - 10*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 10*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 10*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 10*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015);
            - 10*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) - 10*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1) - 10*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1),                                                                                                                                                                             10*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) + 10*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1) + 10*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1);
            x/50,y/50];
    end
    if nargout>3
        H(:,:,1) =[- (1125899906842624*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (281474976710656*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (703687441776640*(x + y)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/3125302502557517 - 20*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015) - 20*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015), - (1125899906842624*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (281474976710656*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (703687441776640*(x + y)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/3125302502557517 - 20*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015) - 20*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015);
            - (1125899906842624*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (281474976710656*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (703687441776640*(x + y)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/3125302502557517 - 20*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015) - 20*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015), - (1125899906842624*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (281474976710656*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/555609333788003 - (703687441776640*(x + y)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1))/3125302502557517 - 20*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015) - 20*(x + y)*((35184372088832*(x + y)^2)/3125302502557517 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((70368744177664*(x + y)^2)/2778046668940015 - 1)*((562949953421312*x)/2778046668940015 + (562949953421312*y)/2778046668940015) - 20*(x + y)*((70368744177664*x)/3125302502557517 + (70368744177664*y)/3125302502557517)*((281474976710656*(x + y)^2)/2778046668940015 - 1)*((140737488355328*x)/2778046668940015 + (140737488355328*y)/2778046668940015)];
        H(:,:,2) = [- (4503599627370496*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/555609333788003 - (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/8681395840437547 - 1))/3125302502557517 - (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/8681395840437547 - 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517) - 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) - 20*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015),   (4503599627370496*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/555609333788003 + (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/8681395840437547 - 1))/3125302502557517 + (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/8681395840437547 + 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517) + 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) + 20*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015);
            (4503599627370496*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/555609333788003 + (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/8681395840437547 - 1))/3125302502557517 + (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/8681395840437547 + 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517) + 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) + 20*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015), - (4503599627370496*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/555609333788003 - (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/8681395840437547 - 1))/3125302502557517 - (2814749767106560*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((140737488355328*(x - y)^2)/3125302502557517 - 1))/8681395840437547 - 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((1125899906842624*(x - y)^2)/2778046668940015 - 1)*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517) - 20*((281474976710656*x)/8681395840437547 - (281474976710656*y)/8681395840437547)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)*((140737488355328*(x - y)^2)/3125302502557517 - 1) - 20*((281474976710656*x)/3125302502557517 - (281474976710656*y)/3125302502557517)*((140737488355328*(x - y)^2)/8681395840437547 - 1)*((2251799813685248*x)/2778046668940015 - (2251799813685248*y)/2778046668940015)];
        H(:,:,3) = [1/50,    0; 0, 1/50];
    end
    if nargin>4
        error("Too many outputs")
    end


end

function [F,R,J,H] = Himmelblau(X,varargin)
x = X(1);   y = X(2);
R(1,1) =  x^2+y -11;
R(2,1) = x+y^2 -7;
F = (sum(R.^2));
if nargout>2
    J(1,1) =  2*x;
    J(2,1) = 1;
    J(1,2) = 1;
    J(2,2) = 2*y;
end
if nargout>3
    H(1,1,1) = 2;
    H(2,2,1) = 0;
    H(1,2,1) = 0;
    H(2,1,1) =  H(1,2,1);
    H(1,1,2) = 0;
    H(2,2,2) = 2;
    H(1,2,2) = 0;
    H(2,1,2) =  H(1,2,2);
end
if nargin>4
    error("Too many outputs")
end
end


