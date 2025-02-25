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

% options.edgecolour = 'white';                 %Changes contour outlines

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
f.Position = [-50 10 20 14];
linestyleorder('mixedstyles')
%
print(f, 'sec2fig.eps', '-depsc')


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
f.Units = 'centimeters';
f.Position = [-50 10 20 8];
print(f, 'sec3fig.eps', '-depsc')


%% Section 3.3(?) figure comparing convergence rates with different epsilon

clear all
obj_fun = @Himmelblau;
obj_fun = @FTrig;
x0=[0;-1];
method = 'Good_GN';
NDeflations = 6;
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0);
[Iterations0,options0] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0.01);
[Iterations01,options01] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0.1);
[Iterations1,options1] = DMin(obj_fun,x0,Opt);
%
f = figure(1);
clf
N = 5;
options.ShowDeflations = 1:N; options.ShowLegend = false;

YLim = [1e-12,1e2];
XLim = [0,8];
Dashed = false;
subplot(2,3,1)
title('\epsilon = 0')
% PlotxConvergence(Iterations0,options)
options0.ShowDeflations = 1:N; options0.ShowLegend = false; options0.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations0,obj_fun,options0)


axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')
% legend('location','eastoutside')

subplot(2,3,2)
options.ShowLegend = true;
title('\epsilon = 0.01')
% PlotxConvergence(Iterations01,options)
options01.ShowDeflations = 1:N; options01.ShowLegend = false;     options01.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations01,obj_fun,options01)

axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')

subplot(2,3,3)
options.ShowLegend = false;
title('\epsilon = 0.1')
% PlotxConvergence(Iterations1,options)
options1.ShowDeflations = 1:N; options1.ShowLegend = false; options1.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations1,obj_fun,options1)

axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')





method = 'Bad_GN';
NDeflations = 6;
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0);
[Iterations0,options0] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0.01);
[Iterations01,options01] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',NDeflations,'Method',method,'epsilon',0.1);
[Iterations1,options1] = DMin(obj_fun,x0,Opt);
YLim = [1e-12,1e2];
XLim = [0,14];
Dashed = false;
subplot(2,3,4)
title('\epsilon = 0')
% PlotxConvergence(Iterations0,options)
options0.ShowDeflations = 1:N; options0.ShowLegend = true; options0.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations0,obj_fun,options0)
legend('location','southwest')
legend('fontsize',7)
axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')
% legend('location','eastoutside')

subplot(2,3,5)
options.ShowLegend = true;
title('\epsilon = 0.01')
% PlotxConvergence(Iterations01,options)
options01.ShowDeflations = 1:N; options01.ShowLegend = false;     options01.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations01,obj_fun,options01)

axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')
% legend('location','south')

subplot(2,3,6)
options.ShowLegend = false;
title('\epsilon = 0.1')
% PlotxConvergence(Iterations1,options)
options1.ShowDeflations = 1:N; options1.ShowLegend = false; options1.Dashed = Dashed;
PlotDeflatedxConvergence(Iterations1,obj_fun,options1)

axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('error')


f.Units = 'centimeters';
f.Position = [-50 10 20 15];
linestyleorder('default')
print(f, 'sec2figEpsilonComparison.eps', '-depsc')

%% Figure to show deflated steps in convergence Good GN
clear all
method = 'Good_GN';
obj_fun = @FTrig;
x0=[1;3];

options.ShowLegend = false;

Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-10,'epsilon',0.0);
[Iterations0,options0] = DMin(obj_fun,x0,Opt);
options0.ShowDeflations = 1:7; options0.ShowLegend = options.ShowLegend ;
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-10,'epsilon',0.01);
[Iterations01,options01] = DMin(obj_fun,x0,Opt);
options01.ShowDeflations = 1:7; options01.ShowLegend = options.ShowLegend ;
%
f = figure(1);
clf

subplot(1,2,1)
title('\epsilon = 0')
PlotDeflatedJFConvergence(Iterations0,obj_fun,options0)
axis on
% xlim([0,15])
% ylim([1e2,120])
grid on
xlabel('k')
ylabel('error')
% legend('location','eastoutside')

subplot(1,2,2)
title('\epsilon = 0.001')
PlotDeflatedJFConvergence(Iterations01,obj_fun,options01)
axis on
% xlim([0,15])
% ylim([1e2,120])
grid on
xlabel('k')
ylabel('error')
% legend('location','south')
linestyleorder('default')

f.Units = 'centimeters';
f.Position = [-50 10 20 8];
%
print(f, 'figShowingDeflatedStepsGoodGN.eps', '-depsc')
%% Figure to show deflated steps and comparison of epsilon convergence for good and bad GN 

clear all
method = 'Good_GN';
obj_fun = @FTrig;
x0=[8;9];
x0=[2;4];

Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.0);
[Iterations0,options0] = DMin(obj_fun,x0,Opt);
options0.ShowDeflations = 1:7; options0.ShowLegend = false ;
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.01);
[Iterations01,options01] = DMin(obj_fun,x0,Opt);
options01.ShowDeflations = 1:7; options01.ShowLegend = false;
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.1);
[Iterations1,options1] = DMin(obj_fun,x0,Opt);
options1.ShowDeflations = 1:7; options1.ShowLegend = false;
%
f = figure(1);
clf
YLim = [1e-15,1e10];
subplot(2,3,1)
title('\epsilon = 0')

PlotDeflatedJFConvergence(Iterations0,obj_fun,options0)
axis on
% xlim([0,15])
ylim(YLim)
grid on
xlabel('k')
ylabel('||\nabla f(x^k)||')

subplot(2,3,2)
title('\epsilon = 0.001')
PlotDeflatedJFConvergence(Iterations01,obj_fun,options01)
axis on
% xlim([0,15])
ylim(YLim)
grid on
xlabel('k')
% legend('location','south')
subplot(2,3,3)
title('\epsilon = 0.01')
PlotDeflatedJFConvergence(Iterations1,obj_fun,options1)
axis on
% xlim([0,15])
ylim(YLim)
grid on
xlabel('k')


method = 'Bad_GN';
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.0);
[Iterations0,options0] = DMin(obj_fun,x0,Opt);
options0.ShowDeflations = 1:7; options0.ShowLegend = true ;
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.01);
[Iterations01,options01] = DMin(obj_fun,x0,Opt);
options01.ShowDeflations = 1:7; options01.ShowLegend = false;
Opt = struct('NDeflations',10,'Method',method,'StepTolerance',1e-9,'epsilon',0.1);
[Iterations1,options1] = DMin(obj_fun,x0,Opt);
options1.ShowDeflations = 1:7; options1.ShowLegend = false;
%



subplot(2,3,4)
title('\epsilon = 0')
XLim = [0,30];
% YLim = [1e-15,1e15];
PlotDeflatedJFConvergence(Iterations0,obj_fun,options0)
axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')
ylabel('||\nabla f(x^k)||')
legend('location','southeast','FontSize',7)


subplot(2,3,5)
title('\epsilon = 0.001')
PlotDeflatedJFConvergence(Iterations01,obj_fun,options01)
axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')

%
subplot(2,3,6)
title('\epsilon = 0.01')
PlotDeflatedJFConvergence(Iterations1,obj_fun,options1)
axis on
xlim(XLim)
ylim(YLim)
grid on
xlabel('k')

f.Units = 'centimeters';
f.Position = [-50 10 25 15];
%
print(f, 'figShowingDeflatedStepsWithEpsilon.eps', '-depsc')




%% Section ManyMinima Deflation vs Multistart
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
f.Position = [-50 10 20 8];
print(f, 'sec4Ftrigfig.eps', '-depsc')
%%
 clear all
 x0=[1;3];
obj_fun1 = @FTrig;
Opt = struct('NDeflations',42,'Method','Good_GN','Linesearch','Armijo','RecordTimes',true);
tic
[GoodIterations,~] = DMin(obj_fun1,x0,Opt);

Opt = struct('NDeflations',42,'Method','Bad_GN','Linesearch','Armijo','RecordTimes',true);
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
%%
f = figure(1);
clf
subplot(1,2,1)
plot(0:length(localSolTableNoJ.FuncCount(sort(IDnoJ)))-1,localSolTableNoJ.FuncCount(sort(IDnoJ)),'linewidth',1)
hold on
plot(0:length(localSolTableJ.FuncCount(sort(IDJ)))-1,localSolTableJ.FuncCount(sort(IDJ)),'linewidth',1)
plot(0:length(localSolTableNoJ20.FuncCount(sort(IDnoJ20)))-1,localSolTableNoJ20.FuncCount(sort(IDnoJ20)),'linewidth',1)
plot(0:length(localSolTableJ20.FuncCount(sort(IDJ20)))-1,localSolTableJ20.FuncCount(sort(IDJ20)),'linewidth',1)
plot(0:length(cumsum([GoodIterations.FuncCount]))-1,cumsum([GoodIterations.FuncCount]),'linewidth',1)
plot(0:length(cumsum([BadIterations.FuncCount]))-1,cumsum([BadIterations.FuncCount]),'linewidth',1)
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
plot(0:length(temp(sort(IDnoJ)))-1,temp(sort(IDnoJ)),'linewidth',1)
hold on
temp = cumsum(timesJ);
plot(0:length(temp(sort(IDJ)))-1,temp(sort(IDJ)),'linewidth',1)
temp = cumsum(timesnoJ20);
plot(0:length(temp(sort(IDnoJ20)))-1,temp(sort(IDnoJ20)),'linewidth',1)
temp = cumsum(timesJ20);
plot(0:length(temp(sort(IDJ20)))-1,temp(sort(IDJ20)),'linewidth',1)

plot(0:length(cumsum([GoodIterations.Times]))-1,cumsum([GoodIterations.Times]),'linewidth',1)
plot(0:length(cumsum([BadIterations.Times]))-1,cumsum([BadIterations.Times]),'linewidth',1)
hold off
setMarkerNumber(f.Children(1),10)
xlabel('Minima')
xlim([0,42])
ylabel('Time in seconds')


f.Units = 'centimeters';
f.Position = [-50 20 30 9];
print(f, 'sec4MultiStartfig.eps', '-depsc')


%% Section 4.3.1 figure comparison
clear all
obj_fun = @FTrig;
x0=[1;3];
Opt = struct('NDeflations',42,'Method','Good_GN');
[GoodIterations,options] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',42,'Method','Bad_GN');
[BadIterations] = DMin(obj_fun,x0,Opt);
Opt = struct('NDeflations',143,'Method','Newton','Regularisation',1e-4,'epsilon',0.01);
[NewtonIterations] = DMin(obj_fun,x0,Opt);
%%
clf
f = figure(1);

subplot(3,1,1)
options.ShowLegend = false; options.ShowDeflations = 1:length(GoodIterations);
PlotxConvergence(GoodIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])
yticks([1e-10,1e-5,1,1e5])
ylim([1e-10,1e5])
rearangelegend(GoodIterations,options)

subplot(3,1,2)
options.ShowDeflations = 1:length(BadIterations);
PlotxConvergence(BadIterations,options)
rearangelegend(BadIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])
yticks([1e-10,1e-5,1,1e5])
ylim([1e-10,1e5])

subplot(3,1,3)
options.ShowDeflations = 1:length(NewtonIterations);
PlotxConvergence(NewtonIterations,options)
rearangelegend(NewtonIterations,options)
grid on
xlabel('k')
ylabel('error')
xlim([0,150])
yticks([1e-10,1e-5,1,1e5])
ylim([1e-10,1e5])


linestyleorder('default')
f.Units = 'centimeters';
f.Position = [-50 10 20 13];
print(f, 'sec4FtrigfigCompare.eps', '-depsc')

%%  Bratu
clear all

n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = (1/m+1)*exp(1i*constants.xi'*(-n:n)*2*pi/T);
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
ylabel('residual')
ylim([1e-15,1e10])
subplot(2,1,1)
rearangelegend(Iterations,options)

subplot(2,1,2)
legend("off")
% ylim([1e-20,1e10])
yticks([1e-20,1e-10,1,1e10])
linesnotmarkers = true;       %Change to switch between line and marker styles
if linesnotmarkers
    linestyleorder('Mixedstyles');
else
    linestyleorder('MixedMarkers');
    setMarkerNumber(f.Children(1),100)
    setMarkerNumber(f.Children(3),21)
end
f.Units = 'centimeters';
f.Position = [-50 10 20 10];
print(f, 'sec4Bratufig.eps', '-depsc')

%%  Carrier Equation
clear all
rng(1)
n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = (1/m+1)*exp(1i*constants.xi'*(-n:n)*2*pi/T);
constants.D2 = spdiags(-((-n:n)'*2*pi/T).^2,0,2*n+ 1,2*n+1);
constants.e0 = ones(1,2*n+1);
constants.et = exp(1i*t*(-n:n)*2*pi/T);
constants.et = (-1).^(-n:n);

% constants.epsilon =  (0.2)       %Carrier equation constant
% % constants.epsilon =  sqrt(0.2)
% constants.epsilon =  sqrt(0.05)

obj_fun=@(x,~)EvaluateCarrier(x,constants);
evalfun = @(c,x) exp(1i*x*(-n:n)*2*pi/T)*c;

x0 = zeros(2*n+1,1); 
% x0 = ones(2*n+1,1).*1e-8; 
 % x0 =@(i) randn(2*n+1,1)*1e-7;
Method = 'Good_GN';
NDeflations = 15;

Opt = struct('NDeflations',NDeflations,'Method',Method,'MaxIter',100,'NormWeighting',constants.A,'ObjectiveTolerance',1e-10,...
    'StepTolerance',1e-5,'LinearSolver','lsqminnorm','Linesearch','Quadratic'...
    ,'MaxNonMinima',NDeflations);
[Iterations,options] = DMin(obj_fun,x0,Opt);
%

f = figure(1);
clf
subplot(2,1,1)
options.ShowNonMinima = false;  options.ShowLegend = false;
PlotFE(xi,Iterations,obj_fun,evalfun,constants,options)
xlabel('x')
ylabel('u(x)')
ylim([-2.5,4.5])
subplot(2,1,2)
options.ShowDeflations = 1:length(Iterations); options.ShowLegend = true;
PlotFConvergence(Iterations,options,obj_fun)

grid on
xlabel('k')
ylabel('residual')
% xlim([0,70])
% ylim([1e-15,1e10])

subplot(2,1,1)
rearangelegend(Iterations,options)
% legend(lgnd)
% legend('location','westoutside')
subplot(2,1,2)
ylim([1e-20,1e10])
yticks([1e-20,1e-10,1,1e10])
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
f.Position = [-50 10 20 10];
print(f, 'sec4Carrierfig.eps', '-depsc')


%% Section Carrier Equation comparison

clear all
n = 100;    
m = 4*n;
t = 1;
T = 2*t;
xi = linspace(0,t,m);
constants.m = m;
constants.xi = xi;
constants.A = (1/m+1)*exp(1i*constants.xi'*(-n:n)*2*pi/T);
constants.D2 = spdiags(-((-n:n)'*2*pi/T).^2,0,2*n+1,2*n+1);
constants.e0 = ones(1,2*n+1);
constants.et = (-1).^(-n:n);
obj_fun=@(x,~)EvaluateCarrier(x,constants);
evalfun = @(c,x) exp(1i*x*(-n:n)*2*pi/T)*c;





x0 = zeros(2*n+1,1);


NDeflations = 8;


% Opt = struct('NDeflations',NDeflations,...
%     'MaxIter',200,'NormWeighting',constants.A,'ObjectiveTolerance',1e-10,...
%     'StepTolerance',1e-4,'LinearSolver','lsqminnorm','linesearch','Quadratic');

Opt = struct('NDeflations',NDeflations,'MaxIter',200,'NormWeighting',constants.A,'ObjectiveTolerance',1e-10,...
    'StepTolerance',1e-5,'LinearSolver','lsqminnorm','Linesearch','Quadratic'...
    ,'MaxNonMinima',NDeflations,'MinAlpha',1e-8);

Opt.Method = "Good_GN";
GoodIterations = DMin(obj_fun,x0,Opt);
%
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

%
f = figure(1);
clf
Xlim = [0,100];
Xlim = 'tickaligned';
subplot(4,1,1)

options.ShowDeflations = 1:length(GoodIterations);
options.ShowNonMinima = false;
options.fontsize = 7;
PlotFConvergence(GoodIterations,options,obj_fun)
grid on
ylabel('residual')
    xlim(Xlim)
ylim([1e-15,1e10])
% setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(GoodIterations,options,loc);end
ylim([1e-20,1e10])
yticks([1e-30,1e-20,1e-10,1,1e10])

subplot(4,1,2)

options.ShowDeflations = 1:length(BadIterations); 
PlotFConvergence(BadIterations,options,obj_fun)
grid on
ylabel('residual')
xlim(Xlim)
ylim([1e-15,1e10])
% setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BadIterations,options,loc);end
ylim([1e-20,1e10])
yticks([1e-30,1e-20,1e-10,1,1e10])

subplot(4,1,3)

options.ShowDeflations = 1:length(BetterIterations); 
PlotFConvergence(BetterIterations,options,obj_fun)
grid on
xlabel('k')
ylabel('residual')
xlim(Xlim)
ylim([1e-15,1e10])
% setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BetterIterations,options,loc);end
ylim([1e-20,1e10])
yticks([1e-30,1e-20,1e-10,1,1e10])

subplot(4,1,4)
%
options.ShowDeflations = 1:length(BadBetterIterations);
PlotFConvergence(BadBetterIterations,options,obj_fun)
grid on
xlabel('k')
ylabel('residual')
xlim(Xlim)
ylim([1e-15,1e10])
% setMarkerNumber(f.Children(2),50)
if options.ShowLegend;rearangelegend(BadBetterIterations,options,loc);end
ylim([1e-20,1e10])
yticks([1e-30,1e-20,1e-10,1,1e10])


f.Units = 'centimeters';
f.Position = [-50 10 20 18];
print(f, 'sec4CarrierfigCompareLines.eps', '-depsc')


%% Mn12 example
clear all
rcm = 29979.2458;    % reciprocal cm to MHz
meV = rcm*8.065;  
B20 = -0.0570*meV/3; %(D = 3*B02)
B40 = (-2.78*10^-6)*meV;
B44 = (-3.2*10^-6)*meV;
B22 = (6.8*10^-4)*meV;
x0 = round([B20;B40;B44;B22;],1,'significant');
Sys.S = 10; 
Sys.B2 = [B22 0 B20 0 0];        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0];  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
H = ham(Sys,[0,0,0]);  [~,E]=eig(H);
EE = diag(E);  Exp.ev=EE-EE(1);

%The Stevens Operators
A{1} = stev(10,[2,0]);
A{2} = stev(10,[4,0]);
A{3} = stev(10,[4,4]);
A{4} = stev(10,[2,2]);


constants.A0=sparse(length(A{1}));
constants.A = A;
constants.ev = EE;
%
Opt = struct('NDeflations',4,'Method','Good_GN','Linesearch','Quadratic',...
'c1',1e-8,'constants',constants,'StepTolerance',1e-5);
obj_fun = @INSEvaulateDifference;
[SysOutGood,options]= DMin(obj_fun,x0(1:end),Opt);


obj_fun = @INSEvaulateDifference;
Opt = struct('NDeflations',5,'Method','Bad_GN','Linesearch','Quadratic',...
    'scaled',false,'c1',1e-5,'constants',constants,'theta',2,'MaxIter',2e4);
[SysOutBad,options]= DMin(obj_fun,x0(1:end),Opt);


obj_fun = @INSEvaulateDifference;
Opt = struct('NDeflations',4,'Method','Newton','Linesearch','Quadratic',...
    'Verbose',false,'scaled',false,'c1',1e-4,'constants',constants);
[SysOutNewton]= DMin(obj_fun,x0(1:end),Opt);
%
clf
f=figure(1);
subplot(3,1,1)
options.ShowLegend = false; options.ShowDeflations = 1:length(SysOutGood);
options.ShowNonMinima = false;
PlotxConvergence(SysOutGood,options)
grid on
xlabel('k')
ylabel('error')
% legend('location','eastoutside')
ylim([1e-10,1e5])

subplot(3,1,2)
options.ShowDeflations = 1:length(SysOutBad);
options.ShowLegend = true;
PlotxConvergence(SysOutBad,options)
legend('location','eastoutside')
grid on
xlabel('k')
ylabel('error')
ylim([1e-10,1e5])
xlim([0,350])

subplot(3,1,3)
options.ShowLegend = false;
options.ShowDeflations = 1:length(SysOutNewton);
PlotxConvergence(SysOutNewton,options)
% legend('location','eastoutside')


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
f.Position = [10 10 20 12];
print(f, 'sec4Mn12.eps', '-depsc')

%%


function [f,r,Jr] = EvaluateBratu(c,constants)
m = constants.m;
A = constants.A;
D2 = constants.D2;
e0 = constants.e0;
et = constants.et;

% Residual for the BVP: u'' + 3*exp(u) = 0, u(0) = 0, u(1) = 0:
r =  [(1/m+1)*(A*(D2*c) + 3*exp(A*c)); ...
     e0*c;
     et*c];

f =  dot(r,r)/2;

if nargout>2
    Jr =  [(1/m+1)*(A*D2 + 3*spdiags(exp(A*c),0,m,m)*A); ...
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
% epsilon = constants.epsilon;

% Residual for the BVP: 0.05u'' + u^2 + 8x(1-x)u = 1, u(0) = 0, u(1) = 0:
r =  [(1/m+1)*(0.05*A*(D2*c) + (A*c).^2 + 8*(xi.*(1-xi))'.*(A*c) - 1); ...
     e0*c;
     et*c];
f =  dot(r,r)/2;

if nargout>2
    Jr =  [(1/m+1)*(0.05*A*D2 + 2*spdiags(A*c,0,m,m)*A+ 8*(xi.*(1-xi))'.*A); ...
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
lgnd = strcat(hLegend(1).String, string(newline), format_errs([Iterations(all([str~="Max Iterations reached";str~="Merit line search terminated with rank deficient Jacobian"],1)).ErrorAtDeflatedPoint]));
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




function [F,R,J,H] = INSEvaulateDifference(x,constants)
A = constants.A; 
    Ad = x(1)*A{1};
    for i = 2:length(x)
        Ad = Ad + x(i)*A{i};
    end
    [Q,D] = eig(full(Ad),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
    R = ((D(2:end) - D(1:end-1)) - (constants.ev(2:end) - constants.ev(1:end-1)));
    F = sqrt(sum((R).^2));
    if nargout>2
        l = length(A);
        m = size(Q,2);
        LJ = zeros(m,l);
        for k = 1:l
            LJ(:,k) =real(sum((Q.'*A{k}).*Q',2)); 
        end
        J = LJ(2:end,:) - LJ(1:end-1,:);
    end
    if nargout>3
        m = length(constants.ev); l=length(A);
        LH = zeros(l,l,m);
        QAQ = cell(1,l);
        for i = 1:l
            QAQ{i} = Q'*A{i}*Q;
            QAQ{i} =QAQ{i}(1:m,1:m);
        end
        DD=D'-D;
        DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                LH(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
        H = LH(:,:,2:end) - LH(:,:,1:end-1);
    end
    
end



function [F,R,J,H] = INSEvaulate(x,constants)
A = constants.A; 
    Ad = x(1)*A{1};
    for i = 2:length(x)
        Ad = Ad + x(i)*A{i};
    end
    [Q,D] = eig(full(Ad),'vector');
    D = D(1:length(constants.ev));
    Q = Q(:,1:length(constants.ev));
    F = sqrt(sum((D-constants.ev).^2));
    if nargout>1
        R = (D-constants.ev);
    end
    if nargout>2
        l = length(A);
        m = size(Q,2);
        J = zeros(m,l);
        for k = 1:l
            J(:,k) =real(sum((Q.'*A{k}).*Q',2)); 
        end
    end
    if nargout>3
        m = length(constants.ev); l=length(A);
        H = zeros(l,l,m);
        QAQ = cell(1,l);
        for i = 1:l
            QAQ{i} = Q'*A{i}*Q;
            QAQ{i} =QAQ{i}(1:m,1:m);
        end
        DD=D'-D;
        DD(abs(DD)<1e-15) = Inf;
        for j=1:l
            for k = 1:l
                H(k,j,:) = real(2*sum(QAQ{k}.*QAQ{j}./DD));
            end
        end
    end
end


