function PlotFContours(problem,options,obj_fun)
% ShowNDeflations,epsilon,xylim,NPoints,plotlines
% minima,Iterates,NIter,Flags,obj_fun
xl = -options.xylim;    xh = options.xylim;     yl=-options.xylim;      yh=options.xylim;
x = linspace(xl,xh,options.NPoints);
y = linspace(yl,yh,options.NPoints);
[X,Y] = meshgrid(x,y);
Z=zeros(length(y),length(x));   
if isempty(problem(options.ShowDeflations))
    DeflatedPoints = [];
else
    DeflatedPoints  = [problem(options.ShowDeflations).DeflatedPoint];
end
for i=1:length(y)
    for j=1:length(x)
%         [Mu,gradMu] = deflation(problem.minima(:,options.ShowNDeflations),[X(i,j);Y(i,j)]);
%         [Z(i,j)] = obj_fun([X(i,j),Y(i,j)]);
          Z(i,j) = deflation(DeflatedPoints,[X(i,j);Y(i,j)])^2*obj_fun([X(i,j);Y(i,j)],options.constants);
    end
end


% N=6;   L=-20;   n=50;
% [~,P] = contourf(X,Y,log10(Z),L:(N/n):N,'edgecolor','none','HandleVisibility','off');
contourf(X,Y,log10(Z),options.levels,'edgecolor','none','HandleVisibility','off');
hold on
set(gca,'ColorOrderIndex',1)

for i = 1:length(problem)
    if (problem(i).ConvergenceFlag == "Objective less than tolerance") ...
            || (problem(i).ConvergenceFlag == "Gradient less than tolerance")...
            || (problem(i).ConvergenceFlag == "Step Size below tolerance")||...
            problem(i).ConvergenceFlag == "Merit line search terminated"
        if any(ismember(i,options.plotlines))
            plot(problem(i).Iterates(1,:),problem(i).Iterates(2,:),'-*','LineWidth',2.5,'HandleVisibility','off'); %Plots convergence paths
        else
            plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'LineWidth',1.5,'HandleVisibility','off'); %Only plots found DeflatedPoint
        end
        xlim([-options.xylim options.xylim]);
        ylim([-options.xylim options.xylim]);
        plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'s','MarkerSize',13,'MarkerFaceColor','black','MarkerEdgeColor','black','HandleVisibility','off')
    else
        if any(ismember(i,options.plotlines))
            plot(problem(i).Iterates(1,:),problem(i).Iterates(2,:),'-*','LineWidth',2.5,'HandleVisibility','off'); %Plots convergence paths
        else
            plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'LineWidth',1.5,'HandleVisibility','off'); %Only plots found DeflatedPoint
        end
        xlim([-options.xylim options.xylim]);
        ylim([-options.xylim options.xylim]);
        plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'s','MarkerSize',7,'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[ 0, 0.4470, 0.7410],'HandleVisibility','off')

    end
end
hold off