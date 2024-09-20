function PlotBetaContours(problem,options,obj_fun)
% ShowNDeflations,epsilon,xylim,NPoints,plotlines
% DeflatedPoint,Iterates,NIter,Flags,obj_fun
xl = -options.xylim;    xh = options.xylim;     yl=-options.xylim;      yh=options.xylim;
x = linspace(xl,xh,options.NPoints);
y = linspace(yl,yh,options.NPoints);
[X,Y] = meshgrid(x,y);
Beta=zeros(length(y),length(x));
for i=1:length(y)
    for j=1:length(x)
        [Mu,gradMu] = deflation([problem(options.ShowDeflations).DeflatedPoint],[X(i,j);Y(i,j)],options.theta,options.sigma,options.SingleShift);
        if options.Method =="Good_GN"||options.Method =="Bad_GN"||options.Method =="LM"
            [~,Rx,Jx] = obj_fun([X(i,j),Y(i,j)],options.constants);
            p =  -lsqminnorm(Jx,Rx);
        elseif options.Method =="Newton"
            [~,Rx,Jx,Hx] = obj_fun([X(i,j);Y(i,j)],options.constants);
            S = zeros(length(gradMu));
            for k = 1:length(Rx)
                S = S+Hx(:,:,k)*Rx(k);
            end
            p = - lsqminnorm(Jx'*Jx+S,Jx'*Rx);
        elseif options.Method =="GradientDescent"
            [~,Rx,Jx] = obj_fun([X(i,j),Y(i,j)],options.constants);
            p = -0.1*Jx'*Rx;
        else
            error("")
        end
        Beta(i,j) = 1-dot((1/Mu)*gradMu,p);
    end
end


% NColours = 6;
% parulas = parula(NColours);
% parulas = parulas(NColours:-1:1,:);
% colours = colours(NColours:1,:)
% colours = turbo(NColours);
% colours = summer(NColours);
% colours = [parulas(2,:);parulas(floor(NColours/2+0),:);parulas(NColours-1,:)];
red = [255,160,128];
orange = [255,223,128];
green = [165,212,106];
% red = [255,160,128];
orange = [240, 225, 149];
green = [222, 245, 208];


colours = [red;orange;green]./255;
hold on;
[~,h1] = contourf(X,Y,Beta,[-inf,-inf],'edgecolor','none','HandleVisibility','off');
[~,h2] = contourf(X,Y,Beta,[0,0],'edgecolor','none','HandleVisibility','off');
[~,h3] = contourf(X,Y,Beta,[1-options.epsilon,1-options.epsilon],'edgecolor','none','HandleVisibility','off');
set(h1,'FaceColor',colours(1,:)); plot(nan, nan,'o','MarkerFaceColor',colours(1,:),'MarkerEdgeColor','none');
set(h2,'FaceColor',colours(2,:)); plot(nan, nan,'o','MarkerFaceColor',colours(2,:),'MarkerEdgeColor','none');
set(h3,'FaceColor',colours(3,:)); plot(nan, nan,'o','MarkerFaceColor',colours(3,:),'MarkerEdgeColor','none' );

if options.ShowLegend
    if options.epsilon ==0
        legend("\beta<0","0<\beta<1","1<\beta",'fontsize',options.FontSize);
    else
        legend("\beta<0","0<\beta<1-\epsilon","1-\epsilon<\beta",'fontsize',options.FontSize);
    end
end

for i = 1:length(problem)
    set(gca,'ColorOrderIndex',i)
    if any(ismember(i,options.plotlines))
        plot(problem(i).Iterates(1,:),problem(i).Iterates(2,:),'-*','LineWidth',2,'HandleVisibility','off'); %Plots convergence paths
    end
    xlim([-options.xylim options.xylim]);
    ylim([-options.xylim options.xylim]);

    if (problem(i).ConvergenceFlag == "Objective less than tolerance") ...
            || (problem(i).ConvergenceFlag == "Gradient less than tolerance")...
            || (problem(i).ConvergenceFlag == "Step Size below tolerance")||...
            problem(i).ConvergenceFlag == "Merit line search terminated"
        plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'s','MarkerSize',7,'MarkerFaceColor','black','MarkerEdgeColor','black','HandleVisibility','off')
    elseif problem(i).ConvergenceFlag == ""
    else
        plot(problem(i).DeflatedPoint(1),problem(i).DeflatedPoint(2),'s','MarkerSize',7,'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[ 0, 0.4470, 0.7410],'HandleVisibility','off')
    end
end
hold off