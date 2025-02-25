function PlotDeflatedJFConvergence(problem,obj_fun,options)

if ~isfield(options,"constants")
    options.constants = [];
end
if ~isfield(options,"fontsize")
    options.fontsize = 9;
end

lgnd = string;
hold on
markers = ["-o","-^","-x","-s","-h","-d","-v"];
for i = options.ShowDeflations
    x = 0:length(problem(i).Iterates)-1;
    y = nan(length(problem(i).Iterates),1);
    Beta = zeros(length(problem(i).Iterates),1);
    for j = 1:size(problem(i).Iterates,2)
        % y(j) = obj_fun(problem(i).Iterates(:,j),options.constants);
        [~,Rx,Jx] = obj_fun(problem(i).Iterates(:,j),options.constants);
        y(j) = norm(Jx'*Rx);
        [Mu,gradMu] = deflation([problem(1:i-1).DeflatedPoint],problem(i).Iterates(:,j),options.theta,options.sigma,options.SingleShift);
        if options.Method =="Good_GN"||options.Method =="Bad_GN"||options.Method =="LM"
            p =  -lsqminnorm(Jx,Rx);
        elseif options.Method =="Newton"
            [~,Rx,Jx,Hx] = obj_fun(problem(i).Iterates(:,j),options.constants);
            S = zeros(length(gradMu));
            for k = 1:length(Rx)
                S = S+Hx(:,:,k)*Rx(k);
            end
            p = - lsqminnorm(Jx'*Jx+S,Jx'*Rx);
        elseif options.Method =="GradientDescent"
            p = -0.1*Jx'*Rx;
        else
            error("")
        end
        Beta(j) = 1-dot((1/Mu)*gradMu,p);
    end
    Deflated = Beta<1-options.epsilon;
    plotDeflated = any([Deflated,[0;Deflated(1:end-1)]],2);
    plotNotDeflated = any([~Deflated,[0;~Deflated(1:end-1)]],2);
    yDeflated = y; yDeflated(~plotDeflated) = NaN;
    y(~plotNotDeflated) = NaN;
    colours = colororder;

    if ~contains(problem(i).ConvergenceFlag,["Max Iterations reached","Merit line search terminated with rank deficient Jacobian"])
        %     msemilogy(x,y,markers(mod(i-1,7)+1),'linewidth',1,'color',[colours(mod(i-1,7)+1,:),0.5],'HandleVisibility','off')
        % semilogy(x,yDeflated,markers(mod(i-1,7)+1),'linewidth',3,'color',[colours(mod(i-1,7)+1,:),1],'HandleVisibility','on')
                   semilogy(x,y,markers(mod(i-1,7)+1),'linewidth',2,'color',[colours(mod(i-1,7)+1,:),1],'HandleVisibility','on',MarkerSize=5)
        semilogy(x,yDeflated,markers(mod(i-1,7)+1),'linewidth',1,'color',[colours(mod(i-1,7)+1,:),0.5],'HandleVisibility','off',MarkerSize=5)
            entry = ['Deflation ', num2str(i-1)];
    lgnd = [lgnd; entry];
    elseif options.ShowNonMinima
        colorOrder = get(gca, 'ColorOrder');
        semilogy(x,y,'linewidth',1,'Color', [colorOrder(mod((get(gca,'ColorOrderIndex'))-1, size(colorOrder, 1))+1, :), 0.2])   
        entry = ['Deflation ', num2str(i-1)];
    lgnd = [lgnd; entry];
        else
        continue
    end

end
hold off
set(gca, 'YScale', 'log')
set(gca,'YMinorGrid','off')
if options.ShowLegend
    lgnd = lgnd(2:end,:);
    if lgnd(1)=="Deflation 0"
        lgnd(1) = ["Undeflated "];
        % % lgnd = lgnd(2:end,:);
    end
    legend(lgnd,'fontsize',options.fontsize)
end
