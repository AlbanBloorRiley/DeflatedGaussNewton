function PlotFConvergence(problem,options,obj_fun)

if ~isfield(options,"constants")
    options.constants = [];
end
if ~isfield(options,"fontsize")
    options.fontsize = 9;
end
% clf
lgnd = ["Undeflated "];
lgnd = string;
hold on

for i = options.ShowDeflations
    x = 0:length(problem(i).Iterates)-1;
    y = nan(length(problem(i).Iterates),1);
    for j = 1:size(problem(i).Iterates,2)
        y(j) = obj_fun(problem(i).Iterates(:,j),options.constants);
    end
    if ~contains(problem(i).ConvergenceFlag,["Max Iterations reached","Merit line search terminated with rank deficient Jacobian"])
        semilogy(x,y,'linewidth',1)
            entry = ['Deflation ', num2str(i-1)];
    lgnd = [lgnd; entry];
    elseif options.ShowNonMinima
        colorOrder = get(gca, 'ColorOrder');
        semilogy(x,y,'linewidth',1,'Color', [colorOrder(mod((get(gca,'ColorOrderIndex'))-1, size(colorOrder, 1))+1, :), 0.2])   
        entry = ['Deflation ', num2str(i-1)];
    lgnd = [lgnd; entry];
        % else
    %     continue
    end

end
hold off
set(gca, 'YScale', 'log')
set(gca,'YMinorGrid','off')
if options.ShowLegend
    lgnd = lgnd(2:end,:);
    if lgnd(1)=="Deflation 0"
        lgnd(1) = ["Undeflated "];
        % lgnd = lgnd(2:end,:);
    end
    legend(lgnd,'fontsize',options.fontsize)
end
