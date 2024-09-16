function PlotFE(xi,Iterations,obj_fun,evalfun,constants,options)
% if isempty(varargin)
%     ShowNonMinima = true;
% else
%     ShowNonMinima = varargin{1};
% end

% clf
hold on;    lgnd = [];
warning("off",'MATLAB:plot:IgnoreImaginaryXYPart')
for i = 1:length(Iterations)
    if ~contains(Iterations(i).ConvergenceFlag,["Max Iterations reached","Merit line search terminated with rank deficient Jacobian"])
        plot(xi,evalfun(Iterations(i).DeflatedPoint,xi'),'linewidth',2)
        lgnd = [lgnd;num2str(obj_fun([Iterations(i).DeflatedPoint],constants),'%0.4e')];
    elseif options.ShowNonMinima
        colorOrder = get(gca, 'ColorOrder');
        plot(xi,evalfun(Iterations(i).DeflatedPoint,xi'),'-','HandleVisibility','on','linewidth',1,...
            'Color', [colorOrder(mod((get(gca,'ColorOrderIndex'))-1, size(colorOrder, 1))+1, :), 0.2])
        set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex'));

    else
        % set(gca, 'ColorOrderIndex', get(gca,'ColorOrderIndex')+1);
    end
end
if options.ShowLegend
    legend(lgnd);
end
hold off
warning("on",'MATLAB:plot:IgnoreImaginaryXYPart')
end