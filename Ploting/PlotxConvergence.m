function PlotxConvergence(problem,varargin)
if nargin ==1
    options.ShowDeflations = 1:length(problem);
    options.ShowLegend = false;
elseif nargin ==2
    options = varargin{1};
else
    error('Too many inputs')
end
if ~isfield(options,"constants")
    options.constants = [];
end
% clf
MaxIter = 0;
for i = options.ShowDeflations
    if MaxIter < length(problem(i).Iterates)
        MaxIter = length(problem(i).Iterates);
    end
end
hold on
lgnd = ["Undeflated "];
for i = options.ShowDeflations
    x = 1:length(problem(i).Iterates);
    xx = nan(length(problem(i).Iterates),1);
    for j = 1:size(problem(i).Iterates,2)
        xx(j) = norm(problem(i).Iterates(:,j)-problem(i).Iterates(:,end));

    end
    if ~contains(problem(i).ConvergenceFlag,["Max Iterations reached","Merit line search terminated with rank deficient Jacobian"])
        semilogy(x,xx,'linewidth',1)
    else
        colorOrder = get(gca, 'ColorOrder');
        semilogy(x,xx,'linewidth',1,'Color', [colorOrder(mod((get(gca,'ColorOrderIndex'))-1, size(colorOrder, 1))+1, :), 0.2])
    end
    if options.ShowLegend
        entry = ['Deflation ', num2str(i)];
        lgnd = [lgnd; entry];
    end
end
if options.ShowLegend
    lgnd = lgnd(1:end-1,:);
    legend(lgnd)
end
set(gca, 'YScale', 'log')
set(gca,'YMinorGrid','off')
%  semilogy(x,xx,'linewidth',2)


hold off