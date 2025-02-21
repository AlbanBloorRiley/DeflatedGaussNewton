function PlotDeflatedFConvergence(problem,obj_fun,options)

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
    x = 1:length(problem(i).Iterates);
    y = nan(length(problem(i).Iterates),1);
    Beta = zeros(length(problem(i).Iterates),1);
    for j = 1:size(problem(i).Iterates,2)
        y(j) = obj_fun(problem(i).Iterates(:,j),options.constants);
        [Mu,gradMu] = deflation([problem(1:i-1).DeflatedPoint],problem(i).Iterates(:,j),options.theta,options.sigma,options.SingleShift);
        if options.Method =="Good_GN"||options.Method =="Bad_GN"||options.Method =="LM"
            [~,Rx,Jx] = obj_fun(problem(i).Iterates(:,j),options.constants);
            p =  -lsqminnorm(Jx,Rx);
        elseif options.Method =="Newton"
            [~,Rx,Jx,Hx] = obj_fun(problem(i).Iterates(:,j),options.constants);
            S = zeros(length(gradMu));
            for k = 1:length(Rx)
                S = S+Hx(:,:,k)*Rx(k);
            end
            p = - lsqminnorm(Jx'*Jx+S,Jx'*Rx);
        elseif options.Method =="GradientDescent"
            [~,Rx,Jx] = obj_fun(x,options.constants);
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
        semilogy(x,y,'linewidth',4,'color',colours(mod(i-1,7)+1,:))
        semilogy(x,yDeflated,'linewidth',1,'color',colours(mod(i-1,7)+1,:),'HandleVisibility','off')
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
% 
% function PlotDeflatedFConvergence(problem,obj_fun,varargin)
% if nargin ==2
%     options.ShowDeflations = 1:length(problem);
%     options.ShowLegend = false;
% elseif nargin ==3
%     options = varargin{1};
% else
%     error('Too many inputs')
% end
% if ~isfield(options,"constants")
%     options.constants = [];
% end
% % clf
% MaxIter = 0;
% for i = options.ShowDeflations
%     if MaxIter < length(problem(i).Iterates)
%         MaxIter = length(problem(i).Iterates);
%     end
% end
% hold on
% lgnd = ["Undeflated "];
% for i = options.ShowDeflations
%     X = 1:length(problem(i).Iterates);
%     xx = nan(length(problem(i).Iterates),1);
%     Beta = zeros(length(problem(i).Iterates),1);
%     for j = 1:size(problem(i).Iterates,2)
%         x = problem(i).Iterates(:,j);
%         xx(j) = norm(x-problem(i).Iterates(:,end));
%         [Mu,gradMu] = deflation([problem(1:i-1).DeflatedPoint],x,options.theta,options.sigma,options.SingleShift);
%         if options.Method =="Good_GN"||options.Method =="Bad_GN"||options.Method =="LM"
%             [~,Rx,Jx] = obj_fun(x,options.constants);
%             p =  -lsqminnorm(Jx,Rx);
%         elseif options.Method =="Newton"
%             [~,Rx,Jx,Hx] = obj_fun(x,options.constants);
%             S = zeros(length(gradMu));
%             for k = 1:length(Rx)
%                 S = S+Hx(:,:,k)*Rx(k);
%             end
%             p = - lsqminnorm(Jx'*Jx+S,Jx'*Rx);
%         elseif options.Method =="GradientDescent"
%             [~,Rx,Jx] = obj_fun(x,options.constants);
%             p = -0.1*Jx'*Rx;
%         else
%             error("")
%         end
%         Beta(j) = 1-dot((1/Mu)*gradMu,p);
%     end
%     Deflated = Beta<1-options.epsilon;
%     plotDeflated = any([Deflated,[0;Deflated(1:end-1)]],2);
%     plotNotDeflated = any([~Deflated,[0;~Deflated(1:end-1)]],2);
%     xxDeflated = xx; xxDeflated(~plotDeflated) = NaN;
%     xx(~plotNotDeflated) = NaN;
%     colours = colororder;
%     if ~contains(problem(i).ConvergenceFlag,["Max Iterations reached","Merit line search terminated with rank deficient Jacobian"])
%         semilogy(X,xx,'linewidth',2,'color',colours(mod(i-1,7)+1,:))
%         semilogy(X,xxDeflated,'--','linewidth',2,'color',colours(mod(i-1,7)+1,:),'HandleVisibility','off')
%     else
%         colorOrder = get(gca, 'ColorOrder');
%         semilogy(X,xx,'linewidth',1,'Color', [colorOrder(mod((get(gca,'ColorOrderIndex'))-1, size(colorOrder, 1))+1, :), 0.2])
%     end
%     if options.ShowLegend
%         entry = ['Deflation ', num2str(i)];
%         lgnd = [lgnd; entry];
%     end
% end
% if options.ShowLegend
%     lgnd = lgnd(1:end-1,:);
%     legend(lgnd)
% end
% set(gca, 'YScale', 'log')
% set(gca,'YMinorGrid','off')
% %  semilogy(x,xx,'linewidth',2)
% 
% 
% hold off
% 
% 
