%% Chebfun Carrier equation
 
 
% clear all
clf
cheb.x

% N = chebop(-1,1); N.lbc = 0; N.rbc = 0;
% N.op = @(x,y) 0.04*diff(y,2) + 2*(1-x^2)*y + y^2; 
x = chebfun(@(x) x,[0,1]);
N = chebop(0,1); N.lbc = 0; N.rbc = 0;
N.op = @(x,y) 0.05*diff(y,2) + 8*x*(1-x)*y + y^2; 
N.init = (2*x-1).^2-1; y1 = N\1; plot(y1), hold on
N.init = 1-(2*x-1).^2; y2 = N\1; plot(y2)
N.init = sin(pi*(2*x-1)); y3 = N\1; plot(y3)
N.init = -sin(pi*(2*x-1)); y4 = N\1; plot(y4)
N.init = 4*sin(pi*(0.5*x-1))+1; y5 = N\1; plot(y5)