function [Mu,gradMu] = deflation(Y,x,varargin)
%    DEFLATION Calculates shifted deflation operators
%    Mu = DEFLATION(y,x) produces the shifted deflation operator Mu at point
%    x, given a row vector of deflated points Y. 
%
%    mu(x) = prod_i (i
%    The default operator uses the
%    square of the 2-norm and a shift of 1.
%
%    [Mu,gradMu] = DEFLATION(y,x) also produces the gradient of the deflation
%    operator.
%
%    [Mu,gradMu] = DEFLATION(y,x,theta) changes the exponent on the norm to
%    theta. Unless theta is 'exp', then it calculates the exponential deflation
%    operator.
%
%    [Mu,gradMu] = DEFLATION(y,x,theta,sigma) Changes the value of the shift to
%    sigma
%
%    [Mu,gradMu] = DEFLATION(y,x,theta,sigma,true) Uses the single shift
%    strategy for multiple deflations.
%
%    [Mu,gradMu] = DEFLATION(y,x,DeflationParameters) An alternative way to
%    pass in the optional variables, given that all options are included.

if isempty(Y)
    Mu = 1;
    gradMu = zeros(1,length(x));
    return
end

if nargin == 3 && isstruct(varargin{1})
    if ~isfield(varargin{1},'epsilon')&&length(fieldnames(varargin{1})) ~= 4 ...
            || isfield(varargin{1},'epsilon')&&length(fieldnames(varargin{1})) ~= 5
        error('Please input all optional parameters as structure')
    end
    theta = varargin{1}.theta;
    sigma = varargin{1}.sigma;
    SingleShift = varargin{1}.singleshift;
    NormWeighting = varargin{1}.NormWeighting;
else
    if nargin < 2
        error("Deflation requires the input of the deflated points, and the point to be evaluated")
    elseif nargin > 6
        error('Too many inputs')
    else
        if nargin > 2
            theta = varargin{1};
        else
            theta = 2; % default
        end
        if nargin > 3
            sigma = varargin{2};
        else
            sigma = 1; % default
        end
        if nargin > 4
            SingleShift = varargin{3};
        else
            SingleShift = false; % default
        end
        if nargin > 5
            NormWeighting = varargin{4};
        else
            NormWeighting = speye(length(x)); % default
        end
    end
end

if length(x)~=size(Y,1)   
    error("The dimension of the current iterate and of deflated points must be the same.")
end

l = length(x);
MUs = NaN(1,size(Y,2));
if strcmp(theta,"exp")  % strcmp needed because theta could be a double
    if SingleShift
        for i = 1:size(Y,2)
            MUs(i) = exp(1 / norm(NormWeighting*(x - Y(:,i)))) - 1;
        end
        Mu = sum(MUs)+sigma;
        if nargout>1
            gradMu = zeros(1,l);
            for i = 1:size(Y,2)
                gradMu = gradMu - ((NormWeighting'*NormWeighting)*(x - Y(:,i))/norm(NormWeighting*(x - Y(:,i)))^3.*(exp(1/norm(NormWeighting*(x - Y(:,i))))-1)).';
            end
        end
    else
        for i = 1:size(Y,2)
            MUs(i) = exp(1/norm(NormWeighting*(x - Y(:,i))))-1+sigma;
        end
        Mu = prod(MUs);
        if nargout > 1
            gradMu = zeros(1,l);
            for i = 1:size(Y,2)
                gradMu = gradMu +prod([MUs(1:i-1),MUs(i+1:size(Y,2))])*(-(NormWeighting'*NormWeighting)*(x - Y(:,i))/norm(NormWeighting*(x - Y(:,i)))^3.*(sigma-1+exp(1/norm(NormWeighting*(x - Y(:,i)))))).';
            end
        end
    end
else
    if SingleShift
        MUs = (1./sum(abs(NormWeighting*(x-Y)).^2,1).^(theta/2));
        Mu = sigma + prod(MUs);
    else
        MUs = sigma+ (1./sum(abs(NormWeighting*(x-Y)).^2,1).^(theta/2));
        Mu = prod(MUs);
    end
        if nargout > 1
            gradMu = ((-theta*(NormWeighting'*NormWeighting))*sum(((Mu./MUs)./sum(abs(NormWeighting*(x-Y)).^2,1).^(1+theta/2)).*(x-Y),2)).';
        end
end
if ~exist('gradMu','var')
    gradMu = [];
end
if isinf(Mu)|| any(isnan(gradMu))
    if any(abs(Y - x)<1e-16)
        warning('The current iterate is a deflated point')
    end

        warning('Deflation operators calculated at the current point are Inf/NaN')
end
end

