function [alpha, FCount]= Linesearch(x, p, DeflatedPts, params, linesearchparams, obj_fun, constants_at_x)
alpha = linesearchparams.alpha0;    FCount = 0;
if linesearchparams.method ~= "No"
    phix = linesearchparams.phi(obj_fun, x, params.method.constants, DeflatedPts, params.deflation);
    newphi = @(alpha)linesearchparams.phi(obj_fun, x + alpha*p, params.method.constants, DeflatedPts, params.deflation);
    % newphix = linesearchparams.phi(obj_fun, x + alpha*p, params.method.constants, DeflatedPts, params.deflation);
    newphix = newphi(alpha);
    FCount = FCount+1;
    switch linesearchparams.method
        case 'Basic'
            while (newphix > phix) && (alpha > linesearchparams.minalpha)
                alpha = alpha*linesearchparams.tau;
                newphix = newphi(alpha);
                FCount = FCount+1;
            end
        case 'Armijo'
            if ~isfield(constants_at_x,'S')
                constants_at_x.S = [];
            end
            gradphix = linesearchparams.gradphi(constants_at_x);
            while (newphix > phix + linesearchparams.c1*alpha*dot(gradphix,p)) ...
                    && (alpha > linesearchparams.minalpha)
                alpha = alpha*linesearchparams.tau;
                newphix = newphi(alpha);
                FCount = FCount+1;
            end
        case 'Quadratic'
            if ~isfield(constants_at_x,'S')
                constants_at_x.S = [];
            end
            gradphix = linesearchparams.gradphi(constants_at_x);
            alphaprev = alpha;
            while (newphix > phix + linesearchparams.c1*alpha*dot(gradphix,p)) ...
                    && (alpha > linesearchparams.minalpha)
                alpha = dot(gradphix,p) / (2*(newphix+phix+dot(gradphix,p)));
                if alpha < 0.3*alphaprev
                    alpha = 0.3*alphaprev;
                elseif alpha > 0.5*alphaprev
                    alpha = 0.5*alphaprev;
                end
                alphaprev = alpha;
                newphix = newphi(alpha);
                FCount = FCount+1;
            end
        otherwise
            error("Please input a valid line search method: No, Basic, Armijo, Quadratic")
    end
end
end % function

