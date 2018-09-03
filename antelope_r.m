%%ANTELOPE_R Residual/Jacobian function for the antelope problem.

function [r,J] = antelope_r(x,t,y)
r = (x(1).*exp(x(2).*t) - y);
if nargout > 1
    % Analytic Jacobian
    %J = [exp(x(2).*t),x(1).*exp(x(2).*t).*t];
    
    % Numeric Jacobian
    J = jacapprox(@antelope_r,x,1e-6,{t,y});
end
