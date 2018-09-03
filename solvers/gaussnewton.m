function [x,FLAG,nIter,X,f,ALPHA] = gaussnewton(fun,x0,params,varargin)
%%GAUSSNEWTON
% solves the non-linear least-square problem.
%
%  [x,nIter] = gaussnewton(fun,x0,params)
%
% finds the optimal variables x that minimizes the residual function fun.
% Maximum number of iterations is set to $50$, and the convergence
% tolerance is set to $10^{-8}$.
%
%INPUT
%
% * fun A function handle to the residual and jacobian of function.
%   [r,J] = fun(x,params{:})
%   will return the return the residuals r
%   and the Jacobian J for the variable x for some funtion. The
%   cell params contains additional parameters required to the
%   function.
%
% * x0 Initial guess for the variable x.
%
% * params A cell with additional parameters required for the residual
%   unction fun.
%
%OUTPUT
%
% * x The optimal variables x that minimizes the funtion residual
%   function.
%
% * FLAG A flag indicating wether or not the algorithm converged.
%   If 0, the algorithm has converged. If -1 the algorithm did not converge.
%   If -2, the algorithm reached maximum number of iterations.
%
% * nIter The number of itertion required for convergence.
%
% * X An array with values for x ar each iteration.
%
% * ALPHA An array with values of step-lengths used by the alogrithm
%   each iteration. See examples.
%
%EXAMPLES
%
% Consider a data-set of number of antelopes y at times t. We want to fit
% an exponetial function to this data-set. We construct the residual
% function:
%
%   function [r,J] = antelope_r(x,t,y)
%   r = (x(1).*exp(x(2).*t) - y);
%   J = [exp(x(2).*t),x(1).*exp(x(2).*t).*t];
%
% [x,FLAG,nIter,X,alpha] = gaussnewton(@antelope_r, x0, params),
% where params{1} = t, and params{2} = y, will find the optimal values of
% x(1) and x(2) that minimizes the residual in the function antelope_r.
%
%  [x,FLAG,nIter,X,alpha] = gaussnewton(@antelope_r, x0, params,
%        'maxIter',30)
%
% will set the maximum number of iterations to 30.
%
%  [x,FLAG,nIter,X,alpha] = gaussnewton(@antelope_r, x0,
%       params,'convTol',1e-6)
%
% will set the convergence tolerance to 1e-6.
%
%  [x,FLAG,nIter,X,alpha] = gaussnewton(@antelope_r, x0,
% params,'Search','Armijo','SearchParams',[alpha_min, c1])
%
% will use Armijo line-search with the specified parameters for alpha_min
% and c1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check input parameters and parse them.
search_methods = {'Armijo','Wolfe','no'};

input_params = inputParser;
addOptional(input_params,'convTol',1e-8,@(x)validateattributes(x, ...
    {'numeric'},{'positive'},{'nonempty'}))
addOptional(input_params,'maxIter',50, @(x)validateattributes(x,...
    {'numeric'},{'nonempty','integer','positive'}))
addOptional(input_params,'Method','Gauss-Newton', ...
    @(x)validateattributes(x,{'char'},{'nonempty'}))
addOptional(input_params,'Search','no', ...
    @(x)any(validatestring(x,search_methods)))
addOptional(input_params,'SearchParams',nan, ...
    @(x)any(validatestring(x,plot_options)))
parse(input_params,varargin{:});

% Store variables and parameters
maxIter = input_params.Results.maxIter;
convTol = input_params.Results.convTol;
search_params = num2cell(input_params.Results.SearchParams);
search_name = input_params.Results.Search;
method_name = ['Gauss-Newton method with ', search_name, ' line-search.'];

% Asign line-search functions
% If no paramteres for Armijo line-search is give from input, use
% $alpha_{min} = 10^{-3}$, and $c_1 = 10^{-4}$.
switch search_name
    case 'no'
        search_fun = @(x,p) 1;
    case 'Armijo'
        errmsg = ['Search parameters for ',search_name,' line-search', ...
            ' must be of length 2, that is [alpha_min, c1]'];
        assert(length(search_params) == 2 || isnan(search_params{1}), ...
            errmsg);
        
        if isnan(search_params{1})
            search_params{1} = 1e-3;
            search_params{2} = 1e-4;
        end
        search_fun = @(x,p) feval(@armijo_search,fun,x,p, ...
            params,search_params{:});
    case 'Wolfe'
        errmsg = ['Search parameters for ',search_name,' line-search', ...
            ' must be of length 4, that is [alpha_min, alpha_max, c1, c2]'];
        assert(length(search_params) == 4 || isnan(search_params{1}),errmsg);

        if isnan(search_params{1})
            search_params{1} = 1e-5; % alpha_min
            search_params{2} = 10; % alpha_max
            search_params{3} = 1e-4; % c1
            search_params{4} = 0.9; % c2
        end

        search_fun = @(x, p) feval(@wolfe_search, fun, x, p, maxIter, params, ...
            search_params{:});
end

% Print method name and parameters used.
printMessage(2,method_name,'Search',search_name, ...
    'SearchParams',search_params);

% Initialize output vector and initial values.
nDim = length(x0);
X = zeros(nDim,maxIter);
f = zeros(1,maxIter);
ALPHA = zeros(1,maxIter);
x = x0;
alpha = 1;

% The big loop, also stores data into vectors for output.
for k = 1:maxIter
    [r,J] = feval(fun,x,params{:});

    X(:,k) = x;
    ALPHA(k) = alpha;
    f(k) = 0.5 * (r' * r);
    % Evaluates the residual/Jacobian function.
    % Recall that our objective function is $f(r) = \frac{1}{2} r^{T} r$,
    % and the gradient $g_f = J^T r$. Its Hessian is approximated as
    % $H \approx J^T J$.
    g_f = J' * r;
    H = makeposdef(J'*J);
    p = -H\(g_f);
    
    
    % Stopping condition
    % If $\|g_f\| \leq \epsilon (1- \|r\|)$, we stop the iterations.
    % Also, if the stepping length $\alpha$ is equal to zero we will stop
    % the iterations.
    if ~alpha
        FLAG = -1;
        nIter = k;
        printMessage(FLAG,method_name,'nIter',nIter);
        return;
    end

    if norm(J*p) <= convTol*(1+norm(r))
        FLAG = 0;
        nIter = k;
        ALPHA = ALPHA(1:nIter);
        f = f(1:nIter);
        X = X(:,1:nIter);
        printMessage(FLAG,method_name,'nIter',nIter);
        return;
    end
    % Compute a step length
    % The function
    %   search_fun(x,p)
    % compute a stepping length $\alpha$. The function depends on the input
    % from the user. It can be either Armijo line-search, Wolfe
    % line-search, or no line search. In the latter case, $\alpha$
    % is always  1, and the program performs the regular
    % Gauss-Newton method.
    alpha = search_fun(x,p);
    x = x+alpha*p;
end
nIter = k;
FLAG = -2;
printMessage(FLAG,method_name,'maxIter',nIter);