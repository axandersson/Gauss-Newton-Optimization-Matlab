function printMessage(flag,method_name,varargin)
flags = inputParser;
search_names = {'no','Armijo','Wolfe'}; % Different search methods
addOptional(flags,'maxIter',50,@isnumeric);
addOptional(flags,'nIter',1,@isnumeric);
addOptional(flags,'SearchParams',nan,@iscell)
addOptional(flags,'Search','no',@(x)any(validatestring(x,search_names)));
parse(flags,varargin{:});

search_params = flags.Results.SearchParams;
search_name = flags.Results.Search;
switch flag
    case 0
        disp([method_name,' converged with number of iterations: ', ...
            num2str(flags.Results.nIter)]);
    case -1
        disp([method_name,' did not converge']);
    case -2
        disp([method_name,' reached maximum number of iterations: maxIter = ', ...
            num2str(flags.Results.maxIter)]);
    case 2
        disp(['Running ', method_name]);
        switch search_name
            case 'no'
            case 'Armijo'
                fprintf('alpha_min = %f, c1 = %f\n', search_params{1}, ...
                    search_params{2});
            case 'Wolfe'
                fprintf('alpha_min = %f, alpha_max = %f, c1 = %f, c2 = %f\n',...
                    search_params{1}, search_params{2}, search_params{3}, ...
                    search_params{4});
        end
end