function alpha_k = wolfe_zoom(fun, xk, pk, maxIter, params, ...
    alpha_low, alpha_high, c1, c2)
%% WOLFE_ZOOM
%   Helps WOLFE_SEARCH to find a suitable step length alpha_k by zooming in
%   a certain interval
%
%   alpha_k = WOLFE_ZOOM(fun, xk, pk, maxIter, params, ...,
%       alpha_min, alpha_max, c1, c2)
%       Returns a suitable step length for the residual function fun at the
%       point of evaluation xk, search direction pk with the parameters params
%       being sent to fun. maxIter iterations are allowed before the function
%       returns. alpha_min and alpha_max are the minimum and maximum allowed
%       step lengths, respectively, and c1 and c2 are constants.
%
%%  Input
%       * fun A function handle to the residual and Jacobian function
%
%       * xk The point of evaluation of fun
%
%       * pk The search direction
%
%       * maxIter The maximum number of iterations before the function returns
%
%       * params Parameters sent to fun
%
%       * alpha_min Minimum step length
%
%       * alpha_max Maximum step length
%
%       * c1 and c2 Constants such that 0 < c1 < c2 < 1
%
%%  Output
%       alpha_k The suitable step length
%
%
alpha_k = 0;

% Returns f = (1/2) * (r') * r and grad_f = (J') * r
phi = @(alpha) get_f_grad(fun, xk, pk, alpha, params);

% Convert to the right form for the conditions
cond_converter = @(dphi) (pk') * dphi;

[phi_0, dphi_0] = phi(0);
dphi_0 = cond_converter(dphi_0);

j = 1;
while j <= maxIter
    % Bisect
    alpha_j = (alpha_low + alpha_high)/2;

    [phi_j, ~] = phi(alpha_j);

    [phi_low, ~] = phi(alpha_low);

    if (phi_j > phi_0 + c1 * alpha_j * dphi_0) || (phi_j >= phi_low)
        alpha_high = alpha_j;

        j = j + 1;
        continue;
    end

    [~, dphi_j] = phi(alpha_j);
    dphi_j = cond_converter(dphi_j);

    if abs(dphi_j) <= -c2 * dphi_0
        alpha_k = alpha_j;
        return;
    end

    if dphi_j * (alpha_high - alpha_low) >= 0
        alpha_high = alpha_low;
    end

    alpha_low = alpha_j;

    j = j + 1;
end
