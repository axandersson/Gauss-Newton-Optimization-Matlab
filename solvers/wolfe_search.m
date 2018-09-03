function alpha_k = wolfe_search(fun, xk, pk, maxIter, params,alpha_min, alpha_max, c1, c2)
%% WOLFE_SEARCH
%   Uses the Wolfe line search algorithm in order to find a suitable step length
%
%   alpha_k = WOLFE_SEARCH(fun, xk, pk, maxIter, params, ...,
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

alpha_k = 0;

% Returns f = (1/2) * (r') * r and grad_f = (J') * r
phi = @(alpha) get_f_grad(fun, xk, pk, alpha, params);

% Convert to the right form for the conditions
cond_converter = @(dphi) (pk') * dphi;

% [alpha_{i - 1}; alpha_{i}] (no need to keep values from all iterations)
alpha = zeros(2, 1);

% alpha(2) \in (0, alpha_max)
alpha(2) = alpha_max/10;

[phi_0, dphi_0] = phi(0);
dphi_0 = cond_converter(dphi_0);

i = 1;
while i <= maxIter
    % Evaluate phi(alpha_{i}) and phi(alpha_{i-1})
    [phi_i, ~] = phi(alpha(2));
    [phi_i_old, ~] = phi(alpha(1));

    if (phi_i > phi_0 + c1 * alpha(2) * dphi_0) || ...
        ((phi_i >= phi_i_old) && (i > 1))
        alpha_k = wolfe_zoom(fun, xk, pk, maxIter, params, ...
            alpha(1), alpha(2), c1, c2);

        return;
    end

    [~, dphi_i] = phi(alpha(2));
    dphi_i = cond_converter(dphi_i);

    if abs(dphi_i) <= -c2 * dphi_0
        alpha_k = alpha(2);
        return;
    end

    if dphi_i >= 0
        alpha_k = wolfe_zoom(fun, xk, pk, maxIter, params, ...
            alpha(1), alpha(2), c1, c2, maxIter);

        return;
    end

    % alpha_{i + 1} \in (alpha_{i}, alpha_max)
    alpha(1) = alpha(2);
    alpha(2) = (alpha(2) + alpha_max)/2;

    i = i + 1;
end
