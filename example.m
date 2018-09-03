close all, clear all;
addpath('solvers');
%% Example for gaussnewton optimization.
%  We have experimental data, y, of a population of antelopes at different
%  times t. We wish to find the coefficients x(1) and x(2) such that the 
%  function f(t) = x(1)exp(x(2)*t)) fits the experimental data of the
%  antelope population. Mathematically, we wish to solve:
%  min_x sum_i^n ||f(t_i) - y_i ||.

%  This can be done using a gauss-newton method.
%  First create a residual function see antelope_r. 
%  The residual function must return the residuals, i.e. 
%  r = x(1)*exp(x(2)*t) - y, and the Jacobian of the residual function,
%  i.e. dr/dx. See antelope_r.m.

% Here we simulate some "experimental data" with noise
alpha = 1; beta = 2; t = 0:0.1:1;
y = alpha .* exp(beta.*t) + rand(size(t));

% Here we make an initial guess 
x0 = [1.2 ; 2.1];

% Parameters for the residual function. See antelope_r.
params{1} = t(:); params{2} = y(:);

% Perform gaussnewton optimization. 
[x,FLAG,nIter,X,f,ALPHA] = gaussnewton(@antelope_r,x0,params);

% We can also perform gauusnewton with line-search. A line-search asserts
% that the solution does not blow up.
% [x,FLAG,nIter,X,f,ALPHA] = gaussnewton(@antelope_r,x0,params, ...
% 'Search','Armijo');

figure(1);
plot(x(1)*exp(x(2)*t)); hold on;
plot(y,'o');
leg = legend('Optimized','Experimental data','location','best');
set(leg,'interpreter','latex');
ylabel('Antelope population','interpreter','latex');
xlabel('Time','interpreter','latex');
