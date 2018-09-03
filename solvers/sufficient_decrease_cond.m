function bool = sufficient_decrease_cond(fun,x,params,p,alpha,c1)
[r,J] = feval(fun,x,params{:});
f = 0.5 * (r' * r);
g_f = J' * r;
RHS = f + c1 * alpha .* (g_f' * p);
r = feval(fun,x+alpha*p,params{:});
LHS = 0.5 * (r' * r);
bool = LHS <= RHS;

