function alpha = armijo_search(fun,x,p,params,alpha_min,c1)  
for alpha = 2.^-(0:52)
    if alpha < alpha_min
        alpha = 0;
        return
    end

    if sufficient_decrease_cond(fun,x,params,p,alpha,c1)
        return
    end
end
end