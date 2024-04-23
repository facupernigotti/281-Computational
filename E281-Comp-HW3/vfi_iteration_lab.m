function [v_new,c] = vfi_iteration_lab(v0,param,num,grid)
% implement the algorithm in the slides.

v_old = v0 ;

% Perform the upwind scheme.
[vp_upwind1] = vp_upwind(v_old,param,num,grid) ;

% Infer consumption from Va_Upwind
c = max(vp_upwind1,1e-08).^(-1);
u = utility(c) ;

% Compute the change between rho*V and the new iteration u + v'*s
Vchange = u +  vp_upwind1 .* (param.y + param.r * grid.a -c) - param.rho * v_old ;

% update the value function with a step parameter Delta
v_new = v_old + num.Delta*Vchange;


end
