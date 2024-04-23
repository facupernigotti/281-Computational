function [v_new,c,fchoice] = vfi_iteration_im(v0,param,num,grid,type)
% implement the algorithm in the slides.

v_old = v0 ;

% Perform the upwind scheme.
[sf1, sb1, vp_upwind1, fchoice] = shifts(v_old,param,num,grid,type) ;

% Infer consumption from Va_Upwind
c = max(vp_upwind1,1e-08).^(-1);
u = utility(c) ;

% Generate the matrix A(vn)
left = -min(sb1,0)/grid.da ;
mid = min(sb1,0)/grid.da - max(sf1, 0)/grid.da ;
right = max(sf1, 0)/grid.da ;
Bin=[left, mid, right] ;
Av1 = spdiags(mid,0,num.a_n,num.a_n)+spdiags(left(2:num.a_n),-1,num.a_n,num.a_n)+spdiags([0;right(1:num.a_n-1)],1,num.a_n,num.a_n);


% update the value function with a step parameter Delta
%v_new = inv((param.rho+(1/num.delta_im)) * speye(num.a_n) - Av1)* (u+v_old/num.delta_im);
v_new = ((param.rho+(1/num.delta_im)) * speye(num.a_n) - Av1) \ (u+v_old/num.delta_im);
fchoice=fchoice ;

end
