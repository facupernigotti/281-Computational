function[EK] = excessk(r)
% create structure with structural parameters
call_parameters;
% create structure with numerical parameters
numerical_parameters ;
% create structure with grids and initial guesses for v
grid = create_grid(par,num) ;

w = (1-par.alpha)*(par.alpha/(r+par.delta))^(par.alpha/(1-par.alpha));

v_old = grid.v0 ;

dist = 1 ;
while dist > num.tol 
    [v_new,~,~] = vfi_iteration(v_old,par,num,grid,w,r) ;
    dist = max(abs((v_new(:) - v_old(:)))) ;
    v_old = v_new ;
   % disp(dist)
end
[v_new,c,A] = vfi_iteration(v_new,par,num,grid,w,r) ;

[gg] = kf_equation(A,grid,num) ;
g = [gg(1:num.a_n),gg(num.a_n+1:2*num.a_n)];
av_wealth= sum(sum(grid.a .*g.*grid.da)) ;

av_labor= sum(sum(par.e .*g.*grid.da)) ;

[fK,fL] = firm_problem(par,num,grid,r,av_labor) ;

EK=(av_wealth-fK)^2;
end
