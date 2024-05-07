%% Part 1: homogeneous returns
clear all ; close all ; clc ;

rstar=fminsearch(@excessk,0.012)

call_parameters;
% create structure with numerical parameters
numerical_parameters ;
% create structure with grids and initial guesses for v
grid = create_grid(par,num) ;

w = (1-par.alpha)*(par.alpha/(rstar+par.delta))^(par.alpha/(1-par.alpha));

v_old = grid.v0 ;

dist = 1 ;
while dist > num.tol 
    [v_new,~,~] = vfi_iteration(v_old,par,num,grid,w,rstar) ;
    dist = max(abs((v_new(:) - v_old(:)))) ;
    v_old = v_new ;
   % disp(dist)
end
[v_new,c,A] = vfi_iteration(v_new,par,num,grid,w,rstar) ;

[gg] = kf_equation(A,grid,num) ;
g = [gg(1:num.a_n),gg(num.a_n+1:2*num.a_n)];

figure()
plot(grid.a,g)
xlim([num.a_min 1])

%% Part 2: heterogeneous returns
%The only difference is on the file vp_upwind_hr, the files
%vfi_iteration_hr and excessk_hr are the same with the exepction they take
%vp_upwind_hr instead of vp_upwind as the input.
clear all ; close all ; clc ;

rstar=fminsearch(@excessk_hr,0.012)

call_parameters;
% create structure with numerical parameters
numerical_parameters ;
% create structure with grids and initial guesses for v
grid = create_grid(par,num) ;

w = (1-par.alpha)*(par.alpha/(rstar+par.delta))^(par.alpha/(1-par.alpha));

v_old = grid.v0 ;

dist = 1 ;
while dist > num.tol 
    [v_new,~,~] = vfi_iteration_hr(v_old,par,num,grid,w,rstar) ;
    dist = max(abs((v_new(:) - v_old(:)))) ;
    v_old = v_new ;
   % disp(dist)
end
[v_new,c,A] = vfi_iteration_hr(v_new,par,num,grid,w,rstar) ;

[gg] = kf_equation(A,grid,num) ;
g = [gg(1:num.a_n),gg(num.a_n+1:2*num.a_n)];

figure()
plot(grid.a,g)
xlim([num.a_min 1])





