%% PROBLEM 1 

clear all ; close all ; clc ;
tic ;
% create structure with:
call_parameters; % structural parameters
numerical_parameters ; %numerical parameters
grid = create_grid(param,num) ; %grids and initial guesses for v

% Iterate v until convergence
v_old = grid.v0 ;
dist = 1 ;
%iteration = 0; 
while dist > num.tol 
    [v_new,~] = vfi_iteration_im(v_old,param,num,grid,'endowment') ;
    dist = max(abs((v_new - v_old))) ;
    v_old = v_new ;
    %iteration = iteration + 1 ;
    %disp(['Iteration number: ', num2str(iteration), 'Error=',num2str(dist) ]);

end

[v_new,c,fchoice] = vfi_iteration_im(v_old,param,num,grid,'endowment') ;
toc ;

figure();
plot(grid.a, c, 'o', 'DisplayName', 'Consumption Policy Function');
hold on;
plot(grid.a, param.y + param.r * grid.a, 'r', 'DisplayName', 'True Consumption');
xlabel('Initial Assets');
ylabel('Optimal Consumption');
legend('Location', 'best');

%% PROBLEM 3 

clear all ; close all ; clc ;
tic ;
% create structure with:
call_parameters; % structural parameters
numerical_parameters ; %numerical parameters
grid = create_grid(param,num) ; %grids and initial guesses for v

% Iterate v until convergence
v_old = grid.v0 ;
dist = 1 ;
%iteration = 0 ;
while dist > num.tol 
    [v_new,~] = vfi_iteration_im(v_old,param,num,grid,'production') ;
    dist = max(abs((v_new - v_old))) ;
    v_old = v_new ;
    %iteration = iteration + 1 ;
    %disp(['Iteration number: ', num2str(iteration), 'Error=',num2str(dist) ]);

end

[v_new,c,fchoice] = vfi_iteration_im(v_old,param,num,grid,'production') ;
toc ;


%Need to define k and pi in order to make true consumption
kstar= ((param.r) / (num.fb_prod * num.fb_alpha))^(1/(num.fb_alpha-1)) ;
pi = num.fb_prod* kstar^num.fb_alpha - param.r * kstar;

figure();
plot(grid.a, c, 'o', 'DisplayName', 'Consumption Policy Function w production');
hold on;
plot(grid.a, pi + param.r * grid.a, 'b', 'DisplayName', 'True Consumption w production');
plot(grid.a, param.y + param.r * grid.a, 'r', 'DisplayName', 'True Consumption w endowment');
xlabel('Initial Assets');
ylabel('Optimal Consumption');
legend('Location', 'best');
ylim([0, 7]);


%% PROBLEM 4 
clear all ; close all ; clc ;
tic ;
% create structure with:
call_parameters; % structural parameters
numerical_parameters ; %numerical parameters
grid = create_grid(param,num) ; %grids and initial guesses for v


%I run the VFI for each value of kappa
kappa_values = [0, 0.1, 0.4];
v_new_kappa = zeros(length(grid.a), numel(kappa_values));
c_kappa = zeros(length(grid.a), numel(kappa_values));

% Loop over each value of kappa
for i = 1:numel(kappa_values)
    % Iterate v until convergence
    num.fg_kappa = kappa_values(i);
    v_old = grid.v0 ;
    dist = 1 ;
    %iteration = 0 ;
    while dist > num.tol 
        [v_new,~] = vfi_iteration_im(v_old,param,num,grid,'production-BG') ;
        dist = max(abs((v_new - v_old))) ;
        v_old = v_new ;
        %iteration = iteration + 1 ;
        %disp(['Iteration number: ', num2str(iteration), 'Error=',num2str(dist) ]);
    
    end
    [v_new_kappa(:,i), c_kappa(:,i),fchoice(:,i)] = vfi_iteration_im(v_old, param, num, grid, 'production-BG');
    toc ;
    end 

% Plot all solutions
figure();
hold on;
colors = lines(numel(kappa_values)); % Define colors for each kappa value
for i = 1:numel(kappa_values)
    % Separate data points based on fchoice
    fchoice_zero_idx = fchoice(:,i) == 0; % Indices where fchoice is 0
    fchoice_one_idx = fchoice(:,i) == 1;  % Indices where fchoice is 1
    
    % Plot data points with fchoice equal to 0 (lines) without labels
    plot(grid.a(fchoice_zero_idx), c_kappa(fchoice_zero_idx,i), '-', 'Color', colors(i,:),'HandleVisibility','off');
    
    % Plot data points with fchoice equal to 1 (thicker lines)
    plot(grid.a(fchoice_one_idx), c_kappa(fchoice_one_idx,i), 'LineWidth', 2, 'Color', colors(i,:), 'DisplayName', ['Consumption Policy Function for \kappa = ', num2str(kappa_values(i))]);
end
xlabel('Initial Assets');
ylabel('Optimal Consumption');
legend('Location', 'best');
annotation('textbox', [0.6 0.0 0.1 0.1], 'String', '{Note: Thinner (thicker) lines represents choice of f=f_B (f = f_G)}', 'LineStyle', 'none');

