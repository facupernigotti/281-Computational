%% Solves the Krusell and Smith (1998)
%{
In order to adapt into our original Aiyagari problem I 
(1) replaced the corresponding parameters in set_parameters_aiyagari
(2) defined a utility function to account for the case where ggamma is 1
(3) changed the initial guess v0 in compute_state_aiyagari
After that compute_steady_state is exactly the same as the old codes put
together
%}

tstart = tic;

%% Setup the toolbox
% Just need to include folders containing the files in the path. (manually)

addpath('C:\Users\fpern\OneDrive\Comp Methods\Topic6/MATLABAutoDiff-master');
addpath('C:\Users\fpern\OneDrive\Comp Methods\Topic6\phact-master');


% initialize shocks for simulation
T = 200;
N = 2000;
vAggregateShock = zeros(1,N);
vAggregateShock(1,1) = 1;

%% Step 0: Set Parameters
% The script sets up parameters relevant for the model
set_parameters_aiyagari;

%% Step 1: Solve for Steady State
tStart = tic;

fprintf('Computing steady state...\n')
global IfSS IbSS I0SS varsSS A

[rSS,wSS,KSS,ASS,uSS,cSS,VSS,gSS,dVUSS,dVfSS,dVbSS,IfSS,IbSS,I0SS] = ...
    compute_steady_state(); %(ESTO HAY QUE TOCAR)

fprintf('Time to compute steady state: %.3g seconds\n\n\n',toc(tStart));

% Store steady state values
varsSS = zeros(nVars,1);
varsSS(1:2*I,1) = reshape(VSS,2*I,1);
ggSS = reshape(gSS,2*I,1);
varsSS(2*I+1:4*I-1,1) = ggSS(1:2*I-1);
varsSS(4*I,1) = 0;
varsSS(4*I+1,1) = KSS;
varsSS(4*I+2,1) = rSS;
varsSS(4*I+3,1) = wSS;
varsSS(4*I+4,1) = (KSS ^ aalpha) * (zAvg ^ (1 - aalpha));
CSS = sum(cSS(:) .* gSS(:) * da);
varsSS(4*I+5,1) = CSS;
varsSS(4*I+6,1) = ddelta * KSS;

%% Step 2: Linearize Model Equations
% For computing derivatives, the codes written for solving for the
%    steady-state can be used almost verbatim using automatic
%    differentiation toolbox as long as only the functions supported by
%    automatic differentation are used. For list of supported functions and
%    documentation of relevant syntax check <<https://github.com/sehyoun/MATLABAutoDiff>>
fprintf('Taking derivatives of equilibrium conditions...\n')
t0 = tic;

% Prepare automatic differentiation
vars = zeros(2*nVars+nEErrors+1,1);
vars = myAD(vars);

% Evaluate derivatives (ESTO HAY QUE TOCAR)
derivativesIntermediate = equilibrium_conditions(vars);

% Extract out derivative values
derivs = getderivs(derivativesIntermediate);

% Unpackage derivatives
mVarsDerivs = derivs(:,1:nVars);
mVarsDotDerivs = derivs(:,nVars+1:2*nVars);
mEErrorsDerivs = derivs(:,2*nVars+1:2*nVars+nEErrors);
mShocksDerivs = derivs(:,2*nVars+nEErrors+1);

% rename derivatives to match notation in paper
g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
c = sparse(nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;


[state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi);
n_g_red = n_g;

%% Step 3: Solve Out Static Constraints and/or Reduce Models
    reduceDist_hor = 5;        % m of K_m(A,b)

% rename derivatives to match notation in paper
g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
c = sparse(nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;

t0 = tic;
fprintf('Model Reduction ...\n')

%%
% *State space reduction using Krylov subspace method*
    % State space reduction
    [state_red,inv_state_red,n_g_red] = krylov_reduction(g0,g1,n_v,n_g,reduceDist_hor);
    [g1,psi,pi,c,g0] = change_basis(state_red,inv_state_red,g1,psi,pi,c,g0);

%%
% *Value function reduction using spline inspired bases*
    % Create knot points for spline (the knot points are not uniformly spaced)
    knots = linspace(amin,amax,n_knots-1)';
    knots = (amax-amin)/(2^c_power-1)*((knots-amin)/(amax-amin)+1).^c_power+amin-(amax-amin)/(2^c_power-1);
    % Function calls to create basis reduction
    [from_spline, to_spline] = oneDquad_spline(x,knots);
    [from_spline, to_spline] = extend_to_nd(from_spline,to_spline,n_prior,n_post);
    n_splined = size(from_spline,2);
    [from_spline, to_spline] = projection_for_subset(from_spline,to_spline,0,n_g_red);
    
    % Reduce the decision vector
    [g1,psi,~,c,g0] = change_basis(to_spline,from_spline,g1,psi,pi,c,g0);
    pi = to_spline * pi * from_spline(1:n_v,1:n_splined);
fprintf('...Done!\n')
fprintf('Time to reduce dimensionality: %2.4f seconds\n\n\n',toc(t0))


%% Step 4: Solve Linear System
t0 = tic;
fprintf('Solving reduced linear system...\n')

[G1,~,impact,eu,F] = schur_solver(g0,g1,c,psi,pi,1,1,1);

fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0))

%% Step 5: Simulate Impulse Response Functions
fprintf('Simulating Model...\n')
t0 = tic;

trans_mat = inv_state_red*from_spline;
[simulated,vTime] = simulate(G1,impact,T,N,vAggregateShock,'implicit',trans_mat);

fprintf('...Done!\n')
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0))

% Add state-states back in to get values in levels
varsSS_small = varsSS(4*I:4*I+6,1);
vAggregateTFP = simulated(400,:) + varsSS(400);
vAggregateOutput = simulated(404,:) + varsSS(404);
vAggregateConsumption = simulated(405,:) + varsSS(405);
vAggregateInvestment = simulated(406,:) + varsSS(406);

% Compute log differences for plotting
vAggregateTFP_reduced = vAggregateTFP;
vAggregateOutput_reduced = log(vAggregateOutput) - log(varsSS(404));
vAggregateConsumption_reduced = log(vAggregateConsumption) - log(varsSS(405));
vAggregateInvestment_reduced = log(vAggregateInvestment) - log(varsSS(406));

%% Step 6: Internal Consistency Check

    %{
    g1 = -mVarsDerivs;
    psi = -mShocksDerivs;
    from_red = inv_state_red * from_spline;
    to_red = to_spline * state_red;
    [epsilon] = internal_consistency_check(G1,impact,n_g_red,from_red,to_red,g1,psi,F,n_v,n_g,600,varsSS,1,0);
    %}

%% Step 7: Plot relevant values
% Plot impulse response functions
figure
subplot(2,2,1)
hold on
plot(vTime,100 * vAggregateTFP_reduced,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('TFP','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation from s.s.','interpreter','latex')
hold off

subplot(2,2,2)
hold on
plot(vTime,100 * vAggregateOutput_reduced,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Output','interpreter','latex','fontsize',14)
hold off

subplot(2,2,3)
hold on
plot(vTime,100 * vAggregateConsumption_reduced,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Consumption','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation from s.s.','interpreter','latex')
xlabel('Quarters','interpreter','latex')
hold off

subplot(2,2,4)
hold on
plot(vTime,100 * vAggregateInvestment_reduced,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Investment','interpreter','latex','fontsize',14)
xlabel('Quarters','interpreter','latex')
hold off

toc(tstart)


