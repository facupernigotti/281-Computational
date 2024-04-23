function [sf, sb, Va_Upwind, fchoice] = shifts(v0,param,num,grid,type)

% unpack the initial guess for v0 and call it V
V = v0; 


% Define the value based on the type
switch type
    case 'endowment'
        y_value = param.y;
        fchoice= zeros(length(grid.a),1);
    case 'production'
        kstar= (param.r / (num.fb_prod * num.fb_alpha))^(1/(num.fb_alpha-1)) ;
        y_value = num.fb_prod* kstar^num.fb_alpha - param.r * kstar;
        fchoice= zeros(length(grid.a),1);        
    case 'production-BG'
        % Choice of k* given fb:
        kbstar= (param.r / (num.fb_prod * num.fb_alpha))^(1/(num.fb_alpha-1)) ;
        kbstar= min(grid.a,kbstar) ;
        pi_b = num.fb_prod* kbstar.^num.fb_alpha - param.r * kbstar;
        % Choice of k* given fg:
        kgstar= num.fg_kappa+(param.r / (num.fg_prod * num.fg_alpha))^(1/(num.fg_alpha-1)) ;
        kgstar= min(grid.a,kgstar) ;
        pi_g = num.fg_prod* (max(kgstar-num.fg_kappa,0)).^num.fg_alpha - param.r * kgstar;
        %Choice of production function:
        y_value=max(pi_b,pi_g) ;
        fchoice = (y_value == pi_g);
    otherwise
        error('Invalid type value. Type must be ''endowment'', ''production'', or ''production-BG''.');
end

%initialize forward and backwards differences with zeros: Vaf Vab
Vaf = zeros(num.a_n,1) ;
Vab = zeros(num.a_n,1) ;

% use V to compute backward and forward differences
Vaf(1:end-1) = (V(2:end)-V(1:end-1))/grid.da;
Vab(2:end) = (V(2:end)-V(1:end-1))/grid.da;

% impose the following boundary conditions. We will talk about them.
Vaf(end) = 0; 
%state constraint boundary condition
Vab(1) = (param.r*grid.a(1) + y_value(1)).^(-1);     


%consumption and savings with forward difference
cf = max(Vaf,1e-08).^(-1);
sf = y_value + param.r * grid.a - cf;

%consumption and savings with backward difference
cb = max(Vab,1e-08).^(-1);
sb = y_value + param.r * grid.a - cb;

%consumption and derivative of value function at steady state
c0 = y_value + param.r * grid.a ;
Va0 = c0.^(-1);

% compute indicator functions that capture the upwind scheme.
If = sf>0; %positive drift --> forward difference
Ib = sb<0; %positive drift --> forward difference
I0 = (1-If-Ib) ; %positive drift --> forward difference

% Compute the upwind scheme
Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0 ;

end
