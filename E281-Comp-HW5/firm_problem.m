function [K,L] = firm_problem(par,num,grid,r,av_labor)

  L = av_labor;
  K = (par.alpha/(r+par.delta))^(1/(1-par.alpha)) *L;

end