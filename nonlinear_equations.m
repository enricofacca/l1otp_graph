function [f] = nonlinear_equations(gfvar,gfvar_old,gradpot,problem,ctrl,solution)	
	tdens=gfvar2tdens(ctrl.study_transformation,ctrl.power_transformation,gfvar);
	[trans_prime,trans_second] = gfvar2f(ctrl.study_transformation,ctrl.power_transformation,gfvar);
	f=zeros(problem.ntdens+problem.npot,1);
	f(1:problem.npot)= problem.matrix_E*(tdens.*gradpot)-problem.forcing;
	f(problem.npot+1:problem.npot+problem.ntdens)= ...
	problem.weight.*(0.5*trans_prime.*(gradpot.^2-1.0))  ...
	- problem.weight.*(gfvar-gfvar_old)/ctrl.deltat;
	if (solution.npot_dirichlet > 0)
		f(solution.pot_dirichlet(1:solution.npot_dirichlet))=0;
	end
	if (solution.ntdens_dirichlet > 0)
		f(solution.npot + solution.tdens_dirichlet(1:solution.ntdens_dirichlet))=0;
	end
end
	
