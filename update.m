function [solution, ierr_update, info_update] = update(problem, solution,ctrl)
	ierr_update = 0;

	% collect info
	info_udpate={};
	info_update.iterations=0;
	info_update.cpu_linalg=0.0;
	
	

	npot=problem.npot;
	ntdens=problem.ntdens;
	if (ctrl.scheme == 'gf')
		% solve gfvar with dumped newton


		% move in the gfvar variable
		gfvar = tdens2gfvar(ctrl.study_transformation, 2.0, solution.tdens);
		gfvar_old = gfvar;

		% init jacobian structure
		jacobian=jacobian_structure;
		jacobian.init(problem.matrix_E,problem.weight);

		% run newton algorithm
		info_newton = 0;
		iter_newton = 0;
		cpu_linalg = 0;

		while info_newton == 0 
			% evalute non-linear functional (including dirichlet bcs)
			gradpot = problem.compute_grad(solution.pot);
			[f]=nonlinear_equations(gfvar,gfvar_old,gradpot,problem,ctrl,solution);
			
			
			% check convergence
			f_norm=norm(f);
			f_norm_pot=norm(f(1:npot));
			f_norm_gfvar=norm(f(npot+1:npot+ntdens));

			print_msg(2,...
								sprintf('%2d : |F| %1.2e (%1.2e, %1.2e)',...
												iter_newton,f_norm, f_norm_pot,f_norm_gfvar),...
								ctrl)
			
			if (f_norm < ctrl.tol_nonlinear)
				info_newton = 0;
				print_msg(1,...
								sprintf('NEWTON CONVERGED in %2d iter. with |F| %1.2e (%1.2e, %1.2e) ',...
												iter_newton,f_norm, f_norm_pot,f_norm_gfvar),...
								ctrl)
				break
			end
			
			
			% check max iterations
			iter_newton = iter_newton+1;
			if (iter_newton > ctrl.max_nonlinear_iterations)
				info_newton = 1;
				break
			end
			
			% update jacobian components
			jacobian.set(gfvar,gradpot,ctrl,solution);

			% print_extreme(gfvar,'gfvar')
			% print_extreme(gradpot,'grad')
			% print_extreme(jacobian.w11,'w11')
			% print_extreme(jacobian.w12,'w12')
			% print_extreme(jacobian.w21,'w21')
			% print_extreme(jacobian.w22,'w22')
			
			% solve J s = -F
			tic;
			[increment,ierr_linear_solver,info_linear_solver] = solvesys(jacobian, -f, ctrl);
			cpu_linalg = cpu_linalg + toc;
			if ~ierr_linear_solver == 0
				info_newton = 2;
			end
					
			% set update lenght
			step = 1.0;
			
			% update solution
			solution.pot=solution.pot+step*increment(1:npot);
			gfvar=gfvar+step*increment(npot+1:npot+ntdens);			
		end

		% convert to tdens
		solution.tdens = gfvar2tdens(ctrl.study_transformation, 2.0, gfvar);

		% strore info
		info_update.iterations=iter_newton;
		info_update.cpu_linalg=cpu_linalg;
		
	end

	
end
