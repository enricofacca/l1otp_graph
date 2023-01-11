function [tdens, pot, ierr]=l1_solver(problem, ctrl, tdens_initial, pot_initial)
	% problem dimension
	ntdens=size(tdens_initial,1);
	npot=size(pot_initial,1);
	
	% create a structure to store the solution [pot,tdens]
	% including the node where the solution is freezed
	solution={};	
	solution.ntdens = ntdens;
	solution.tdens = tdens_initial; 

	solution.npot = npot;
	solution.pot   = pot_initial;
	
	solution.ntdens_dirichlet = 0;
	solution.tdens_dirichlet = zeros(solution.ntdens,1);
	solution.npot_dirichlet = 0;
	solution.pot_dirichlet = zeros(solution.npot,1);

	% start iterate
	iter = 0;
	ierr = 0;

	% keep updating the solution until an error occurs.
	% Update restart are possible
	while ( ierr == 0)
		% update with restart
		info_restart = 0; 
		nrestarts = 0;

		% store solution before update
		solution_old = solution;
		
		print_msg(1,sprintf('UPDATE %d DELTA %1.2e',iter,ctrl.deltat),ctrl)
		while info_restart == 0;
			% udpate
			[solution, ierr_update, info_update] = update(problem, solution,ctrl);
			if (ierr_update == 0)
				% break cycle
				break
			else
				nrestarts = nrestarts + 1
				if ( nrestarts == ctrl.max_restarts)
					info_restart = 1;
					ierr = 1;
					break
				end
				% set controls for next update
				ctrl.deltat = ctrl.deltat * ctrl.contraction_deltat;

				% restore solution
				solution = solution_old;
			end
			print_msg(1,sprintf('UPDATE %d RESTART %d DELTA %1.2e',iter,nrestart,ctrl.deltat),ctrl)
		end


		% check errors
		if (info_restart ~= 0 )
			ierr = 1;
		else
			%fprintf('UPDATE %d SUCCEED \n',iter)
			
			% print info system
			gradpot = problem.compute_grad(solution.pot);
			fprintf('%3d : max (gradpot)-1=%1.2e\n',iter,max(abs(gradpot))-1.0)
			fprintf('%3d : %1.2e<= TDENS<=%1.2e\n',iter,min(solution.tdens),max(solution.tdens))
			if (problem.opt_pot_exists)
				[m,imax]=max(problem.forcing);
				err_pot=norm(problem.opt_pot-solution.pot-solution.pot(imax));
				fprintf('%3d : err_pot=%1.2e\n',iter,err_pot)
			end

			if (problem.opt_tdens_exists)
				err_tdens = problem.weighted_norm(problem.opt_tdens-solution.tdens,1);
				fprintf('%3d : err_tdens = %1.2e\n',iter,err_tdens)
			end			
			
			% check if convergence is achived
			var_tdens = problem.weighted_norm(solution.tdens-solution_old.tdens,1);
			fprintf('%3d : var_tdens=%1.2e\n',iter,var_tdens)
			if var_tdens < ctrl.tol
				fprintf('CONVEGENCE ACHIEVED \n')
				ierr = 0;
				break
			end
			% print info update

			% save data	
			iter = iter + 1;
			if (iter > ctrl.max_iterations)
				ierr = 2;
			end

			%
			% set controls for next update
			%
			if (info_restart == 0)
				% enalrge time step
				ctrl.deltat = ctrl.deltat * ctrl.expansion_deltat;
				% set new dirichlet regions
				if (ctrl.selection > 0 )
					solution = threshold_solution(solution, problem, ctrl);
				end
			end

			
		end
	end

	% return the solution
	pot=solution.pot;
	tdens=solution.tdens;
end
