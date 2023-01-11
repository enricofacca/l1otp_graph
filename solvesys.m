function[sol,ierr,info_linear_solver] = solvesys(jacobian, rhs, ctrl);

ntdens=jacobian.ntdens;
npot=jacobian.npot;

%if (ctrl.approach_linear_solver == 'reduced'):
	% form components
	[A,B1T,B2,C]=jacobian.form();

	% form C^{-1} , the primal schur complement, and its inverse
  inv_C=1.0./diag(C);
	inv_C = sparse(1:ntdens,1:ntdens,(inv_C)',ntdens,ntdens);

	%S_primal = jacobian.matrix_E *W_S*jacobian.matrix_E';
	S_primal = A + B1T * inv_C * B2;

	% inverse of primal schur complement
	ctrl_S = ctrl_solver;
	ctrl_S.init('agmg',... % method
							1e-1,... % tolerance
							20,...% iter
							0.0,...% omega (not used)
							0,... % verbose
							'invSp',...% label
							1);% agmg_seed
						 
	inv_S = sparse_inverse;
	inv_S.init(S_primal+ctrl.relax_prec11*speye(npot,npot), ctrl_S);

	
	% block precondtiner
	prec = @(z) PrimalSchur_prec(z, ...
															 @(x) inv_S.apply(x),...
															 @(y) inv_C*y,...
															 @(x) B1T*x,...
															 @(y) B2*y,...
															 npot,ntdens,...
															 'full');

	% set controls outer solver
	ctrl_outer = ctrl_solver;
	ctrl_outer.init('fgmres',1e-4,200); 
	[sol,info_linear_solver]= apply_iterative_solver(@(z) jacobian.apply(z), rhs, ...
																		 ctrl_outer, ...
																		 [],prec);

	ierr=0;

%end

end
