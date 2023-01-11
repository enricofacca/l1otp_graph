label_problem=strcat("rect_0064")
label_graph=strcat("grid_0064_graph")
loader;

% set controls
ctrl={};
%
% globals
%
ctrl.study_transformation="square";
ctrl.power_transformation=2.0;

ctrl.max_iterations=50;
ctrl.tol=1e-3;

ctrl.verbose=1;

ctrl.selection = 1;
ctrl.threshold_selection = 1e-8;


% relaxation controls
ctrl.deltat=0.1;
ctrl.min_deltat=1e-6;
ctrl.max_deltat=1e6;
ctrl.expansion_deltat=1.5;

% restart controls
ctrl.max_restarts=10;
ctrl.contraction_deltat=2.0;


% non-linear solver 
ctrl.scheme="gf";
ctrl.max_nonlinear_iterations=10;
ctrl.tol_nonlinear=1e-8;

% linear solver
ctrl.relax_prec11=1.0e-11;


%
% init. problem structure
%
prob = problem_structure;
prob.init(topol, rhs, sizecell);
optpot;
prob.set_opt_pot(optpot);
prob.set_opt_tdens(opttdens);


%
% run solver 
%
[tdens, pot, ierr] = l1_solver(prob, ctrl, tdens, pot);


%
% save solution to fiel
%
write_td(label_problem+'tdens.dat',tdens);
write_td(label_problem+'pot.dat',pot);
