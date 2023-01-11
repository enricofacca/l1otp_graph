classdef jacobian_structure <handle
	% structure describing the jacobian matrix
	% given by
	% (A   )(W11 W12)(A^T )
	% (  I )(W21 W22)(   I)
  properties
		initalized=0;
    ntdens=0;
		npot=0;
		matrix_E;
		weight;
		w11;w12;w21;w22;
		mat11;mat12;mat21;mat22;
    name='empty';
		ntdens_dirichlet;
		tdens_dirichlet;
		tdens_penalty; % sparse matrix with penalty in the tdens-edges
		npot_dirichlet;
		pot_dirichlet;
		pot_penalty; % sparse matrix with penalty in the pot nodes
  end
  methods
    function obj = init(obj,matrix_E,weight)
			obj.initalized=1;
      obj.ntdens=size(matrix_E,2);
			obj.npot=size(matrix_E,1);
      obj.matrix_E=matrix_E;
      obj.weight=weight;
			obj.w11=zeros(obj.ntdens,1);
			obj.w12=zeros(obj.ntdens,1);
			obj.w21=zeros(obj.ntdens,1);
			obj.w21=zeros(obj.ntdens,1);
    end

		% set the vectors w11,w12,w21,w22
		function y = set(obj,gfvar,gradpot,ctrl,solution)
			% short hands
			ntdens=obj.ntdens;
			npot=obj.npot;
			
			tdens=gfvar2tdens(ctrl.study_transformation,ctrl.power_transformation,gfvar);
			[trans_prime,trans_second] = gfvar2f(ctrl.study_transformation,ctrl.power_transformation,gfvar);
			% update jacobian components
			obj.w11 = tdens./obj.weight;
			obj.w12 = trans_prime.*gradpot;
			obj.w21 = trans_prime.*gradpot;
			obj.w22 = ( obj.weight.*(0.5*trans_second.*(gradpot.^2-1.0))-...
									obj.weight/ctrl.deltat);

			% set dirichlet matrix for pot and tdens (may be a null matrix)
			penalty=zeros(npot,1);
			penalty(obj.pot_dirichlet)=1e30;
			obj.pot_penalty = sparse(1:npot,1:npot,penalty,npot,npot);

			penalty=zeros(ntdens,1);
			penalty(obj.tdens_dirichlet)=1e30;
			obj.tdens_penalty = sparse(1:ntdens,1:ntdens,penalty,ntdens,ntdens);

			
			
			% update matrix W
			obj.mat11=sparse(1:ntdens,1:ntdens,obj.w11',ntdens,ntdens);
			obj.mat12=sparse(1:ntdens,1:ntdens,obj.w12',ntdens,ntdens);
			obj.mat21=sparse(1:ntdens,1:ntdens,obj.w21',ntdens,ntdens);
			obj.mat22=sparse(1:ntdens,1:ntdens,obj.w22',ntdens,ntdens);
		end
		
    % y = J * x
    function y = apply(obj,x)
			npot=obj.npot; ntdens=obj.ntdens;
			y=zeros(npot+ntdens,1);
			y(1:npot) = obj.matrix_E * (obj.w11 .* (obj.matrix_E'*x(1:npot)))+...
									obj.matrix_E * (obj.w12 .* x(npot+1:npot+ntdens))+...
									obj.pot_penalty*x(1:npot);
			
			y(npot+1:npot+ntdens) = obj.w21 .* (obj.matrix_E'*x(1:npot))+...
															obj.w22 .* x(npot+1:npot+ntdens)+...
															obj.tdens_penalty * x(npot+1:npot+ntdens);
		end

		% form explicitely the component of J
		function [A,B1T,B2,C] = form(obj)
			A = obj.matrix_E*obj.mat11*obj.matrix_E'+obj.pot_penalty;
			B1T = obj.matrix_E*obj.mat12;
			B2 = obj.mat21*obj.matrix_E';
			C = -obj.mat22-obj.tdens_penalty;
		end
		
		
		% destructor
		function obj = kill(obj)
			if (obj.initalized == 1)					
				clear obj.ntdens;
				clear obj.npot;
				clear obj.matrixA;
				clear obj.weight;
				clear obj.w11;
				clear obj.w12;
				clear obj.w21;
				clear obj.w21;
				obj.initalized=0;
			end
		end

		% info
		function obj = info(obj,fid)
      if (~exist('fid','var') )
				fid=1;
			end
			obj.ctrl.info(fid);
			if ( strcmp(obj.ctrl.approach ,'krylov'))
				if ( obj.is_symmetric) 
					fprintf(fid,'nnz(IL)/nnz(A),nnz(A) %f3 %d \n',...
									nnz(obj.IL)/nnz(obj.matrix),nnz(obj.matrix));
				else
					fprintf(fid,'(nnz(IL),nnz(IU))/nnz(A)| nnz(A) %f3 %f3 %d\n',...
									nnz(obj.IL)/nnz(obj.matrix),...
									nnz(obj.IU)/nnz(obj.matrix),...
									nnz(obj.matrix));
				end
			end
		end
	end
end
