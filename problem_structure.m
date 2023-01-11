classdef problem_structure <handle
	% structure describing the jacobian matrix given by (A
	% )(W11 W12)(A^T ) (  I )(W21 W22)(   I)
  properties
		% size of the problem
		% ntdens = number of edges
		% npot = nnode
	  ntdens=0; npot=0;
		% graph topology
		topol;
		% signed incidence matrix
		matrix_E;
		% forcing term
		forcing;
		% weights and derivated quantities
		weight;inv_weight;sqrt_weight;
    name='empty';
		% reference solutions (if the exist)
		opt_pot;
		opt_tdens;
		opt_pot_exists=0;
		opt_tdens_exists=0;
  end
  methods
		function obj = init(obj,topol, forcing, weight)
			% constructor
			% args:
			% topol(int): (2,nedge) array with graph topology.
			%             [n1,n2]=topol(:,i) are the node indeces
			%             in i^th-edge
			% forcing(float) : array (len=number of nodes) with the forcing term.
			%                  MUST SUM TO ZERO (tol=1e-10)
			% weight(float) : array (len=number of edges) with edge weight 
			%               Optional (DEFAULT = 1.0)
			% returns:
			% obj : instatiated object
			obj.ntdens = size(topol,2);
			obj.npot = cast(max(topol(:))-min(topol(:))+1,"like",1); %handle 0 or
			%1-based numbering

			obj.forcing=forcing;
			if sum(forcing) > 1e-10
				error('Forcing term passed does not sum to zero')
			end

			% define sparse matrix with signed incidence matrix
			obj.topol = topol;
			obj.matrix_E = myincidence(topol);

			% set edge-weight
			if (exist('weight','var'))
				obj.weight = weight;				
			else
				obj.weight = ones(ntdens,1);
			end
			% variable usefull in certain computation
			obj.inv_weight = 1.0./weight;
			obj.sqrt_weight = sqrt(obj.weight);
		end

		function gradpot = compute_grad(obj,pot)
			% compute the gradient of the potential
			gradpot = obj.inv_weight.*(obj.matrix_E'*pot);
		end

		function obj = set_opt_pot(obj,opt_pot)
			% compute the gradient of the potential
			obj.opt_pot_exists = 1;
			obj.opt_pot = opt_pot;
		end

		function obj = set_opt_tdens(obj,opt_tdens)
			% compute the gradient of the potential
			obj.opt_tdens_exists = 1;
			obj.opt_tdens = opt_tdens;
		end

		function w_norm = weighted_norm(obj,vec,exponent)
			% compute the weighted norm of edge-based vectors
			if (~exist('exp','var'))
				exponent = 2;
			end

			if (exponent == 2)
				w_norm = norm(vec .* obj.sqrt_weight);
			elseif (exponent == 1)
				w_norm = norm(vec .* obj.weight,1);
			else
				error('Only 1 and 2 norm are computed in problem_structure')
			end	
		end

	

		function select_topol = select_graph(topol,edge_values,threeshold)
			% create a direct topology according to pot
			select_topol = topol(:,edge_values > threshold)
		end

		function [gamma] = build_plan(obj, topol,coord,vel, pot)
			approx_zero = 1e-13;

			% create container for plan
			obj.npot
			gamma = sparse(obj.npot,obj.npot);
			source_nodes = find(obj.forcing > approx_zero);

			% get flow on directed graph
			flow = abs(vel);
			rhs = obj.forcing;
			direct_topol = direct_graph(topol,pot);
			
			% cycle all source point
			for i = 1:size(source_nodes,1)
				isource = source_nodes(i);
				
				% take the mass 
				% if all mass was descarted we can move to the next
				% source point
				f_source = rhs(isource);
				inode = isource;
				fprintf('node %d source = %.1e (x,y)=%f,%f\n',inode,f_source,coord(1,inode),coord(2,inode));
				path_length = 0;
				while (f_source > approx_zero)
					%fprintf('node %d (x,y)=%f,%f\n',inode,coord(1,inode),coord(2,inode))
					% get the nodes reached by the flux
					% 
					edges = find(direct_topol(1,:) == inode);
					other_edges = find(direct_topol(2,:) == inode);
					nodes = direct_topol(2,edges);
					flow(edges);
					ind = find(flow(edges) > approx_zero);
					edges = edges(ind);
					flow(edges)';
					[sorted_flow, perm]=sort(flow(edges),'descend');
					perm;
					sorted_flow;
					edges=edges(perm);
					flow(edges)';
					nodes = direct_topol(2,edges);
					nedges = size(edges,1);

					if ( nedges == 0)
						% go back to previous biforcation
						edge = other_edge;
						fprintf('  went back to edge %d\n',edge)
					else
						% select the first edge and the corresponding node
						edge = edges(1);
						path_length = 	path_length + 1;
					end
					
					node = direct_topol(2,edge);
					fprintf('nedge %d - path %d- edge %d node %d flow %.1e\n',nedges,	path_length, edge, node,flow(edge))
					fprintf('node %d source = %.1e (x,y)=%f,%f\n',node,f_source,coord(1,node),coord(2,node))
					
					if (size(edges,2) >= 2 )
						% store biforcation point
						other_edge = edges(2);
					end
					
					% get the value of forcing at the new node
					f_node = obj.forcing(node);

					% reduce the flow by the mass injected
					flow(edge) = flow(edge) - f_source;
					
					% reduce the rhs by the flow passing through the edge
					if (inode == isource)
						rhs(inode) = rhs(inode) - flow(edge);
					end

					% if we get to a sink node reduce f_source and assign gamma
					if (f_node < - approx_zero)
						fprintf('sink node f=%.1e \n',f_node)
						f_source = f_source + f_node;
						gamma(isource,inode) = -f_node;
						fprintf('transport f=%.1e \n',f_source)
					end
					inode = node;
				end
			end

		end

	end
end

