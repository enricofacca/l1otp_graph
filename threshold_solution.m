function [solution] = threshold_solution(solution, problem, ctrl)

	ntdens = solution.ntdens;
	npot   = solution.npot;
	
  if (ctrl.selection == 1)
		% remove edges below threshold and disconneted nodes
		solution.ntdens_dirichlet=0;
		nactive=0; tdens_active=zeros(ntdens,1);
    for iedge = 1:ntdens
			if (solution.tdens(iedge) < ctrl.threshold_selection )
				solution.ntdens_dirichlet = solution.ntdens_dirichlet+1;
				solution.tdens_dirichlet(solution.ntdens_dirichlet) = iedge;
			else
				nactive = nactive + 1;
				tdens_active(nactive) = iedge;
			end
    end
		
		% add one if node is within an active edge
		selection_node = zeros(npot);
	  for iactive = 1 : nactive
      iedge = tdens_active(iactive);
			[nodes, lines, values] = find( problem.matrix_E(:,iedge) );
      selection_node(nodes)=selection_node(nodes)+1;
    end
    
		solution.npot_dirichlet=0;
    for inode = 1 : npot
			% if no edge gave contribution is a Dirichlet node
      if (selection_node(inode) == 0 )
				solution.npot_dirichlet = solution.npot_dirichlet + 1;
				solution.pot_dirichlet(solution.npot_dirichlet)=inode;
			end
		end
	end
end
