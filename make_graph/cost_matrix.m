function [D] = cost_matrix(coord)
	% compute the distance matrix (the cost of the LP optimal transport problem)
	% for the quadrilater graph
	nnode = size(coord,2);
	D=zeros(nnode,nnode);
	for j = 1: nnode
		for i = 1: nnode
			D(i,j) = abs(coord(1,i)-coord(1,j)) + abs(coord(2,i)-coord(2,j));
		end
	end
	
