function [topol, coord, weight] = quadrilater_graph(ndiv)
	% create a graph of a regualr quadrilateral grid of the unit square
	% quadrilater_graph(3) will create
	% 
	% 7 - 8 - 9
	% |   |   |
	% 4 - 5 - 6
	% |   |   |
	% 1 - 2 - 3
	nnode = (ndiv)^2;
	nedge = 2*(ndiv-1)*ndiv;

	coord = zeros(2,nnode);
	topol = zeros(2,nedge);
	weight = zeros(nedge,1);
	h=1/(ndiv);
	weight(:) = h;
	
	x = linspace(h/2,1-h/2,ndiv);
	y = linspace(h/2,1-h/2,ndiv);
	[x,y] = meshgrid(x,y);

	coord(1,:) = reshape(x,1,[]);
	coord(2,:) = reshape(y,1,[]);

	sequence_h = zeros(2,ndiv-1);
	for i = 1 : ndiv-1;
		sequence_h(1,i) = i;
		sequence_h(2,i) = i+1;
	end
	sequence_h;

	sequence_v = zeros(2,ndiv);
	for i = 1:ndiv
		sequence_v(1,i) = i;
		sequence_v(2,i) = i+ndiv;
	end
	sequence_v;

	
	for j = 1:ndiv-1;
		start  = (j-1)*(ndiv-1)+(j-1)*ndiv+1;
		finish = (j-1)*ndiv+j*(ndiv-1);
		topol(:,start:finish) = (j-1)*(ndiv)+sequence_h;
		
		start  = (j-1)*ndiv + j*(ndiv-1) + 1;
		finish = (j-1)*ndiv + j*(ndiv-1) + ndiv;
		topol(:,start:finish) = (j-1)*(ndiv)+sequence_v;
	end

	j=ndiv;
	start  = (j-1)*(ndiv-1)+(j-1)*ndiv+1;
	finish = (j-1)*ndiv+j*(ndiv-1);
	topol(:,start:finish) = (j-1)*ndiv+sequence_h;
	topol(:,:);
