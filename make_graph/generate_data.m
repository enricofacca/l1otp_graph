function [cost_mat, initial_density, final_density] = generate_data(ndiv)
	[topol, coord, weight] = quadrilater_graph(ndiv);
	cost_mat = cost_matrix(coord);
	forcing = rect_forcing(topol,coord); 
	initial_density = forcing;
	initial_density(forcing < 0) = 0;
	final_density = forcing;
	final_density(forcing > 0) = 0;
