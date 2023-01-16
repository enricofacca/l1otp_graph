function [forcing] = rect_forcing(topol,coord)
	nnode = size(coord,2);
	forcing = zeros(nnode,1);

	for i = 1:nnode
		x = coord(1,i);
		y = coord(2,i);

		if ( y>0.25 && y<0.75)
			if ( x>1/8 && x<3/8)
				forcing(i) = 1;
			end
			if ( x>5/8 && x<7/8)
				forcing(i) = -1;
			end
		end
	end
	mass = sum(abs(forcing))/2;
	forcing = forcing/mass;	
end


