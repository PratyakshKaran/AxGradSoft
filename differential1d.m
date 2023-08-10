%		function computing derivative of variable defined over only one coordinate
%		arguments are: var - variable of which derivative is computed; diff - differentiation type 'x', 'xx', 'y', 'yy'; dirr - 
%		direction of derivative 'fd', 'cd', 'bd', with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy'; ix - 
%		index of first coordinate x where derivative is to be obtained; iy - index of second coordinate y where derivative is to be 
%		obtained; dimn - coordinate for which the variable is defined; x - x-grid; y - y-grid

function differential1d = differential1d (var,diff,dir,ix,iy,dimn,x,y)
	differential1d = 0.0;
	if (	((strcmp(diff,'x') || strcmp(diff,'xx')) && strcmp(dimn,'x')) || ...
			((strcmp(diff,'y') || strcmp(diff,'yy')) && strcmp(dimn,'y'))		)
		acoeff = a(diff,dir,ix,iy,x,y);
		if (dimn == 'x') 
			for jx = -2+3:2+3
				jy = 0+3;
				if (acoeff(jx,jy) ~= 0.0)
					differential1d = differential1d + acoeff(jx,jy)*var(ix+jx-3);
				end
			end
		else
			jx = 0+3;
			for jy = -2+3:2+3
				if (acoeff(jx,jy) ~= 0.0) 
					differential1d = differential1d + acoeff(jx,jy)*var(iy+jy-3);
				end
			end
		end
	end
end
