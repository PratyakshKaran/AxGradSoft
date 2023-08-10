%		function computing derivative of variable defined over both coordinates
%		arguments are: var - variable of which derivative is computed; diff - differentiation type 'x', 'xx', 'y', 'yy', 'xy'; 
%		dirr - direction of derivative 'fd', 'cd', 'bd', with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy'
%		and 'xy'; ix - index of first coordinate x where derivative is to be obtained; iy - index of second coordinate y where 
%		derivative is to be obtained; dimn - dummy variable with 'xy' the only applicable option; x - x-grid; y - y-grid
function differential2d = differential2d(var,diff,dir,ix,iy,~,x,y)
	differential2d = 0.0;
	acoeff = a(diff,dir,ix,iy,x,y);
	for jx = -2+3:2+3
		for jy = -2+3:2+3
			if (acoeff(jx,jy) ~= 0.0)
				differential2d = differential2d + acoeff(jx,jy)*var(ix+jx-3,iy+jy-3);
			end
		end 
	end
end
