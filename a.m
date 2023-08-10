%		function defining coefficients of the appropriate boundaries for different derivatives
%		arguments are: diff - differentiation type 'x', 'xx', 'y', 'yy', 'xy'; dirr - direction of derivative 'fd', 'cd', 'bd', 
%		with only 'cd' being applicable for 2nd order derivatives i.e. 'xx', 'yy', 'xy'; ix - index of first coordinate x where 
%		derivative is to be obtained; iy - index of second coordinate y where derivative is to be obtained; x - x-grid; y - y-grid

function a = a(diff,dir,ix,iy,x,y)
	a = zeros(5,5);
	switch (diff)
	case ('')
		a(0+3,0+3) =	1.0;
	case ('x')
		switch (dir)
		case ('fd')
			dx =			x(ix+1)-x(ix);
			dxi =			x(ix+2)-x(ix+1);
			a(2+3,0+3) =    -(dx^2.0)/(dxi*dx*(dx+dxi));
			a(1+3,0+3) =    ((dx+dxi)^2.0)/(dxi*dx*(dx+dxi));
			a(0+3,0+3) =    -((dx+dxi)^2.0-dx^2.0)/(dxi*dx*(dx+dxi));
		case ('bd')
			idx =			x(ix)-x(ix-1);
			iidx =			x(ix-1)-x(ix-2);
			a(0+3,0+3) =    ((idx+iidx)^2.0-idx^2.0)/(iidx*idx*(idx+iidx));
			a(-1+3,0+3) =   -((idx+iidx)^2.0)/(iidx*idx*(idx+iidx));
			a(-2+3,0+3) =   (idx^2.0)/(iidx*idx*(idx+iidx));
		case ('cd')
			dx =			x(ix+1)-x(ix);
			idx =		    x(ix)-x(ix-1);
			a(1+3,0+3) =    idx/(2.0*dx*idx);
			a(0+3,0+3) =    (dx-idx)/(2.0*dx*idx);
			a(-1+3,0+3) =   -dx/(2.0*dx*idx);
		end
	case ('y')
		switch (dir)
		case ('fd')
			dy =			y(iy+1)-y(iy);
			dyj =		    y(iy+2)-y(iy+1);
			a(0+3,2+3) =    -(dy^2.0)/(dyj*dy*(dy+dyj));
			a(0+3,1+3) =    ((dy+dyj)^2.0)/(dyj*dy*(dy+dyj));
			a(0+3,0+3) =    -((dy+dyj)^2.0-dy^2.0)/(dyj*dy*(dy+dyj));
		case ('bd')
			jdy =			y(iy)-y(iy-1);
			jjdy =		    y(iy-1)-y(iy-2);
			a(0+3,0+3) =    ((jdy+jjdy)^2.0-jdy^2.0)/(jjdy*jdy*(jdy+jjdy));
			a(0+3,-1+3) =   -((jdy+jjdy)^2.0)/(jjdy*jdy*(jdy+jjdy));
			a(0+3,-2+3) =   (jdy^2.0)/(jjdy*jdy*(jdy+jjdy));
		case ('cd')
			dy =			y(iy+1)-y(iy);
			jdy =			y(iy)-y(iy-1);
			a(0+3,1+3) =    jdy/(2.0*dy*jdy);
			a(0+3,0+3) =    (dy-jdy)/(2.0*dy*jdy);
			a(0+3,-1+3) =   -dy/(2.0*dy*jdy);
		end
	case ('xx')
		a(-1+3,0+3) = 		2.0/((x(ix)-x(ix-1))*(x(ix+1)-x(ix-1)));
		a(1+3,0+3) = 		2.0/((x(ix+1)-x(ix))*(x(ix+1)-x(ix-1)));
		a(0+3,0+3) = 		-2.0/((x(ix+1)-x(ix))*(x(ix)-x(ix-1)));
	case ('yy')
		a(0+3,-1+3) = 		2.0/((y(iy)-y(iy-1))*(y(iy+1)-y(iy-1)));
		a(0+3,1+3) = 		2.0/((y(iy+1)-y(iy))*(y(iy+1)-y(iy-1)));
		a(0+3,0+3) = 		-2.0/((y(iy+1)-y(iy))*(y(iy)-y(iy-1)));
	case ('xy')
		dx =		      	x(ix+1)-x(ix);
		idx =				x(ix)-x(ix-1);
		dy =				y(iy+1)-y(iy);
		jdy =				y(iy)-y(iy-1);
		a(1+3,1+3) =    	(idx*jdy)/(4.0*idx*dx*jdy*dy);
		a(1+3,0+3) =    	(idx*(dy-jdy))/(4.0*idx*dx*jdy*dy);
		a(0+3,1+3) =    	((dx-idx)*jdy)/(4.0*idx*dx*jdy*dy);
		a(1+3,-1+3) =   	(-idx*dy)/(4.0*idx*dx*jdy*dy);
		a(0+3,0+3) =    	((dx-idx)*(dy-jdy))/(4.0*idx*dx*jdy*dy);
		a(-1+3,1+3) =   	(-dx*jdy)/(4.0*idx*dx*jdy*dy);
		a(0+3,-1+3) =   	(-(dx-idx)*dy)/(4.0*idx*dx*jdy*dy);
		a(-1+3,0+3) =   	(-dx*(dy-jdy))/(4.0*idx*dx*jdy*dy);
		a(-1+3,-1+3) =  	(dx*dy)/(4.0*idx*dx*jdy*dy);
	end
end
