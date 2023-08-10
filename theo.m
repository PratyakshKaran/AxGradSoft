%-----------------------------------------------------------------------------------------------------------------------------------
% obtains the theoretical solution of deformation characteristics and flow dynamics for a circular microchannel in a thick medium 
% with axially-graded elasticity  
%-----------------------------------------------------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------------------------------------------------
function [rd, xd, zd, pfd, pfc, pwd, dltd, vrd, vzd, srd, uxd, Ey, nupois] = theo(soln,solid,fluid,geom)

%	theoretical solution

%	reading value from input
%	solution parameters
	nr =								soln('nr');						% grid points on r-axis
	nx =								soln('nx');						% grid points on x-axis
	nz =								soln('nz');						% grid points of z-axis
	scheme =							soln('scheme');					% scheme for two-way coupling solution (1=NR, 2=STL)
	tol =								soln('tol');					% error tolerance for iterative solvers
	maxiter =							soln('maxiter');				% maximum iterations for iterative solvers
	adapt =								soln('adapt');					% relaxation parameter for iterative solvers
	calcflow =							soln('calcflow');				% calculate the flow-field? (1=yes, 0=no)
	calcdfrm =							soln('calcdfrm');				% calculate the substrate-bulk deformation? (1=yes, 0=no)
	tolpert =							soln('tolpert');				% maximum value of small-parameter for perturbation solution
	pert =								soln('pert');					% solve with perturbation or not? (1=yes, 0=no)
	head =								soln('head');					% whether to calculate pressure head or match scaled 
																		% deflectionat outlet (1=yes, 0=no)
%	solid material properties
	cl =								solid('cl');					% cross-linking ratio (array)
	Ey_0 =								solid('Ey_0');					% first Young's modulus variation parameter
	Ey_1 =								solid('Ey_1');					% second Young's modulus variation parameter
	Ey_k =								solid('Ey_k');					% third Young's modulus variation parameter
	nupois_0 =							solid('nupois_0');				% first Poisson's ratio variation parameter
	nupois_1 =							solid('nupois_1');				% second Poisson's ratio variation parameter
	nupois_k =							solid('nupois_k');				% third Poisson's ratio variation parameter
%	fluid material properties
	mu_visc =							fluid('mu_visc');				% power law flow consistency index
	n_visc =							fluid('n_visc');				% power law flow behavior index
%	geometric parameters
	a =									geom('a');						% microchannel inner radius (when solving with head)
	A =									geom('A');						% microchannel outer radius (when solving with head)
	L =									geom('L');						% microchannel length (when solving with head)
	L_tygon =							geom('L_tygon');				% outlet drainage tube length (when solving with head)
	a_tygon =							geom('a_tygon');				% outlet drainage tube inner radius (when solving with head)
	Q =									geom('Q');						% flow-rate
	dlt_sc_0 =							geom('dlt_sc_0');				% outlet scaled deflection (when solving without head)
	inf_bulk =							geom('inf_bulk');				% is the substrate bulk practically infinite radially
%	calculating preliminary values
%	calculating axially varying substrate compliance values as per power law
	Ey =								Ey_0-Ey_1*cl.^-Ey_k;			
	nupois =							nupois_0-nupois_1*cl.^-nupois_k;
	lambd =								(Ey.*nupois)./((1.0+nupois).*(1.0-2*nupois));
	G =									Ey./(2*(1.0+nupois));
	G0 =								min(G);
	lambd0 =							min(lambd);
%	calculating simulation parameter values
	epsi =								a/L;
	gamm =								A/L;
	p0 =								8*(L_tygon/L)*((a/a_tygon)^4);
	M0 =								(-2*((3.0+(1.0/n_visc))^n_visc))/(1+p0);
	
	thet =								((gamm*mu_visc*(1+p0))/(epsi^2*(lambd0+2*G0)))*((Q/(pi*epsi^3*L^3))^n_visc);
	if (inf_bulk == 0)
		K =								(((epsi/(2*gamm))./(1.0+(1.0+lambd./G)*(epsi^2/gamm^2)))* ...
										(1-(epsi^2/gamm^2))*(lambd0+2*G0))./G;
	else
		K =								(epsi/(2*gamm))*((lambd0+2*G0)./G);
	end
	if (head == 1)
		p_0 =							(p0/(1+p0));
		dlt_sc_0 =						thet*(p_0*K(nz));
	end
	if (head ~= 1)
		disp('CAUTION: PRESSURE SCALING NOT FIXED FOR SCALED DEFLECTION PINNING MODE OF COMPUTATION');
	end
%	grid generation (non-dimensional)
	z =									linspace(0,1,nz);
	dz =								z(2)-z(1);

%	pre-allocation
	pw =								zeros( 1,nz);
	dlt =								zeros( 1,nz);

%	solution
	if ((pert == 1) && (thet <= tolpert))
%		asymoptotic leading order solution
		for iz = 1:nz
			pw(iz) =					(p0-2*((3+(1/n_visc))^n_visc)*(z(iz)-1))/(1+p0);
		end
		for iz = 1:nz
			dlt(iz) =					K(iz)*pw(iz);
		end
			disp(						['Asymptotic leading order solution with thet ',num2str(thet)]);
	else
%		non-linear solution
		Jac =							zeros(2*nz,2*nz);
		Res =                           zeros(2*nz,1);
		if (scheme == 1)
			Res(1:nz-1) =				M0;
			iz = nz;
			Res(iz) =					dlt_sc_0/(thet*K(iz));
			Res(nz+iz) =				dlt_sc_0/thet;
		end
		iz = nz;
		Jac(iz,iz) =					1.0;
		Jac(nz+iz,nz+iz) =				1.0;
		for iz = 1:nz-1
			Jac(nz+iz,nz+iz) =			1.0;
			if (scheme == 0)
				Jac(nz+iz,iz) =			-K(iz);
			end
		end
		iter = 0;
		errit = 1.0e1;
		while((errit > tol) && (iter < maxiter))
			pprev =						pw;
			dltprev =					dlt;
			if (scheme == 0)
				iz = 1;
				Res(iz) =				((1.0+thet*dlt(iz))^(1+3*n_visc))*((-3*pw(iz)+4*pw(iz+1)-pw(iz+2))/(2*dz))-M0;
				for iz = 2:nz-1
					Res(iz) =			((1.0+thet*dlt(iz))^(1+3*n_visc))*((-pw(iz-1)+pw(iz+1))/(2*dz))-M0;
				end
				iz = nz;
				Res(iz) =				pw(iz)-dlt_sc_0/(thet*K(iz));
			end
			iz = 1;
			Jac(iz,iz) =				(-3/(2*dz))*((1.0+thet*dlt(iz))^(1+3*n_visc));
			Jac(iz,iz+1) =				(4/(2*dz))*((1.0+thet*dlt(iz))^(1+3*n_visc));
			Jac(iz,iz+2) =				(-1/(2*dz))*((1.0+thet*dlt(iz))^(1+3*n_visc));
			for iz = 2:nz-1
				Jac(iz,iz+1) =			(1/(2*dz))*((1.0+thet*dlt(iz))^(1+3*n_visc));
				Jac(iz,iz-1) =			(-1/(2*dz))*((1.0+thet*dlt(iz))^(1+3*n_visc));
			end
			if (scheme == 0)
				iz = 1;
				Jac(iz,nz+iz) =			(1.0+3*n_visc)*thet*((1.0+thet*dlt(iz))^(1+3*n_visc))* ...
										((-3*pw(iz)+4*pw(iz+1)-pw(iz+2))/(2*dz));
				for iz = 2:nz-1
					Jac(iz,nz+iz) =		(1.0+3*n_visc)*thet*((1.0+thet*dlt(iz))^(1+3*n_visc))*((-pw(iz-1)+pw(iz+1))/(2*dz));
				end
			end
			if (scheme == 0)
				for iz = 1:nz-1
					Res(nz+iz) =		dlt(iz)-K(iz)*pw(iz);
				end
				iz = nz;
				Res(nz+iz) =			dlt(iz)-dlt_sc_0/thet;
			else
				for iz = 1:nz-1
					Res(nz+iz) =		K(iz)*pw(iz);
				end
			end
			Sol = Jac\Res;
			if (scheme == 0)
				for iz = 1:nz
					pw(iz) =			pw(iz)-adapt*Sol(iz);
					dlt(iz) =			dlt(iz)-adapt*Sol(nz+iz);
				end
			else
				for iz = 1:nz
					pw(iz) =			(1.0-adapt)*pw(iz)+adapt*Sol(iz);
					dlt(iz) =			(1.0-adapt)*dlt(iz)+adapt*Sol(nz+iz);
				end
			end
			if (scheme == 0)
				errit =					sqrt(sum(Sol.^2))/(2*nz);
			else
				errit =					sqrt(sum((pw-pprev).^2)+sum((dlt-dltprev).^2))/(2*nz);
			end
			disp(						['Scheme ',num2str(scheme),' iteration ',num2str(iter),' with error ',num2str(errit)]);
		end
	end

%	field variables
	pf =								(1+p0)*pw-p0;
	if (calcflow >= 1)
%		pre-allocating arrays
		r =								zeros(nr,nz);
		for iz = 1:nz
			r(:,iz) =					linspace(0,1+thet*dlt(iz),nr)';
		end
		vr =							zeros(nr,nz);
		vz =							zeros(nr,nz);
		sr =							zeros(nr,nz);
%		obtaining derivatives of pressure
		dpfdz =							zeros(3,nz);
		ddltdz =						zeros(2,nz);
		iz = 1;
		dpfdz(1,iz) =					((-3*pf(iz)+4*pf(iz+1)-pf(iz+2))/(2*dz));
		for iz = 2:nz-1
			dpfdz(1,iz) =				((-pf(iz-1)+pf(iz+1))/(2*dz));
		end
		iz = nz;
		dpfdz(1,iz) =					((3*pf(iz)-4*pf(iz-1)+pf(iz-2))/(2*dz));
		iz = 1;
		ddltdz(1,iz) =					((-3*dlt(iz)+4*dlt(iz+1)-dlt(iz+2))/(2*dz));
		for iz = 2:nz-1
			ddltdz(1,iz) =				((-dlt(iz-1)+dlt(iz+1))/(2*dz));
		end
		iz = nz;
		ddltdz(1,iz) =					((3*dlt(iz)-4*dlt(iz-1)+dlt(iz-2))/(2*dz));
		iz = 1;
		dpfdz(2,iz) =					((-3*dpfdz(1,iz)+4*dpfdz(1,iz+1)-dpfdz(1,iz+2))/(2*dz));
		for iz = 2:nz-1
			dpfdz(2,iz) =				((-dpfdz(1,iz-1)+dpfdz(1,iz+1))/(2*dz));
		end
		iz = nz;
		dpfdz(2,iz) =					((3*dpfdz(1,iz)-4*dpfdz(1,iz-1)+dpfdz(1,iz-2))/(2*dz));
		if (calcflow == 2)
			iz = 1;
			ddltdz(2,iz) =				((-3*ddltdz(1,iz)+4*ddltdz(1,iz+1)-ddltdz(1,iz+2))/(2*dz));
			for iz = 2:nz-1
				ddltdz(2,iz) =			((-ddltdz(1,iz-1)+ddltdz(1,iz+1))/(2*dz));
			end
			iz = nz;
			ddltdz(2,iz) =				((3*ddltdz(1,iz)-4*ddltdz(1,iz-1)+ddltdz(1,iz-2))/(2*dz));
			iz = 1;
			dpfdz(3,iz) =				((-3*dpfdz(2,iz)+4*dpfdz(2,iz+1)-dpfdz(2,iz+2))/(2*dz));
			for iz = 2:nz-1
				dpfdz(3,iz) =			((-dpfdz(2,iz-1)+dpfdz(2,iz+1))/(2*dz));
			end
			iz = nz;
			dpfdz(3,iz) =				((3*dpfdz(2,iz)-4*dpfdz(2,iz-1)+dpfdz(2,iz-2))/(2*dz));
		end	
%		obtaining axial velocity
		for iz = 1:nz
			for ir = 1:nr
				vz(ir,iz) =				(n_visc/(2^(1/n_visc)*(n_visc+1)))*((abs(dpfdz(1,iz)))^(1/n_visc))* ...
										((1+thet*dlt(iz))^((1+n_visc)/n_visc)-r(ir,iz)^((1+n_visc)/n_visc));
			end
		end
%		obtaining radial velocity
		for iz = 1:nz
			for ir = 1:nr
				vr(ir,iz) =				-(1.0/(2^((1+n_visc)/n_visc)))* ...
										(abs(dpfdz(1,nz)))^((1-n_visc)/n_visc)* ...
										((1.0/(1+n_visc))*dpfdz(2,nz)* ...
										(r(ir,iz)*(1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										((2*n_visc)/(1+3*n_visc))*r(ir,iz)^((1+2*n_visc)/n_visc)) + ...
										thet*r(ir,iz)*(1+thet*dlt(iz))^((1+n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz));
			end
		end
%		obtaining the strain-rate components and the aggregate second invariant strain rate	(dimensional values calculated here)	
		if (calcflow == 2)
			sr =						zeros(7,nr,nz);
			srd =						zeros(7,nr,nz);
			for ir = 1:nr
				for iz = 1:nz
					sr(1,ir,iz) =		-(1.0/(2^((1+n_visc)/n_visc)))* ...
										(abs(dpfdz(1,nz)))^((1-n_visc)/n_visc)* ...
										((1.0/(1+n_visc))*dpfdz(2,nz)* ...
										((1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										((2*(1+2*n_visc))/(1+3*n_visc))*r(ir,iz)^((1+n_visc)/n_visc)) + ...
										thet*(1+thet*dlt(iz))^((1+n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz));
					sr(2,ir,iz) =		-(1.0/(2^((1+n_visc)/n_visc)))* ...
										(abs(dpfdz(1,nz)))^((1-n_visc)/n_visc)* ...
										((1.0/(1+n_visc))*dpfdz(2,nz)* ...
										((1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										((2*n_visc)/(1+3*n_visc))*r(ir,iz)^((1+n_visc)/n_visc)) + ...
										thet*(1+thet*dlt(iz))^((1+n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz));
					sr(3,ir,iz) =		(1.0/(2^(1/n_visc)))* ...
										(abs(dpfdz(1,nz)))^((1-n_visc)/n_visc)* ...
										((1.0/(1+n_visc))*dpfdz(2,nz)* ...
										((1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										r(ir,iz)^((1+n_visc)/n_visc)) + ...
										thet*(1+thet*dlt(iz))^((1+n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz));
					if (n_visc ~= 1)
						sr(6,ir,iz) =	-((1-n_visc)/((2^((1+n_visc)/n_visc))*n_visc))*((dpfdz(1,iz))^((1-2*n_visc)/n_visc))* ...
										((1.0/(1+n_visc))*dpfdz(2,nz)* ...
										(r(ir,iz)*(1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										((2*n_visc)/(1+3*n_visc))*r(ir,iz)^((1+2*n_visc)/n_visc)) + ...
										thet*r(ir,iz)*(1+thet*dlt(iz))^((1+n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz))- ...
										((1-n_visc)/((2^((1+n_visc)/n_visc))*n_visc))*((dpfdz(1,iz))^((1-n_visc)/n_visc))* ...
										((1.0/(1+n_visc))*dpfdz(3,nz)* ...
										(r(ir,iz)*(1+thet*dlt(iz))^((1+n_visc)/n_visc)- ...
										((2*n_visc)/(1+3*n_visc))*r(ir,iz)^((1+2*n_visc)/n_visc)) + ...
										thet*r(ir,iz)*(1+thet*dlt(iz))^(1.0/n_visc)* ...
										(((1+n_visc)/n_visc)*dpfdz(2,iz)*ddltdz(1,iz)+dpfdz(1,iz)*ddltdz(2,iz))+ ...
										(((thet^2)*r(ir,iz))/n_visc)* ...
										(1+thet*dlt(iz))^((1-n_visc)/n_visc)*dpfdz(1,iz)*ddltdz(1,iz)^2);
					else
						sr(6,ir,iz) =	-(1/16)*(dpfdz(3,nz)*(2*r(ir,iz)*(1+thet*dlt(iz))^2-r(ir,iz)^3) + ...
										4*thet*r(ir,iz)*(1+thet*dlt(iz))*(2*dpfdz(2,iz)*ddltdz(1,iz)+dpfdz(1,iz)*ddltdz(2,iz))+ ...
										4*(thet^2)*r(ir,iz)*dpfdz(1,iz)*ddltdz(1,iz)^2);
					end
					sr(7,ir,iz) =		-(((r(ir,iz)/2)*dpfdz(1,iz))^(1.0/n_visc));
					sr(4,ir,iz) =		epsi^2*sr(6,ir,iz)+sr(7,ir,iz);
					sr(5,ir,iz) =		abs(sr(4,ir,iz));
					srd(1,ir,iz) =		epsi*(Q/(pi*a^3))*sr(1,ir,iz);
					srd(2,ir,iz) =		epsi*(Q/(pi*a^3))*sr(2,ir,iz);
					srd(3,ir,iz) =		epsi*(Q/(pi*a^3))*sr(3,ir,iz);
					srd(6,ir,iz) =		epsi^2*(Q/(pi*a^3))*sr(6,ir,iz);
					srd(7,ir,iz) =		(Q/(pi*a^3))*sr(7,ir,iz);
					srd(4,ir,iz) =		srd(6,ir,iz)+srd(7,ir,iz);
					srd(5,ir,iz) =		sqrt(abs(srd(1,ir,iz)*srd(2,ir,iz)+srd(2,ir,iz)*srd(3,ir,iz)+srd(3,ir,iz)*srd(1,ir,iz)- ...
										srd(4,ir,iz)^2));
				end
			end
		end
	if (calcdfrm == 1)
		x =								linspace(epsi/gamm,1,nx);
		ux =							zeros(nx,nz);
		if (inf_bulk == 0)
			for iz = 1:nz
				for ix = 1:nx
					ux(ix,iz) =			(((epsi^2)/(2*gamm^2))/(1+(1+(lambd(iz)/G(iz)))*(epsi^2/gamm^2)))* ...
										((1.0/x(ix))-x(ix))*((lambd0+2*G0)/G(iz))*pw(iz);
				end
			end
		else
			for iz = 1:nz
				for ix = 1:nx
					ux(ix,iz) =			(((epsi^2)*(lambd0+2*G0))/(2*(gamm^2)*G(iz)))*((1.0/x(ix))-x(ix))*pw(iz);
				end
			end
		end
	end

%   calculating dimensional values (dimensional values of strain rates calculated earlier)
	rd =								r*(epsi*L);
	xd =								x*(gamm*L);
	zd =								z*L;
	pfd =								(mu_visc/epsi)*((Q/(pi*epsi^3*L^3))^n_visc)*pf;
	pwd =								(mu_visc/epsi)*((Q/(pi*epsi^3*L^3))^n_visc)*(1+p0)*pw;
	pfc =								pfd/(((8*mu_visc*L)/a)*((Q/(pi*a^3))^n_visc));
	dltd =								thet*epsi*L*dlt;
	if (calcflow ~= 2)
		srd =							zeros(7,nr,nz);
	end
	if ((calcflow ~= 2) && (calcflow ~= 1))
		vr =							zeros(nr,nz);
		vz =							zeros(nr,nz);
	end
	if (calcdfrm == 0)
		ux =							zeros(nx,nz);
	end
	if ((calcflow == 1) || (calcflow == 2))
		vzd =							(Q/(pi*epsi^2*L^2))*vz;
		vrd =							(Q/(pi*epsi*L^2))*vr;
	else
		vzd =							0.0;
		vrd =							0.0;
	end
	if (calcdfrm == 1)
		uxd =							thet*epsi*L*ux;
	else
		uxd =							0.0;
	end
	
end
