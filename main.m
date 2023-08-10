fclose all;
close all;
clear;
clc;

%	input parameters for theoretical computations (constant)
%	solution parameters
soln =										containers.Map('UniformValues',false);
soln('nr') =								501;
soln('nx') =								501;
soln('nz') =								501;
soln('scheme') =							0;
soln('tol') =								1e-10;
soln('maxiter') =							2500;
soln('adapt') =								0.75;
soln('calcflow') =							2;
soln('calcdfrm') =							1;
soln('tolpert') =							0.1;
soln('pert') =								0;
soln('head') =								1;	
soln('plotall') =							0;
soln('vistog') =							'on';
soln('flowisinit') =						0;
soln('centralsmooth') =						0;
soln('plotsrhm') =							1;
%	solid material properties
solid =										containers.Map('UniformValues',false);
solid('Ey_0') =								-0.00e6;
solid('Ey_1') =								-410.0e6;
solid('Ey_k') =								2.235;
solid('nupois_0') =							0.49;
solid('nupois_1') =							0.00;
solid('nupois_k') =							0.00;
%	fluid material properties
fluid =										containers.Map('UniformValues',false);
fluid('rh_dens') =							1.0e3;
fluid('mu_visc') =							1.0e-3;
fluid('n_visc') =							1.0;
%	geometric parameters
geom =										containers.Map('UniformValues',false);
geom('A') =									0.5e-2;
geom('L') =									5.0e-2;
geom('L_tygon') =							85.0e-2;
geom('a_tygon') =							300.0e-6;
geom('inf_bulk') =							1;
% constants for experimental data
zmax =										geom('L');
batch =										'consol';
type =										'glassmidthin';
npart =										1*1+18*0;
nseg =										npart+2;
midseg =									'hard';

% pre-processing
if (strcmp(type,'glass') == 1)
	nchamber =								1;
	solid('ktanh') =						0.325;
	solid('clfit') =						2;
elseif (strcmp(type,'PMMA') == 1)
	nchamber =								3;
	solid('ktanh') =						1.00;
	solid('clfit') =						1;
elseif (strcmp(type,'glassmidthin') == 1)
	nchamber =								2;
	solid('ktanh') =						0.325;
	solid('clfit') =						2;
end
outfolder =									[batch,'/',type,'/'];
infolder =									['../res/',outfolder];

% reading and analysing experimental data
if ((strcmp(type,'glass') == 1) || (strcmp(type,'PMMA') == 1))
	nz =									5;
	z =										linspace(0,zmax,nz);
elseif (strcmp(type,'glassmidthin') == 1)
	nz =									7;
	z =										zeros(1,nz);
	z(1:5) =								linspace(0,zmax,nz-2);
	z(6:7) =								z(4:5);
	z(4) =									z(3);
	z(3) =									z(4)-0.25e-2;
	z(5) =									z(4)+0.25e-2;
end
for ichamber = 0:0
	nflow =									5;
	if (strcmp(type,'glass'))
		if (ichamber ~= 0)
			nrun =							3;
		else
			nrun =							6;
		end
		if (npart == 1)
			cl_seg =						[5,15,25];
		else
			cl_seg =						[5,10,15,20];
		end
	elseif (strcmp(type,'PMMA'))
		if (ichamber ~= 0)
			nrun =							5;
		else
			nrun =							15;
		end
		if (npart == 1)
			cl_seg =						[5,15,25];
		else
			cl_seg =						[5,10,15,20];
		end
	elseif (strcmp(type,'glassmidthin'))
		if (ichamber ~= 0)
			nrun =							3;
		else
			nrun =							6;
		end
		cl_seg =							zeros(1,20);
		if (strcmp(midseg,'hard'))
			cl_seg(1:8) =					25;
			cl_seg(9:10) =					5;
			cl_seg(11:20) =					25;
		elseif (strcmp(midseg,'soft'))
			cl_seg(1:8) =					5;
			cl_seg(9:10) =					25;
			cl_seg(12:20) =					5;
		end
	end
%	setting cross-linking ratios along channel length
	z_seg =									linspace(0,zmax,nseg);
	z_theo_cl =								linspace(0,zmax,soln('nz'));
	if ((strcmp(type,'PMMA')) || (strcmp(type,'glass')))
		if (solid('clfit') == 1)
			pol =							polyfit(z_seg,cl_seg,nseg-1);
			cl =							pchip(z_seg,cl_seg,z_theo_cl);	
		elseif (solid('clfit') == 2)
			l =								zeros(1,soln('nz'));
			npin =							2*nseg+1;
			cl(1:floor(soln('nz')/npin)) =	cl_seg(1);
			cl(floor(soln('nz')/npin)+ceil((nseg-1)*ceil(soln('nz')/nseg)):soln('nz')) =	cl_seg(nseg);
			for ipartn = 1:nseg-1
				ileft =						floor(soln('nz')/npin)+ceil((ipartn-1)*ceil(soln('nz')/nseg));
				iright =					floor(soln('nz')/npin)+ceil((ipartn)*ceil(soln('nz')/nseg));
				for iz = ileft:iright
					cl(iz) =				cl_seg(ipartn)+(cl_seg(ipartn+1)-cl_seg(ipartn))*0.5*(1+ ...
											tanh(((2*z_theo_cl(iz)-z_theo_cl(iright)-z_theo_cl(ileft)))/ ...
											(solid('ktanh')*(z_theo_cl(iright)-z_theo_cl(ileft)))));
				end
			end
		end
	elseif (strcmp(type,'glassmidthin'))
		ileft =								((soln('nz')-1)/2)+1-((soln('nz')-1)/10);
		iright =							((soln('nz')-1)/2)+1+((soln('nz')-1)/10);
		imid =								((soln('nz')-1)/2)+1;
		nsmooth =							((soln('nz')-1)/10);
		if (strcmp(midseg,'hard'))
			clmid =							5;
			clout =							25;
		elseif (strcmp(midseg,'soft'))
			clmid =							25;
			clout =							5;
		end
		cl =								clout*ones(1,soln('nz'));
		for iz = 1:imid
			cl(iz) =						clout+(clmid-clout)*0.5*(1+ ...
											tanh((z_theo_cl(iz)-z_theo_cl(ileft))/ ...
											(solid('ktanh')*0.5*(z_theo_cl(iright)-z_theo_cl(ileft)))));
		end
		for iz = imid:soln('nz')
			cl(iz) =						clmid+(clout-clmid)*0.5*(1+ ...
											tanh((z_theo_cl(iz)-z_theo_cl(iright))/ ...
											(solid('ktanh')*0.5*(z_theo_cl(iright)-z_theo_cl(ileft)))));
		end
		ksmooth =							(1.0/(2*(z_theo_cl(imid+nsmooth)-z_theo_cl(imid))))* ...
											mean([ ...
											((cl(imid+nsmooth+1)-cl(imid+nsmooth))/ ...
											(z_theo_cl(imid+nsmooth+1)-z_theo_cl(imid+nsmooth))), ...
											((cl(imid+nsmooth+2)-cl(imid+nsmooth+1))/ ...
											(z_theo_cl(imid+nsmooth+2)-z_theo_cl(imid+nsmooth+1))), ...
											((cl(imid+nsmooth+3)-cl(imid+nsmooth+2))/ ...
											(z_theo_cl(imid+nsmooth+3)-z_theo_cl(imid+nsmooth+2))), ...
											]);
		clmid =								cl(imid+nsmooth)-ksmooth*((z_theo_cl(imid+nsmooth)-z_theo_cl(imid))^2);
		if (soln('centralsmooth') == 1)
			for ismooth=imid-nsmooth:imid+nsmooth
				cl(ismooth) =				clmid+ksmooth*((z_theo_cl(ismooth)-z_theo_cl(imid))^2);
			end
		end
	end
	solid('cl') =							cl;
%	loading experimental data from file
	if (((strcmp(type,'glass')) == 1) || ((strcmp(type,'PMMA')) == 1))
		inpt =								load([infolder,num2str(nseg),'part_',num2str(ichamber),'.csv']);
	elseif (strcmp(type,'glassmidthin'))
		inpt =								load([infolder,midseg,'mid_',num2str(ichamber),'.csv']);
	end
%	converting diameter to radius and scaling experimental data to SI units
	inpt(:,2:nz+1) =						((1.0e-6)/2.0)*inpt(:,2:nz+1);
	inpt(:,1) =								(((1.0e-6)*(1.0e-3))/60.0)*inpt(:,1);
%	obtaining initial undeformed radius
%	pre-allocating
	flow =									zeros(1,nflow);
	init_r =								zeros(nz,nrun);
	r =										zeros(nz,nflow,nrun);		
	r_avg =									zeros(nz,nflow);
	r_err =									zeros(nz,nflow);
	r_err1 =								zeros(nz,nflow);
	r_sc =									zeros(nz,nflow,nrun);		
	r_sc_avg =								zeros(nz,nflow);
	r_sc_err =								zeros(nz,nflow);
	r_sc_err1 =								zeros(nz,nflow);
	dlt =									zeros(nz,nflow,nrun);		
	dlt_avg =								zeros(nz,nflow);	
	dlt_err =								zeros(nz,nflow);
	dlt_err1 =								zeros(nz,nflow);
	dlt_sc =								zeros(nz,nflow,nrun);		
	dlt_sc_avg =							zeros(nz,nflow);	
	dlt_sc_err =							zeros(nz,nflow);
	dlt_sc_err1 =							zeros(nz,nflow);
	run =									zeros(nflow,nrun);
	crun =									nrun*ones(1,nflow);
	r_theo =								zeros(soln('nr'),soln('nz'),nflow);
	x_theo =								zeros(soln('nx'),nflow);
	z_theo =								zeros(soln('nz'),nflow);
	Ey_theo =								zeros(soln('nz'),nflow);
	nuPois_theo =							zeros(soln('nz'),nflow);
	p_theo =								zeros(soln('nz'),nflow);
	p_f_theo =								zeros(soln('nz'),nflow);
	p_f_c_theo =							zeros(soln('nz'),nflow);
	dlt_theo =								zeros(soln('nz'),nflow);
	vr_theo =								zeros(soln('nr'),soln('nz'),nflow);
	vz_theo =								zeros(soln('nr'),soln('nz'),nflow);
	sr_theo =								zeros(7,soln('nr'),soln('nz'),nflow);
	ux_theo =								zeros(soln('nx'),soln('nz'),nflow);
%	analysing experimental data
	for irun = 1:nrun		
		for iz = 1:nz
			init_r(iz,irun) =				inpt(irun,iz+1);
		end
	end
	if (soln('flowisinit') == 1)
		for iflow = 1:2
			for irun = 1:nrun	
				for iz = 1:nz
					init_r(iz,irun) =		init_r(iz,irun)+inpt(iflow*nrun+irun,iz+1);
				end
			end
		end
		init_r = init_r/3;
	end
	for iflow = 1:nflow
		flow(iflow) =						inpt(iflow*nrun+1,1);
		for irun = 1:nrun	
			for iz = 1:nz
				r(iz,iflow,irun) =			inpt(iflow*nrun+irun,1+iz);
				r_sc(iz,iflow,irun) =		r(iz,iflow,irun)/init_r(iz);
				r_avg(iz,iflow) =			r_avg(iz,iflow)+r(iz,iflow,irun)/nrun;
				r_sc_avg(iz,iflow) =		r_sc_avg(iz,iflow)+r_sc(iz,iflow,irun)/nrun;
				dlt(iz,iflow,irun) =		r(iz,iflow,irun)-init_r(iz);
				dlt_sc(iz,iflow,irun) =		r_sc(iz,iflow,irun)-1.0;
				dlt_avg(iz,iflow) =			dlt_avg(iz,iflow)+dlt(iz,iflow,irun)/nrun;
				dlt_sc_avg(iz,iflow) =		dlt_sc_avg(iz,iflow)+dlt_sc(iz,iflow,irun)/nrun;
				run(iflow,irun) =			inpt(1+(iflow-1)*nrun+irun,nz+2);
			end
		end
		for iz = 1:nz
			r_err(iz,iflow) =				sum((r(iz,iflow,:)-r_avg(iz,iflow)).^2)/nflow;
			r_sc_err(iz,iflow) =			sum((r_sc(iz,iflow,:)-r_sc_avg(iz,iflow)).^2)/nflow;
			dlt_err(iz,iflow) =				sum((dlt(iz,iflow,:)-dlt_avg(iz,iflow)).^2)/nflow;
			dlt_sc_err(iz,iflow) =			sum((dlt_sc(iz,iflow,:)-dlt_sc_avg(iz,iflow)).^2)/nflow;
			r_err1(iz,iflow) =				max(abs(r(iz,iflow,:)-r_avg(iz,iflow)));
			r_sc_err1(iz,iflow) =			max(abs(r_sc(iz,iflow,:)-r_sc_avg(iz,iflow)));
			dlt_err1(iz,iflow) =			max(abs(dlt(iz,iflow,:)-dlt_avg(iz,iflow)));
			dlt_sc_err1(iz,iflow) =			max(abs(dlt_sc(iz,iflow,:)-dlt_sc_avg(iz,iflow)));
		end
		disp([	'Expr: batch, type, seg, chamber - ',batch,', ', type,', ', ...
				num2str(nseg),', ', num2str(ichamber),', and flow rate ',num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0))]);
%		geometric constants for the particular setup for theoretical calculations
		geom('a') =							mean(mean(init_r));
		geom('Q') =							flow(iflow);
		geom('dlt_sc_0') =					dlt_sc_avg(1,iflow);	
		[r_theo(:,:,iflow), x_theo(:,iflow), z_theo(:,iflow), p_f_theo(:,iflow), p_f_c_theo(:,iflow), p_theo(:,iflow), ...
		dlt_theo(:,iflow), vr_theo(:,:,iflow), vz_theo(:,:,iflow), sr_theo(:,:,:,iflow), ux_theo(:,:,iflow), Ey_theo(:,iflow), ...
		nuPois_theo(:,iflow)] =				theo(soln,solid,fluid,geom);
		disp([	'Theo: batch, type, seg, chamber - ',batch,', ', type,', ', ...
				num2str(nseg),', ', num2str(ichamber),', and flow rate ',num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0))]);
	end
%	noting how many runs were done for the specific flow-rate
	for iflow = 1:nflow
		for irun = 1:nrun
			if (ismember(irun,run(iflow,:)) == 0)
				crun(iflow) =				crun(iflow)-1;
			end
		end
	end
	if ((strcmp(type,'PMMA') == 1) || (strcmp(type,'glass') == 1))
		setup =	[num2str(nseg),'_',num2str(ichamber)];
	else
		setup =	[midseg,'_',num2str(ichamber)];
	end
%	plotting deflected radius with the original for reference
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,(dlt_theo(:,iflow)+geom('a'))*1e6,'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(r_avg(:,iflow))*1e6,flip(r_err(:,iflow))*1e6,'o','color',col1(1:3,iflow),'LineWidth',2);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(r(:,iflow,irun))*1e6,'s','LineWidth',1);
			end
		end
	end
	plot(ax1,z*1e2,flip(init_r)*1e6,'-','Color',[0.5,0.5,0.5],'LineWidth',2);
	plot(ax1,z_theo(:,1)*1e2,geom('a')*ones(1,length(z_theo(:,1)))*1e6,':','Color',[0.5,0.5,0.5],'LineWidth',1);
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/((((1.0e-6)*(1.0e-3)))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deformed radius ($\mu$m)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_r'],'fig');
	close(fig1);
%	plotting deflection
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,dlt_theo(:,iflow)*1e6,'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(dlt_avg(:,iflow))*1e6,flip(dlt_err(:,iflow))*1e6,'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(dlt(:,iflow,irun))*1e6,'o','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deflection ($\mu$m)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_dlt'],'fig');
	close(fig1);
%	plotting scaled deflected radius
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,(dlt_theo(:,iflow)+geom('a'))/geom('a'),'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(r_sc_avg(:,iflow)),flip(r_sc_err(:,iflow)),'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(r_sc(:,iflow,irun)),'s','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deformed radius (scaled)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_r_sc'],'fig');
	close(fig1);
%	plotting scaled deflection
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,dlt_theo(:,iflow)/geom('a'),'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(dlt_sc_avg(:,iflow)),flip(dlt_sc_err(:,iflow)),'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(dlt_sc(:,iflow,irun)),'s','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deflection (scaled)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_dlt_sc'],'fig');
	close(fig1);
%	plotting deflected radius with the original for reference
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,(dlt_theo(:,iflow)+geom('a'))*1e6,'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(r_avg(:,iflow))*1e6,flip(r_err1(:,iflow))*1e6,'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(r(:,iflow,irun))*1e6,'s','LineWidth',1);
			end
		end
	end
	plot(ax1,z*1e2,flip(init_r)*1e6,'-','Color',[0.5,0.5,0.5],'LineWidth',2);
	plot(ax1,z_theo(:,1)*1e2,geom('a')*ones(1,length(z_theo(:,1)))*1e6,':','Color',[0.5,0.5,0.5],'LineWidth',1);
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/((((1.0e-6)*(1.0e-3)))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deformed radius ($\mu$m)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_r1'],'fig');
	close(fig1);
%	plotting deflection
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,dlt_theo(:,iflow)*1e6,'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(dlt_avg(:,iflow))*1e6,flip(dlt_err1(:,iflow))*1e6,'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(dlt(:,iflow,irun))*1e6,'o','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deflection ($\mu$m)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_dlt1'],'fig');
	close(fig1);
%	plotting scaled deflected radius
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,(dlt_theo(:,iflow)+geom('a'))/geom('a'),'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(r_sc_avg(:,iflow)),flip(r_sc_err1(:,iflow)),'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(r_sc(:,iflow,irun)),'s','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deformed radius (scaled)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_r_sc1'],'fig');
	close(fig1);
	
%	plotting scaled deflection
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	col1 = zeros(3,nflow);
	for iflow = 1:nflow
		plt = plot(ax1,z_theo(:,iflow)*1e2,dlt_theo(:,iflow)/geom('a'),'-','LineWidth',3);
		col1(1:3,iflow) = plt.Color();
	end
	for iflow = 1:nflow
		errorbar(ax1,z*1e2,flip(dlt_sc_avg(:,iflow)),flip(dlt_sc_err1(:,iflow)),'o','color',col1(1:3,iflow),'LineWidth',3);
	end
	if (soln('plotall') == 1)
		for iflow = 1:nflow
			for irun = 1:nrun
				plot(ax1,z*1e2,flip(dlt_sc(:,iflow,irun)),'s','LineWidth',1);
			end
		end
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf deflection (scaled)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_dlt_sc1'],'fig');
	close(fig1);
	
%	plotting outlet pressure
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	plot(ax1,flow(:)/(((1.0e-6)*(1.0e-3))/60.0),p_theo(end,:)/1e6,'-','LineWidth',3);
	hold on;
	plot(ax1,	flow(:)/(((1.0e-6)*(1.0e-3))/60.0),(((2.0*fluid('mu_visc')*geom('L_tygon'))/geom('a_tygon'))* ...
				((3.0+(1.0/fluid('n_visc')))^fluid('n_visc'))*((flow(:)/(pi*geom('a_tygon')^3)).^fluid('n_visc')))/1e6, ...
				'o','LineWidth',3);
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf flow rate ($\mu$L/min)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf outlet pressure (MPa)','interpreter','latex','FontSize',25);
	saveas(fig1,[outfolder,setup,'_p_out'],'fig');
	close(fig1);
	
%	plotting scaled pressure
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	for iflow = 1:nflow
		plot(ax1,z_theo(:,iflow)*1e2,p_theo(:,iflow)/p_theo(end,iflow),'-','LineWidth',3);
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf pressure (scaled)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_p_sc'],'fig');
	close(fig1);
	
%	plotting absolute pressure (gauged with atmospheric)
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	for iflow = 1:nflow
		plot(ax1,z_theo(:,iflow)*1e2,p_theo(:,iflow)*1e-6,'-','LineWidth',3);
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf pressure (guaged with atmospheric) (MPa)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_pw'],'fig');
	close(fig1);
	
%	plotting absolute pressure (gauged with outlet)
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	for iflow = 1:nflow
		plot(ax1,z_theo(:,iflow)*1e2,p_f_theo(:,iflow)*1e-6,'-','LineWidth',3);
	end
	ent_lgnd = {};
	for iflow = 1:nflow
		ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
	end
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf pressure (guaged with outlet) (MPa)','interpreter','latex','FontSize',25);
	lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
	saveas(fig1,[outfolder,setup,'_pf'],'fig');
	close(fig1);
	
%	plotting channel pressure drop compared to that for rigid channel
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	plot(ax1,flow(:)/(((1.0e-6)*(1.0e-3))/60.0),p_f_c_theo(1,:),'-','LineWidth',3);
	hold on;
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf flow rate ($\mu$L/min)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'$\Delta p/\Delta p_{\rm rigid}$','interpreter','latex','FontSize',25);
	saveas(fig1,[outfolder,setup,'_p_f_c_out'],'fig');
	close(fig1);
%	plotting channel pressure drop compared to that for rigid channel
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	plot(ax1,flow(:)/(((1.0e-6)*(1.0e-3))/60.0),p_f_theo(1,:)./p_theo(end,:),'-','LineWidth',3);
	hold on;
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf flow rate ($\mu$L/min)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'$\Delta p/p_{\rm out}$','interpreter','latex','FontSize',25);
	saveas(fig1,[outfolder,setup,'_p_f_by_p_w'],'fig');
	close(fig1);
	if (soln('calcflow') ~= 0)
%		plotting extensional rate at centerline
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,z_theo(:,iflow)*1e2,permute(sr_theo(3,soln('nr'),:,iflow),[3,1,2]),'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf extensional rate at centerline (m/m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_exrcl'],'fig');
		close(fig1);
%		plotting shear rate at wall
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,z_theo(:,iflow)*1e2,permute(sr_theo(7,soln('nr'),:,iflow),[3,1,2]),'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf shear rate at wall (m/m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_shrw'],'fig');
		close(fig1);
%		plotting shear stress at wall
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,z_theo(:,iflow)*1e2,fluid('mu_visc')*permute(sr_theo(7,soln('nr'),:,iflow),[3,1,2])*1e-3,'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf shear stress at wall kPa','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_shsw'],'fig');
		close(fig1);
%		plotting axial velocity field at channel center
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,vz_theo(:,(soln('nz')-1)/2+1,iflow)*1e6,r_theo(:,(soln('nz')-1)/2+1,iflow)*1e6,'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xl1 = xlabel(	ax1,'\bf axial velocity ($\mu$m/s)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf distance along r ($\mu$m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_vz'],'fig');
		close(fig1);
%		plotting radial velocity field at channel center
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,vr_theo(:,(soln('nz')-1)/2+1,iflow)*1e6,r_theo(:,(soln('nz')-1)/2+1,iflow)*1e6,'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		xl1 = xlabel(	ax1,'\bf radial velocity ($\mu$m/s)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf distance along r ($\mu$m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_vr'],'fig');
		close(fig1);
	end
	if (soln('calcdfrm') ~= 0)
%		plotting deformation field at channel middle
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,ux_theo(:,(soln('nz')-1)/2+1,iflow)*1e6,x_theo(:,iflow)*1e6,'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'box','on');
		set(ax1,'xscale','log');
		set(ax1,'yscale','log');
		xl1 = xlabel(	ax1,'\bf radial deformation ($\mu$m)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf distance along x ($\mu$m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_ux'],'fig');
		close(fig1);
%		plotting radial strain rate at channel middle
		fig1 = figure('visible',soln('vistog'));
		set(fig1,'Position',[33.33,33.33,745,745]);
		ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
		hold on;
		for iflow = 1:nflow
			plot(ax1,ux_theo(:,(soln('nz')-1)/2+1,iflow)./x_theo(:,iflow),x_theo(:,iflow)*1e6,'-','LineWidth',3);
		end
		ent_lgnd = {};
		for iflow = 1:nflow
			ent_lgnd{iflow} =				num2str(flow(iflow)/(((1.0e-6)*(1.0e-3))/60.0));
		end
		set(ax1,'TickDir','out');
		set(ax1,'XMinorTick','on','YMinorTick','on');
		set(ax1,'FontSize',25);	
		set(ax1,'TickLabelInterpreter','latex');
		set(ax1,'xscale','log');
		set(ax1,'yscale','log');
		xl1 = xlabel(	ax1,'\bf radial strain rate (m/m)','interpreter','latex','FontSize',25);
		yl1 = ylabel(	ax1,'\bf distance along x ($\mu$m)','interpreter','latex','FontSize',25);
		lgnd = legend(ent_lgnd,'interpreter','latex','box','off','FontSize',25,'Location','best');
		saveas(fig1,[outfolder,setup,'_ux_x'],'fig');
		close(fig1);
	end
%	plotting channel Young's modulus
	fig1 = figure('visible',soln('vistog'));
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	iflow = 1;
	plt = plot(ax1,z_theo(:,iflow)*1e2,Ey_theo(:,iflow)*1e-6,'-','LineWidth',3);
	ent_lgnd = {};
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf Youngs modulus (MPa)','interpreter','latex','FontSize',25);
	saveas(fig1,[outfolder,setup,'_Ey'],'fig');
	close(fig1);
%	plotting channel cross-linking ratio
	fig1 = figure('visible',soln('vistog'));
	cl = solid('cl');
	set(fig1,'Position',[33.33,33.33,745,745]);
	ax1 =   axes(fig1,'Position',[0.2,0.15,0.775,0.825]);
	hold on;
	iflow = 1;
	plt = plot(ax1,z_theo(:,iflow)*1e2,cl,'-','LineWidth',3);
	ent_lgnd = {};
	set(ax1,'TickDir','out');
	set(ax1,'XMinorTick','on','YMinorTick','on');
	set(ax1,'FontSize',25);	
	set(ax1,'TickLabelInterpreter','latex');
	set(ax1,'box','on');
	xl1 = xlabel(	ax1,'\bf distance along z (cm)','interpreter','latex','FontSize',25);
	yl1 = ylabel(	ax1,'\bf Cross-Linking Ratio','interpreter','latex','FontSize',25);
	saveas(fig1,[outfolder,setup,'_cl'],'fig');
	close(fig1);
%	plotting heatmaps for the term values modulus
	if (soln('plotsrhm') == 1)
		for iterm = 1:7
			for iflow = 1:nflow
				fig1 = figure('visible',soln('vistog'));
				set(fig1,'Position',[33.33,33.33,745,745]);
				clf
				heatmap(flip(permute(sr_theo(iterm,:,:,iflow),[2,3,1,4]),1),'FontSize',28);
				set(gca,'grid','off'); set(gcf,'Position',[10.00,90.000,1500,250]); set(gca,'Position',[0.003,0.045,0.890,0.755]);
				saveas(fig1,[outfolder,'/heatmap_',setup,'_term_',num2str(iterm),'_flow_',num2str(iflow)],'fig');
				close(fig1);
			end
		end
	end
end
