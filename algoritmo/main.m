
time0=clock();
format long;

out2 = fopen('time_energ.out','w')

Nx = 64;
Ny = 64;
NxNy= Nx*Ny;
dx = 1.0;
dy = 1.0; 


%--- Time integration parameters:
nstep =   5000;
nprint=    50;
dtime = 1.0e-2;
ttime = 0.0;

%--- Material specific Parameters:
cr0=0.35;


Rconst=8.314472;
tempr=723.0;
RT=Rconst*tempr;

a = 0.2866e-9;
nn = (a/2)*(sqrt(3));
int_par = 18600 + 0.1*tempr;
%grcoef_cr = (fz/2)*(a^2)*ddGexddcr
%grcoef_cr = 4e-16; %Coeficiente gradiente de energia - Segundo 2.1 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
Hmix = cr0*(1-cr0)*int_par;
grcoef_cr = (2/3)*Hmix*nn;



%--Parâmetros de difusividade para mobilidade:
%--Segundo 2.3 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
QCr=3.08e5; %Valor dummy - falta alterar com os dados do Thermo-Calc
D0Cr=0.2;
DCr=(D0Cr*exp(-QCr/RT));
%--- Mobilidade
mcoef_cr=cr0*(1.0-cr0)*(DCr/RT);

%---inicialização da microestrutura com concentração c0
[cr] = micro_ch_pre(Nx,Ny,cr0);

%--- Evolução do sistema
for istep =1:nstep
	ttime = ttime +dtime;
	

	for i=1:Nx
		for j=1:Ny
			
			jp=j+1;
			jm=j-1;
			
			ip=i+1;
			im=i-1;
			
			jp=j+1;
			jm=j-1;
			
			ip=i+1;
			im=i-1;
			
			if(im == 0)
			 im=Nx;
			end
			if(ip == (Nx+1))
			  ip=1;
			end
			
			if(jm == 0) 
			  jm = Ny;
			end
			
			if(jp == (Ny+1))
			  jp=1;
			end
			
			hne=cr(ip,j);
			hnw=cr(im,j);
			hns=cr(i,jm);
			hnn=cr(i,jp);
			hnc=cr(i,j);
			
			lap_cr(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
			
		end 
	end 

	%--- Derivada da energia livre
	[dgdcr] = Fe_Cr_free_energy(Nx,Ny,cr0,cr,Rconst,tempr);

	%--- Derivada parcial da energia livre em relação a concentração Cr
	delFdelcrgrid = dgdcr - grcoef_cr*lap_cr;

	%laplaciano da derivada parcial da energia livre em relação a concentração Cr
	for i=1:Nx
		for j=1:Ny
			
			jp=j+1;
			jm=j-1;
			
			ip=i+1;
			im=i-1;
			
			jp=j+1;
			jm=j-1;
			
			ip=i+1;
			im=i-1;
			
			if(im == 0)
			 im=Nx;
			end
			if(ip == (Nx+1))
			  ip=1;
			end
			
			if(jm == 0) 
			  jm = Ny;
			end
			
			if(jp == (Ny+1))
			  jp=1;
			end
			
			hne=delFdelcrgrid(ip,j);
			hnw=delFdelcrgrid(im,j);
			hns=delFdelcrgrid(i,jm);
			hnn=delFdelcrgrid(i,jp);
			hnc=delFdelcrgrid(i,j);
			
			lap_delFdelcr(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
			
		end 
	end 
	

	%--- integração discreta no tempo
	for i=1:Nx
		for j=1:Ny
		cr(i,j) = cr(i,j) + dtime*mcoef_cr*(dgdcr(i,j) - grcoef_cr*lap_delFdelcr(i,j));
		
		if(cr(i,j) >= 0.9999);
			cr(i,j)= 0.9999;
		end
		if(cr(i,j) < 0.00001);
			cr(i,j) = 0.00001;
		end
		
		end 
	end 
	
	
	if((mod(istep,nprint) == 0) || (istep == 1) )
	
	fprintf('done step: %5d\n',istep);
	
	
	%--- write vtk file:
	
	%-- calculate total energy
	
	[energ] = calculate_energ(Nx,Ny,cr,grcoef_cr);
	
	fprintf(out2,'%14.6e %14.6e\n',ttime,energ);
	
	write_vtk_grid_values(Nx,Ny,dx,dy,istep,cr);  
	
	end %if 
end %istep

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);
