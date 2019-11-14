time0=clock();
format long;

out2 = fopen('time_energ.out','w')

Nx = 8;
Ny = 8;
NxNy= Nx*Ny;
dx = 1.0;
dy = 1.0; 


%--- Time integration parameters:
nstep =   100*60*60*24*2;
nprint=    50;
dtime = 1.0e-3;
ttime = 0.0;

%--- Material specific Parameters:
cr0=0.40;


Rconst=8.314472;
tempr=740;
RT=Rconst*tempr;

a = 0.2866e-9;
nn = (a/2)*(sqrt(3));
int_par = 18600 + 0.1*tempr;
%grcoef_cr = (fz/2)*(a^2)*ddGexddcr
%grcoef_cr = 4e-16; % Coeficiente gradiente de energia - Segundo 2.1 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
Hmix = cr0*(1-cr0)*int_par;
grcoef_cr = (2/3)*Hmix*nn;
L = 20500 - 9.68*tempr;
grcoef_cr = (1/6)*(a*a)*L; %O nde r é a distância interatomica livre de stress e mudança de composição: doi:10.3390/ma10121431
%grcoef_cr = (1/6)*(r*r)*L; %O nde r é a distância interatomica livre de stress e mudança de composição: doi:10.3390/ma10121431
%grcoef_cr = 0.5;
%grcoef_cr = 4e-16;

%--ParÃ¢metros de difusividade para mobilidade:
%--Segundo 2.3 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
QCr=3.08e5;
D0Cr=2e-5;
%DCr=(D0Cr*exp(-QCr/RT));
DCr=D0Cr;

QFe=2.94e5; 
D0Fe=1e-4;
%DFe=(D0Fe*exp(-QFe/RT));
DFe = D0Fe
%--- Mobilidade -> ameniza a velocidade segregação a medida que diminui. *Decomposição depende da difusão
MFe = DFe/RT;
MCr = DCr/RT;
mcoef_cr=cr0*(1.0-cr0)*(DCr/RT);
mcoef_cr=cr0*(1.0-cr0)*(cr0*MFe + (1-cr0)*MCr);

%equação de Darken: mobilidade quimica
%MFe = DFe/RT;
%MCr = DCr/RT;
%mcoef_cr = cr0*(1.0-cr0)*(cr0*MFe + (1-cr0)*MCr) %O nde r é a distância interatomica livre de stress e mudança de composição: doi:10.3390/ma10121431
%mcoef_cr = mcoef_cr*1e26;

%mcoef_cr = 1;

%---inicializaÃ§Ã£o da microestrutura com concentraÃ§Ã£o c0
[cr] = micro_ch_pre(Nx,Ny,cr0);

[kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy);

%--- EvoluÃ§Ã£o do sistema
for istep =1:nstep
	ttime = ttime +dtime;
	
  cr
input('press return to continue');
  crk = fft2(cr);
crk
input('press return to continue');
##	for i=1:Nx
##		for j=1:Ny
##			
##			jp=j+1;
##			jm=j-1;
##			
##			ip=i+1;
##			im=i-1;
##			
##			jp=j+1;
##			jm=j-1;
##			

##			ip=i+1;
##			im=i-1;
##			
##			if(im == 0)
##			 im=Nx;
##			end
##			if(ip == (Nx+1))
##			  ip=1;
##			end
##			
##			if(jm == 0) 
##			  jm = Ny;
##			end
##			
##			if(jp == (Ny+1))
##			  jp=1;
##			end
##			
##			hne=cr(ip,j);
##			hnw=cr(im,j);
##			hns=cr(i,jm);
##			hnn=cr(i,jp);
##			hnc=cr(i,j);
##			
##			lap_cr(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
##			
##		end 
##	end 

	%--- Derivada da energia livre
	[dgdcr] = Fe_Cr_free_energy(Nx,Ny,cr0,cr,Rconst,tempr);

  
	dgdcrk = fft2(dgdcr);
dgdcrk
input('press return to continue');
  
  
##	%--- Derivada parcial da energia livre em relaÃ§Ã£o a concentraÃ§Ã£o Cr
##	delFdelcrgrid = dgdcr - grcoef_cr*lap_cr;
##
##	%laplaciano da derivada parcial da energia livre em relaÃ§Ã£o a concentraÃ§Ã£o Cr
##	for i=1:Nx
##		for j=1:Ny
##			
##			jp=j+1;
##			jm=j-1;
##			
##			ip=i+1;
##			im=i-1;
##			
##			jp=j+1;
##			jm=j-1;
##			
##			ip=i+1;
##			im=i-1;
##			
##			if(im == 0)
##			 im=Nx;
##			end
##			if(ip == (Nx+1))
##			  ip=1;
##			end
##			
##			if(jm == 0) 
##			  jm = Ny;
##			end
##			
##			if(jp == (Ny+1))
##			  jp=1;
##			end
##			
##			hne=delFdelcrgrid(ip,j);
##			hnw=delFdelcrgrid(im,j);
##			hns=delFdelcrgrid(i,jm);
##			hnn=delFdelcrgrid(i,jp);
##			hnc=delFdelcrgrid(i,j);
##			
##			lap_delFdelcr(i,j) =(hnw + hne + hns + hnn -4.0*hnc)/(dx*dy);
##			
##		end 
##	end 
	

	%--- integraÃ§Ã£o discreta no tempo
  crk = (crk - dtime*mcoef_cr.*dgdcrk.*k2) ./ ...
      (1.0 + dtime*grcoef_cr*k4.*mcoef_cr);
      
crk
input('press return to continue');
  
##	for i=1:Nx
##		for j=1:Ny
##  crk(i,j) = (crk(i,j) - dtime*mcoef_cr*dgdcrk(i,j)*k2(i,j)) / (1.0 + dtime*grcoef_cr*k4(i,j)*mcoef_cr);
##		
##		if(cr(i,j) >= 0.9999);
##			cr(i,j)= 0.9999;
##		end
##		if(cr(i,j) < 0.00001);
##			cr(i,j) = 0.00001;
##		end
##		
##		end 
##	end 
	
  cr = real(ifft2(crk));
cr
input('press return to continue');
 inrange =(cr >= 0.9999);
cr(inrange) = 0.9999;
inrange =(cr < 0.00001);
cr(inrange) = 0.00001;

##
## 	for i=1:Nx
##		for j=1:Ny
##
##		
##		if(cr(i,j) >= 0.9999);
##			cr(i,j)= 0.9999;
##		end
##		if(cr(i,j) < 0.00001);
##			cr(i,j) = 0.00001;
##		end
##		
##		end 
##	end 
	
	
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
