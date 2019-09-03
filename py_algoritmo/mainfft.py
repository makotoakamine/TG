import numpy as np
import math

def micro_ch_pre(Nx,Ny,c0):
	noise = 0.001
	return c0*np.ones((Nx,Ny)) + noise*(0.5*np.ones((Nx,Ny))-np.random.rand(Nx,Ny))

def prepare_fft(Nx,Ny,dx,dy):
  
	Nx21 = math.ceil(Nx/2)
	Ny21 = math.ceil(Ny/2)
	
	Nx2=Nx+2
	Ny2=Ny+2
	
	#--
	
	delkx=(2.0*pi)/(Nx*dx)
	delky=(2.0*pi)/(Ny*dy)
	
	#--
	fk1=np.arange(0,Nx
	fk1=(i-1)*delkx

	for i in range(0,Nx21-1):
		fk1=(i-1)*delkx
		kx(i)=fk1
		kx(Nx2-i)=-fk1
	
	for j in range(1,Ny21):
		fk2=(j-1)*delky
		ky(j)=fk2
		ky(Ny2-j)=-fk2
	
	#---
	
	for i in (1,Nx):
		for j in (1,Ny):
	
			k2(i,j)=kx(i)**2+ky(j)**2
	
	
	#%--
	
	k4 = k2**2
	return [kx,ky,k2,k4]
	
	
	out2 = open("tim_energ.out","w")
	
	Nx = 64
	Ny = 64
	NxNy = Nx*Ny
	dx = 1
	dy = 1

#Parametros para integracao no tempo:
nstep = 100*60*60*24; #24h
nprint = 50
dtime = 1e-2

#Parametros do material:

cr0 = 0.40

R = 8.314472
T = 723
RT = R*T

a = 0.2866e-9

nn = (a/2)*(math.sqrt(3))
int_par = 18600 + 0.1*T
#grcoef_cr = (fz/2)*(a^2)*ddGexddcr
#grcoef_cr = 4e-16; %Coeficiente gradiente de energia - Segundo 2.1 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
Hmix = cr0*(1-cr0)*int_par
grcoef_cr = (2/3)*Hmix*nn
#grcoef_cr = Hmix*((nn^2)/2)
grcoef_cr = grcoef_cr*1e10


#Parametros de difusibilidade para mobilidade:
#--Segundo 2.3 BARKAR et al - Effect of concentration dependent gradient energy coefficient on spinodal decomposition in the Fe-Cr system
QCr=3.08e5 # Valor dummy - falta alterar com os dados do Thermo-Calc
D0Cr=2e-5
DCr=(D0Cr*math.exp(-QCr/RT))

QFe=2.94e5 # Valor dummy - falta alterar com os dados do Thermo-Calc
D0Fe=1e-4
DFe=(D0Fe*math.exp(-QFe/RT))

MFe = DFe/RT
MCr = DCr/RT

M=cr0*(1.0-cr0)*(DCr/RT);
M=cr0*(1.0-cr0)*(cr0*MFe + (1-cr0)*MCr);
#M = M*1e26;
#inicialização da microestrutura com concentração c0
cr = micro_ch_pre(Nx,Ny,cr0);
#[kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy);
