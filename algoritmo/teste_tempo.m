time0=clock();
format long;

Nx= 128;
Ny= 128;
NxNy =Nx*Ny;

dx=0.5;
dy=0.5;

nstep =   100;
nprint=    50;
dtime = 1.0e-2;
ttime = 0.0;
coefA = 1.0;

cr = 0.4;

[cr] = micro_ch_pre(Nx,Ny,cr0);


for istep =1:nstep
	ttime = ttime +dtime;

##	for i=1:Nx
##  	for j=1:Ny
##      cr(i,j) = (cr(i,j)*cr(i,j) + 1)/cr(i,j) ;
##    end
##  end
	cr = (cr.*cr + 1)./cr ;
end %istep

compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);