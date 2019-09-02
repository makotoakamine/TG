function [kx,ky,k2,k4] = prepare_fft(Nx,Ny,dx,dy)
  
format long;

Nx21 = Nx/2 + 1;
Ny21 = Ny/2 + 1;

Nx2=Nx+2;
Ny2=Ny+2;

%--

delkx=(2.0*pi)/(Nx*dx);
delky=(2.0*pi)/(Ny*dy);

%--

for i=1:Nx21
fk1=(i-1)*delkx;
kx(i)=fk1;
kx(Nx2-i)=-fk1;
end

for j=1:Ny21
fk2=(j-1)*delky;
ky(j)=fk2;
ky(Ny2-j)=-fk2;
end

%---

for i=1:Nx
for j=1:Ny

k2(i,j)=kx(i)^2+ky(j)^2;

end
end

%%--

k4 = k2.^2;

end %endfunction
