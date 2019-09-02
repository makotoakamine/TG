function [con] =micro_ch_pre(Nx,Ny,c0)

format long;

noise = 0.001;

for i=1:Nx
	for j=1:Ny
		con(i,j)=c0 + noise*(0.5-rand);
	end
end


end %endfunction
