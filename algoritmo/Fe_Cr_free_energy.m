function [dgdcr] = Fe_Cr_free_energy(Nx,Ny,cr0,cr,R,tempr)

format long;

RT=R*tempr;

%Coef=17000.0;
L = 20500 - 9.68*tempr;

for i=1:Nx
	for j=1:Ny
	
		c2=cr(i,j);
		c1=1.0-c2;

    			%Dados SGTE da energia livre em função da temperatura
			%G do Fe em relação a BCC_A2 paramagnetico
			p = 0.4;
			D = 1.55828482; % D = 510/1125 + (11692/15975)*(1/p-1)
			B0 =  2.22;
			Tc = tempr/1043;
			tal = 1/Tc;
			if(tal <= 1)
				gtal = 1-((79/140)*(1/p)*(tal^(-1)) + (474/497)*(1/p-1)*((tal^3)/6 + (tal^9)/135 + (tal^15)/600))/D;
			else
				gtal = -( ((tal^(-5))/10 + (tal^(-15))/315 + (tal^(-25))/1500 )/D );
			end %if
			Gmag = RT*log(B0+1)*gtal;

			a0 = 2.3987e-5;
			a1 = 2.569e-8;
			P = 10e5;
			A = 7.042095e-6;
			Gpres = A*P*(1+a0*tempr + (a1/2)*(tempr^2));

			Gfe = 0 + Gmag  + Gpres; % 298.15 < T < 6000.00

			%G do Cr em relação a BCC_A2 paramagnetico 

			a0 = 1.5e-5;
			a1 = 1.84e-8;
			P = 10e5;
			A = 7.188e-6;
			Gpres = A*P*(1+a0*tempr + (a1/2)*(tempr^2));
			Gcr = 0 + Gpres; % 311.50 < T < 6000.00


			%== Equação 5.27 BINER segundo SGTE
      % gcr=(1-c2)*Gfe+c2*Gcr+L*(1-c2)*c2-RT*(c2*log(c2)-(1-c2)*log(1-c2);
			%derivada da função G em relação a concentração de Cro
			%falta levar em consideração G_excesso
			dgcr=(-Gfe+Gcr+L*(1-2*c2)+RT*log(c2/c1))/RT; % Gex
%dgcua=1.4613878411949395E-4*(36076.894*c1+6842.810456*(log(c2)-log(c1))-36076.894*c2+2984.135);
			dgdcr(i,j)=dgcr;
    
##    
##		if(c1 < 0.9995 && c1 > 0.0005)
##			%Dados SGTE da energia livre em função da temperatura
##			%G do Fe em relação a BCC_A2 paramagnetico
##			p = 0.4;
##			D = 1.55828482; % D = 510/1125 + (11692/15975)*(1/p-1)
##			B0 =  2.22;
##			Tc = tempr/1043;
##			tal = 1/Tc;
##			if(tal <= 1)
##				gtal = 1-((79/140)*(1/p)*(tal^(-1)) + (474/497)*(1/p-1)*((tal^3)/6 + (tal^9)/135 + (tal^15)/600))/D;
##			else
##				gtal = -( ((tal^(-5))/10 + (tal^(-15))/315 + (tal^(-25))/1500 )/D );
##			end %if
##			Gmag = RT*log(B0+1)*gtal;
##
##			a0 = 2.3987e-5;
##			a1 = 2.569e-8;
##			P = 10e5;
##			A = 7.042095e-6;
##			Gpres = A*P*(1+a0*tempr + (a1/2)*(tempr^2));
##
##			Gfe = 0 + Gmag  + Gpres; % 298.15 < T < 6000.00
##
##			%G do Cr em relação a BCC_A2 paramagnetico 
##
##			a0 = 1.5e-5;
##			a1 = 1.84e-8;
##			P = 10e5;
##			A = 7.188e-6;
##			Gpres = A*P*(1+a0*tempr + (a1/2)*(tempr^2));
##			Gcr = 0 + Gpres; % 311.50 < T < 6000.00
##
##
##			%== Equação 5.27 BINER segundo SGTE
##			%derivada da função G em relação a concentração de Cro
##			%falta levar em consideração G_excesso
##			dgcr=-Gfe+Gcr+L*(1-2*c2)-RT*log(c2/c1); % Gex
##%dgcua=1.4613878411949395E-4*(36076.894*c1+6842.810456*(log(c2)-log(c1))-36076.894*c2+2984.135);
##			dgdcr(i,j)=dgcr;
##		end % if 
##		 
##		if(c1 >=0.9995)
##			dgcr=-2.0*L*((1.0-c2)-0.9995);
##			dgdcr(i,j)=dgcr;
##		end %if
##		
##		if( c1 <= 0.0005)
##			dgcr=-2.0*L*((1.0-c2)-0.0005);
##			dgdcu(i,j)=dgcr;
##		end %if
		
	end %i
end %j
