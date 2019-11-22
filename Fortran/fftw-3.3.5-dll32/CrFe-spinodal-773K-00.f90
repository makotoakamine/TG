!------------------------------------------------------------------------------
! Cahn-Hilliard program
! Tutorial from MOOSE
! https://mooseframework.inl.gov/old/wiki/MooseTutorials/IronChromiumDecomposition/
!------------------------------------------------------------------------------

program teste
    use fftw3
    implicit none

!------------------------------------------------------------------------------
 
    integer,parameter           :: nx=64,ny=64,nz=1
    integer,parameter           :: freq=100
    integer,parameter           :: tmax= 600000
    integer,parameter           :: zero=0
    integer                     :: nxy,t,ierror,i,j
    !integer iret
	real :: start,finish

    real*16,parameter  :: l=25.0
    real*16,parameter  :: calloy = 0.40
    real*16,parameter  :: dt=0.50
    real*16            :: dx,dy
    real*16            :: R,Temp,RT,QCr,QFe,D0Cr,D0Fe,DCr,DFe,MCr,MFe

real(C_DOUBLE), pointer            :: c(:,:)
type(C_PTR)                        :: pc
complex(C_DOUBLE_COMPLEX), pointer :: ck(:,:)
type(C_PTR)                        :: pck
complex(C_DOUBLE_COMPLEX), pointer :: cktemp(:,:)
type(C_PTR)                        :: pcktemp
real(C_DOUBLE), pointer            :: c0(:,:)
type(C_PTR)                        :: pc0
real(C_DOUBLE), pointer            :: gr2c(:,:)
type(C_PTR)                        :: pgr2c
real(C_DOUBLE), pointer            :: grcoef_cr(:,:)
type(C_PTR)                        :: pgrcoef_cr
real(C_DOUBLE), pointer            :: mcoef_cr(:,:)
type(C_PTR)                        :: pmcoef_cr
real(C_DOUBLE), pointer            :: mu1(:,:)
type(C_PTR)                        :: pmu1
real(C_DOUBLE), pointer            :: mu2(:,:)
type(C_PTR)                        :: pmu2
real(C_DOUBLE), pointer            :: df(:,:)
type(C_PTR)                        :: pdf
complex(C_DOUBLE_COMPLEX), pointer :: dfk(:,:)
type(C_PTR)                        :: pdfk
real(C_DOUBLE), pointer            :: g(:,:)
type(C_PTR)                        :: pg
real(C_DOUBLE), pointer            :: gr1mu1(:,:)
type(C_PTR)                        :: pgr1mu1
real(C_DOUBLE), pointer            :: gr1mu2(:,:)
type(C_PTR)                        :: pgr1mu2
real(C_DOUBLE), pointer            :: mob(:,:)
type(C_PTR)                        :: pmob
real(C_DOUBLE), pointer            :: k2(:,:)
type(C_PTR)                        :: pk2
real(C_DOUBLE), pointer            :: k4(:,:)
type(C_PTR)                        :: pk4
complex(C_DOUBLE_COMPLEX), pointer :: in(:,:)
type(C_PTR)                        :: pin
complex(C_DOUBLE_COMPLEX), pointer :: out(:,:)
type(C_PTR)                        :: pout

real*16                            :: k
type(C_PTR)                        :: plan_c2ck,plan_df2dfk,plan_ck2c

real*16                            :: gmed,cmed

character                   :: vtkfile*30

    !--- total number of grid points
nxy = nx*ny

!Alocacao de memoria para as matrizes
pc = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pc,c,[nx,ny])

pck = fftw_alloc_complex(int(nxy,C_SIZE_T))
call c_f_pointer(pck,ck,[nx,ny])

pcktemp = fftw_alloc_complex(int(nxy,C_SIZE_T))
call c_f_pointer(pcktemp,cktemp,[nx,ny])

pc0 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pc0,c0,[nx,ny])

pgr2c = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pgr2c,gr2c,[nx,ny])

pgrcoef_cr = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pgrcoef_cr,grcoef_cr,[nx,ny])

pmcoef_cr = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pmcoef_cr,mcoef_cr,[nx,ny])

pmu1 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pmu1,mu1,[nx,ny])

pmu2 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pmu2,mu2,[nx,ny])

pdf = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pdf,df,[nx,ny])

pdfk = fftw_alloc_complex(int(nxy,C_SIZE_T))
call c_f_pointer(pdfk,dfk,[nx,ny])

pg = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pg,g,[nx,ny])

pgr1mu1 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pgr1mu1,gr1mu1,[nx,ny])

pgr1mu2 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pgr1mu2,gr1mu2,[nx,ny])

pmob = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pmob,mob,[nx,ny])

pk2 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pk2,k2,[nx,ny])

pk4 = fftw_alloc_real(int(nxy,C_SIZE_T))
call c_f_pointer(pk4,k4,[nx,ny])

pin = fftw_alloc_complex(int(nxy,C_SIZE_T))
call c_f_pointer(pin,in,[nx,ny])

pout = fftw_alloc_complex(int(nxy,C_SIZE_T))
call c_f_pointer(pout,out,[nx,ny])

!mult thread
!call dfftw_init_threads(iret)
!call dfftw_plan_with_nthreads(8)



!-------------------------------------------------------------------------------

    write(*,*)
    write(*,"(A)") "! ! ! ! ! ! ! ! ! Starting simulation ! ! ! ! ! ! ! ! !"

	plan_c2ck = fftw_plan_dft_2d( nx, ny, in, ck, FFTW_FORWARD, FFTW_MEASURE )
	plan_df2dfk = fftw_plan_dft_2d( nx, ny, in, dfk, FFTW_FORWARD, FFTW_MEASURE )
	plan_ck2c = fftw_plan_dft_2d( nx, ny, cktemp, out, FFTW_BACKWARD, FFTW_MEASURE )

    !--- open file to save data
    open (20,FILE='data.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)

    !--- space discretization
    dx = l/nx; dy = dx


    R = 8.314472
	T = 700.00
	RT = R*T
    grcoef_cr = 8.125E-16
    mcoef_cr =  2.2841E-26
    k = 8.125E-16

    QCr=3.08e5
    QFe=2.94e5 
    D0Cr=2.00E-5
    D0Fe=1.00E-4

    !DCr=(D0Cr*exp(-QCr/RT));
    !DFe=(D0Fe*exp(-QFe/RT));
    !DCr = D0Cr
    !DFe = D0Fe
    !MFe = DFe/RT
    !MCr = DCr/RT

    !--- set initial c field
    call initial(calloy,nx,ny,c0)
	c = c0
	call prepare_fft(k2,k4,nx,ny,dx,dy)

    !--- save initial c field to vtk file
    vtkfile = 'vtk/time_0.vtk'
    call savevtk(c0,nx,ny,nz,dx,dy,vtkfile)



    write(*,*)

    write(*,"(A)") "! ! ! ! ! ! ! ! ! Initialization done ! ! ! ! ! ! ! ! !"
    write(*,*)
    write(*,"(A)") "! ! ! ! ! ! ! ! ! Starting time loop  ! ! ! ! ! ! ! ! !"

    write(*,*)

    !--- start time loop
    do t=1,tmax
		!call cpu_time(start)

        !mob = c*(1.0-c)*(c*MFe + (1.0-c)*MCr)
        !call mobility(c,nx,ny,mob)
		!print *, mob
!----------------------
!-------fft c -> ck
		in = c
		call fftw_execute_dft ( plan_c2ck, in, ck)
!----------------------
!----------------------

        !call d1g(c,nx,ny,df)
        call d1g(c,nx,ny,df,Temp,RT)

!----------------------
!-------fft df -> dfk
		in = df
		call fftw_execute_dft (plan_df2dfk,in, dfk)
!----------------------
!----------------------

        !cktemp = ((ck - dt*(1E18*mob)*dfk*k2) / (1.0 + dt*(1E18*k)*k4*(1E18*mob)))/nxy;
        cktemp = ((ck - dt*(1E18*mcoef_cr)*dfk*k2) / (1.0 + dt*(1E18*grcoef_cr)*k4*(1E18*mcoef_cr)))/nxy;
        !cktemp = ((ck - dt*(mcoef_cr)*dfk*k2) / (1.0 + dt*(grcoef_cr)*k4*(mcoef_cr)))/nxy;

!----------------------
!-------fft ck -> c
		call fftw_execute_dft ( plan_ck2c, cktemp, out)
		c = out
		!call cpu_time(finish)
    	!print '("Time = ",f6.3," seconds.")',finish-start
!----------------------
!----------------------

!		do i=1,nx
!			do j=1,ny
!				if (c(i,j) >= 0.9999) then
!					c(i,j) = 0.9999
!				end if
!				if (c(i,j) <= 0.00001) then
!					c(i,j) = 0.00001
!				end if
!			enddo
!		enddo


        !--- output
        if(mod(t,freq)==0) then

            !--- write c field to vtk file
            write(vtkfile,"(A9,I0,A4)") "vtk/time_",t,".vtk"
            call savevtk(c,nx,ny,nz,dx,dy,vtkfile)

            !--- calculate sum of free energy and c 
            call gchem(c,nx,ny,g)
            gmed = sum(g(1:nx,1:ny))/nxy
            cmed = sum(c(1:nx,1:ny))/nxy

            !--- write c and f sum to file
            write(20,*) t,cmed,gmed
            
            !--- write to screen
            write(*,5) "Timestep     : ",t
            5 format (a,i0)

            write(*,6) "Concentration: ",cmed
            write(*,6) "Free energy  : ",gmed
            6 format (a,1PE13.6)

            write(*,*)

            !write(*,7)"-------------------------------------------------------"
            !7 format (a)

        end if

    end do ! end of time loop

    !call fftw_free(p)

    write(*,"(A)") "! ! ! ! ! ! ! ! ! ! Simulation end ! ! ! ! ! ! ! ! ! ! "
 
end program teste

!-------------------------------------------------------------------------------
! Returns a normally distributed deviate with mean=0.0 and variance=1.0
! uses ran1(idum) as the source of uniform deviates.
! from: Numerical Recipes in Fortran 77 (2nd ed) pg I-280

function gasdev(idum)

    integer idum
    real*16    gasdev
    integer iset
    real*16    fac,gset,rsq,v1,v2,ran1

    save    iset,gset
    data    iset/0/

    if (idum.lt.0) iset=0                 ! reinitialize

    if (iset.eq.0) then 
     1  v1=2.0*ran1(idum)-1.0             ! pick 02 uniform numbers
        v2=2.0*ran1(idum)-1.0             ! from -1 to +1

        rsq=v1**2+v2**2                   ! see if v1,v2 are in the unit circle

        if(rsq.ge.1..or.rsq.eq.0.)goto 1  ! if not try again

            fac=sqrt(-2.*log(rsq)/rsq)    ! Box-Muller transformation to get
            gset=v1*fac                   ! two normal deviates. Return one
            gasdev=v2*fac                 ! and save the other for next time

            iset=1                        ! set flag

        else
            gasdev=gset                   ! return 2nd deviate
            iset=0                        ! unset flag
    end if

return
 
end function gasdev

!------------------------------------------------------------------------------
! Returns a uniform random deviate between 0.0 and 1.0 (excluding endpoints)
! Call with idum a negative integer to initialize;
! thereafter, do not alter idum between successive deviates in a sequence.
! RNMX should approximate the largest floating value that is less than 1
! from: Numerical Recipes in Fortran 77 (2nd ed) pg I-271

function ran1(idum)

    integer   idum,IA,IM,IQ,IR,NTAB,NDIV
    real*16      ran1,AM,EPS,RNMX
    parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
    parameter (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    integer   j,k,iv(NTAB),iy

    save iv,iy
    data iv /NTAB*0/, iy /0/

    if (idum.le.0.or.iy.eq.0) then       ! be sure to prevent idum=0
        idum=max(-idum,1)

        do 11 j=NTAB+8,1,-1              ! load shuffle table after 8 warm-ups
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
            11 continue

        iy=iv(1)
    endif

    k=idum/IQ                           ! start here when not initializing          

    idum=IA*(idum-k*IQ)-IR*k            ! compute idum=mod(IA*idum,IM) without
                                        ! overflows by Schrage’s method
    if (idum.lt.0) idum=idum+IM
    j=1+iy/NDIV                         ! Will be in the range 1:NTAB

    iy=iv(j)                            ! Output previously stored value and 
    iv(j)=idum                          ! refill the shuffle table

    ran1=min(AM*iy,RNMX)                ! users don’t expect endpoint values
    return

    end function ran1

!-------------------------------------------------------------------------------
! set initial concentration
! gasdev: normally distributed random variables with mean=0.0, var=1.0 -> N(0,1)
! 
! if x is N(0,1) then y=a+v*x is N(a,v^2) i.e
! 
! normally distributed random variables with mean=a, var=v^2

subroutine initial(calloy,nx,ny,m)

    integer                     :: nx,ny,i,j,idum
    real*16,parameter           :: a = -0.002,b = 0.002
    real*16,parameter           :: v=0.01
    !double precision            :: m(0:nx+1,0:ny+1),calloy
    real*8,dimension(nx,ny)    :: m
    !real(C_DOUBLE), intent(out)    :: m(:,:)
    real*16                     :: calloy
    real*16,external            :: gasdev

    idum = 92923932

    do i=1,nx,1
            do j=1,ny,1

                !m(i,j) = calloy + (ran1(idum)*(b-a))+a  ! uniform distribution
                m(i,j) = calloy + gasdev(idum)*v        ! gaussian distribution

             end do
        end do

    return

end subroutine initial

!------------------------------------------------------------------------------
! chemical free energy

!subroutine gchem(c,nx,ny,g,T,RT)
subroutine gchem(c,nx,ny,g)

implicit none
    integer, intent(in) :: nx,ny
    real*8, intent(in) :: c(1:nx,1:ny)
    real*8, intent(out) :: g(1:nx,1:ny)
    real*16 :: A1,A2,A3,A4,A5,A6,A7
    !double precision :: g(0:nx+1,0:ny+1)
    !double precision :: c(0:nx+1,0:ny+1)
    !double precision, dimension(nx,ny) :: g,c

    A1 =-24468.31
    A2 =-28275.33
    A3 = 4167.994
    A4 = 7052.907
    A5 = 12089.93
    A6 = 2568.625
    A7 =-2345.293

    g = A1*c + A2*(1-c) + A3*c*log(c) + A4*(1-c)*log(1-c) + A5*c*(1-c) &
        + A6*c*(1-c)*(2*c-1) + A7*c*(1-c)*(2*c-1)*(2*c-1)

!-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/
!-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/

!    real*16 :: GFe,GCr,p,D,B0,Tc,tal,a0,a1,P,A,L
!    real*16 :: g(0:nx+1,0:ny+1)
!    real*16 :: c(0:nx+1,0:ny+1)
!
!!Dados SGTE da energia livre em funcao da Temperatura para Fe BCC_A2 paramagnetico
!
!    p = 0.4
!    D = 1.55828482; ! D = 510/1125 + (11692/15975)*(1/p-1)
!    B0 =  2.22;
!    Tc = T/1043;
!    tal = 1/Tc;
!    if (tal <= 1) then
!        gtal = 1-((79/140)*(1/p)*(tal^(-1)) + (474/497)*(1/p-1)*((tal^3)/6 + (tal^9)/135 + (tal^15)/600))/D;
!    else
!        gtal = -( ((tal^(-5))/10 + (tal^(-15))/315 + (tal^(-25))/1500 )/D );
!    end if
!    Gmag = RT*log(B0+1)*gtal;
!
!    a0 = 2.3987E-5;
!    a1 = 2.569E-8;
!    P = 10E5;
!    A = 7.042095E-6;
!    Gpres = A*P*(1+a0*T + (a1/2)*(T^2));
!
!    Gfe = 0 + Gmag  + Gpres; ! 298.15 < T < 6000.00
!
!
!!	Dados SGTE da energia livre em funcao da Temperatura para Cr BCC_A2 paramagnetico
!    a0 = 1.5E-5;
!    a1 = 1.84E-8;
!    P = 10E5;
!    A = 7.188E-6;
!    Gpres = A*P*(1+a0*T + (a1/2)*(T^2));
!
!    Gcr = 0 + Gpres; ! 311.50 < T < 6000.00
!
!    L = L = 20500 - 9.68*T;
!
!    g =(1-c)*GFe+c*GCr+L*(1-c)*c+RT*(c*log(c)-(1-c)*log(1-c);

    return

end subroutine gchem

!------------------------------------------------------------------------------
! 1st derivative of free energy

subroutine d1g(c,nx,ny,dg,Temp,RT)
!subroutine d1g(c,nx,ny,dg)
implicit none

    integer, intent(in) :: nx,ny
    real*8, intent(in) :: c(1:nx,1:ny)
    real*16, intent(in) :: Temp
    real*16, intent(in) :: RT
    real*8, intent(out) :: dg(1:nx,1:ny)
    real*16 :: GFe,GCr,Gpres,Gmag,gtal,p,D,B0,Tc,tal,a0,a1,Pres,A,L
!    real*16 :: A1,A2,A3,A4,A5,A6,A7
!    !double precision :: g(0:nx+1,0:ny+1)
!    !double precision :: c(0:nx+1,0:ny+1)
!    !double precision, dimension(nx,ny) :: g,c
!	
!    A1 =-24468.31
!    A2 =-28275.33
!    A3 = 4167.994
!    A4 = 7052.907
!    A5 = 12089.93
!    A6 = 2568.625
!    A7 =-2345.293
!
!    dg = A1-A2 + A3*(1+log(c)) - A4*(1+log(1-c)) + A5*(1-2*c) + &
!         A6*(-6*c**2 + 6*c -1) + A7*(-16*c**3 + 24*c**2 - 10*c + 1)

	
!	Dados SGTE da energia livre em funcao da Temperatura para Fe BCC_A2 paramagnetico

    p = 0.4
    D = 1.55828482; ! D = 510/1125 + (11692/15975)*(1/p-1)
    B0 =  2.22;
    Tc = Temp/1043;
    tal = 1/Tc;
    if (tal <= 1) then
        gtal = 1-((79/140)*(1/p)*(tal**(-1)) + (474/497)*(1/p-1)*((tal**3)/6 + (tal**9)/135 + (tal**15)/600))/D;
    else
        gtal = -( ((tal**(-5))/10 + (tal**(-15))/315 + (tal**(-25))/1500 )/D );
    end if
    Gmag = RT*log(B0+1)*gtal;

    a0 = 2.3987E-5;
    a1 = 2.569E-8;
    Pres = 10E5;
    A = 7.042095E-6;
    Gpres = A*Pres*(1+a0*Temp + (a1/2)*(Temp**2));

    Gfe = 0 + Gmag  + Gpres; ! 298.15 < Temp < 6000.00


!	Dados SGTE da energia livre em funcao da Temperatura para Cr BCC_A2 paramagnetico
    a0 = 1.5E-5;
    a1 = 1.84E-8;
    Pres = 10E5;
    A = 7.188E-6;
    Gpres = A*Pres*(1+a0*Temp + (a1/2)*(Temp**2));

    Gcr = 0 + Gpres; ! 311.50 < Temp < 6000.00

    L = 20500 - 9.68*Temp;

    dg =(-GFe+GCr+L*(1-2*c)+RT*log(c/(1-c)));

    return

end subroutine d1g

!------------------------------------------------------------------------------
! mobilities

subroutine mobility(c,nx,ny,mob)
implicit none

    real*8, intent(in)  :: c(1:nx,1:ny)
	integer, intent(in)           :: nx,ny
    real*8, intent(out) :: mob(1:nx,1:ny)
    real*16              :: A1,A2,A3,A4,A5,A6,A7
    real*16              :: B1,B2,B3,B4,B5,B6,B7
!    real*16 :: c(0:nx+1,0:ny+1)
!    real*16 :: gcr(0:nx+1,0:ny+1)
!    real*16 :: gfe(0:nx+1,0:ny+1)
!    real*16 :: mob(0:nx+1,0:ny+1)
    real*16 :: gcr(1:nx,1:ny)
    real*16 :: gfe(1:nx,1:ny)

    A1 = -32.770969
    A2 = -25.8186669
    A3 = -3.29612744
    A4 = 17.669757
    A5 = 37.6197853
    A6 = 20.6941796
    A7 = 10.8095813

    B1 = -31.687117
    B2 = -26.0291774
    B3 = 0.2286581
    B4 = 24.3633544
    B5 = 44.3334237
    B6 = 8.72990497
    B7 = 20.956768

    gcr = A1*c + A2*(1-c) + A3*c*log(c) + A4*(1-c)*log(1-c) + A5*c*(1-c) &
         + A6*c*(1-c)*(2*c-1) + A7*c*(1-c)*(2*c-1)*(2*c-1)

    gfe = B1*c + B2*(1-c) + B3*c*log(c) + B4*(1-c)*log(1-c) + B5*c*(1-c) &
         + B6*c*(1-c)*(2*c-1) + B7*c*(1-c)*(2*c-1)*(2*c-1)

    mob = c*(1-c)*(1-c)*(10**gcr) + (1-c)*c*c*(10**gfe)
        
    return

end subroutine mobility

!------------------------------------------------------------------------------
! finite difference calculation of gradient of A; output in gr1

subroutine nabla1(a,nx,ny,dx,dy,gr1)
implicit none

    integer, intent(in)           :: nx,ny
!    real*16 :: gr1(0:nx+1,0:ny+1)
!    real*16 :: a(0:nx+1,0:ny+1)
    real*16 :: gr1(1:nx,1:ny)
    real*16 :: a(1:nx,1:ny)
	integer :: i,j
    !double precision, dimension(nx,ny) :: gr1g
    !double precision, dimension(nx,ny) :: a
    real*16 :: dx,dy

    do i=1,nx
        do j=1,ny

            gr1(i,j) = (a(i+1,j)-a(i-1,j)+a(i,j+1)-a(i,j-1))/(2.0*dx)
            
        end do
    end do

end subroutine nabla1

!------------------------------------------------------------------------------
! finite difference calculation of Laplacian of A; output in gr2
! from:  Nikolas Provatas and Ken Elder
!        Phase-Field Methods in Materials Science and Engineering
!        (appendix -- eq B6)

subroutine nabla2(a,nx,ny,dx,dy,gr2)
implicit none

    integer, intent(in)           :: nx,ny
    real*16 :: gr2(0:nx+1,0:ny+1)
    real*16 :: a(0:nx+1,0:ny+1)
    real*16 :: dx,dy,dx2_in
	integer :: i,j

    dx2_in = 1.0/dx**2

    do i=1,nx
        do j=1,ny

            gr2(i,j) = ((a(i+1,j)+a(i-1,j)+a(i,j+1)+a(i,j-1))/2.0d0+ &
                        (a(i+1,j+1)+a(i-1,j+1)+a(i+1,j-1)+a(i-1,j-1))/4.0d0- &
                         3.0d0*a(i,j) )*dx2_in

        end do
    end do

end subroutine nabla2

!------------------------------------------------------------------------------
! periodic boundary conditions - pbc
! apply the following modification to c field before each time step commences

subroutine pbc(nx,ny,a)
implicit none

    real*16  :: a(0:nx+1,0:ny+1)
    integer, intent(in)           :: nx,ny

    a(nx+1,:) = a(1,:)
    a(0,:)    = a(nx,:)
    a(:,ny+1) = a(:,1)
    a(:,0)    = a(:,ny)

end subroutine pbc

!------------------------------------------------------------------------------
! write data to file in vtk format

subroutine savevtk(c,nx,ny,nz,dx,dy,outfile)
implicit none

    integer, intent(in)           :: nx,ny,nz
    real*16, intent(in)  :: dx,dy
    real*8  :: c(1:nx,1:ny)
    real*8  :: x,y,z

    integer           :: ierror,unit,npoin,i,j

    character         :: outfile*30

    npoin = nx*ny*nz

    !--- open vtk file

    unit = 25

    open (UNIT=unit, FILE=outfile, STATUS='REPLACE', ACTION='WRITE', &
          IOSTAT=ierror)

    !--- vtk file header

    write(25,1) '# vtk DataFile Version 2.0'
    write(25,1) 'title'
    write(25,1) 'ASCII'
    write(25,1) 'DATASET STRUCTURED_GRID'

    1 format (a)

    !--- coords of grid points

    write(25,2) 'DIMENSIONS', nx,ny,nz
    2 format(a,i5,i5,i5)

    write(25,3) 'POINTS ',npoin,'  float'
    3 format(a,i6,a)

    do i=1,nx
        do j=1,ny

            x = (i-1)*dx
            y = (j-1)*dy
            z = 0.0

            write(25,4) x,y,z
            4 format(1PE14.6,1PE14.6,1PE14.6)

        end do
    end do

    !--- write grid point values

    write(25,5) 'POINT_DATA ',npoin
    5 format(a,i6)

    write(25,1) 'SCALARS CON  float  1'
    write(25,1) 'LOOKUP_TABLE default'

    do i=1,nx
        do j=1,ny
            write(25,6) c(i,j)
            6 format(1PE14.6)
        end do
    end do

    return

end subroutine savevtk

!------------------------------------------------------------------------------

!-------------------
subroutine prepare_fft(k2,k4,nx,ny,dx,dy)
  integer, intent(in)           :: nx,ny
  real*16, intent(in)  :: dx,dy
  real*16 , parameter :: PI=4.D0*DATAN(1.D0)
  integer :: i ,j
  real*16 :: w
  integer :: Nx2
  integer :: Ny2
  !double precision, dimension(nx) :: kx
  real*8 :: k2(1:nx,1:ny)
  real*8 :: k4(1:nx,1:ny)
  real*8 :: kx(1:nx)
  !double precision, dimension(ny) :: ky
  real*8 :: ky(1:ny)
  Nx2 = nx/2
  Ny2 = ny/2

  w = 2*PI/(nx*dx)
  do i= 2,Nx2
    kx(i) = (i-1)*w
    kx(i+Nx2) = (i-1-Nx2)*w
  enddo

  w = 2*PI/(ny*dy)
  do i= 2,Ny2
    ky(i) = (i-1)*w
    ky(i+Ny2) = (i-1-Ny2)*w
  enddo

  do i=1,nx
    do j=1,ny
      k2(i,j) = kx(i)**2 + ky(j)**2
    enddo
  enddo
  k4 = k2**2
end subroutine prepare_fft
