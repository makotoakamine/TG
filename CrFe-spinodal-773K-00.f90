!------------------------------------------------------------------------------
! Cahn-Hilliard program
! Tutorial from MOOSE
! https://mooseframework.inl.gov/old/wiki/MooseTutorials/IronChromiumDecomposition/
!------------------------------------------------------------------------------

program teste

    implicit none

!------------------------------------------------------------------------------
 
    integer,parameter           :: nx=64,ny=64,nz=1
    integer,parameter           :: freq=1000
    integer,parameter           :: tmax= 6000000
    integer,parameter           :: zero=0
    integer                     :: nxy,t,ierror,i,j

    double precision,parameter  :: l=25.0
    double precision,parameter  :: calloy = 0.46774
    double precision,parameter  :: dt=0.05
    double precision            :: dx,dy
    double precision            ::  c(0:nx+1,0:ny+1)
    double precision            :: c0(0:nx+1,0:ny+1)
    double precision            :: gr2c(0:nx+1,0:ny+1)
    double precision            :: mu1(0:nx+1,0:ny+1)
    double precision            :: mu2(0:nx+1,0:ny+1)
    double precision            :: df(0:nx+1,0:ny+1)
    double precision            :: g(0:nx+1,0:ny+1)
    double precision            :: gr1mu1(0:nx+1,0:ny+1)
    double precision            :: gr1mu2(0:nx+1,0:ny+1)
    double precision            :: mob(0:nx+1,0:ny+1)
    double precision            :: k

    double precision            :: gmed,cmed

    character                   :: vtkfile*30

!-------------------------------------------------------------------------------

    write(*,*)
    write(*,"(A)") "! ! ! ! ! ! ! ! ! Starting simulation ! ! ! ! ! ! ! ! !"

    !--- open file to save data
    open (20,FILE='data.txt',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)

    !--- space discretization
    dx = l/nx; dy = dx

    !--- total number of grid points
    nxy = nx*ny

    k = (1E18)*(8.125E-16)

    !--- set initial c field
    call initial(calloy,nx,ny,c0)

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

        !--- apply pbc to c field
        call pbc(nx,ny,c0)

        !--- calculate Laplacian of c
        call nabla2(c0,nx,ny,dx,dy,gr2c)

        !--- calculate derivative of free energy
        call d1g(c0,nx,ny,df)

        !--- calculate mu1
        mu1 = df-k*gr2c

        !--- apply pbc to mu1
        call pbc(nx,ny,mu1)

        !--- calculate gradient of mu1
        call nabla1(mu1,nx,ny,dx,dy,gr1mu1)

        !--- calculate mobility
        call mobility(c0,nx,ny,mob)

        !--- calculate mu2
        mu2 = (1E18)*mob*gr1mu1
        
        !--- apply pbc to mu2
        call pbc(nx,ny,mu2)

        !--- calculate gradient of mu2
        call nabla1(mu2,nx,ny,dx,dy,gr1mu2)

        !--- update c field
        c  = c0 + dt*gr1mu2
        c0 = c
        
        !--- output
        if(mod(t,freq)==0) then

            !--- write c field to vtk file
            write(vtkfile,"(A9,I0,A4)") "vtk/time_",t,".vtk"
            call savevtk(c0,nx,ny,nz,dx,dy,vtkfile)

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

    write(*,"(A)") "! ! ! ! ! ! ! ! ! ! Simulation end ! ! ! ! ! ! ! ! ! ! "
 
end program teste

!-------------------------------------------------------------------------------
! Returns a normally distributed deviate with mean=0.0 and variance=1.0
! uses ran1(idum) as the source of uniform deviates.
! from: Numerical Recipes in Fortran 77 (2nd ed) pg I-280

function gasdev(idum)

    integer idum
    real    gasdev
    integer iset
    real    fac,gset,rsq,v1,v2,ran1

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
    real      ran1,AM,EPS,RNMX
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
    double precision,parameter  :: a = -0.002,b = 0.002
    double precision,parameter  :: v=0.01
    double precision            :: m(0:nx+1,0:ny+1),calloy

    idum = 92923932

    do i=0,nx+1
            do j=0,ny+1

                !m(i,j) = calloy + (ran1(idum)*(b-a))+a  ! uniform distribution
                m(i,j) = calloy + gasdev(idum)*v        ! gaussian distribution

             end do
        end do

    return

end subroutine initial

!------------------------------------------------------------------------------
! chemical free energy

subroutine gchem(c,nx,ny,g)

    double precision :: A1,A2,A3,A4,A5,A6,A7
    double precision :: g(0:nx+1,0:ny+1)
    double precision :: c(0:nx+1,0:ny+1)

    A1 =-24468.31
    A2 =-28275.33
    A3 = 4167.994
    A4 = 7052.907
    A5 = 12089.93
    A6 = 2568.625
    A7 =-2345.293

    g = A1*c + A2*(1-c) + A3*c*log(c) + A4*(1-c)*log(1-c) + A5*c*(1-c) &
        + A6*c*(1-c)*(2*c-1) + A7*c*(1-c)*(2*c-1)*(2*c-1)

    return

end subroutine gchem

!------------------------------------------------------------------------------
! 1st derivative of free energy

subroutine d1g(c,nx,ny,dg)

    double precision :: A1,A2,A3,A4,A5,A6,A7
    double precision :: c(0:nx+1,0:ny+1)
    double precision :: dg(0:nx+1,0:ny+1)

    A1 =-24468.31
    A2 =-28275.33
    A3 = 4167.994
    A4 = 7052.907
    A5 = 12089.93
    A6 = 2568.625
    A7 =-2345.293

    dg = A1-A2 + A3*(1+log(c)) - A4*(1+log(1-c)) + A5*(1-2*c) + &
         A6*(-6*c**2 + 6*c -1) + A7*(-16*c**3 + 24*c**2 - 10*c + 1)

    return

end subroutine d1g

!------------------------------------------------------------------------------
! mobilities

subroutine mobility(c,nx,ny,mob)

    double precision :: A1,A2,A3,A4,A5,A6,A7
    double precision :: B1,B2,B3,B4,B5,B6,B7
    double precision :: c(0:nx+1,0:ny+1)
    double precision :: gcr(0:nx+1,0:ny+1)
    double precision :: gfe(0:nx+1,0:ny+1)
    double precision :: mob(0:nx+1,0:ny+1)

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

    double precision :: gr1(0:nx+1,0:ny+1)
    double precision :: a(0:nx+1,0:ny+1)
    double precision :: dx,dy

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

    double precision :: gr2(0:nx+1,0:ny+1)
    double precision :: a(0:nx+1,0:ny+1)
    double precision :: dx,dy,dx2_in

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

    double precision  :: a(0:nx+1,0:ny+1)
    integer           :: nx,ny

    a(nx+1,:) = a(1,:)
    a(0,:)    = a(nx,:)
    a(:,ny+1) = a(:,1)
    a(:,0)    = a(:,ny)

end subroutine pbc

!------------------------------------------------------------------------------
! write data to file in vtk format

subroutine savevtk(c,nx,ny,nz,dx,dy,outfile)

    double precision  :: c(0:nx+1,0:ny+1)
    double precision  :: dx,dy,x,y,z

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

