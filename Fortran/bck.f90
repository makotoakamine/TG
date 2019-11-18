program test
use, intrinsic :: iso_c_binding 
include 'fftw3.f03'

implicit none

real*16 , parameter :: PI=4.D0*DATAN(1.D0)
integer, parameter :: N=16
real*16, parameter :: a=0
real*16, parameter :: b=4*PI
real*16, parameter :: L=(b-a)
real*16,parameter :: dx=L/N
integer :: i
integer, parameter :: N2 = N/2
real*16, dimension(N) :: x = (/(i*dx, i=0,(N-1),1)/)
real*16, dimension(N) :: y
complex(C_DOUBLE_COMPLEX), dimension(N) :: yt
real*16, dimension(N) :: dy
complex(C_DOUBLE_COMPLEX), dimension(N) :: idyt
complex(C_DOUBLE_COMPLEX), dimension(N) :: k
double precision :: w

type(C_PTR) plan
  !double  complex :: in, out
 complex(C_DOUBLE_COMPLEX), dimension(N) :: in
 complex(C_DOUBLE_COMPLEX), dimension(N) :: out
!        double precision y
!        dimension y(N)
!        double complex out
!        dimension out(N/2+1)

  y = sin(x)
  dy = cos(x)
  k(1) = 0
  k(N2+1) = 0

  do i= 2,N2
    w = 2*PI/(b-a)
	k(i) = dcmplx(float(0),(i-1)*w)
	k(i+N2) = dcmplx(float(0),(i-1-N2)*w)
  enddo

!  do i = 1,N,1
!    in(i) = dcmplx(y(i),float(0))
!    write(*,*) '    in(',i,') = ',in(i)
!  enddo

  in = y
  plan = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_ESTIMATE )

  call fftw_execute_dft( plan, in,out )

  call fftw_destroy_plan( plan )

    yt = out
  idyt = yt*k
  out = idyt


  plan = fftw_plan_dft_1d( N, out, in, FFTW_FORWARD, FFTW_ESTIMATE )

  call fftw_execute_dft( plan, out , in )

  write(*,*) 'derivada espectral'
  do i = 1,N,1
    write(*,*) '    ',N,' * in(',i,') = ',in(i)/N
  enddo

  write(*,*) 'derivada real'
  print *, dy

  call fftw_destroy_plan ( plan )

  end
