program test
use, intrinsic :: iso_c_binding 
use fftw3
!include fftw3
implicit none
real :: start, finish
real :: xteste
complex(C_DOUBLE_COMPLEX) :: xcomplex

real*16 , parameter :: PI=4.D0*DATAN(1.D0)
integer, parameter :: N=256
real*16, parameter :: a=0
real*16, parameter :: b=4*PI
real*16, parameter :: L=(b-a)
real*16,parameter :: dx=L/N
integer :: i
real*16 ::dev
integer, parameter :: N2 = N/2
real(C_DOUBLE), pointer :: x(:)
type(C_PTR) px
!real*16, dimension(N) :: x = (/(i*dx, i=0,(N-1),1)/)
real(C_DOUBLE), pointer :: y(:)
type(C_PTR) py
!real*16, dimension(N) :: y

complex(C_DOUBLE_COMPLEX), pointer :: yt(:)
type(C_PTR) pyt

!real*16, dimension(N) :: dy
real(C_DOUBLE), pointer :: dy(:)
type(C_PTR) pdy

complex(C_DOUBLE_COMPLEX), pointer :: idyt(:)
type(C_PTR) :: pidyt

complex(C_DOUBLE_COMPLEX), pointer :: k(:)
type(C_PTR) :: pk

double precision :: w

type(C_PTR) plan
type(C_PTR) plan2
type(C_PTR) planinv
type(C_PTR) plan2inv
  !double  complex :: in, out
complex(C_DOUBLE_COMPLEX), pointer :: in(:)
type(C_PTR) :: pin

complex(C_DOUBLE_COMPLEX), pointer :: out(:)
type(C_PTR) :: pout

complex(C_DOUBLE_COMPLEX), pointer :: teste(:,:)
type(C_PTR) :: pteste
!        double precision y
!        dimension y(N)
!        double complex out
!        dimension out(N/2+1)
pteste = fftw_alloc_complex(int(N*N,C_SIZE_T))
call c_f_pointer(pteste, teste, [N,N])

px = fftw_alloc_real(int(N,C_SIZE_T))
call c_f_pointer(px, x, [N])

py = fftw_alloc_real(int(N,C_SIZE_T))
call c_f_pointer(py, y, [N])

pyt = fftw_alloc_complex(int(N,C_SIZE_T))
call c_f_pointer(pyt, yt, [N])

pdy = fftw_alloc_real(int(N,C_SIZE_T))
call c_f_pointer(pdy, dy, [N])

pk = fftw_alloc_complex(int(N,C_SIZE_T))
call c_f_pointer(pk, k, [N])

pidyt = fftw_alloc_complex(int(N,C_SIZE_T))
call c_f_pointer(pidyt, idyt, [N])

pin = fftw_alloc_complex(int(N,C_SIZE_T))
call c_f_pointer(pin, in, [N])

pout = fftw_alloc_complex(int(N,C_SIZE_T))
call c_f_pointer(pout, out, [N])

  plan = fftw_plan_dft_1d( N, in, out, FFTW_FORWARD, FFTW_MEASURE )
  plan2 = fftw_plan_dft_r2c_1d( N, y, yt, FFTW_MEASURE )
  planinv = fftw_plan_dft_1d( N, in, out, FFTW_BACKWARD, FFTW_MEASURE )
  plan2inv = fftw_plan_dft_c2r_1d( N, in, y, FFTW_MEASURE )


x = (/(i*dx, i=0,(N-1),1)/)

  y = sin(x)
  dy = cos(x)

  k(1) = 0
  k(N2+1) = 0

  do i= 2,N2
    w = 2*PI/(b-a)
	k(i) = dcmplx(float(0),(i-1)*w)
	k(i+N2) = dcmplx(float(0),(i-1-N2)*w)
  enddo
  xteste = -1.0
  xteste = sqrt(xteste)
  k(4) = xteste
  xcomplex = cmplx(float(0),xteste)
  print *, imagpart(xcomplex)
  if (isnan(imagpart(xcomplex))) stop '"x" is a NaN'

!  yt = cmplx(float(0),float(0))
!  call fftw_execute_dft_r2c( plan2, y,yt )
!  in = (yt*k)/N
!  y = 0
!  call fftw_execute_dft_c2r( plan2inv, in, y )
!
!
!!  write(*,*) 'derivada espectral | derivada real'
!!  do i = 1,N,1
!!    write(*,*) '    ',N,' * y(',i,') = ',y(i) , ' | ' , dy(i)
!!  enddo
!
!  !y = out
!
!  dev = sum(abs(dy-y))/N
!  print *, "dev= ", dev


  in = y
  call fftw_execute_dft( plan, in ,out )
  in = (out*k)/N
  call fftw_execute_dft( planinv, in,out )
  y=out
  dev = sum(abs(dy-y))/N
  print *, "dev= ", dev

  call fftw_destroy_plan( plan )
  call fftw_destroy_plan( plan2 )
  call fftw_destroy_plan ( planinv )
  call fftw_destroy_plan ( plan2inv )
  end
