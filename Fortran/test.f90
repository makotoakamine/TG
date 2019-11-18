program test
use, intrinsic :: iso_c_binding 
use fftw3
!include fftw3
implicit none

real*16 , parameter :: PI=4.D0*DATAN(1.D0)
integer, parameter :: N=16
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


  plan = fftw_plan_dft_1d( N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE )

  call fftw_execute_dft( plan, out , in )

  write(*,*) 'derivada espectral'
  do i = 1,N,1
    write(*,*) '    ',N,' * in(',i,') = ',in(i)/N
  enddo

  write(*,*) 'derivada real'
  print *, dy

  call fftw_destroy_plan ( plan )
  dev = sum(dy-(in(i)/N))
  print *, "dev= ", dev

  end
