program arraycons
  implicit none
  integer, parameter :: N=8
  real, parameter :: a=0
  real, parameter :: b=4
  real, parameter :: L=(b-a)
  real, parameter :: dx=L/N
  integer :: i
  real :: x(N) = (/(i*dx, i=0,(N-1),1)/)
  print *, x
end program arraycons
