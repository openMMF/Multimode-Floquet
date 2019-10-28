!!$program main
!!$  implicit none
!!$  double precision :: pi,a,b,alpha
!!$  COMPLEX*16, DIMENSION(:), ALLOCATABLE:: in,out
!!$  integer i,n,q,Ly
!!$
!!$  n = 512
!!$  allocate(in(n))
!!$  allocate(out(n))
!!$
!!$  call timestamp ( )
!!$  
!!$ ! write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) '#FFTW3_PRB'
!!$  write ( *, '(a)' ) '#  FORTRAN77 version'
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#  Demonstrate the use of the FFTW3 library.'
!!$  alpha = 1.0/7.0
!!$  pi = 4.0*atan(1.0)
!!$  q  = 12
!!$  Ly = 256
!!$  do i = 1, n
!!$     a =  cos(cos(2*pi*i*alpha + 2*pi*q/Ly))!r8_uniform_01 ( seed )
!!$     b = -sin(cos(2*pi*i*alpha + 2*pi*q/Ly))!r8_uniform_01 ( seed )
!!$     in(i) = dcmplx ( a, b )
!!$  end do
!!$  
!!$  call FFTWinFORTRAN ( in , out,  n )
!!$
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#FFTW3_PRB'
!!$  write ( *, '(a)' ) '#  Normal end of execution.'
!!$  
!!$  write ( *, '(a)' ) '# '
!!$  call timestamp ( )
!!$  
!!$  stop
!!$end program main

subroutine FFTWinFORTRAN( in, out,n)
  
!n : vector length
!in: COMPLEX VECTOR IN
!out: fourier components


  implicit none
  
  integer,intent(in) :: n  !
  complex*16,dimension(n),intent(in)  :: in
  complex*16,dimension(n),intent(out) :: out
  !parameter ( n = 1024 )
  
  include "fftw3.f"
  
  double precision a
  double precision b
  double precision r8_uniform_01,pi
  integer i
  complex*16 in2(n)

  integer*8 plan_backward
  integer*8 plan_forward
  integer seed
  
  seed = 123456789
  
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '# TEST01'
!!$  write ( *, '(a)' ) '#  Demonstrate FFTW3 on a single vector '
!!$  write ( *, '(a)' ) '#  of complex data.'
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#  Transform data to FFT coefficients.'
!!$  write ( *, '(a)' ) '#  Backtransform FFT coefficients to recover'
!!$  write ( *, '(a)' ) '#  the data.'
!!$  write ( *, '(a)' ) '#  Compare recovered data to original data.'

!!$  pi = 4.0*atan(1.0)
!!$  do i = 1, n
!!$     a = cos(2*pi*i/n) + 0.5*cos(4*pi*i/n) + 0.1*sin(3*pi*i/n) + 2.0 + 2.0*cos((n-1)*2*pi*i/(4*n))!r8_uniform_01 ( seed )
!!$     b = 0.0 !r8_uniform_01 ( seed )
!!$     in(i) = dcmplx ( a, b )
!!$  end do
  
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#  Input Data:'
!!$  write ( *, '(a)' ) '# '
!!$  
!!$  do i = 1, n
!!$     write ( *, '(2x,i4,2x,2g14.6)' ) i, in(i)
!!$  end do
!!$  write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) ' '

  call dfftw_plan_dft_1d ( plan_forward, n, in, out, &
       &  FFTW_FORWARD, FFTW_ESTIMATE )
  
  call dfftw_execute ( plan_forward )
  out = out/n

!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#  Output FFT Coefficients:'
!!$  write ( *, '(a)' ) '# '
!!$  
!!$  do i = 1, n
!!$     write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)/n
!!$  end do
!!$  write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) ' '

!!$  call dfftw_plan_dft_1d ( plan_backward, n, out, in2, &
!!$       &  FFTW_BACKWARD, FFTW_ESTIMATE )
!!$  
!!$  call dfftw_execute ( plan_backward )
!!$  
!!$  write ( *, '(a)' ) '# '
!!$  write ( *, '(a)' ) '#  Recovered input data divided by N:'
!!$  write ( *, '(a)' ) '# '
!!$  
!!$  do i = 1, n
!!$     write ( *, '(2x,i4,2x,2g14.6)' ) i, in2(i) / dble ( n )
!!$  end do
!!$  write ( *, '(a)' ) ' '
!!$  write ( *, '(a)' ) ' '

  call dfftw_destroy_plan ( plan_forward )
!!$  call dfftw_destroy_plan ( plan_backward )
  
  return
end subroutine FFTWinFORTRAN

function r8_uniform_01 ( seed )

  implicit none
  
  double precision r8_uniform_01
  integer k
  integer seed
  
  k = seed / 127773
  
  seed = 16807 * ( seed - k * 127773 ) - k * 2836
  
  if ( seed < 0 ) then
     seed = seed + 2147483647
  end if
  
  r8_uniform_01 = dble ( seed ) * 4.656612875D-10
  
  return
end function r8_uniform_01

subroutine timestamp ( )

  implicit none
  
  character ( len = 8 ) date
  character ( len = 10 ) time
  
  call date_and_time ( date, time )
  
!  write ( *, '(a8,2x,a10)' ) date, time
  
  return
end subroutine timestamp
