subroutine comp_next ( n, k, iarray, more )
!
!*******************************************************************************
!
!! COMP_NEXT computes the compositions of the integer N into K parts.
!
!
!  Discussion:
!
!    The routine computes one composition on each call until there are no more.
!    For instance, one composition of 6 into 3 parts is
!    3+2+1, another would be 6+0+0.
!
!  Example:
!
!    The 28 compositions of 6 into three parts are:
!
!      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
!      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
!      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
!      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
!      0 3 3,  0 2 4,  0 1 5,  0 0 6.
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the integer whose compositions are desired.
!
!    Input, integer K, the number of parts in the composition.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th part of
!    the composition.
!
!    Input/output, logical MORE.
!    Set MORE = .FALSE. on first call.  It will be reset to .TRUE. on return
!    with a new composition.  Each new call returns another composition until
!    MORE is set to .FALSE. when the last composition has been computed
!    and returned.
!
  implicit none
!
  integer k
!
  integer i
  integer iarray(k)
  integer, save :: ih = 0
  integer, save :: it = 0
  logical more
  integer n
!
  if ( .not. more ) then

    it = n
    ih = 0
    iarray(1) = n
    iarray(2:k) = 0

  else

    if ( it > 1 ) then
      ih = 0
    end if

    ih = ih + 1
    it = iarray(ih)
    iarray(ih) = 0
    iarray(1) = it - 1
    iarray(ih+1) = iarray(ih+1) + 1

  end if

  more = iarray(k) /= n

  return
end subroutine comp_next
