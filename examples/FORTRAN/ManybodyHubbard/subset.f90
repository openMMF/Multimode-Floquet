subroutine asm_enum ( n, asm_num )
!
!*******************************************************************************
!
!! ASM_ENUM returns the number of alternating sign matrices of a given order.
!
!
!  Discussion:
!
!    N     ASM_NUM
!
!    0       1
!    1       1
!    2       2
!    3       7
!    4      42
!    5     429
!    6    7436
!    7  218348
!
!    A direct formula is
!
!      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrices.
!
!    Output, integer ASM_NUM, the number of alternating sign matrices of
!    order N.
!
  implicit none
!
  integer n
!
  integer a(n+1)
  integer asm_num
  integer b(n+1)
  integer c(n+1)
  integer i
  integer nn
!
  asm_num = 0

  nn = n + 1

  if ( n+1 <= 0 ) then
    return
  end if
!
!  Row 1
!
  a(1) = 1

  if ( n+1 == 1 ) then
    asm_num = a(1)
    return
  end if
!
!  Row 2
!
  a(1) = 1

  if ( n+1 == 2 ) then
    asm_num = a(1)
    return
  end if

  b(1) = 2
  c(1) = 2
  a(2) = 1
!
!  Row 3 and on.
!
  do nn = 3, n

    b(nn-1) = nn
    do i = nn-2, 2, -1
      b(i) = b(i) + b(i-1)
    end do
    b(1) = 2

    c(nn-1) = 2
    do i = nn-2, 2, -1
      c(i) = c(i) + c(i-1)
    end do
    c(1) = nn

    a(1) = sum ( a(1:nn-1) )
    do i = 2, nn
      a(i) = a(i-1) * c(i-1) / b(i-1)
    end do

  end do

  asm_num = sum ( a(1:n) )

  return
end
subroutine asm_triangle ( n, a )
!
!*******************************************************************************
!
!! ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
!
!
!  Discussion:
!
!    The first seven rows of the triangle are as follows:
!
!          1      2      3      4      5      6     7
!
!    1     1
!    2     1      1
!    3     2      3      2
!    4     7     14     14      7
!    5    42    105    135    105     42
!    6   429   1287   2002   2002   1287    429
!    7  7436  26026  47320  56784  47320  26026  7436
!
!    For a given N, the value of A(J) represents entry A(I,J) of
!    the triangular matrix, and gives the number of alternating sign matrices
!    of order N in which the (unique) 1 in row 1 occurs in column J.
!
!    Thus, of alternating sign matrices of order 3, there are
!    2 with a leading 1 in column 1:
!
!      1 0 0  1 0 0
!      0 1 0  0 0 1
!      0 0 1  0 1 0
!
!    3 with a leading 1 in column 2, and
!
!      0 1 0  0 1 0  0 1 0
!      1 0 0  0 0 1  1-1 1
!      0 0 1  1 0 0  0 1 0
!
!    2 with a leading 1 in column 3:
!
!      0 0 1  0 0 1
!      1 0 0  0 1 0
!      0 1 0  1 0 0
!
!  Modified:
!
!    01 January 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the desired row.
!
!    Output, integer A(N), the entries of the row.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer b(n)
  integer c(n)
  integer i
  integer nn
!
  if ( n <= 0 ) then
    return
  end if
!
!  Row 1
!
  a(1) = 1

  if ( n == 1 ) then
    return
  end if
!
!  Row 2
!
  nn = 2
  b(1) = 2
  c(1) = nn

  a(1) = sum ( a(1:nn-1) )
  do i = 2, nn
    a(i) = a(i-1) * c(i-1) / b(i-1)
  end do

  if ( n == 2 ) then
    return
  end if
!
!  Row 3 and on.
!
  do nn = 3, n

    b(nn-1) = nn
    do i = nn-2, 2, -1
      b(i) = b(i) + b(i-1)
    end do
    b(1) = 2

    c(nn-1) = 2
    do i = nn-2, 2, -1
      c(i) = c(i) + c(i-1)
    end do
    c(1) = nn

    a(1) = sum ( a(1:nn-1) )
    do i = 2, nn
      a(i) = a(i-1) * c(i-1) / b(i-1)
    end do

  end do

  return
end
subroutine bell ( b, n )
!
!*******************************************************************************
!
!! BELL returns the Bell numbers from 0 to N.
!
!
!  Discussion:
!
!    The Bell number B(N) is the number of restricted growth functions
!    on N.
!
!    Note that the Stirling numbers of the second kind, S^m_n, count the
!    number of partitions of N objects into M classes, and so it is
!    true that
!
!      B(N) = S^1_N + S^2_N + ... + S^N_N.
!
!  Definition:
!
!    The Bell number B(N) is defined as the number of partitions (of
!    any size) of a set of N distinguishable objects.
!
!    A partition of a set is a division of the objects of the set into
!    subsets.
!
!  Examples:
!
!    There are 15 partitions of a set of 4 objects:
!
!      (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
!      (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
!      (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
!
!    and so B(4) = 15.
!
!  First values:
!
!     N         B(N)
!     0           1
!     1           1
!     2           2
!     3           5
!     4          15
!     5          52
!     6         203
!     7         877
!     8        4140
!     9       21147
!    10      115975
!
!  Recursion:
!
!    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer B(0:N), the Bell numbers from 0 to N.
!
!    Input, integer N, the number of Bell numbers desired.
!
  implicit none
!
  integer n
!
  integer b(0:n)
  integer combo
  integer i
  integer j
!
  b(0) = 1

  do i = 1, n
    b(i) = 0
    do j = 1, i
      call combin2 ( i-1, j-1, combo )
      b(i) = b(i) + combo * b(i-j)
    end do
  end do

  return
end
subroutine catalan ( n, c )
!
!*******************************************************************************
!
!! CATALAN computes the Catalan numbers, from C(0) to C(N).
!
!
!  First values:
!
!     C(0)     1
!     C(1)     1
!     C(2)     2
!     C(3)     5
!     C(4)    14
!     C(5)    42
!     C(6)   132
!     C(7)   429
!     C(8)  1430
!     C(9)  4862
!    C(10) 16796
!
!  Formula:
!
!    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
!         = 1 / (N+1) * COMB ( 2N, N )
!         = 1 / (2N+1) * COMB ( 2N+1, N+1).
!
!  Recursion:
!
!    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
!    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
!
!  Comments:
!
!    The Catalan number C(N) counts:
!
!    1) the number of binary trees on N vertices;
!    2) the number of ordered trees on N+1 vertices;
!    3) the number of full binary trees on 2N+1 vertices;
!    4) the number of well formed sequences of 2N parentheses;
!    5) number of ways 2N ballots can be counted, in order,
!       with N positive and N negative, so that the running sum
!       is never negative;
!    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
!    7) the number of monotone functions from [1..N} to [1..N} which
!       satisfy f(i) <= i for all i;
!    8) the number of ways to triangulate a polygon with N+2 vertices.
!
!  Example:
!
!    N = 3
!
!    ()()()
!    ()(())
!    (()())
!    (())()
!    ((()))
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    14 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of Catalan numbers desired.
!
!    Output, integer C(0:N), the Catalan numbers from C(0) to C(N).
!
  implicit none
!
  integer n
!
  integer c(0:n)
  integer i
!
  c(0) = 1
!
!  The extra parentheses ensure that the integer division is
!  done AFTER the integer multiplication.
!
  do i = 1, n
    c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
  end do

  return
end
subroutine cfrac_to_rat ( a, n, p, q )
!
!*******************************************************************************
!
!! CFRAC_TO_RAT converts a monic continued fraction to an ordinary fraction.
!
!
!  Discussion:
!
!    The routine is given the monic or "simple" continued fraction with
!    integer coefficients:
!
!      A(1) + 1 / ( A(2) + 1 / ( A(3) ... + 1 / A(N) ) )
!
!    and returns the N successive approximants P(I)/Q(I)
!    to the value of the rational number represented by the continued
!    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    01 October 2000
!
!  Parameters:
!
!    Input, integer A(N), the continued fraction coefficients.
!
!    Input, integer N, the number of continued fraction coefficients.
!
!    Output, integer P(N), Q(N), the N successive approximations
!    to the value of the continued fraction.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  integer p(n)
  integer q(n)
!
  do i = 1, n

    if ( i == 1 ) then
      p(i) = a(i) * 1 + 0
      q(i) = a(i) * 0 + 1
    else if ( i == 2 ) then
      p(i) = a(i) * p(i-1) + 1
      q(i) = a(i) * q(i-1) + 0
    else
      p(i) = a(i) * p(i-1) + p(i-2)
      q(i) = a(i) * q(i-1) + q(i-2)
    end if

  end do

  return
end
subroutine cfrac_to_rfrac ( g, h, m, p, q )
!
!*******************************************************************************
!
!! CFRAC_TO_RFRAC converts a polynomial fraction from continued to rational form.
!
!
!  Discussion:
!
!    The routine accepts a continued polynomial fraction:
!
!      G(1)     / ( H(1) +
!      G(2) * X / ( H(2) +
!      G(3) * X / ( H(3) + ...
!      G(M) * X / ( H(M) )...) ) )
!
!    and returns the equivalent rational polynomial fraction:
!
!      P(1) + P(2) * X + ... + P(L1) * X**(L1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(L2) * X**(L2-1)
!
!    where
!
!      L1 = (M+1)/2
!      L2 = (M+2)/2.
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    17 April 2000
!
!  Parameters:
!
!    Input, real G(M), H(M), the continued polynomial fraction coefficients.
!
!    Input, integer M, the number of continued fraction polynomial coefficients.
!
!    Output, real P((M+1)/2), Q((M+2)/2), the rational polynomial fraction
!    coefficients.
!
  implicit none
!
  integer m
!
  real a(m,(m+2)/2)
  real g(m)
  real h(m)
  integer i
  integer j
  real p((m+1)/2)
  real q((m+2)/2)
!
  if ( m == 1 ) then
    p(1) = g(1)
    q(1) = h(1)
    return
  end if

  do i = 1, m
    do j = 1, (m+2)/2
      a(i,j) = 0.0E+00
    end do
  end do
!
!  Solve for P's
!
  a(1,1) = g(1)
  a(2,1) = g(1) * h(2)

  do i = 3, m
    a(i,1) = h(i) * a(i-1,1)
    do j = 2, (i+1)/2
      a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
    end do
  end do

  do i = 1, (m+1)/2
    p(i) = a(m,i)
  end do

  a(1,1) = h(1)
  a(2,1) = h(1) * h(2)
  a(2,2) = g(2)

  do i = 3, m
    a(i,1) = h(i) * a(i-1,1)
    do j = 2, (i+2) / 2
      a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
    end do
  end do

  do i = 1, (m+2)/2
    q(i) = a(m,i)
  end do

  return
end
subroutine ch_cap ( c )
!
!*******************************************************************************
!
!! CH_CAP capitalizes a single character.
!
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none
!
  character c
  integer itemp
!
  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine chinese_check ( n, m, ierror )
!
!*******************************************************************************
!
!! CHINESE_CHECK checks the Chinese remainder moduluses.
!
!
!  Discussion:
!
!    For a Chinese remainder representation, the moduluses M(I) must
!    be positive and pairwise prime.  Also, in case this is not obvious,
!    no more than one of the moduluses may be 1.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of moduluses.
!
!    Input, integer M(N), the moduluses.  These should be positive
!    and pairwise prime.
!
!    Output, integer IERROR, an error flag.
!    0, no error was detected.
!    nonzero, an error was detected.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  logical ivec_pairwise_prime
  integer j
  integer m(n)
!
  ierror = 0
!
!  Do not allow nonpositive entries.
!
  if ( any ( m(1:n) <= 0 ) ) then
    ierror = 1
    return
  end if
!
!  Allow one entry to be 1, but not two entries.
!
  do i = 1, n
    do j = i+1, n
      if ( m(i) == 1 .and. m(j) == 1 ) then
        ierror = 2
        return
      end if
    end do
  end do
!
!  Now check pairwise primeness.
!
  if ( .not. ivec_pairwise_prime ( n, m ) ) then
    ierror = 3
    return
  end if

  return
end
subroutine chinese_to_i ( n, m, r, j )
!
!*******************************************************************************
!
!! CHINESE_TO_I converts a set of Chinese remainders to an equivalent integer.
!
!
!  Discussion:
!
!    Given a set of N pairwise prime, positive moduluses M(I), and
!    a corresponding set of remainders M(I), this routine finds an
!    integer J such that, for all I,
!
!      J = R(I) mod M(I)
!
!  Modified:
!
!    10 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of moduluses.
!
!    Input, integer M(N), the moduluses.  These should be positive
!    and pairwise prime.
!
!    Input, integer R(N), the Chinese remainder representation of the integer.
!
!    Output, integer J, the corresponding integer.
!
  implicit none
!
  integer n
!
  integer a
  integer b(n)
  integer big_m
  integer c
  integer i
  integer ierror
  integer j
  integer m(n)
  integer r(n)
!
  call chinese_check ( n, m, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHINESE_TO_I - Fatal error!'
    write ( *, '(a)' ) '  The moduluses are not legal.'
    stop
  end if
!
!  Set BIG_M.
!
  big_m = product ( m )
!
!  Solve BIG_M / M(I) * B(I) = 1, mod M(I)
!
  do i = 1, n
    a = big_m / m(i)
    c = 1
    call congruence ( a, m(i), c, ierror, b(i) )
  end do
!
!  Set J = sum ( 1 <= I <= N ) ( R(I) * B(I) * BIG_M / M(I) ) mod M
!
  j = 0
  do i = 1, n
    j = mod ( j + r(i) * b(i) * ( big_m / m(i) ), big_m )
  end do

  return
end
subroutine comb ( m, n, l, iarray )
!
!*******************************************************************************
!
!! COMB returns the L-th combination of N things out of M.
!
!
!  Discussion:
!
!    The combinations are ordered lexically.
!
!    Lexical order can be illustrated for the general case of N and M as
!    follows:
!
!    1:       1,     2,     3,     ..., N-2, N-1, N
!    2:       1,     2,     3,     ..., N-2, N-1, N+1
!    3:       1,     2,     3,     ..., N-2, N-1, N+2
!    ...
!    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
!    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
!    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
!    ...
!    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
!    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
!    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
!
!    There are a total of M!/(N!*(M-N)!) combinations of M
!    things taken N at a time.
!
!  Reference:
!
!    B P Buckles and M Lybanon,
!    Algorithm 515,
!    Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, pages 180-182, June 1977.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer M, the size of the set.
!
!    Input, integer N, the number of things in the combination.
!    N must be greater than 0, and no greater than M.
!
!    Input, integer L, the lexicographical index of combination
!    sought.  L must be at least 1, and no greater than M!/(N!*(M-N)!).
!
!    Output, integer IARRAY(N), array containing the combination set.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer j
  integer k
  integer l
  integer m
!
!  Initialize lower bound index at zero.
!
  k = 0
!
!  Loop to select elements in ascending order.
!
  do i = 1, n - 1
!
!  Set lower bound element number for next element value.
!
    iarray(i) = 0

    if ( i /= 1 ) then
      iarray(i) = iarray(i-1)
    end if
!
!  Check each element value.
!
    do

      iarray(i) = iarray(i) + 1
      call combin2 ( m-iarray(i), n-i, j )
      k = k + j

      if ( k >= l ) then
        exit
      end if

    end do

    k = k - j

  end do

  iarray(n) = iarray(n-1) + l - k

  return
end
subroutine comb_next ( done, n, k, iarray )
!
!*******************************************************************************
!
!! COMB_NEXT computes combinations of K things out of N.
!
!
!  Discussion:
!
!    The combinations are computed one at a time, in lexicographical order.
!
!  Reference:
!
!    Charles Mifsud,
!    Combination in Lexicographic Order,
!    ACM algorithm 154,
!    March 1963.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input/output, logical DONE, indicator.
!    On input, if this is the first call, the user should set
!    DONE to FALSE.  On each subsequent call, the value of
!    DONE should be left at its output value from the previous
!    call.
!
!    On output, if DONE is TRUE, then a new combination was
!    computed and returned.  If DONE is FALSE, then the list
!    of combinations was finished on the previous call.
!
!    Input, integer N, the total number of things.
!
!    Input, integer K, the number of things in each combination.
!
!    Output, integer IARRAY(K), contains the list of elements in
!    the current combination.
!
  implicit none
!
  integer k
!
  logical done
  integer i
  integer iarray(k)
  integer j
  integer n
!
  if ( done ) then

    call ivec_identity ( k, iarray )

    if ( k > 1 ) then
      done = .false.
    else
      done = .true.
    end if

  else

    if ( iarray(k) < n ) then
      iarray(k) = iarray(k) + 1
      return
    end if

    do i = k, 2, -1

      if ( iarray(i-1) < n-k+i-1 ) then

        iarray(i-1) = iarray(i-1) + 1

        do j = i, k
          iarray(j) = iarray(i-1) + j - ( i-1 )
        end do

        return

      end if

    end do

    done = .true.

  end if

  return
end
subroutine comb_row ( ido, n, irow )
!
!*******************************************************************************
!
!! COMB_ROW computes a row of the combinatorial coefficients.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IDO, startup information.
!
!    0, build up the table from scratch.
!    1, N = 0, or IROW contains the previous row of combinatorial coefficients.
!
!    Input, integer N, the desired row.  N must be 0 or greater.
!
!    Output, integer IROW(N+1), the combinatorial coefficients.
!
  implicit none
!
  integer n
!
  integer i
  integer ido
  integer irow(n+1)
  integer k
!
  if ( n < 0 ) then
    return
  end if

  if ( ido == 1 ) then

    do i = n, 2, -1
      irow(i) = irow(i) + irow(i-1)
    end do

    irow(n+1) = 1

  else

    irow(1) = 1
    irow(2:n+1) = 0

    do k = 1, n
      do i = k+1, 2, -1
        irow(i) = irow(i) + irow(i-1)
      end do
    end do

  end if

  return
end
subroutine combin ( n, k, cnk )
!
!*******************************************************************************
!
!! COMBIN computes the combinatorial coefficient C(N,K).
!
!
!  Method:
!
!    Real arithmetic is used, and C(N,K) is computed directly, via
!    Gamma functions, rather than recursively.
!
!  Definition:
!
!    C(N,K) is the number of distinct combinations of K objects
!    chosen from a set of N distinct objects.  A combination is
!    like a set, in that order does not matter.
!
!  Examples:
!
!    The number of combinations of 2 things chosen from 5 is 10.
!
!    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
!
!    The actual combinations may be represented as:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3),
!      (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Formula:
!
!    C(N,K) = N! / ( (N-K)! * K! )
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value of N.
!
!    Input, integer K, the value of K.
!
!    Output, real CNK, the value of C(N,K)
!
  implicit none
!
  real arg
  real cnk
  real fack
  real facn
  real facnmk
  real gamma_log
  integer k
  integer n
!
  if ( n < 0 ) then

    cnk = 0.0E+00

  else if ( k == 0 ) then

    cnk = 1.0E+00

  else if ( k == 1 ) then

    cnk = real ( n )

  else if ( k > 1 .and. k < n-1 ) then

    arg = real ( n + 1 )
    facn = gamma_log ( arg )

    arg = real ( k + 1 )
    fack = gamma_log ( arg )

    arg = real ( n - k + 1 )
    facnmk = gamma_log ( arg )

    cnk = anint ( exp ( facn - fack - facnmk ) )

  else if ( k == n-1 ) then

    cnk = real ( n )

  else if ( k == n ) then

    cnk = 1.0E+00

  else

    cnk = 0.0E+00

  end if

  return
end
subroutine combin2 ( n, k, icnk )
!
!*******************************************************************************
!
!! COMBIN2 computes the binomial coefficient C(N,K).
!
!
!  Discussion:
!
!    The value is calculated in such a way as to avoid overflow and
!    roundoff.  The calculation is done in integer arithmetic.
!
!  Formula:
!
!    ICNK = C(N,K) = N! / ( K! * (N-K)! )
!
!  Reference:
!
!    M L Wolfson and H V Wright,
!    Combinatorial of M Things Taken N at a Time,
!    ACM algorithm 160,
!    Communications of the ACM,
!    April, 1963.
!
!  Modified:
!
!    17 January 1999
!
!  Parameters:
!
!    Input, integer N, K, are the values of N and K.
!
!    Output, integer ICNK, the number of combinations of N
!    things taken K at a time.
!
  implicit none
!
  integer i
  integer icnk
  integer k
  integer mn
  integer mx
  integer n
!
  mn = min ( k, n-k )

  if ( mn < 0 ) then

    icnk = 0

  else if ( mn == 0 ) then

    icnk = 1

  else

    mx = max ( k, n-k )
    icnk = mx + 1

    do i = 2, mn
      icnk = ( icnk * ( mx + i ) ) / i
    end do

  end if

  return
end
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
end
subroutine comp_random ( n, k, iarray )
!
!*******************************************************************************
!
!! COMP_RANDOM selects a random composition of the integer N into K parts.
!
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
!    Input, integer N, the integer to be decomposed.
!
!    Input, integer K, the number of parts in the composition.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th part of the
!    composition.
!
  implicit none
!
  integer k
!
  integer i
  integer iarray(k)
  integer l
  integer m
  integer n
!
  call ksub_random ( n+k-1, k-1, iarray )

  iarray(k) = n + k
  l = 0

  do i = 1, k
    m = iarray(i)
    iarray(i) = iarray(i) - l - 1
    l = m
  end do

  return
end
subroutine congruence ( a, b, c, ierror, x )
!
!*******************************************************************************
!
!! CONGRUENCE solves a congruence of the form A * X = C mod ( B ).
!
!
!  Discussion:
!
!    A, B and C are given integers.  The equation is solvable if and only
!    if the greatest common divisor of A and B also divides C.
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Modified:
!
!    09 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, C, the coefficients of the Diophantine equation.
!
!    Output, integer IERROR, error flag.
!    0, no error, X was computed.
!    1, A = B = 0, C is nonzero.
!    2, A = 0, B and C nonzero, but C is not a multiple of B.
!    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
!    4, A, B, C nonzero, but GCD of A and B does not divide C.
!    5, algorithm ran out of internal space.
!
!    Output, integer X, the solution of the Diophantine equation.
!    X will be between 0 and B-1.
!
  implicit none
!
  integer, parameter :: nmax = 100
!
  integer a
  integer a_copy
  integer a_mag
  integer a_sign
  integer b
  integer b_copy
  integer b_mag
  integer b_sign
  integer c
  integer c_copy
  integer g
  integer i_gcd
  integer ierror
  integer k
  integer n
  real norm_new
  real norm_old
  integer q(nmax)
  logical swap
  integer temp
  integer x
  integer xnew
  integer y
  integer ynew
!
!  Defaults for output parameters.
!
  ierror = 0
  x = 0
  y = 0
!
!  Special cases.
!
  if ( a == 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b == 0 .and. c /= 0 ) then
    ierror = 1
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c /= 0 ) then
    x = 0
    if ( mod ( c, b ) /= 0 ) then
      ierror = 2
    end if
    return
  else if ( a /= 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    return
  else if ( a /= 0 .and. b == 0 .and. c /= 0 ) then
    x = c / a
    if ( mod ( c, a ) /= 0 ) then
      ierror = 3
    end if
    return
  else if ( a /= 0 .and. b /= 0 .and. c == 0 ) then
    g = i_gcd ( a, b )
    x = b / g
    return
  end if
!
!  Now handle the "general" case: A, B and C are nonzero.
!
!  Step 1: Compute the GCD of A and B, which must also divide C.
!
  g = i_gcd ( a, b )

  if ( mod ( c, g ) /= 0 ) then
    ierror = 4
    return
  end if

  a_copy = a / g
  b_copy = b / g
  c_copy = c / g
!
!  Step 2: Split A and B into sign and magnitude.
!
  a_mag = abs ( a_copy )
  a_sign = sign ( 1, a_copy )
  b_mag = abs ( b_copy )
  b_sign = sign ( 1, b_copy )
!
!  Another special case, A_MAG = 1 or B_MAG = 1.
!
  if ( a_mag == 1 ) then
    x = a_sign * c_copy
    return
  else if ( b_mag == 1 ) then
    x = 0
    return
  end if
!
!  Step 3: Produce the Euclidean remainder sequence.
!
  if ( a_mag >= b_mag ) then

    swap = .false.
    q(1) = a_mag
    q(2) = b_mag

  else

    swap = .true.
    q(1) = b_mag
    q(2) = a_mag

  end if

  n = 3

  do

    q(n) = mod ( q(n-2), q(n-1) )

    if ( q(n) == 1 ) then
      exit
    end if

    n = n + 1

    if ( n > nmax ) then
      ierror = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CONGRUENCE - Fatal error!'
      write ( *, '(a)' ) '  Exceeded number of iterations.'
      stop
    end if

  end do
!
!  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
!
  y = 0
  do k = n, 2, -1
    x = y
    y = ( 1 - x * q(k-1) ) / q(k)
  end do
!
!  Step 5: Undo the swapping.
!
  if ( swap ) then
    call i_swap ( x, y )
  end if
!
!  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
!
  x = x * a_sign
!
!  Step 7: Multiply by C, so that X * A + Y * B = C.
!
  x = x * c_copy
!
!  Step 8: Now force 0 <= X < B.
!
  x = mod ( x, b )

  return
end
subroutine count_pose ( blocks, goal )
!
!*******************************************************************************
!
!! COUNT_POSE poses a problem for the game "The Count is Good"
!
!
!  Discussion:
!
!    The French television show "The Count is Good" has a game that goes
!    as follows:
!
!      A number is chosen at random between 100 and 999.  This is the GOAL.
!
!      Six numbers are randomly chosen from the set 1, 2, 3, 4, 5, 6, 7, 8,
!      9, 10, 25, 50, 75, 100.  These numbers are the BLOCKS.
!
!      The player must construct a formula, using some or all of the blocks,
!      (but not more than once), and the operations of addition, subtraction,
!      multiplication and division.  Parentheses should be used to remove
!      all ambiguity.  However, it is forbidden to use subtraction in a
!      way that produces a negative result, and all division must come out
!      exactly, with no remainder.
!
!    This routine poses a sample problem from the show.  The point is,
!    to determine how to write a program that can solve such a problem.
!
!  Reference:
!
!    Raymond Seroul,
!    Programming for Mathematicians,
!    Springer Verlag, 2000, page 355-357.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer BLOCKS(6), the six numbers available for the formula.
!
!    Output, integer GOAL, the goal number.
!
  implicit none
!
  integer blocks(6)
  integer goal
  integer ind(6)
  integer, dimension ( 14 ) :: stuff = &
    (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100 /)
!
  call i_random ( 100, 999, goal )

  call ksub_random ( 14, 6, ind )

  blocks(1:6) = stuff(ind(1:6))

  return
end
function d_pi ( )
!
!*******************************************************************************
!
!! D_PI returns the value of pi as a double precision quantity.
!
!
!  Modified:
!
!    19 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision D_PI, the value of pi.
!
  implicit none
!
  double precision d_pi
!
  d_pi = 3.14159265358979323846264338327950288419716939937510D+00

  return
end
subroutine d_to_cfrac ( r, n, a, p, q )
!
!*******************************************************************************
!
!! D_TO_CFRAC converts a double precision value to a continued fraction.
!
!
!  Discussion:
!
!    The routine is given a real number R.  It computes a sequence of
!    continued fraction approximations to R, returning the results as
!    simple fractions of the form P(I) / Q(I).
!
!  Example:
!
!    X = 2 * PI
!    N = 7
!
!    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
!    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
!    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
!
!    (This ignores roundoff error, which will cause later terms to differ).
!
!  Reference:
!
!    Norman Richert,
!    Strang's Strange Figures,
!    American Mathematical Monthly,
!    Volume 99, Number 2, February 1992, pages 101-107.
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision R, the real value.
!
!    Input, integer N, the number of convergents to compute.
!
!    Output, integer A(0:N), the partial quotients.
!
!    Output, integer P(-1:N), Q(-1:N), the numerators and denominators
!    of the continued fraction approximations.
!
  implicit none
!
  integer n
!
  integer a(0:n)
  integer i
  integer p(-1:n)
  integer q(-1:n)
  double precision r
  double precision r_copy
  double precision x(0:n)
!
  if ( r == 0.0D+00 ) then
    a(0:n) = 0
    p(-1:n) = 0
    q(-1:n) = 1
    return
  end if

  r_copy = abs ( r )

  p(-1) = 1
  q(-1) = 0

  p(0) = int ( r_copy )
  q(0) = 1
  x(0) = r_copy
  a(0) = int ( x(0) )

  do i = 1, n
    x(i) = 1.0E+00 / ( x(i-1) - dble ( a(i-1) ) )
    a(i) = int ( x(i) )
    p(i) = a(i) * p(i-1) + p(i-2)
    q(i) = a(i) * q(i-1) + q(i-2)
  end do

  if ( r < 0.0D+00 ) then
    p(-1:n) = - p(-1:n)
  end if

  return
end
subroutine d_to_dec ( dval, itop, ibot )
!
!*******************************************************************************
!
!! D_TO_DEC converts a double precision quantity to a decimal representation.
!
!
!  Discussion:
!
!    Given the double precision value DVAL, the routine computes integers
!    ITOP and IBOT so that it is approximately true that:
!
!      DVAL = ITOP * 10 ** IBOT
!
!    In particular, only DEC_DIGIT digits of DVAL are used in constructing the
!    representation.
!
!  Modified:
!
!    25 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision DVAL, the value whose decimal representation
!    is desired.
!
!    Output, integer ITOP, IBOT, the approximate decimal representation of DVAL.
!    ITOP is an integer, strictly between -10**DEC_DIGIT and 10**DEC_DIGIT.
!    IBOT is an integer exponent of 10.
!
  implicit none
!
  integer dec_digit
  double precision dtop
  double precision dval
  integer ibot
  integer itop
  real ten1
  real ten2
!
!  Special cases.
!
  if ( dval == 0.0D+00 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  Factor DVAL = DTOP * 10**IBOT
!
  dtop = dval
  ibot = 0
!
!  Get the number of decimal digits.
!
  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )
!
!  Now normalize so that 10**(DEC_DIGIT-1) <= ABS(DTOP) < 10**(DEC_DIGIT)
!
  ten1 = 10.0E+00**( dec_digit - 1 )
  ten2 = 10.0E+00**dec_digit

  do while ( abs ( dtop ) < ten1 )
    dtop = dtop * 10.0D+00
    ibot = ibot - 1
  end do

  do while ( abs ( dtop ) >= ten2 )
    dtop = dtop / 10.0D+00
    ibot = ibot + 1
  end do
!
!  ITOP is the integer part of DTOP, rounded.
!
  itop = nint ( dtop )
!
!  Now divide out any factors of ten from ITOP.
!
  if ( itop /= 0 ) then
    do while ( 10 * ( itop / 10 ) == itop )
      itop = itop / 10
      ibot = ibot + 1
    end do
  end if

  return
end
subroutine debruijn ( m, n, string )
!
!*******************************************************************************
!
!! DEBRUIJN constructs a de Bruijn sequence.
!
!
!  Description:
!
!    Suppose we have an alphabet of M letters, and we are interested in
!    all possible strings of length N.  If M = 2 and N = 3, then we are
!    interested in the M**N strings:
!
!      000
!      001
!      010
!      011
!      100
!      101
!      110
!      111
!
!    Now, instead of making a list like this, we prefer, if possible, to
!    write a string of letters, such that every consecutive sequence of
!    N letters is one of the strings, and every string occurs once, if
!    we allow wraparound.
!
!    For the above example, a suitable sequence would be the 8 characters:
!
!      00011101(00...
!
!    where we have suggested the wraparound feature by repeating the first
!    two characters at the end.
!
!    Such a sequence is called a de Bruijn sequence.  It can easily be
!    constructed by considering a directed graph, whose nodes are all
!    M**(N-1) strings of length N-1.  A node I has a directed edge to
!    node J (labeled with character K) if the string at node J can
!    be constructed by beheading the string at node I and adding character K.
!
!    In this setting, a de Bruijn sequence is simply an Eulerian circuit
!    of the directed graph, with the edge labels being the entries of the
!    sequence.  In general, there are many distinct de Bruijn sequences
!    for the same parameter M and N.  This program will only find one
!    of them.
!
!  Modified:
!
!    06 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of letters in the alphabet.
!
!    Input, integer N, the number of letters in a codeword.
!
!    Output, integer STRING(M**N), a deBruijn string.
!
  implicit none
!
  integer m
  integer n
!
  integer i
  integer iedge
  integer inode(m**n)
  integer ivec(n-1)
  integer j
  integer jnode(m**n)
  integer jvec(n-1)
  integer k
  integer knode(m**n)
  integer nedge
  integer nnode
  logical success
  integer string(m**n)
  integer trail(m**n)
!
!  Construct the adjacency information.
!
  nnode = m**(n-1)
  nedge = m**n

  iedge = 0

  do i = 1, nnode

    call index_unrank0 ( n-1, m, i, ivec )

    do k = 1, m
      jvec(1:n-2) = ivec(2:n-1)
      jvec(n-1) = k
      call index_rank0 ( n-1, m, jvec, j )
      iedge = iedge + 1
      inode(iedge) = i
      jnode(iedge) = j
      knode(iedge) = k
    end do

  end do
!
!  Determine a circuit.
!
  call digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )
!
!  The string is constructed from the labels of the edges in the circuit.
!
  string(1:nedge) = knode(trail(1:nedge))

  return
end
subroutine dec_add ( itop1, ibot1, itop2, ibot2, itop, ibot )
!
!*******************************************************************************
!
!! DEC_ADD adds two decimal quantities.
!
!
!  Discussion:
!
!    The routine computes
!
!      ITOP * 10**IBOT = ITOP1 * 10**IBOT1 + ITOP2 * 10**IBOT2
!
!    while trying to avoid integer overflow.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the first number to be added.
!
!    Input, integer ITOP2, IBOT2, the second number to be added.
!
!    Output, integer ITOP, IBOT, the sum.
!
  implicit none
!
  integer ibot
  integer ibot1
  integer ibot2
  integer itop
  integer itop1
  integer itop2
  integer jtop1
  integer jtop2
!
  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  else if ( ibot1 == ibot2 ) then
    itop = itop1 + itop2
    ibot = ibot1
    call dec_round ( itop, ibot )
    return
  end if
!
!  Line up the exponents.
!
  jtop1 = itop1
  jtop2 = itop2

  if ( ibot1 < ibot2 ) then
    jtop2 = jtop2 * 10**( ibot2 - ibot1 )
  else
    jtop1 = jtop1 * 10**( ibot1 - ibot2 )
  end if
!
!  Add the coefficients.
!
  itop = jtop1 + jtop2
  ibot = min ( ibot1, ibot2 )
!
!  Clean up the result.
!
  call dec_round ( itop, ibot )

  return
end
subroutine dec_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
!
!*******************************************************************************
!
!! DEC_DIV divides two decimal values.
!
!
!  Discussion:
!
!    A decimal quantity is stored as
!
!      (ITOP,IBOT)
!
!    representing the value
!
!      ITOP * 10 ** IBOT.
!
!    The routine computes
!
!      ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) / (ITOP2 * 10**IBOT2)
!
!                      = (ITOP1/ITOP2) * 10**(IBOT1-IBOT2)
!
!    while avoiding integer overflow.
!
!  Modified:
!
!    23 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the numerator.
!
!    Input, integer ITOP2, IBOT2, the denominator.
!
!    Output, integer ITOP, IBOT, the result.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.
!
  implicit none
!
  double precision dval
  integer ibot
  integer ibot1
  integer ibot2
  integer ibot3
  integer ierror
  integer itop
  integer itop1
  integer itop2
  integer itop3
!
!  First special case, top fraction is 0.
!
  if ( itop1 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  First error, bottom of fraction is 0.
!
  if ( itop2 == 0 ) then
    ierror = 1
    itop = 0
    ibot = 0
    return
  end if
!
!  Second special case, result is 1.
!
  if ( itop1 == itop2 .and. ibot1 == ibot2 ) then
    itop = 1
    ibot = 0
    return
  end if
!
!  Third special case, result is power of 10.
!
  if ( itop1 == itop2 ) then
    itop = 1
    ibot = ibot1 - ibot2
    return
  end if
!
!  Fourth special case: ITOP1/ITOP2 is exact.
!
  if ( ( itop1 / itop2 ) * itop2 == itop1 ) then
    itop = itop1 / itop2
    ibot = ibot1 - ibot2
    return
  end if
!
!  General case.
!
  dval = dble ( itop1 ) / dble ( itop2 )

  call d_to_dec ( dval, itop3, ibot3 )

  itop = itop3
  ibot = ibot3 + ibot1 - ibot2

  return
end
subroutine dec_mul ( itop1, ibot1, itop2, ibot2, itop, ibot )
!
!*******************************************************************************
!
!! DEC_MUL multiplies two decimals.
!
!
!  Discussion:
!
!    The routine computes
!
!      ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) * (ITOP2 * 10**IBOT2)
!                      = (ITOP1*ITOP2) * 10**(IBOT1+IBOT2)
!
!    while avoiding integer overflow.
!
!  Modified:
!
!    25 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  Input, integer ITOP1, IBOT1, the first multiplier.
!
!  Input, integer ITOP2, IBOT2, the second multiplier.
!
!  Output, integer ITOP, IBOT, the product.
!
  implicit none
!
  double precision dval
  integer i_max
  integer ibot
  integer ibot1
  integer ibot2
  integer ibot3
  integer itop
  integer itop1
  integer itop2
  integer itop3
  real temp
!
  i_max = huge ( i_max )
!
!  The result is zero if either ITOP1 or ITOP2 is zero.
!
  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  The result is simple if either ITOP1 or ITOP2 is one.
!
  if ( abs ( itop1 ) == 1 .or. abs ( itop2 ) == 1 ) then
    itop = itop1 * itop2
    ibot = ibot1 + ibot2
    return
  end if

  temp = log ( real ( abs ( itop1 ) ) ) + log ( real ( abs ( itop2 ) ) )

  if ( temp < log ( real ( i_max ) ) ) then

    itop = itop1 * itop2
    ibot = ibot1 + ibot2

  else

    dval = dble ( itop1 ) * dble ( itop2 )

    call d_to_dec ( dval, itop3, ibot3 )

    itop = itop3
    ibot = ibot3 + ( ibot1 + ibot2 )

  end if

  call dec_round ( itop, ibot )

  return
end
subroutine dec_round ( itop, ibot )
!
!*******************************************************************************
!
!! DEC_ROUND rounds a decimal fraction to a given number of digits.
!
!
!  Discussion:
!
!    The routine takes an arbitrary decimal fraction represented by
!
!      ITOP * 10**IBOT
!
!    and makes sure that ITOP has no more than DEC_DIGIT digits.
!    DEC_DIGIT can be set or reported by calling I_DATA.
!
!  Modified:
!
!    25 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ITOP, IBOT, the coefficient and exponent
!    of a decimal fraction.  On return, ITOP has no more than
!    DEC_DIGIT decimal digits.
!
  implicit none
!
  integer dec_digit
  integer ibot
  integer itop
!
  if ( itop == 0 ) then
    ibot = 0
    return
  end if

  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )

  do while ( abs ( itop ) >= 10**dec_digit )
    itop = nint ( real ( itop ) / 10.0E+00 )
    ibot = ibot + 1
  end do
!
!  Absorb trailing 0's into the exponent.
!
  do while ( ( itop / 10 ) * 10 == itop )
    itop = itop / 10
    ibot = ibot + 1
  end do

  return
end
subroutine dec_to_r ( a, itop, ibot )
!
!*****************************************************************************
!
!! DEC_TO_R converts a decimal ITOP * 10**IBOT to a real value.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real A, the equivalent real value.
!
!    Input, integer ITOP, IBOT, the coefficient and exponent
!    of the decimal value.
!
  implicit none
!
  real a
  integer ibot
  integer itop
!
  a = itop * 10.0E+00**ibot

  return
end
subroutine dec_to_rat ( iatop, iabot )
!
!*******************************************************************************
!
!! DEC_TO_RAT converts a decimal to a rational representation.
!
!
!  Discussion:
!
!    On input, a value is represented as IATOP * 10**IABOT.
!
!    On output, approximately the same value is represented as IATOP / IABOT.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer IATOP, IABOT.
!    On input, these quantities represent the value IATOP * 10 ** IABOT.
!    On output, these quantities represent the value IATOP / IABOT.
!
  integer iabot
  integer iatop
  integer i_gcd
  integer itmp
!
  if ( iabot >= 0 ) then
    iatop = iatop * 10**iabot
    iabot = 1
  else
    iabot = 10**( - iabot )
    itmp = i_gcd ( iatop, iabot )
    iatop = iatop / itmp
    iabot = iabot / itmp
  end if

  return
end
subroutine dec_to_s_left ( ival, jval, s )
!
!*******************************************************************************
!
!! DEC_TO_S_LEFT returns a left-justified representation of IVAL * 10**JVAL.
!
!
!  Examples:
!
!    IVAL     JVAL       S
!    ----     ----       ------
!       0        0       0
!      21        3       21000
!      -3        0       -3
!     147       -2       14.7
!      16       -5       0.00016
!      34       30       Inf
!     123      -21       0.0000000000000000012
!      34      -30       0.0
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, JVAL, integers which represent the decimal.
!
!    Output, character ( len = * ) S, the representation of the value.
!    The string is 'Inf' or '0.0' if the value was too large
!    or small to represent with a fixed point format.
!
  implicit none
!
  character ( len = 22 ) chrrep
  integer i
  integer iget1
  integer iget2
  integer iput1
  integer iput2
  integer ival
  integer jval
  integer maxdigit
  integer ndigit
  integer nleft
  character ( len = * ) s
!
  s = ' '

  if ( ival == 0 ) then
    s = '0'
    return
  end if

  maxdigit = len ( s )
!
!  Store a representation of IVAL in CHRREP.
!
  write ( chrrep, '(i22)' ) ival
  call s_blank_delete ( chrrep )
  ndigit = len_trim ( chrrep )
!
!  Overflow if JVAL is positive, and NDIGIT + JVAL > MAXDIGIT.
!
  if ( jval > 0 ) then
    if ( ndigit + jval > maxdigit ) then
      s = 'Inf'
      return
    end if
  end if
!
!  Underflow if JVAL is negative, and 3 + NDIGIT - JVAL > MAXDIGIT.
!
  if ( jval < 0 ) then
    if ( ival > 0 ) then
      if ( 3 - ndigit - jval > maxdigit ) then
        s = '0.0'
        return
      end if
    else
      if ( 5 - ndigit - jval > maxdigit ) then
        s = '0.0'
        return
      end if
    end if
  end if
!
!  If JVAL is nonnegative, insert trailing zeros.
!
  if ( jval >= 0 ) then

    s(1:ndigit) = chrrep(1:ndigit)

    do i = ndigit+1, ndigit+jval
      s(i:i) = '0'
    end do

  else if ( jval < 0 ) then

    iput2 = 0
    iget2 = 0
!
!  Sign.
!
    if ( ival < 0 ) then
      iput1 = 1
      iput2 = 1
      iget2 = 1
      s(iput1:iput2) = '-'
      ndigit = ndigit - 1
    end if
!
!  Digits of the integral part.
!
    if ( ndigit + jval > 0 ) then
      iput1 = iput2 + 1
      iput2 = iput1 + ndigit + jval -1
      iget1 = iget2 + 1
      iget2 = iget1 + ndigit+jval - 1
      s(iput1:iput2) = chrrep(iget1:iget2)
    else
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end if
!
!  Decimal point.
!
    iput1 = iput2 + 1
    iput2 = iput1
    s(iput1:iput2) = '.'
!
!  Leading zeroes.
!
    do i = 1, - jval - ndigit
      iput1 = iput2 + 1
      iput2 = iput1
      s(iput1:iput2) = '0'
    end do

    nleft = min ( -jval, ndigit )
    nleft = min ( nleft, maxdigit - iput2 )
    iput1 = iput2 + 1
    iput2 = iput1 + nleft - 1
    iget1 = iget2 + 1
    iget2 = iget1 + nleft - 1
    s(iput1:iput2) = chrrep(iget1:iget2)

  end if

  return
end
subroutine decmat_det ( n, iatop, iabot, idtop, idbot, ierror )
!
!*******************************************************************************
!
!! DECMAT_DET finds the determinant of an N by N matrix of decimal entries.
!
!
!  Warning:
!
!    The brute force method is used.  The routine should only be used for
!    small matrices, since this calculation requires the summation of N!
!    products of N numbers.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of A.
!
!    Input, integer IATOP(N,N), IABOT(N,N), the decimal
!    representation of the matrix.
!
!    Output, integer IDTOP, IDBOT, the decimal determinant of the matrix.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred (probably overflow).
!
  implicit none
!
  integer n
!
  logical even
  integer i
  integer iabot(n,n)
  integer iatop(n,n)
  integer iarray(n)
  integer ibot
  integer ibot1
  integer ibot2
  integer idbot
  integer idtop
  integer ierror
  integer itop
  integer itop1
  integer itop2
  logical more
!
  ierror = 0
  more = .false.
  idtop = 0
  idbot = 1
!
!  Compute the next permutation.
!
  do

    call perm_next ( n, iarray, more, even )
!
!  The sign of this term depends on the sign of the permutation.
!
    if ( even ) then
      itop = 1
    else
      itop = - 1
    end if
!
!  Choose one item from each row, as specified by the permutation,
!  and multiply them
!
    ibot = 0

    do i = 1, n

      itop1 = itop
      ibot1 = ibot
      itop2 = iatop(i,iarray(i))
      ibot2 = iabot(i,iarray(i))

      call dec_mul ( itop1, ibot1, itop2, ibot2, itop, ibot )

    end do
!
!  Add this term to the total.
!
    itop1 = itop
    ibot1 = ibot

    itop2 = idtop
    ibot2 = idbot

    call dec_add ( itop1, ibot1, itop2, ibot2, itop, ibot )

    idtop = itop
    idbot = ibot

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine decmat_print ( m, n, a, b, title )
!
!*******************************************************************************
!
!! DECMAT_PRINT prints out decimal vectors and matrices.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!
!    Input, integer A(M,N), B(M,N), the decimal matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer b(m,n)
  character ( len = 22 ) chrtmp
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format2
  integer i
  integer ichi
  integer iclo
  integer imax
  integer imin
  integer izhi
  integer izlo
  integer j
  integer jmax
  integer jmin
  integer khi
  integer klo
  integer kmax
  integer lenc
  integer, parameter :: ncolum = 80
  integer npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how wide we must make each column.
!
  imax = 0
  jmax = 0

  do i = 1, m
    do j = 1, n

      call dec_to_s_left ( a(i,j), b(i,j), chrtmp )
      lenc = len_trim ( chrtmp )
      jmax = max ( jmax, lenc )

    end do
  end do

  kmax = 2 + imax + 1 + jmax
  npline = ncolum / kmax
!
!  Set up the format for the heading.
!
    call i_to_s_left ( npline, chrtmp2 )
    call i_to_s_left ( kmax, chrtmp3 )
    format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format2 )

  do jmin = 1, n, npline

    jmax = min ( jmin+npline-1, n )


    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '

    if ( jmin > 1 .or. jmax < n ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if

    do i = 1, m

      output = ' '

      do j = jmin, jmax
        klo = 4 + (j-jmin) * kmax + 1
        khi = 4 + (j-jmin) * kmax + kmax
        call dec_to_s_left ( a(i,j), b(i,j), chrtmp )
        output(klo:khi) = adjustr ( chrtmp(1:kmax) )
      end do

      write ( *, '(a)' ) trim ( output )

    end do

  end do

  return
end
subroutine derange_back_candidate ( n, iarray, k, nstack, stack, ncan )
!
!*******************************************************************************
!
!! DERANGE_BACK_CANDIDATE finds possible values for the K-th entry of a derangement.
!
!
!  Modified:
!
!    10 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the derangement.
!
!    Input, integer IARRAY(N).  The first K-1 entries of IARRAY
!    record the currently set values of the derangement.
!
!    Input, integer K, the entry of the derangement for which candidates
!    are to be found.
!
!    Input/output, integer NSTACK, the length of the stack.
!
!    Input/output, integer STACK((N*(N+1))/2).  On output, we have added
!    the candidates for entry K to the end of the stack.
!
!    Input/output, integer NCAN(N), the number of candidates for each level.
!
  implicit none
!
  integer n
!
  integer iarray(n)
  integer ican
  integer ifree(n)
  integer k
  integer ncan(n)
  integer nfree
  integer nstack
  integer stack((n*(n+1)/2))
!
!  Consider all the integers from 1 through N that have not been used yet.
!
  nfree = n - k + 1

  call perm_free ( iarray, k-1, ifree, nfree )
!
!  Everything but K is a legitimate candidate for the K-th entry.
!
  ncan(k) = 0

  do ican = 1, nfree

    if ( ifree(ican) /= k ) then
      ncan(k) = ncan(k) + 1
      nstack = nstack + 1
      stack(nstack) = ifree(ican)
    end if

  end do

  return
end
subroutine derange_back_next ( n, iarray, more )
!
!*******************************************************************************
!
!! DERANGE_BACK_NEXT returns the next derangement of N items.
!
!
!  Discussion:
!
!    This routine uses backtracking.
!
!  Modified:
!
!    26 November 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items to be deranged.
!
!    Input/output, integer IARRAY(N).
!    On the first call, the input value of IARRAY is not important.
!    On return with MORE = .TRUE., IARRAY contains the next derangement.
!    On subsequent input, IARRAY should not be changed.
!
!    Input/output, logical MORE.
!    On first call, set MORE to .FALSE, and do not alter it after.
!    On return, MORE is TRUE if another derangement was found, and FALSE
!    if there are no more.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer, save :: indx = 0
  integer, save :: k = 0
  integer, save :: maxstack = 0
  logical more
  integer, save, allocatable, dimension (:) :: ncan
  integer, save :: nstack = 0
  integer, save, allocatable, dimension (:) :: stack
!
  if ( .not. more ) then

    indx = 0
    k = 0
    maxstack = ( n * ( n + 1 ) ) / 2
    nstack = 0

    if ( allocated ( stack ) ) then
      deallocate ( stack )
    end if

    allocate ( stack(1:(n*(n+1))/2) )
    stack(1:maxstack) = 0

    if ( allocated ( ncan ) ) then
      deallocate ( ncan )
    end if

    allocate ( ncan(1:n) )
    ncan(1:n) = 0

    more = .true.

  end if

  do

    call ivec_backtrack ( n, iarray, indx, k, nstack, stack, maxstack, ncan )

    if ( indx == 1 ) then

      exit

    else if ( indx == 2 ) then

      call derange_back_candidate ( n, iarray, k, nstack, stack, ncan )

    else

      more = .false.
      exit

    end if

  end do

  return
end
subroutine derange_check ( n, iarray, deranged )
!
!*******************************************************************************
!
!! DERANGE_CHECK determines whether a permutation is a derangement.
!
!
!  Definition:
!
!    A derangement of the integers 1 through N is a permutation of the
!    integers such that the first value is not 1, the second is not 2,
!    and so on.
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer IARRAY(N), a permutation of the integers 1 through N.
!
!    Output, logical DERANGED, is TRUE if IARRAY is a derangement, and
!    FALSE otherwise.
!
  implicit none
!
  integer n
!
  logical deranged
  integer i
  integer iarray(n)
!
  do i = 1, n
    if ( iarray(i) == i ) then
      deranged = .false.
      return
    end if
  end do

  deranged = .true.

  return
end
function derange_enum ( n )
!
!*******************************************************************************
!
!! DERANGE_ENUM returns the number of derangements of N objects.
!
!
!  Definition:
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!  Comments:
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer DERANGE_ENUM, the number of derangements of N objects.
!
  implicit none
!
  integer derange_enum
  integer dn
  integer dnm1
  integer dnm2
  integer i
  integer n
!
  if ( n < 0 ) then

    dn = 0

  else if ( n == 0 ) then

    dn = 1

  else if ( n == 1 ) then

    dn = 0

  else if ( n == 2 ) then

    dn = 1

  else

    dnm1 = 0
    dn = 1

    do i = 3, n
      dnm2 = dnm1
      dnm1 = dn
      dn = ( i - 1 ) * ( dnm1 + dnm2 )
    end do

  end if

  derange_enum = dn

  return
end
subroutine derange_enum2 ( n, d )
!
!*******************************************************************************
!
!! DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
!
!
!  Definition:
!
!    A derangement of N objects is a permutation with no fixed
!    points.  If we symbolize the permutation operation by "P",
!    then for a derangment, P(I) is never equal to I.
!
!  Recursion:
!
!      D(0) = 1
!      D(1) = 0
!      D(2) = 1
!      D(N) = (N-1) * ( D(N-1) + D(N-2) )
!
!    or
!
!      D(0) = 1
!      D(1) = 0
!      D(N) = N * D(N-1) + (-1)**N
!
!  Formula:
!
!    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
!
!    Based on the inclusion/exclusion law.
!
!  Comments:
!
!    D(N) is the number of ways of placing N non-attacking rooks on
!    an N by N chessboard with one diagonal deleted.
!
!    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
!
!    The number of permutations with exactly K items in the right
!    place is COMB(N,K) * D(N-K).
!
!  First values:
!
!     N         D(N)
!     0           1
!     1           0
!     2           1
!     3           2
!     4           9
!     5          44
!     6         265
!     7        1854
!     8       14833
!     9      133496
!    10     1334961
!
!  Modified:
!
!    28 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the maximum number of objects to be permuted.
!
!    Output, integer D(0:N); D(I) is the number of derangements of
!    I objects.
!
  implicit none
!
  integer n
!
  integer d(0:n)
  integer i
!
  d(0) = 1
  d(1) = 0

  do i = 2, n
    d(i) = ( i - 1 ) * ( d(i-1) + d(i-2) )
  end do

  return
end
subroutine derange_weed_next ( n, iarray, more )
!
!*******************************************************************************
!
!! DERANGE_WEED_NEXT computes all of the derangements of N objects, one at a time.
!
!
!  Definition:
!
!    A derangement of the integers 1 through N is a permutation of the
!    integers such that the first value is not 1, the second is not 2,
!    and so on.
!
!  Discussion:
!
!    This routine simply generates all permutations, one at a time,
!    and weeds out those that are not derangements.
!
!  Examples:
!
!    Here are the derangements when N = 4:
!
!    2143  3142  4123
!    2341  3412  4312
!    2413  3421  4321
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer IARRAY(N).
!
!    On first call, the input contents of IARRAY are unimportant.  But
!    on the second and later calls, the input value of IARRAY should be
!    the output value returned on the previous call.
!
!    On output, IARRAY contains the next derangement.
!
!    Input/output, logical MORE.
!
!    Set MORE = .FALSE. before the first call.
!
!    MORE will be reset to .TRUE. and a derangement will be returned.
!
!    Each new call produces a new derangement until MORE is returned .FALSE.
!
  implicit none
!
  integer n
!
  logical deranged
  integer derange_enum
  integer iarray(n)
  integer, save :: maxder = 0
  logical more
  integer, save :: numder = 0
!
!  Initialization on call with MORE = .FALSE.
!
  if ( .not. more ) then

    maxder = derange_enum ( n )
    numder = 0

  end if
!
!  Watch out for cases where there are no derangements.
!
  if ( maxder == 0 ) then
    more = .false.
    return
  end if
!
!  Get the next permutation.
!
  do

    call perm_lex ( n, iarray, more )
!
!  See if it is a derangment.
!
    call derange_check ( n, iarray, deranged )

    if ( deranged ) then
      exit
    end if

  end do

  numder = numder + 1

  if ( numder >= maxder ) then
    more = .false.
  end if

  return
end
subroutine digit_to_ch ( digit, c )
!
!*******************************************************************************
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
  implicit none
!
  character c
  integer digit
!
  if ( 0 <= digit .and. digit <= 9 ) then

    c = char ( digit + 48 )

  else

    c = '*'

  end if

  return
end
subroutine digraph_arc_euler ( nnode, nedge, inode, jnode, success, trail )
!
!*******************************************************************************
!
!! DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
!
!
!  Description:
!
!    An Euler circuit of a digraph is a path which starts and ends at
!    the same node and uses each directed edge exactly once.  A digraph is
!    eulerian if it has an Euler circuit.  The problem is to decide whether
!    a given digraph is eulerian and to find an Euler circuit if the
!    answer is affirmative.
!
!  Method:
!
!    A digraph has an Euler circuit if and only if the number of incoming
!    edges is equal to the number of outgoing edges at each node.
!
!    This characterization gives a straightforward procedure to decide whether
!    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian digraph
!    G of NEDGE edges can be determined by the following method:
!
!      STEP 1: Choose any node U as the starting node, and traverse any edge
!      ( U, V ) incident to node U, and than traverse any unused edge incident
!      to node U.  Repeat this process of traversing unused edges until the
!      starting node U is reached.  Let P be the resulting walk consisting of
!      all used edges.  If all edges of G are in P, than stop.
!
!      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
!      in P and Y is not in P.  Use node X as the starting node and
!      find another walk Q using all unused edges as in step 1.
!
!      STEP 3: Walk P and walk Q share a common node X, they can be merged
!      to form a walk R by starting at any node S of P and to traverse P
!      until node X is reached; than, detour and traverse all edges of Q
!      until node X is reached and continue to traverse the edges of P until
!      the starting node S is reached.  Set P = R.
!
!      STEP 4: Repeat steps 2 and 3 until all edges are used.
!
!    The running time of the algorithm is O ( NEDGE ).
!
!  Note:
!
!    The digraph is assumed to be connected.
!
!  Reference:
!
!    Hang Tong Lau,
!    Algorithms on Graphs,
!    Tab Books, 1989.
!
!  Modified:
!
!    25 July 2000
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE); the I-th edge starts at node
!    INODE(I) and ends at node JNODE(I).
!
!    Output, logical SUCCESS, is TRUE if an Euler circuit was found,
!    and FALSE otherwise.
!
!    Output, integer TRAIL(NEDGE).  TRAIL(I) is the edge number of the I-th
!    edge in the Euler circuit.
!
  implicit none
!
  integer nedge
!
  logical candid(nedge)
  integer endnod(nedge)
  integer i
  integer inode(nedge)
  integer istak
  integer j
  integer jnode(nedge)
  integer k
  integer l
  integer len
  integer lensol
  integer lenstk
  integer nnode
  integer stack(2*nedge)
  logical success
  integer trail(nedge)
!
!  Check if the digraph is eulerian.
!
  trail(1:nedge) = 0
  endnod(1:nedge) = 0

  do i = 1, nedge
    j = inode(i)
    trail(j) = trail(j) + 1
    j = jnode(i)
    endnod(j) = endnod(j) + 1
  end do

  do i = 1, nnode
    if ( trail(i) /= endnod(i) ) then
      success = .false.
      return
    end if
  end do
!
!  The digraph is eulerian; find an Euler circuit.
!
  success = .true.
  lensol = 1
  lenstk = 0
!
!  Find the next edge.
!
  do

    if ( lensol == 1 ) then

      endnod(1) = inode(1)
      stack(1) = 1
      stack(2) = 1
      lenstk = 2

    else

      l = lensol - 1

      if ( lensol /= 2 ) then
        endnod(l) = inode(trail(l)) + jnode(trail(l)) - endnod(l-1)
      end if

      k = endnod(l)

      do i = 1, nedge
        candid(i) = ( k == jnode(i) )
      end do

      do i = 1, l
        candid(trail(i)) = .false.
      end do

      len = lenstk

      do i = 1, nedge

        if ( candid(i) ) then
          len = len + 1
          stack(len) = i
        end if

      end do

      stack(len+1) = len - lenstk
      lenstk = len + 1

    end if

    do

      istak = stack(lenstk)
      lenstk = lenstk - 1

      if ( istak /= 0 ) then
        exit
      end if

      lensol = lensol - 1

      if ( lensol == 0 ) then
        call ivec_reverse ( nedge, trail )
        return
      end if

    end do

    trail(lensol) = stack(lenstk)
    stack(lenstk) = istak - 1

    if ( lensol == nedge ) then
      exit
    end if

    lensol = lensol + 1

  end do

  call ivec_reverse ( nedge, trail )

  return
end
subroutine digraph_arc_print ( nedge, inode, jnode, title )
!
!*******************************************************************************
!
!! DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
!
!
!  Modified:
!
!    04 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NEDGE, the number of edges.
!
!    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning and end
!    nodes of the edges.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none
!
  integer nedge
!
  integer i
  integer inode(nedge)
  integer jnode(nedge)
  character ( len = * ) title
!
  if ( len_trim ( title ) /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  do i = 1, nedge
    write ( *, '(i6,4x,2i6)' ) i, inode(i), jnode(i)
  end do

  return
end
subroutine diophantine ( a, b, c, ierror, x, y )
!
!*******************************************************************************
!
!! DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
!
!
!  Discussion:
!
!    Given integers A, B and C, produce X and Y so that
!
!      A * X + B * Y = C.
!
!    In general, the equation is solvable if and only if the
!    greatest common divisor of A and B also divides C.
!
!    A solution (X,Y) of the Diophantine equation also gives the solution
!    X to the congruence equation:
!
!      A * X = C mod ( B ).
!
!    Generally, if there is one nontrivial solution, there are an infinite
!    number of solutions to a Diophantine problem.
!    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
!    T is any integer, and D is the greatest common divisor of A and B.
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Modified:
!
!    02 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, C, the coefficients of the Diophantine equation.
!
!    Output, integer IERROR, error flag.
!    0, no error, X and Y were computed.
!    1, A = B = 0, C is nonzero.
!    2, A = 0, B and C nonzero, but C is not a multiple of B.
!    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
!    4, A, B, C nonzero, but GCD of A and B does not divide C.
!    5, the algorithm ran out of internal space.
!
!    Output, integer X, Y, the solution of the Diophantine equation.
!    Note that the algorithm will attempt to return a solution with
!    smallest Euclidean norm.
!
  implicit none
!
  integer, parameter :: nmax = 100
!
  integer a
  integer a_copy
  integer a_mag
  integer a_sign
  integer b
  integer b_copy
  integer b_mag
  integer b_sign
  integer c
  integer c_copy
  integer g
  integer i_gcd
  integer ierror
  integer k
  integer n
  real norm_new
  real norm_old
  integer q(nmax)
  logical swap
  integer temp
  integer x
  integer xnew
  integer y
  integer ynew
!
!  Defaults for output parameters.
!
  ierror = 0
  x = 0
  y = 0
!
!  Special cases.
!
  if ( a == 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b == 0 .and. c /= 0 ) then
    ierror = 1
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a == 0 .and. b /= 0 .and. c /= 0 ) then
    x = 0
    y = c / b
    if ( mod ( c, b ) /= 0 ) then
      ierror = 2
    end if
    return
  else if ( a /= 0 .and. b == 0 .and. c == 0 ) then
    x = 0
    y = 0
    return
  else if ( a /= 0 .and. b == 0 .and. c /= 0 ) then
    x = c / a
    y = 0
    if ( mod ( c, a ) /= 0 ) then
      ierror = 3
    end if
    return
  else if ( a /= 0 .and. b /= 0 .and. c == 0 ) then
    g = i_gcd ( a, b )
    x = b / g
    y = - a / g
    return
  end if
!
!  Now handle the "general" case: A, B and C are nonzero.
!
!  Step 1: Compute the GCD of A and B, which must also divide C.
!
  g = i_gcd ( a, b )

  if ( mod ( c, g ) /= 0 ) then
    ierror = 4
    return
  end if

  a_copy = a / g
  b_copy = b / g
  c_copy = c / g
!
!  Step 2: Split A and B into sign and magnitude.
!
  a_mag = abs ( a_copy )
  a_sign = sign ( 1, a_copy )
  b_mag = abs ( b_copy )
  b_sign = sign ( 1, b_copy )
!
!  Another special case, A_MAG = 1 or B_MAG = 1.
!
  if ( a_mag == 1 ) then
    x = a_sign * c_copy
    y = 0
    return
  else if ( b_mag == 1 ) then
    x = 0
    y = b_sign * c_copy
    return
  end if
!
!  Step 3: Produce the Euclidean remainder sequence.
!
  if ( a_mag >= b_mag ) then

    swap = .false.
    q(1) = a_mag
    q(2) = b_mag

  else

    swap = .true.
    q(1) = b_mag
    q(2) = a_mag

  end if

  n = 3

  do

    q(n) = mod ( q(n-2), q(n-1) )

    if ( q(n) == 1 ) then
      exit
    end if

    n = n + 1

    if ( n > nmax ) then
      ierror = 5
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIOPHANTINE - Fatal error!'
      write ( *, '(a)' ) '  Exceeded number of iterations.'
      stop
    end if

  end do
!
!  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
!
  y = 0
  do k = n, 2, -1
    x = y
    y = ( 1 - x * q(k-1) ) / q(k)
  end do
!
!  Step 5: Undo the swapping.
!
  if ( swap ) then
    call i_swap ( x, y )
  end if
!
!  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
!
  x = x * a_sign
  y = y * b_sign
!
!  Step 7: Multiply by C, so that X * A + Y * B = C.
!
  x = x * c_copy
  y = y * c_copy
!
!  Step 8: Given a solution (X,Y), try to find the solution of
!  minimal magnitude.
!
  call diophantine_solution_minimize ( a_copy, b_copy, x, y )

  return
end
subroutine diophantine_solution_minimize ( a, b, x, y )
!
!*******************************************************************************
!
!! DIOPHANTINE_SOLUTION_MINIMIZE seeks a minimal solution of a Diophantine equation.
!
!
!  Discussion:
!
!    Given a solution (X,Y) of a Diophantine equation:
!
!      A * X + B * Y = C.
!
!    then there are an infinite family of solutions of the form
!
!      ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
!
!    An integral solution of minimal Euclidean norm can be found by
!    tentatively moving along the vectors (B,-A) and (-B,A) one step
!    at a time.
!
!  Reference:
!
!    Eric Weisstein, editor,
!    CRC Concise Encylopedia of Mathematics,
!    CRC Press, 1998, page 446.
!
!  Modified:
!
!    03 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the coefficients of the Diophantine equation.
!    A and B are assumed to be relatively prime.
!
!    Input/output, integer X, Y, on input, a solution of the Diophantine
!    equation.  On output, a solution of minimal Euclidean norm.
!
  implicit none
!
  integer a
  integer b
  real norm
  real norm_new
  real t
  integer x
  integer xnew
  integer y
  integer ynew
!
!  Compute the minimum for T real, and then look nearby.
!
  t = ( - real ( b ) * real ( x ) + real ( a ) * real ( y ) ) &
    / ( real ( a )**2 + real ( b )**2 )

  x = x + nint ( t ) * b
  y = y - nint ( t ) * a
!
!  Now look nearby.
!
  norm = ( real ( x ) )**2 + ( real ( y ) )**2

  do

    xnew = x + b
    ynew = y - a

    norm_new = ( real ( xnew ) )**2 + ( real ( ynew ) )**2

    if ( norm_new >= norm ) then
      exit
    end if

    x = xnew
    y = ynew
    norm = norm_new

  end do

  do

    xnew = x - b
    ynew = y + a

    norm_new = ( real ( xnew ) )**2 + ( real ( ynew ) )**2

    if ( norm_new >= norm ) then
      exit
    end if

    x = xnew
    y = ynew
    norm = norm_new

  end do

  return
end
subroutine equiv_next ( n, npart, jarray, iarray, more )
!
!*******************************************************************************
!
!! EQUIV_NEXT computes the partitions of a set one at a time.
!
!
!  Definition:
!
!    A partition of a set assigns each element to exactly one subset.
!
!  Comments:
!
!    The number of partitions of a set of size N is the Bell number B(N).
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
!    18 July 2000
!
!  Parameters:
!
!    Input, integer N, number of elements in the set to be partitioned.
!
!    Output, integer NPART, number of subsets in the partition.
!
!    Output, integer JARRAY(N).  JARRAY(I) is the number of elements
!    in the I-th subset of the partition.
!
!    Output, integer IARRAY(N).  IARRAY(I) is the class to which
!    element I belongs.
!
!    Input/output, logical MORE.  Set MORE = .FALSE. before first call.
!    It is reset and held at .TRUE. as long as
!    the partition returned is not the last one.
!    When MORE is returned .FALSE., all the partitions
!    have been computed and returned.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer jarray(n)
  integer l
  integer m
  logical more
  integer npart
!
  if ( .not. more ) then

    npart = 1
    iarray(1:n) = 1
    jarray(1) = n

  else

    m = n

    do while ( jarray(iarray(m)) == 1 )
      iarray(m) = 1
      m = m - 1
    end do

    l = iarray(m)
    npart = npart + m - n
    jarray(1) = jarray(1) + n - m

    if ( l == npart ) then
      npart = npart + 1
      jarray(npart) = 0
    end if

    iarray(m) = l + 1
    jarray(l) = jarray(l) - 1
    jarray(l+1) = jarray(l+1) + 1

  end if

  more = npart /= n

  return
end
subroutine equiv_next2 ( done, iarray, n )
!
!*******************************************************************************
!
!! EQUIV_NEXT2 computes, one at a time, the partitions of a set.
!
!
!  Definition:
!
!    A partition of a set assigns each element to exactly one subset.
!
!  Comments:
!
!    The number of partitions of a set of size N is the Bell number B(N).
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input/output, logical DONE.  Before the very first call, the
!    user should set DONE to .TRUE., which prompts the program
!    to initialize its data, and return the first partition.
!
!    Thereafter, the user should call again, for the next
!    partition, and so on, until the routine returns with DONE
!    equal to .TRUE., at which point there are no more partitions
!    to compute.
!
!    Input/output, integer IARRAY(N), contains the information
!    defining the current partition.  The user should not alter
!    IARRAY between calls.  Except for the very first
!    call, the routine uses the previous output value of IARRAY to compute
!    the next value.
!
!    The entries of IARRAY are the partition subset to which each
!    element of the original set belongs.  If there are NPART distinct
!    parts of the partition, then each entry of IARRAY will be a
!    number between 1 and NPART.  Every number from 1 to NPART will
!    occur somewhere in the list.  If the entries of IARRAY are
!    examined in order, then each time a new partition subset occurs,
!    it will be the next unused integer.
!
!    For instance, for N = 4, the program will describe the set
!    where each element is in a separate subset as 1, 2, 3, 4,
!    even though such a partition might also be described as
!    4, 3, 2, 1 or even 1, 5, 8, 19.
!
!    Input, integer N, the number of elements in the set.
!
  implicit none
!
  integer n
!
  logical done
  integer i
  integer iarray(n)
  integer imax
  integer j
  integer jmax
!
  if ( done ) then

    done = .false.

    iarray(1:n) = 1

  else
!
!  Find the last element J that can be increased by 1.
!  This is the element that is not equal to its maximum possible value,
!  which is the maximum value of all preceding elements +1.
!
    jmax = iarray(1)
    imax = 1

    do j = 2, n

      if ( iarray(j) > jmax ) then
        jmax = iarray(j)
      else
        imax = j
      end if

    end do
!
!  If no element can be increased by 1, we are done.
!
    if ( imax == 1 ) then
      done = .true.
      return
    end if
!
!  Increase the value of the IMAX-th element by 1, set its successors to 1.
!
    done = .false.
    iarray(imax) = iarray(imax) + 1
    iarray(imax+1:n) = 1

  end if

  return
end
subroutine equiv_print ( n, iarray )
!
!*******************************************************************************
!
!! EQUIV_PRINT prints a partition of a set.
!
!
!  Modified:
!
!    11 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, number of elements in set to be partitioned.
!
!    Input, integer IARRAY(N), defines the partition or set of equivalence
!    classes.  Element I belongs to subset IARRAY(I).
!
  implicit none
!
  integer n
!
  integer iarray(n)
  integer karray(n)
  integer j
  integer k
  integer s
  integer s_max
  integer s_min
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Set partition:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Set  Size'

  s_min = minval ( iarray(1:n) )
  s_max = maxval ( iarray(1:n) )

  do s = s_min, s_max

    k = 0

    do j = 1, n

      if ( iarray(j) == s ) then
        k = k + 1
        karray(k) = j
      end if

    end do

    if ( k > 0 ) then
      write ( *, '(2i6,4x,(10i4))' ) s, k, karray(1:k)
    end if

  end do

  return
end
subroutine equiv_random ( n, npart, iarray, b )
!
!*******************************************************************************
!
!! EQUIV_RANDOM selects a random partition of a set.
!
!
!  Discussion:
!
!    The user does not control the number of parts in the partition.
!
!    The equivalence classes are numbered in no particular order.
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
!    12 May 2002
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set to be partitioned.
!
!    Output, integer NPART, the number of classes or parts in the 
!    partition.  NPART will be between 1 and N.
!
!    Output, integer IARRAY(N), indicates the class to which each element
!    is assigned.
!
!    Output, real B(N).  B(K) = A(K)/(K!), where
!    A(K) = number of partitions of a set of K objects.
!
  implicit none
!
  integer n
!
  real b(n)
  integer iarray(n)
  integer k
  integer l
  integer m
  integer npart
  real sum1
  real z
  real, parameter :: zhi = 1.0E+00
  real, parameter :: zlo = 0.0E+00
!
  b(1) = 1.0E+00

  do l = 1, n-1

    sum1 = 1.0E+00 / real ( l )
    do k = 1, l-1
      sum1 = ( sum1 + b(k) ) / real ( l - k )
    end do

    b(l+1) = ( sum1 + b(l) ) / real ( l + 1 )

  end do

  m = n
  npart = 0

  do

    call r_random ( zlo, zhi, z )
    z = real ( m  ) * b(m) * z
    k = 0
    npart = npart + 1

    do while ( z >= 0.0E+00 )

      iarray(m) = npart
      m = m - 1

      if ( m == 0 ) then
        exit
      end if

      z = z - b(m)
      k = k + 1
      z = z * k

    end do

    if ( m == 0 ) then
      exit
    end if

  end do
!
!  Randomly permute the assignments.
!
  call perm_random2 ( n, iarray )

  return
end
subroutine euler ( ieuler, n )
!
!*******************************************************************************
!
!! EULER returns the N-th row of Euler's triangle.
!
!
!  Definition:
!
!    E(N,K) counts the number of permutations of the N digits that have
!    exactly K "ascents", that is, K places where the Ith digit is
!    less than the (I+1)th digit.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IEULER(0:N), the N-th row of Euler's
!    triangle, IEULER(K) contains the value of E(N,K).  Note
!    that IEULER(0) should be 1 and IEULER(N) should be 0.
!
!    Input, integer N, the row of Euler's triangle desired.
!
  implicit none
!
  integer n
!
  integer ieuler(0:n)
  integer irow
  integer k
!
  ieuler(0) = 1

  if ( n > 0 ) then

    ieuler(1) = 0

    do irow = 2, n

      ieuler(irow) = 0

      do k = irow-1, 1, -1

        ieuler(k) = ( k + 1 ) * ieuler(k) + ( irow - k ) * ieuler(k-1)

      end do

      ieuler(0) = 1

    end do

  end if

  return
end
function fall ( x, n )
!
!*******************************************************************************
!
!! FALL computes the falling factorial function [X]_N.
!
!
!  Discussion:
!
!    Note that the number of "injections" or 1-to-1 mappings from
!    a set of N elements to a set of M elements is [M]_N.
!
!    The number of permutations of N objects out of M is [M}_N.
!
!    Moreover, the Stirling numbers of the first kind can be used
!    to convert a falling factorial into a polynomial, as follows:
!
!      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
!
!  Formula:
!
!    [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the falling factorial function.
!
!    Input, integer N, the order of the falling factorial function.
!    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
!    negative, a "rising" factorial will be computed.
!
!    Output, real FALL, the value of the falling factorial function.
!
  implicit none
!
  real arg
  real fall
  integer i
  integer n
  real x
!
  fall = 1.0E+00

  arg = x

  if ( n > 0 ) then

    do i = 1, n
      fall = fall * arg
      arg = arg - 1.0E+00
    end do

  else if ( n < 0 ) then

    do i = -1, n
      fall = fall * arg
      arg = arg + 1.0E+00
    end do

  end if

  return
end
function gamma_log ( x )
!
!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
!    approximation for X > 12 is from reference 3, while approximations
!    for X < 12.0 are similar to those in reference 1, but are unpublished.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  References:
!
!    # 1)
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Math. Comp.
!    Volume 21, 1967, pages 198-203.
!
!    # 2)
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    # 3)
!    Hart, Et. Al.,
!    Computer Approximations,
!    Wiley and sons, New York, 1968.
!
!  Parameters:
!
!    Input, real X, the argument of the Gamma function.  X must be positive.
!
!    Output, real GAMMA_LOG, the logarithm of the Gamma function of X.
!    If X <= 0.0, or if overflow would occur, the program returns the
!    value XINF, the largest representable floating point number.
!
!*******************************************************************************
!
!  Explanation of machine-dependent constants
!
!  BETA   - radix for the floating-point representation.
!
!  MAXEXP - the smallest positive power of BETA that overflows.
!
!  XBIG   - largest argument for which LN(GAMMA(X)) is representable
!           in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!  FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!  Approximate values for some important machines are:
!
!                            BETA      MAXEXP         XBIG
!
!  CRAY-1        (S.P.)        2        8191       9.62E+2461
!  Cyber 180/855
!    under NOS   (S.P.)        2        1070       1.72E+319
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)        2         128       4.08E+36
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!  IBM 3033      (D.P.)       16          63       4.29D+73
!  VAX D-Format  (D.P.)        2         127       2.05D+36
!  VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           FRTBIG
!
!  CRAY-1        (S.P.)   3.13E+615
!  Cyber 180/855
!    under NOS   (S.P.)   6.44E+79
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   1.42E+9
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)   2.25D+76
!  IBM 3033      (D.P.)   2.56D+18
!  VAX D-Format  (D.P.)   1.20D+9
!  VAX G-Format  (D.P.)   1.89D+76
!
  implicit none
!
  real, parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728E-03, &
     8.4171387781295E-04, &
    -5.952379913043012E-04, &
     7.93650793500350248E-04, &
    -2.777777777777681622553E-03, &
     8.333333333333333331554247E-02, &
     5.7083835261E-03 /)
  real corr
  real, parameter :: d1 = - 5.772156649015328605195174E-01
  real, parameter :: d2 =   4.227843350984671393993777E-01
  real, parameter :: d4 =   1.791759469228055000094023E+00
  real, parameter :: frtbig = 1.42E+09
  integer i
  real gamma_log
  real, parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888E+00, &
    2.018112620856775083915565E+02, &
    2.290838373831346393026739E+03, &
    1.131967205903380828685045E+04, &
    2.855724635671635335736389E+04, &
    3.848496228443793359990269E+04, &
    2.637748787624195437963534E+04, &
    7.225813979700288197698961E+03 /)
  real, parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064E+00, &
    5.424138599891070494101986E+02, &
    1.550693864978364947665077E+04, &
    1.847932904445632425417223E+05, &
    1.088204769468828767498470E+06, &
    3.338152967987029735917223E+06, &
    5.106661678927352456275255E+06, &
    3.074109054850539556250927E+06 /)
  real, parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062E+04, &
    2.426813369486704502836312E+06, &
    1.214755574045093227939592E+08, &
    2.663432449630976949898078E+09, &
    2.940378956634553899906876E+010, &
    1.702665737765398868392998E+011, &
    4.926125793377430887588120E+011, &
    5.606251856223951465078242E+011 /)
  real, parameter :: pnt68 = 0.6796875E+00
  real, parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036E+01, &
    1.113332393857199323513008E+03, &
    7.738757056935398733233834E+03, &
    2.763987074403340708898585E+04, &
    5.499310206226157329794414E+04, &
    6.161122180066002127833352E+04, &
    3.635127591501940507276287E+04, &
    8.785536302431013170870835E+03 /)
  real, parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942E+02, &
    7.765049321445005871323047E+03, &
    1.331903827966074194402448E+05, &
    1.136705821321969608938755E+06, &
    5.267964117437946917577538E+06, &
    1.346701454311101692290052E+07, &
    1.782736530353274213975932E+07, &
    9.533095591844353613395747E+06 /)
  real, parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843E+03, &
    6.393885654300092398984238E+05, &
    4.135599930241388052042842E+07, &
    1.120872109616147941376570E+09, &
    1.488613728678813811542398E+010, &
    1.016803586272438228077304E+011, &
    3.417476345507377132798597E+011, &
    4.463158187419713286462081E+011 /)
  real res
  real, parameter :: sqrtpi = 0.9189385332046727417803297E+00
  real x
  real, parameter :: xbig = 4.08E+36
  real xden
  real xm1
  real xm2
  real xm4
  real xnum
  real xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0E+00 .or. x > xbig ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  if ( x <= epsilon ( x ) ) then

    res = - log ( x )

  else if ( x <= 1.5E+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0E+00
      xm1 = ( x - 0.5E+00 ) - 0.5E+00
    end if

    if ( x <= 0.5E+00 .or. x >= pnt68 ) then

      xden = 1.0E+00
      xnum = 0.0E+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5E+00 ) - 0.5E+00
      xden = 1.0E+00
      xnum = 0.0E+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0E+00 ) then

    xm2 = x - 2.0E+00
    xden = 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0E+00 ) then

    xm4 = x - 4.0E+00
    xden = - 1.0E+00
    xnum = 0.0E+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0E+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5E+00 * corr
    res = res + x * ( corr - 1.0E+00 )

  end if

  gamma_log = res

  return
end
subroutine gray_rank ( gray, rank )
!
!*******************************************************************************
!
!! GRAY_RANK ranks a Gray code.
!
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer GRAY, the Gray code to be ranked.
!
!    Output, integer RANK, the rank of GRAY, and the integer whose Gray
!    code is GRAY.
!
  implicit none
!
  integer gray
  integer i
  integer ii
  integer j
  integer, parameter :: nbits = 32
  integer rank
!
  j = 0

  if ( btest ( gray, nbits-1 ) ) then
    j = ibset ( j, nbits-1 )
  end if

  do i = 2, nbits

    ii = nbits - i

    if ( btest ( j, ii+1 ) .and. ( .not. btest ( gray, ii ) ) ) then
      j = ibset ( j, ii )
    end if

    if ( ( .not. btest ( j, ii+1 ) ) .and. btest ( gray, ii ) ) then
      j = ibset ( j, ii )
    end if

  end do

  rank = j

  return
end
subroutine gray_unrank ( rank, gray )
!
!*******************************************************************************
!
!! GRAY_UNRANK unranks a Gray code.
!
!
!  Comments:
!
!    The binary values of the Gray codes of successive integers differ in
!    just one bit.
!
!    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
!    Hamiltonian cycle on a graph of the cube in N dimensions.
!
!  Examples:
!
!    RANK IGRAY BINARY
!     0     0       0
!     1     1       1
!     2     3      11
!     3     2      10
!     4     6     110
!     5     7     111
!     6     5     101
!     7     4     100
!     8    12    1100
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer RANK, the integer whose Gray code is desired.
!
!    Output, integer GRAY, the Gray code of the given rank.
!
  implicit none
!
  integer gray
  integer i
  integer ii
  integer j
  integer, parameter :: nbits = 32
  integer rank
!
  j = 0

  if ( btest ( rank, nbits-1 ) ) then
    j = ibset ( j, nbits-1 )
  end if

  do i = 2, nbits

    ii = nbits - i

    if ( btest ( rank, ii+1 ) .neqv. btest ( rank, ii ) ) then
      j = ibset ( j, ii )
    end if

  end do

  gray = j

  return
end
subroutine i_data ( op, var, ival )
!
!*******************************************************************************
!
!! I_DATA stores and retrieves common data items.
!
!
!  Discussion:
!
!    This routine works like a sort of COMMON block.  It stores or returns
!    the values of certain variables.  Thus, it allows routines
!    to "communicate" without having to have data passed up and
!    down the calling tree in argument lists.
!
!    The variables stored by this version of the routine are:
!
!    'DEC_DIGIT', the number of digits stored for decimals;
!    'DEC_DIGIT_MAX', the maximum number of digits stored for decimals;
!    'I_BIG', the biggest integer to use in calculations.
!
!  Modified:
!
!    20 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) OP, describes the operation to be done.
!    'SET' means set a value.
!    'INC' means increment a value (and return its new value)
!    'GET' means get a value.
!
!    Input, character ( len = * ) VAR, the name of the variable.
!
!    Input/output, integer IVAL.
!    If OP is 'SET', then the variable named in VAR is set to the
!    value IVAL.
!    If OP is 'GET', then the value of IVAL is set to the value of
!    the variable named in VAR.
!    If OP is 'INC', then the value of IVAL is incremented by 1,
!    and its new value is returned in VAR.
!
  implicit none
!
  integer, save :: dec_digit = 7
  integer, save :: dec_digit_max = 7
  integer, save :: i_big = huge ( i_big )
  integer ival
  character ( len = * ) op
  logical s_eqi
  character ( len = * ) var
!
  if ( s_eqi ( op, 'SET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = ival
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = ival
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      i_big = ival
    end if

  else if ( s_eqi ( op, 'GET' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      ival = i_big
    end if

  else if ( s_eqi ( op, 'INC' ) ) then

    if ( s_eqi ( var, 'DEC_DIGIT' ) ) then
      dec_digit = dec_digit + 1
      ival = dec_digit
    else if ( s_eqi ( var, 'DEC_DIGIT_MAX' ) ) then
      dec_digit_max = dec_digit_max + 1
      ival = dec_digit_max
    else if ( s_eqi ( var, 'I_BIG' ) ) then
      i_big = i_big + 1
      ival = i_big
    end if

  end if

  return
end
subroutine i_factor ( n, maxfactor, nfactor, factor, power, nleft )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer into prime factors.
!
!
!  Formula:
!
!    N = NLEFT * Product ( I = 1 to NFACTOR ) FACTOR(I)**POWER(I).
!
!  Modified:
!
!    20 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be factored.  N may be positive,
!    negative, or 0.
!
!    Input, integer MAXFACTOR, the maximum number of prime factors for
!    which storage has been allocated.
!
!    Output, integer NFACTOR, the number of prime factors of N discovered
!    by the routine.
!
!    Output, integer FACTOR(MAXFACTOR), the prime factors of N.
!
!    Output, integer POWER(MAXFACTOR).  POWER(I) is the power of
!    the FACTOR(I) in the representation of N.
!
!    Output, integer NLEFT, the factor of N that the routine could not
!    divide out.  If NLEFT is 1, then N has been completely factored.
!    Otherwise, NLEFT represents factors of N involving large primes.
!
  implicit none
!
  integer maxfactor
!
  integer factor(maxfactor)
  integer i
  integer maxprime
  integer n
  integer nleft
  integer nfactor
  integer p
  integer power(maxfactor)
  integer prime
!
  nfactor = 0

  factor(1:maxfactor) = 0
  power(1:maxfactor) = 0

  nleft = n

  if ( n == 0 ) then
    return
  end if

  if ( abs ( n ) == 1 ) then
    nfactor = 1
    factor(1) = 1
    power(1) = 1
    return
  end if
!
!  Find out how many primes we stored.
!
  maxprime = prime ( -1 )
!
!  Try dividing the remainder by each prime.
!
  do i = 1, maxprime

    p = prime ( i )

    if ( mod ( abs ( nleft ), p ) == 0 ) then

      if ( nfactor < maxfactor ) then

        nfactor = nfactor + 1
        factor(nfactor) = p

        do

          power(nfactor) = power(nfactor) + 1
          nleft = nleft / p

          if ( mod ( abs ( nleft ), p ) /= 0 ) then
            exit
          end if

        end do

        if ( abs ( nleft ) == 1 ) then
          exit
        end if

      end if

    end if

  end do

  return
end
function i_gcd ( i, j )
!
!*******************************************************************************
!
!! I_GCD finds the greatest common divisor of I and J.
!
!
!  Modified:
!
!    03 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, two numbers whose greatest common divisor
!    is desired.
!
!    Output, integer I_GCD, the greatest common divisor of I and J.
!
!    Note that only the absolute values of I and J are
!    considered, so that the result is always nonnegative.
!
!    If I or J is 0, I_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
!
!    If I and J have no common factor, I_GCD is returned as 1.
!
!    Otherwise, using the Euclidean algorithm, I_GCD is the
!    largest common factor of I and J.
!
  implicit none
!
  integer i
  integer i_gcd
  integer ip
  integer iq
  integer ir
  integer j
!
  i_gcd = 1
!
!  Return immediately if either I or J is zero.
!
  if ( i == 0 ) then
    i_gcd = max ( 1, abs ( j ) )
    return
  else if ( j == 0 ) then
    i_gcd = max ( 1, abs ( i ) )
    return
  end if
!
!  Set IP to the larger of I and J, IQ to the smaller.
!  This way, we can alter IP and IQ as we go.
!
  ip = max ( abs ( i ), abs ( j ) )
  iq = min ( abs ( i ), abs ( j ) )
!
!  Carry out the Euclidean algorithm.
!
  do

    ir = mod ( ip, iq )

    if ( ir == 0 ) then
      exit
    end if

    ip = iq
    iq = ir

  end do

  i_gcd = iq

  return
end
function i_modp ( i, j )
!
!*******************************************************************************
!
!! I_MODP returns the nonnegative remainder of integer division.
!
!
!  Formula:
!
!    If
!      NREM = I_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I_MODP(A,360) is between 0 and 360, always.
!
!  Examples:
!
!        I         J     MOD  I_MODP   I_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I_MODP, the nonnegative remainder when I is divided by J.
!
  implicit none
!
  integer i
  integer i_modp
  integer j
!
  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_MODP - Fatal error!'
    write ( *, '(a,i6)' ) '  I_MODP ( I, J ) called with J = ', j
    stop
  end if

  i_modp = mod ( i, j )

  if ( i_modp < 0 ) then
    i_modp = i_modp + abs ( j )
  end if

  return
end
subroutine i_moebius ( n, mu )
!
!*******************************************************************************
!
!! I_MOEBIUS returns the value of MU(N), the Moebius function of N.
!
!
!  Definition:
!
!    MU(N) is defined as follows:
!
!      MU(N) = 1 if N = 1;
!              0 if N is divisible by the square of a prime;
!              (-1)**K, if N is the product of K distinct primes.
!
!  First values:
!
!     N  MU(N)
!
!     1    1
!     2   -1
!     3   -1
!     4    0
!     5   -1
!     6    1
!     7   -1
!     8    0
!     9    0
!    10    1
!    11   -1
!    12    0
!    13   -1
!    14    1
!    15    1
!    16    0
!    17   -1
!    18    0
!    19   -1
!    20    0
!
!  Note:
!
!    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
!    if N is a square, cube, etc.
!
!  Formula:
!
!    The Moebius function is related to Euler's totient function:
!
!  PHI(N) = sum ( D divides N ) MU(D) * ( N / D ).
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the value to be analyzed.
!
!    Output, integer MU, the value of MU(N).
!    If N is less than or equal to 0, MU will be returned as -2.
!    If there was not enough internal space for factoring, MU
!    is returned as -3.
!
  implicit none
!
  integer, parameter :: maxfactor = 20
!
  integer factor(maxfactor)
  integer i
  integer mu
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
!
  if ( n <= 0 ) then
    mu = - 2
    return
  end if

  if ( n == 1 ) then
    mu = 1
    return
  end if
!
!  Factor N.
!
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )

  if ( nleft /= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_MOEBIUS - Fatal error!'
    write ( *, '(a)' ) '  Not enough factorization space.'
    mu = - 3
    return
  end if

  mu = 1

  do i = 1, nfactor

    mu = - mu

    if ( power(i) > 1 ) then
      mu = 0
      return
    end if

  end do

  return
end
subroutine i_random ( ilo, ihi, i )
!
!*******************************************************************************
!
!! I_RANDOM returns a random integer in a given range.
!
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ILO, IHI, the minimum and maximum acceptable values.
!
!    Output, integer I, the randomly chosen integer.
!
  implicit none
!
  integer i
  integer ihi
  integer ilo
  real r
  real rhi
  real rlo
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = r )
!
!  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
!  each with a "neighborhood" of width 1.
!
  rlo = real ( ilo ) - 0.5E+00
  rhi = real ( ihi ) + 0.5E+00
!
!  Set I to the integer that is nearest the scaled value of R.
!
  i = nint ( ( 1.0E+00 - r ) * rlo + r * rhi )
!
!  In case of oddball events at the boundary, enforce the limits.
!
  i = max ( i, ilo )
  i = min ( i, ihi )

  return
end
subroutine i_sqrt ( n, q, r )
!
!*******************************************************************************
!
!! I_SQRT finds the integer square root of N by solving N = Q**2 + R.
!
!
!  Discussion:
!
!    The integer square root of N is an integer Q such that
!    Q**2 <= N but (Q+1)**2 > N.
!
!    A simpler calculation would be something like
!
!      Q = INT ( SQRT ( REAL ( N ) ) )
!
!    but this calculation has the virtue of using only integer arithmetic.
!
!    To avoid the tedium of worrying about negative arguments, the routine
!    automatically considers the absolute value of the argument.
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999, pages 294-307.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number whose integer square root is desired.
!    Actually, only the absolute value of N is considered.
!
!    Output, integer Q, R, the integer square root, and positive remainder,
!    of N.
!
  implicit none
!
  integer n
  integer n_abs
  integer q
  integer r
!
  n_abs = abs ( n )

  q = n_abs

  if ( n_abs > 0 ) then

    do while ( q > n_abs / q )
      q = ( q + ( n_abs / q ) ) / 2
    end do

  end if

  r = n_abs - q**2

  return
end
subroutine i_sqrt_cf ( n, max_term, n_term, b )
!
!*******************************************************************************
!
!! I_SQRT_CF finds the continued fraction representation of a square root of an integer.
!
!
!  Discussion:
!
!    The continued fraction representation of the square root of an integer
!    has the form
!
!      [ B0, (B1, B2, B3, ..., BM), ... ]
!
!    where
!
!      B0 = int ( sqrt ( real ( N ) ) )
!      BM = 2 * B0
!      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
!      the value M is termed the period of the representation.
!
!  Examples:
!
!     N  Period  Continued Fraction
!
!     2       1  [ 1, 2, 2, 2, ... ]
!     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
!     4       0  [ 2 ]
!     5       1  [ 2, 4, 4, 4, ... ]
!     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
!     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
!     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
!     9       0  [ 3 ]
!    10       1  [ 3, 6, 6, 6, ... ]
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999, pages 294-307.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number whose continued fraction square root
!    is desired.
!
!    Input, integer MAX_TERM, the maximum number of terms that may
!    be computed.
!
!    Output, integer N_TERM, the number of terms computed beyond the
!    0 term.  The routine should stop if it detects that the period
!    has been reached.
!
!    Output, integer B(0:MAX_TERM), contains the continued fraction
!    coefficients in entries B(0:N_TERM).
!
  implicit none
!
  integer max_term
!
  integer b(0:max_term)
  integer k
  integer n
  integer n_term
  integer p
  integer q
  integer r
  integer s
!
  n_term = 0

  call i_sqrt ( n, s, r )
  b(0) = s

  if ( r > 0 ) then

    p = 0
    q = 1

    do

      p = b(n_term) * q - p
      q = ( n - p**2 ) / q

      if ( n_term >= max_term ) then
        return
      end if

      n_term = n_term + 1
      b(n_term) = ( p + s ) / q

      if ( q == 1 ) then
        exit
      end if

    end do

  end if

  return
end
subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP switches two integer values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none
!
  integer i
  integer j
  integer k
!
  k = i
  i = j
  j = k

  return
end
subroutine i_to_chinese ( j, n, m, r )
!
!*******************************************************************************
!
!! I_TO_CHINESE converts an integer to its Chinese remainder form.
!
!
!  Discussion:
!
!  Modified:
!
!    05 July 2000
!
!  Parameters:
!
!    Input, integer J, the integer to be converted.
!
!    Input, integer N, the number of moduluses.
!
!    Input, integer M(N), the moduluses.  These should be positive
!    and pairwise prime.
!
!    Output, integer R(N), the Chinese remainder representation of the integer.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer i_modp
  integer j
  integer m(n)
  integer r(n)
!
  call chinese_check ( n, m, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I_TO_CHINESE - Fatal error!'
    write ( *, '(a)' ) '  The moduluses are not legal.'
    stop
  end if

  do i = 1, n
    r(i) = i_modp ( j, m(i) )
  end do

  return
end
subroutine i_to_halton ( n, p, r )
!
!*******************************************************************************
!
!! I_TO_HALTON computes an element of a Halton sequence.
!
!
!  Discussion:
!
!    The Halton sequence is often used to generate a "subrandom"
!    sequence of points which have a better covering property
!    than pseudorandom points.
!
!    The Halton sequence generates a sequence of points in [0,1]
!    which (theoretically) never repeats.
!
!    To generate a sequence of points in a 2 dimensional space,
!    it is typical practice to generate a pair of Halton sequences,
!    the first with prime base 2, the second with prime base 3.
!    Similarly, by using the first K primes, a suitable sequence
!    in K-dimensional space can be generated.
!
!    The generation is quite simple.  Given an integer N, the expansion
!    of N in base P is generated.  Then, essentially, the result R
!    is generated by writing a decimal point followed by the digits of
!    the expansion of N, in reverse order.  This decimal value is actually
!    still in base P, so it must be properly interpreted to generate
!    a usable value.
!
!  Example:
!
!    P = 2
!
!    N        N         Halton   Halton
!    decimal  binary    binary   decimal
!    -------  ------    ------   -------
!        1  =     1  =>  .1     = 0.5
!        2  =    10  =>  .01    = 0.25
!        3  =    11  =>  .11    = 0.75
!        4  =   100  =>  .001   = 0.125
!        5  =   101  =>  .101   = 0.625
!        6  =   110  =>  .011   = 0.375
!        7  =   111  =>  .111   = 0.875
!        8  =  1000  =>  .0001  = 0.0625
!
!  Reference:
!
!    J H Halton,
!    Numerische Mathematik,
!    Volume 2, pages 84-90.
!
!  Modified:
!
!    17 December 2000
!
!  Author:
!
!    Lili Ju
!
!  Parameters:
!
!    Input, integer N, the index of the desired element.
!
!    Input, integer P, the Halton base, which should be a prime number.
!
!    Output, real R, the N-th element of the Halton sequence
!    for base P.
!
  implicit none
!
  integer digit
  integer n
  integer n2
  integer p
  real pp
  real r
!
  n2 = abs ( n )
  r = 0.0E+00
  pp = 1.0E+00 / real ( p )

  do while ( n2 /= 0 )
    digit = mod ( n2, p )
    r = r + digit * pp
    pp = pp / real ( p )
    n2 = n2 / p
  end do

  if ( n < 0 ) then
    r = - r
  end if

  return
end
subroutine i_to_ivec_binary ( i, n, ivec )
!
!*******************************************************************************
!
!! I_TO_IVEC_BINARY makes a vector binary representation of an integer.
!
!
!  Example:
!
!     I       IVEC
!    --  --------------
!     1  0,0,...0,0,0,1
!     2  0,0,...0,0,1,0
!     3  0,0,...0,0,1,1
!     4  0,0,...0,1,0,0
!     9  0,0,...1,0,0,1
!    -9  1,1,...0,1,1,1
!
!  Discussion:
!
!    The vector is assumed to have dimension 32.  Negative values have
!    a two's complement operation applied.
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, an integer to be represented.
!
!    Input, integer N, the dimension of the vector, which should be 32.
!
!    Output, integer IVEC(N), the binary representation.
!
  implicit none
!
  integer n
!
  integer i
  integer i_copy
  integer ivec(n)
  integer j
!
  i_copy = abs ( i )

  do j = n, 2, -1

    ivec(j) = mod ( i_copy, 2 )

    i_copy = i_copy / 2

  end do

  ivec(1) = 0

  if ( i < 0 ) then
    call ivec_binary_complement2 ( n, ivec )
  end if

  return
end
subroutine i_to_s_left ( intval, s )
!
!*******************************************************************************
!
!! I_TO_S_LEFT converts an integer to a left-justified string.
!
!
!  Examples:
!
!    Assume that S is 6 characters long:
!
!    INTVAL  S
!
!         1  1
!        -1  -1
!         0  0
!      1952  1952
!    123456  123456
!   1234567  ******  <-- Not enough room!
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer INTVAL, an integer to be converted.
!
!    Output, character ( len = * ) S, the representation of the integer.
!    The integer will be left-justified.  If there is not enough space,
!    the string will be filled with stars.
!
  implicit none
!
  character c
  integer i
  integer idig
  integer ihi
  integer ilo
  integer intval
  integer ipos
  integer ival
  character ( len = * ) s
!
  s = ' '

  ilo = 1
  ihi = len ( s )

  if ( ihi <= 0 ) then
    return
  end if
!
!  Make a copy of the integer.
!
  ival = intval
!
!  Handle the negative sign.
!
  if ( ival < 0 ) then

    if ( ihi <= 1 ) then
      s(1:1) = '*'
      return
    end if

    ival = - ival
    s(1:1) = '-'
    ilo = 2

  end if
!
!  The absolute value of the integer goes into S(ILO:IHI).
!
  ipos = ihi
!
!  Find the last digit of IVAL, strip it off, and stick it into the string.
!
  do

    idig = mod ( ival, 10 )
    ival = ival / 10

    if ( ipos < ilo ) then
      do i = 1, ihi
        s(i:i) = '*'
      end do
      return
    end if

    call digit_to_ch ( idig, c )

    s(ipos:ipos) = c
    ipos = ipos - 1

    if ( ival == 0 ) then
      exit
    end if

  end do
!
!  Shift the string to the left.
!
  s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
  s(ilo+ihi-ipos:ihi) = ' '

  return
end
subroutine imat_01_rowcolsum ( m, n, r, c, a, ierror )
!
!*******************************************************************************
!
!! IMAT_01_ROWCOLSUM creates a 0/1 matrix with given row and column sums.
!
!
!  Discussion:
!
!    Given an M vector R and N vector C, there may exist one or more
!    M by N matrices with entries that are 0 or 1, whose row sums are R
!    and column sums are C.
!
!    For convenience, this routine requires that the entries of R and C
!    be given in nonincreasing order.
!
!    There are several requirements on R and C.  The simple requirements
!    are that the entries of R and C must be nonnegative, that the entries
!    of R must each be no greater than N, and those of C no greater than M,
!    and that the sum of the entries of R must equal the sum of the entries 
!    of C.
!
!    The final technical requirement is that if we form R*, the conjugate
!    partition of R, then C is majorized by R*, that is, that every partial
!    sum from 1 to K of the entries of C is no bigger than the sum of the same
!    entries of R*, for every K from 1 to N.
!
!    Given these conditions on R and C, there is at least one 0/1 matrix
!    with the given row and column sums.
!
!    The conjugate partition of R is constructed as follows:
!      R*(1) is the number of entries of R that are 1 or greater.
!      R*(2) is the number of entries of R that are 2 or greater.
!      ...
!      R*(N) is the number of entries of R that are N (can't be greater).
!
!  Example:
!
!    M = N = 5
!    R = ( 3, 2, 2, 1, 1 )
!    C = ( 2, 2, 2, 2, 1 )
!
!    A =
!      1 0 1 0 1
!      1 0 0 1 0
!      0 1 0 1 0
!      0 1 0 0 0
!      0 0 1 0 0
!
!  Reference:
!
!    J H van Lint and R M Wilson,
!    A Course in Combinatorics,
!    Oxford, 1992, pages 148-156.
!
!    James Sandeson,
!    Testing Ecological Patterns,
!    American Scientist,
!    Volume 88, July-August 2000, pages 332-339.
!
!    Ian Saunders,
!    Algorithm AS 205,
!    Enumeration of R x C Tables with Repeated Row Totals,
!    Applied Statistics,
!    Volume 33, Number 3, pages 340-352, 1984.
!
!  Modified:
!
!    26 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the array.
!
!    Input, integer R(M), C(N), the row and column sums desired for the array.
!    Both vectors must be arranged in descending order.
!    The elements of R must be between 0 and N.
!    The elements of C must be between 0 and M.
!
!    Output, integer A(M,N), the M by N matrix with the given row and
!    column sums.
!    Each entry of A is 0 or 1.
!
!    Output, integer IERROR, an error flag.
!    0, no error was encountered, and the array was computed.
!    1, R and C do not have the same total.
!    2, R is not monotone decreasing, or has illegal entries.
!    3, C is not monotone decreasing, or has illegal entries.
!    4, R and C are not a possible set of row and column sums.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer c(n)
  integer c_sum
  integer i
  integer ierror
  logical ivec_descends
  integer ivec_imax_last
  integer j
  integer k
  integer r(m)
  integer r_conj(n)
  integer r_sum
  integer r2(m)
!
  a(1:m,1:n) = 0
!
!  Check conditions.
!
  ierror = 0

  if ( sum ( r(1:m) ) /= sum ( c(1:n) ) ) then
    ierror = 1
    return
  end if

  if ( .not. ivec_descends ( m, r ) ) then
    ierror = 2
    return
  end if

  if ( r(1) > n .or. r(m) < 0 ) then
    ierror = 2
    return
  end if

  if ( .not. ivec_descends ( n, c ) ) then
    ierror = 3
    return
  end if

  if ( c(1) > m .or. c(n) < 0 ) then
    ierror = 3
    return
  end if
!
!  Compute the conjugate of R.
!
  r_conj(1:n) = 0

  do i = 1, m
    do j = 1, r(i)
      r_conj(j) = r_conj(j) + 1
    end do
  end do
!
!  C must be majorized by R_CONJ.
!
  r_sum = 0
  c_sum = 0
  do i = 1, n
    r_sum = r_sum + r_conj(i)
    c_sum = c_sum + c(i)
    if ( c_sum > r_sum ) then
      ierror = 4
      return
    end if
  end do
!
!  We need a temporary copy of R that we can decrement.
!
  r2(1:m) = r(1:m)

  do j = n, 1, -1

    i = ivec_imax_last ( m, r2 )

    do k = 1, c(j)
!
!  By adding 1 rather than setting A(I,J) to 1, we were able to spot
!  an error where the index was "sticking".
!
      a(i,j) = a(i,j) + 1

      r2(i) = r2(i) - 1
!
!  There's a special case you have to watch out for.
!  If I was 1, and when you decrement R2(1), I is going to be 1 again,
!  and you're staying in the same column, that's not good.
!
      if ( i > 1 ) then
        i = i - 1
      else
        i = ivec_imax_last ( m, r2 )
        if ( i == 1 .and. k < c(j) ) then
          i = 1 + ivec_imax_last ( m-1, r2(2:m) )
        end if
      end if

    end do

  end do

  return
end
subroutine imat_01_rowcolsum2 ( m, n, r, c, a, ierror )
!
!*******************************************************************************
!
!! IMAT_01_ROWCOLSUM2 creates a 0/1 matrix with given row and column sums.
!
!
!  Discussion:
!
!    This routine uses network flow optimization to compute the results.
!
!    Given an M vector R and N vector C, there may exist one or more
!    M by N matrices with entries that are 0 or 1, whose row sums are R
!    and column sums are C.
!
!    For convenience, this routine requires that the entries of R and C
!    be given in nonincreasing order.
!
!    There are several requirements on R and C.  The simple requirements
!    are that the entries of R and C must be nonnegative, that the entries
!    of R must each no greater than N, and those of C no greater than M,
!    and that the sum of the entries of R must equal the sum of the 
!    entries of C.
!
!    The final technical requirement is that if we form R*, the conjugate
!    partition of R, then C is majorized by R*, that is, that every partial
!    sum from 1 to K of the entries of C is no bigger than the sum of the same
!    entries of R*, for every K from 1 to N.
!
!    Given these conditions on R and C, there is at least one 0/1 matrix
!    with the given row and column sums.
!
!    The conjugate partition of R is constructed as follows:
!      R*(1) is the number of entries of R that are 1 or greater.
!      R*(2) is the number of entries of R that are 2 or greater.
!      ...
!      R*(N) is the number of entries of R that are N (can't be greater).
!
!  Example:
!
!    M = N = 5
!    R = ( 3, 2, 2, 1, 1 )
!    C = ( 2, 2, 2, 2, 1 )
!
!    A =
!      1 0 1 0 1
!      1 0 0 1 0
!      0 1 0 1 0
!      0 1 0 0 0
!      0 0 1 0 0
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!    J H van Lint and R M Wilson,
!    A Course in Combinatorics,
!    Oxford, 1992, pages 148-156.
!
!    James Sandeson,
!    Testing Ecological Patterns,
!    American Scientist,
!    Volume 88, July-August 2000, pages 332-339.
!
!  Modified:
!
!    27 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the array.
!    These values do not have to be equal.
!
!    Input, integer R(M), C(N), the row and column sums desired for the array.
!    Both vectors must be arranged in descending order.
!    The elements of R must be between 0 and N.
!    The elements of C must be between 0 and M.
!    One of the conditions for a solution to exist is that the sum of the
!    elements in R equal the sum of the elements in C.
!
!    Output, integer A(M,N), the matrix with the given row and column sums.
!    Each entry of A is 0 or 1.
!
!    Output, integer IERROR, an error flag.
!    0, no error was encountered, and the array was computed.
!    1, R and C are not consistent.  A partial solution may be constructed.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer c(n)
  integer capflo(2,2*(m+m*n+n))
  integer flow(2*(m+m*n+n))
  integer i
  integer icut(m+n+2)
  integer iendpt(2,2*(m+m*n+n))
  integer ierror
  integer isink
  integer isorce
  integer j
  integer k
  integer nedge
  integer nnode
  integer node_flow(m+n+2)
  integer r(m)
!
  ierror = 0
!
!  There are M + N + 2 nodes.  The last two are the special source and sink.
!
  isorce = m + n + 1
  isink = m + n + 2
  nnode = m + n + 2
!
!  The source is connected to each of the R nodes.
!
  k = 0

  do i = 1, m

    k = k + 1
    iendpt(1,k) = isorce
    iendpt(2,k) = i
    capflo(1,k) = r(i)
    capflo(2,k) = 0

    k = k + 1
    iendpt(1,k) = i
    iendpt(2,k) = isorce
    capflo(1,k) = r(i)
    capflo(2,k) = 0

  end do
!
!  Every R node is connected to every C node, with capacity 1.
!
  do i = 1, m
    do j = 1, n

      k = k + 1
      iendpt(1,k) = i
      iendpt(2,k) = j+m
      capflo(1,k) = 1
      capflo(2,k) = 0

      k = k + 1
      iendpt(1,k) = j+m
      iendpt(2,k) = i
      capflo(1,k) = 1
      capflo(2,k) = 0

    end do
  end do
!
!  Every C node is connected to the sink.
!
  do j = 1, n

    k = k + 1
    iendpt(1,k) = j+m
    iendpt(2,k) = isink
    capflo(1,k) = c(j)
    capflo(2,k) = 0

    k = k + 1
    iendpt(1,k) = isink
    iendpt(2,k) = j+m
    capflo(1,k) = c(j)
    capflo(2,k) = 0

  end do
!
!  Determine the maximum flow on the network.
!
  nedge = k

  call network_flow_max ( nnode, nedge, iendpt, capflo, isorce, isink, &
    icut, node_flow )
!
!  We have a perfect solution if, and only if, the edges leading from the
!  source, and the edges leading to the sink, are all saturated.
!
  do k = 1, nedge

    i = iendpt(1,k)
    j = iendpt(2,k) - m

    if ( i <= m .and. 1 <= j .and. j <= n ) then
      if ( capflo(2,k) /= 0 .and. capflo(2,k) /= 1 ) then
        ierror = 1
      end if
    end if

    if ( iendpt(1,k) == isorce ) then
      if ( capflo(1,k) /= capflo(2,k) ) then
        ierror = 1
      end if
    end if

    if ( iendpt(2,k) == isink ) then
      if ( capflo(1,k) /= capflo(2,k) ) then
        ierror = 1
      end if
    end if

  end do
!
!  If we have a solution, then A(I,J) = the flow on the edge from
!  R node I to C node J.
!
  a(1:m,1:n) = 0

  do k = 1, nedge

    i = iendpt(1,k)
    j = iendpt(2,k) - m

    if ( i <= m .and. 1 <= j .and. j <= n ) then
      a(i,j) = capflo(2,k)
    end if

  end do

  return
end
subroutine imat_inverse ( n, a, b )
!
!*******************************************************************************
!
!! IMAT_INVERSE inverts a unit upper triangular matrix.
!
!
!  Discussion:
!
!    A unit upper triangular matrix is a matrix with only 1's on the main
!    diagonal, and only 0's below the main diagonal.  Above the main
!    diagonal, the entries may be assigned any value.
!
!    It may be surprising to note that the inverse of an integer unit upper
!    triangular matrix is also an integer unit upper triangular matrix.
!
!    Note that this routine can invert a matrix in place, that is, with no
!    extra storage.  If the matrix is stored in A, then the call
!
!      call imat_inverse ( n, a, a )
!
!    will result in A being overwritten by its inverse, which can
!    save storage if the original value of A is not needed later.
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
!    05 March 1999
!
!  Parameters:
!
!    Input, integer N, number of rows and columns in matrix.
!
!    Input, integer A(N,N).  Unit upper triangular matrix
!    to be inverted.
!
!    Output, integer B(N,N), the inverse matrix.
!
  implicit none
!
  integer n
!
  integer a(n,n)
  integer b(n,n)
  integer i
  integer isum
  integer j
  integer k
!
  do j = n, 1, -1

    do i = n, 1, -1

      if ( i == j ) then
        isum = 1
      else
        isum = 0
      end if

      do k = i + 1, j
        isum = isum - a(i,k) * b(k,j)
      end do

      b(i,j) = isum

    end do
  end do

  return
end
subroutine imat_perm ( n, a, p )
!
!*******************************************************************************
!
!! IMAT_PERM permutes the rows and columns of a square integer matrix.
!
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
!    27 July 2000
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, integer A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(N), the permutation.  P(I) is the new number of row
!    and column I.
!
  implicit none
!
  integer n
!
  integer a(n,n)
  integer i
  integer i1
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(n)
!
  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call i_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine imat_perm2 ( m, n, a, p, q )
!
!*******************************************************************************
!
!! IMAT_PERM2 permutes the rows and columns of a rectangular integer matrix.
!
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
!    28 October 1999
!
!  Parameters:
!
!    Input, integer M, number of rows in the matrix.
!
!    Input, integer N, number of columns in the matrix.
!
!    Input/output, integer A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(M), the row permutation.  P(I) is the new number of row I.
!
!    Input, integer Q(N).  The column permutation.  Q(I) is the new number
!    of column I.  Note that this routine allows you to pass a single array
!    as both P and Q.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer i
  integer i1
  integer is
  integer it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(m)
  integer q(n)
!
  call perm_cycle ( m, p, is, nc, 1 )

  if ( q(1) > 0 ) then
    call perm_cycle ( n, q, is, nc, 1 )
  end if

  do i = 1, m

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call i_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then
    q(1:n) = abs ( q(1:n) )
  end if

  return
end
subroutine imat_print ( m, n, a, title )
!
!*******************************************************************************
!
!! IMAT_PRINT prints an integer matrix.
!
!
!  Modified:
!
!    08 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, integer A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  do jlo = 1, n, 10
    jhi = min ( jlo + 9, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,10(i7))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,10i7)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
!
!
!  Discussion:
!
!    The box is has center at (IC,JC), and has half-widths N1 and N2.
!    The indices are exactly those which are between (IC-N1,JC-N1) and
!    (IC+N1,JC+N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the half-widths of the box, that is, the
!    maximum distance allowed between (IC,JC) and (I,J).
!
!    Input, integer IC, JC, the central cell of the box.
!
!    Input/output, integer I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer ic
  integer j
  integer jc
  logical more
  integer n1
  integer n2
!
  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. ( i == ic - n1 .or. i == ic + n1 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, more )
!
!*******************************************************************************
!
!! INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
!
!
!  Discussion:
!
!    The box has a central cell of (IC,JC,KC), with a half widths of
!    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
!    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
!    maximum distances from the central cell allowed for I, J and K.
!
!    Input, integer IC, JC, KC, the central cell of the box.
!
!    Input/output, integer I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer ic
  integer j
  integer jc
  integer k
  integer kc
  logical more
  integer n1
  integer n2
  integer n3
!
  if ( .not. more ) then
    more = .true.
    i = ic - n1
    j = jc - n2
    k = kc - n3
    return
  end if

  if ( i == ic + n1 .and. j == jc + n2 .and. k == kc + n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( k > kc + n3 ) then
    k = kc - n3
    j = j + 1
  else if ( k < kc + n3 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      j == jc - n2 .or. j == jc + n2 ) ) then
    return
  else
    k = kc + n3
    return
  end if
!
!  Check J.
!
  if ( j > jc + n2 ) then
    j = jc - n2
    i = i + 1
  else if ( j < jc + n2 .and. &
    ( i == ic - n1 .or. i == ic + n1 .or. &
      k == kc - n3 .or. k == kc + n3 ) ) then
    return
  else
    j = jc + n2
    return
  end if

  return
end
subroutine index_box_next_2d ( n1, n2, i, j, more )
!
!*******************************************************************************
!
!! INDEX_BOX_NEXT_2D produces indices on the surface of a box in 2D.
!
!
!  Discussion:
!
!    The indices are exactly those which are between (1,1) and
!    (N1,N2) with the property that at least one of I and J
!    is an "extreme" value.
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, the "dimensions" of the box, that is, the
!    maximum values allowed for I and J.  The minimum values are
!    assumed to be 1.
!
!    Input/output, integer I, J.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I and J are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer j
  logical more
  integer n1
  integer n2
!
  if ( .not. more ) then
    more = .true.
    i = 1
    j = 1
    return
  end if

  if ( i == n1 .and. j == n2 ) then
    more = .false.
    return
  end if
!
!  Increment J.
!
  j = j + 1
!
!  Check J.
!
  if ( j > n2 ) then
    j = 1
    i = i + 1
  else if ( j < n2 .and. ( i == 1 .or. i == n1 ) ) then
    return
  else
    j = n2
    return
  end if

  return
end
subroutine index_box_next_3d ( n1, n2, n3, i, j, k, more )
!
!*******************************************************************************
!
!! INDEX_BOX_NEXT_3D produces indices on the surface of a box in 3D.
!
!
!  Discussion:
!
!    The indices are exactly those which are between (1,1,1) and
!    (N1,N2,N3) with the property that at least one of I, J, and K
!    is an "extreme" value.
!
!  Modified:
!
!    13 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N1, N2, N3, the "dimensions" of the box, that is, the
!    maximum values allowed for I, J and K.  The minimum values are
!    assumed to be 1.
!
!    Input/output, integer I, J, K.  On input, the previous index set.
!    On output, the next index set.  On the first call, MORE should
!    be set to FALSE, and the input values of I, J, and K are ignored.
!
!    Input/output, logical MORE.
!    On the first call for a given box, the user should set MORE to FALSE.
!    On return, the routine sets MORE to TRUE.
!    When there are no more indices, the routine sets MORE to FALSE.
!
  implicit none
!
  integer i
  integer j
  integer k
  logical more
  integer n1
  integer n2
  integer n3
!
  if ( .not. more ) then
    more = .true.
    i = 1
    j = 1
    k = 1
    return
  end if

  if ( i == n1 .and. j == n2 .and. k == n3 ) then
    more = .false.
    return
  end if
!
!  Increment K.
!
  k = k + 1
!
!  Check K.
!
  if ( k > n3 ) then
    k = 1
    j = j + 1
  else if ( k < n3 .and. &
    ( i == 1 .or. i == n1 .or. j == 1 .or. j == n2 ) ) then
    return
  else
    k = n3
    return
  end if
!
!  Check J.
!
  if ( j > n2 ) then
    j = 1
    i = i + 1
  else if ( j < n2 .and. &
    ( i == 1 .or. i == n1 .or. k == 1 .or. k == n3 ) ) then
    return
  else
    j = n2
    return
  end if

  return
end
subroutine index_next0 ( n, hi, iarray, more )
!
!*******************************************************************************
!
!! INDEX_NEXT0 generates all indices within given upper limits.
!
!
!  Discussion:
!
!    The indices are generated in such a way that the reversed
!    sequences are produced in lexicographic order.
!
!  Example:
!
!    N = 3,
!    HI = 3
!
!    1   2   3
!    ---------
!    1   1   1
!    2   1   1
!    3   1   1
!    1   2   1
!    2   2   1
!    3   2   1
!    1   3   1
!    2   3   1
!    3   3   1
!    1   1   2
!    2   1   2
!    3   1   2
!    1   2   2
!    2   2   2
!    3   2   2
!    1   3   2
!    2   3   2
!    3   3   2
!    1   1   3
!    2   1   3
!    3   1   3
!    1   2   3
!    2   2   3
!    3   2   3
!    1   3   3
!    2   3   3
!    3   3   3
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI, the upper limit for the array indices.
!    The lower limit is implicitly 1 and HI must be at least 1.
!
!    Input/output, integer IARRAY(N).
!    On startup calls, with MORE = .FALSE., the input value of IARRAY
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = .TRUE., the input value of IARRAY must be
!    the output value of IARRAY from the previous call.  (In other words,
!    just leave it alone!).
!    On output, IARRAY contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
  implicit none
!
  integer n
!
  integer hi
  integer i
  integer iarray(n)
  integer inc
  logical more
!
  if ( .not. more ) then

    iarray(1:n) = 1

    if ( hi < 1 ) then
      more = .false.
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INDEX_NEXT0 - Fatal error!'
      write ( *, '(a,i6)' ) '  HI is ', hi
      write ( *, '(a)' ) '  but HI must be at least 1.'
      return
    end if

  else

    inc = 1

    do while ( iarray(inc) >= hi )
      iarray(inc) = 1
      inc = inc + 1
    end do

    iarray(inc) = iarray(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( iarray(i) < hi ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_next1 ( n, hi, iarray, more )
!
!*******************************************************************************
!
!! INDEX_NEXT1 generates all indices within given upper limits.
!
!
!  Discussion:
!
!    The indices are generated in such a way that the reversed
!    sequences are produced in lexicographic order.
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!
!    1   2   3
!    ---------
!    1   1   1
!    2   1   1
!    3   1   1
!    4   1   1
!    1   2   1
!    2   2   1
!    3   2   1
!    4   2   1
!    1   1   2
!    2   1   2
!    3   1   2
!    4   1   2
!    1   2   2
!    2   2   2
!    3   2   2
!    4   2   2
!    1   1   3
!    2   1   3
!    3   1   3
!    4   1   3
!    1   2   3
!    2   2   3
!    3   2   3
!    4   2   3
!
!  Modified:
!
!    27 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input/output, integer IARRAY(N).
!    On startup calls, with MORE = .FALSE., the input value of IARRAY
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = .TRUE., the input value of IARRAY must be
!    the output value of IARRAY from the previous call.  (In other words,
!    just leave it alone!).
!    On output, IARRAY contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer inc
  logical more
!
  if ( .not. more ) then

    iarray(1:n) = 1

    do i = 1, n
      if ( hi(i) < 1 ) then
        more = .false.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INDEX_NEXT1 - Fatal error!'
        write ( *, '(a,i6,a,i6)' ) '  Entry ', i, ' of HI is ', hi(i)
        write ( *, '(a)' ) '  but all entries must be at least 1.'
        return
      end if
    end do

  else

    inc = 1

    do while ( iarray(inc) >= hi(inc) )
      iarray(inc) = 1
      inc = inc + 1
    end do

    iarray(inc) = iarray(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( iarray(i) < hi(i) ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_next2 ( n, lo, hi, iarray, more )
!
!*******************************************************************************
!
!! INDEX_NEXT2 generates all indices within given lower and upper limits.
!
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!
!    1   2   3
!    ---------
!    1  10   4
!    2  10   4
!    1  11   4
!    2  11   4
!    1  10   5
!    2  10   5
!    1  11   5
!    2  11   5
!    1  10   6
!    2  10   6
!    1  11   6
!    2  11   6
!
!  Modified:
!
!    21 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.  The rank of
!    the object being indexed.
!
!    Input, integer LO(N), HI(N), the lower and upper limits for the array
!    indices.  LO(I) should be less than or equal to HI(I), for each I.
!
!    Input/output, integer IARRAY(N).
!    On startup calls, with MORE = .FALSE., the input value of IARRAY
!    doesn't matter, because the routine initializes it.
!    On calls with MORE = .TRUE., the input value of IARRAY must be
!    the output value of IARRAY from the previous call.  (In other words,
!    just leave it alone!).
!    On output, IARRAY contains the successor set of indices to the input
!    value.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer inc
  integer lo(n)
  logical more
!
  if ( .not. more ) then

    iarray(1:n) = lo(1:n)

    do i = 1, n
      if ( hi(i) < lo(i) ) then
        more = .false.
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INDEX_NEXT2 - Fatal error!'
        write ( *, '(a,i6,a,i6)' ) '  Entry ', i, ' of HI is ', hi(i)
        write ( *, '(a,i6,a,i6)' ) '  Entry ', i, ' of LO is ', lo(i)
        write ( *, '(a)' ) '  but LO(I) <= HI(I) is required.'
        return
      end if
    end do

  else

    inc = 1

    do while ( iarray(inc) >= hi(inc) )
      iarray(inc) = lo(inc)
      inc = inc + 1
    end do

    iarray(inc) = iarray(inc) + 1

  end if
!
!  See if there are more entries to compute.
!
  more = .false.

  do i = 1, n
    if ( iarray(i) < hi(i) ) then
      more = .true.
    end if
  end do

  return
end
subroutine index_rank0 ( n, hi, iarray, rank )
!
!*******************************************************************************
!
!! INDEX_RANK0 ranks an index within given upper limits.
!
!
!  Example:
!
!    N = 3,
!    HI = 3
!    IARRAY = ( 3, 1, 2 )
!
!    RANK = 12
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI, the upper limit for the array indices.
!    The lower limit is implicitly 1, and HI should be at least 1.
!
!    Input, integer IARRAY(N), the index to be ranked.
!
!    Output, integer RANK, the rank of the index, or -1 if IARRAY
!    is not a legal index.
!
  implicit none
!
  integer n
!
  integer hi
  integer i
  integer iarray(n)
  integer range
  integer rank
!
  rank = - 1
  do i = 1, n
    if ( iarray(i) < 1 .or. iarray(i) > hi ) then
      return
    end if
  end do

  rank = 0
  do i = n, 1, -1
    rank = hi * rank + iarray(i)
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( iarray(i) - 1 ) * range
    range = range * hi
  end do

  return
end
subroutine index_rank1 ( n, hi, iarray, rank )
!
!*******************************************************************************
!
!! INDEX_RANK1 ranks an index within given upper limits.
!
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!    IARRAY = ( 4, 1, 2 )
!
!    RANK = 12
!
!  Modified:
!
!    22 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input, integer IARRAY(N), the index to be ranked.
!
!    Output, integer RANK, the rank of the index, or -1 if IARRAY
!    is not a legal index.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer range
  integer rank
!
  rank = - 1
  do i = 1, n
    if ( iarray(i) < 1 .or. iarray(i) > hi(i) ) then
      return
    end if
  end do

  rank = 0
  do i = n, 1, -1
    rank = hi(i) * rank + iarray(i)
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( iarray(i) - 1 ) * range
    range = range * hi(i)
  end do

  return
end
subroutine index_rank2 ( n, lo, hi, iarray, rank )
!
!*******************************************************************************
!
!! INDEX_RANK2 ranks an index within given lower and upper limits.
!
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!    IARRAY = ( 1, 11, 5 )
!
!    RANK = 7
!
!  Modified:
!
!    22 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer LO(N), HI(N), the lower and upper limits for the array
!    indices.  LO(I) should be less than or equal to HI(I), for each I.
!
!    Input, integer IARRAY(N), the index to be ranked.
!
!    Output, integer RANK, the rank of the index, or -1 if IARRAY
!    is not a legal index.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer lo(n)
  integer range
  integer rank
!
  rank = - 1
  do i = 1, n
    if ( iarray(i) < lo(i) .or. iarray(i) > hi(i) ) then
      return
    end if
  end do

  rank = 1
  range = 1
  do i = 1, n
    rank = rank + ( iarray(i) - lo(i) ) * range
    range = range * ( hi(i) + 1 - lo(i) )
  end do

  return
end
subroutine index_unrank0 ( n, hi, rank, iarray )
!
!*******************************************************************************
!
!! INDEX_UNRANK0 unranks an index within given upper limits.
!
!
!  Example:
!
!    N = 3,
!    HI = 3
!    RANK = 12
!
!    IARRAY = ( 3, 1, 2 )
!
!  Modified:
!
!    05 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI, the upper limit for the array indices.
!    The lower limit is implicitly 1, and HI should be at least 1.
!
!    Input, integer RANK, the rank of the desired index.
!
!    Output, integer IARRAY(N), the index of the given rank.
!
  implicit none
!
  integer n
!
  integer hi
  integer i
  integer iarray(n)
  integer j
  integer k
  integer range
  integer rank
!
  iarray(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = hi**n
!
!  The rank might be too large.
!
  if ( rank > range ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / hi
    j = k / range
    iarray(i) = j + 1
    k = k - j * range
  end do

  return
end
subroutine index_unrank1 ( n, hi, rank, iarray )
!
!*******************************************************************************
!
!! INDEX_UNRANK1 unranks an index within given upper limits.
!
!
!  Example:
!
!    N = 3,
!    HI(1) = 4, HI(2) = 2, HI(3) = 3
!    RANK = 11
!
!    IARRAY = ( 3, 1, 2 )
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer HI(N), the upper limits for the array indices.
!    The lower limit is implicitly 1, and each HI(I) should be at least 1.
!
!    Input, integer RANK, the rank of the desired index.
!
!    Output, integer IARRAY(N), the index of the given rank.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer j
  integer k
  integer range
  integer rank
!
  iarray(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = product ( hi )
!
!  The rank might be too large.
!
  if ( rank > range ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / hi(i)
    j = k / range
    iarray(i) = j + 1
    k = k - j * range
  end do

  return
end
subroutine index_unrank2 ( n, lo, hi, rank, iarray )
!
!*******************************************************************************
!
!! INDEX_UNRANK2 unranks an index within given lower and upper limits.
!
!
!  Example:
!
!    N = 3,
!    LO(1) = 1, LO(2) = 10, LO(3) = 4
!    HI(1) = 2, HI(2) = 11, HI(3) = 6
!    RANK = 7
!
!    IARRAY = ( 1, 11, 5 )
!
!  Modified:
!
!    24 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in IARRAY.
!
!    Input, integer LO(N), HI(N), the lower and upper limits for the array
!    indices.  It should be the case that LO(I) <= HI(I) for each I.
!
!    Input, integer RANK, the rank of the desired index.
!
!    Output, integer IARRAY(N), the index of the given rank.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer iarray(n)
  integer j
  integer k
  integer lo(n)
  integer range
  integer rank
!
  iarray(1:n) = 0
!
!  The rank might be too small.
!
  if ( rank < 1 ) then
    return
  end if

  range = 1
  do i = 1, n
    range = range * ( hi(i) + 1 - lo(i) )
  end do
!
!  The rank might be too large.
!
  if ( rank > range ) then
    return
  end if

  k = rank - 1
  do i = n, 1, -1
    range = range / ( hi(i) + 1 - lo(i) )
    j = k / range
    iarray(i) = j + lo(i)
    k = k - j * range
  end do

  return
end
subroutine ins_perm ( n, ins, p )
!
!*******************************************************************************
!
!! INS_PERM computes a permutation from its inversion sequence.
!
!
!  Definition:
!
!    For a given permutation P acting on objects 1 through N, the
!    inversion sequence INS is defined as:
!
!      INS(1) = 0
!      INS(I) = number of values J < I for which P(J) > P(I).
!
!  Example:
!
!    Input:
!
!      ( 0, 0, 2, 1, 3 )
!
!    Output:
!
!      ( 3, 5, 1, 4, 2 )
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    17 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input, integer INS(N), the inversion sequence of a permutation.
!    It must be the case that 0 <= INS(I) < I for I = 1 to N.
!
!    Output, integer P(N), the permutation.
!
  implicit none
!
  integer n
!
  integer i
  integer ins(n)
  integer itemp
  integer j
  integer p(n)
!
  call ivec_identity ( n, p )

  do i = n, 2, -1

    itemp = p(i-ins(i))

    do j = i-ins(i), i-1
      p(j) = p(j+1)
    end do

    p(i) = itemp

  end do

  return
end
subroutine intp_conj ( iarray1, mult1, npart1, iarray2, mult2, npart2, n )
!
!*******************************************************************************
!
!! INTP_CONJ computes the conjugate of a partition.
!
!
!  Definition:
!
!    A partition of an integer N is a set of positive integers which
!    add up to N.  The conjugate of a partition P1 of N is another partition
!    P2 of N obtained in the following way:
!
!      The first element of P2 is the number of parts of P1 greater than
!      or equal to 1.
!
!      The K-th element of P2 is the number of parts of P1 greater than
!      or equal to K.
!
!    Clearly, P2 will have no more than N elements; it may be surprising
!    to find that P2 is guaranteed to be a partition of N.  However, if
!    we symbolize the initial partition P1 by rows of X's, then we can
!    see that P2 is simply produced by grouping by columns:
!
!        6 3 2 2 1
!      5 X X X X X
!      4 X X X X
!      2 X X
!      1 X
!      1 X
!      1 X
!
!  Example:
!
!    14 = 5 + 4 + 2 + 1 + 1 + 1
!
!    The conjugate partition is:
!
!    14 = 6 + 3 + 2 + 2 + 1
!
!  Modified:
!
!    20 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer IARRAY1(NPART1).  IARRAY1 contains the parts of
!    the partition.  The value of N is represented by
!
!      sum ( 1 <= I <= NPART1 ) MULT1(I) * IARRAY1(I).
!
!    Intput, integer MULT1(NPART1).  MULT1 counts the multiplicity of
!    the parts of the partition.  MULT1(I) is the multiplicity
!    of the part IARRAY1(I), for I = 1 to NPART1.
!
!    Input, integer N, the integer to be partitioned.
!
!    Input, integer NPART1, the number of "parts" in the partition.
!
!    Output, integer IARRAY2(N).  IARRAY contains the parts of
!    the conjugate partition in entries 1 through NPART2.
!
!    Output, integer MULT2(N).  MULT2 counts the multiplicity of
!    the parts of the conjugate partition in entries 1 through NPART2.
!
!    Output, integer NPART2, the number of "parts" in the conjugate partition.
!
!    Input, integer N, the integer to be partitioned.
!
  implicit none
!
  integer n
  integer npart1
!
  integer i
  integer iarray1(npart1)
  integer iarray2(n)
  integer itemp
  integer itest
  integer mult1(npart1)
  integer mult2(n)
  integer npart2
!
  iarray2(1:n) = 0
  mult2(1:n) = 0
  npart2 = 0

  itest = 0

  do

    itest = itest + 1

    itemp = 0

    do i = 1, npart1
      if ( iarray1(i) >= itest ) then
        itemp = itemp + mult1(i)
      end if
    end do

    if ( itemp <= 0 ) then
      exit
    end if

    if ( npart2 > 0 ) then
      if ( itemp == iarray2(npart2) ) then
        mult2(npart2) = mult2(npart2) + 1
      else
        npart2 = npart2 + 1
        iarray2(npart2) = itemp
        mult2(npart2) = 1
      end if
    else
      npart2 = npart2 + 1
      iarray2(npart2) = itemp
      mult2(npart2) = 1
    end if

  end do

  return
end
subroutine intp_enum ( n, p )
!
!*******************************************************************************
!
!! INTP_ENUM computes the number of partitions of an integer.
!
!
!  Method:
!
!    Partition numbers are difficult to compute.  This routine uses
!    Euler's method, which observes that:
!
!      P(0) = 1
!      P(N) =   P(N-1)  + P(N-2)
!             - P(N-5)  - P(N-7)
!             + P(N-12) + P(N-15)
!             - ...
!      where the numbers subtracted from N are the pentagonal numbers,
!      (with both positive and negative indices) with the computation
!      stopping when a negative index is reached.
!
!  First values:
!
!    N   P
!
!    0   1
!    1   1
!    2   2
!    3   3
!    4   5
!    5   7
!    6  11
!    7  15
!    8  22
!    9  30
!   10  42
!
!  Reference:
!
!    John Conway and Richard Guy,
!    The Book of Numbers,
!    Springer Verlag, 1996, page 95.
!
!  Modified:
!
!    18 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the highest partition number desired.
!
!    Output, integer P(0:N), the partition numbers.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  integer p(0:n)
  integer pj
  integer sgn
!
  p(0) = 1

  do i = 1, n

    p(i) = 0

    j = 0
    sgn = 1

    do

      j = j + 1
      call pent_enum ( j, pj )

      if ( pj > i ) then
        exit
      end if

      p(i) = p(i) + sgn * p(i-pj)
      sgn = - sgn

    end do

    j = 0
    sgn = 1

    do

      j = j - 1
      call pent_enum ( j, pj )

      if ( pj > i ) then
        exit
      end if

      p(i) = p(i) + sgn * p(i-pj)
      sgn = - sgn

    end do

  end do

  return
end
subroutine intp_next ( done, iarray, mult, n, npart )
!
!*******************************************************************************
!
!! INTP_NEXT generates the partitions of an integer, one at a time.
!
!
!  Comments:
!
!    The number of partitions of N is:
!
!      1    1
!      2    2
!      3    3
!      4    5
!      5    7
!      6   11
!      7   15
!      8   22
!      9   30
!     10   42
!     11   56
!     12   77
!     13  101
!     14  135
!     15  176
!     16  231
!     17  297
!     18  385
!     19  490
!     20  627
!     21  792
!     22 1002
!     23 1255
!     24 1575
!     25 1958
!     26 2436
!     27 3010
!     28 3718
!     29 4565
!     30 5604
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
!    Input/output, logical DONE.
!
!    On first call, the user should set DONE to .TRUE. to signal
!    that the program should initialize data.
!
!    On each return, the programs sets DONE to .FALSE. if it
!    has another partition to return.  If the program returns
!    with DONE .TRUE., then there are no more partitions.
!
!    Output, integer IARRAY(N).  IARRAY contains the parts of
!    the partition.  The value of N is represented by
!
!      sum ( 1 <= I <= NPART ) MULT(I) * IARRAY(I).
!
!    Output, integer MULT(N).  MULT counts the multiplicity of
!    the parts of the partition.  MULT(I) is the multiplicity
!    of the part IARRAY(I), for I = 1 to NPART.
!
!    Input, integer N, the integer to be partitioned.
!
!    Output, integer NPART, the number of "parts" in the partition.
!
  implicit none
!
  integer n
!
  logical done
  integer i
  integer iarray(n)
  integer is
  integer iu
  integer iv
  integer iw
  integer k
  integer k1
  integer mult(n)
  integer npart
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'INTP_NEXT - Fatal error!'
    write ( *, '(a)' ) '  N must be positive.'
    write ( *, '(a,i6)' ) '  The input value of N was ', n
    stop
  end if

  if ( done ) then

    iarray(1) = n
    iarray(2:n) = 0

    mult(1) = 1
    mult(2:n) = 0

    npart = 1
    done = .false.

  else

    if ( iarray(npart) > 1 .or. npart > 1 ) then

      done = .false.

      if ( iarray(npart) == 1 ) then
        is = iarray(npart-1) + mult(npart)
        k = npart - 1
      else
        is = iarray(npart)
        k = npart
      end if

      iw = iarray(k) - 1
      iu = is / iw
      iv = mod ( is, iw )
      mult(k) = mult(k) - 1

      if ( mult(k) == 0 ) then
        k1 = k
      else
        k1 = k + 1
      end if

      mult(k1) = iu
      iarray(k1) = iw

      if ( iv == 0 ) then
        npart = k1
      else
        mult(k1+1) = 1
        iarray(k1+1) = iv
        npart = k1 + 1
      end if

    else
      done = .true.
    end if

  end if

  return
end
subroutine intp_next2 ( n, iarray, mult, npart, more )
!
!*******************************************************************************
!
!! INTP_NEXT2 computes the partitions of the integer N one at a time.
!
!
!  Discussion:
!
!    Unlike compositions, order is not important in a partition.  Thus
!    the sequences 3+2+1 and 1+2+3 represent distinct compositions, but
!    not distinct partitions.  Also 0 is never returned as one of the
!    elements of the partition.
!
!  Examples:
!
!    Sample partitions of 6 include:
!
!      6 = 4+1+1 = 3+2+1 = 2+2+2
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the integer whose partitions are desired.
!
!    Output, integer IARRAY(N).  IARRAY(I) is the I-th distinct part
!    of the partition, for I = 1, NPART.  Note that if a certain number
!    shows up several times in the partition, it is listed only
!    once in IARRAY, and its multiplicity is counted in MULT.
!
!    Output, integer MULT(N).  MULT(I) is the multiplicity of IARRAY(I)
!    in the partition, for I = 1, NPART; that is, the number of repeated
!    times that IARRAY(I) is used in the partition.
!
!    Output, integer NPART, the number of distinct, nonzero parts in the
!    output partition.
!
!    Input/output, logical MORE.  Set MORE = .FALSE. on first call.  It
!    will be reset .TRUE. on return with the first partition.
!    Keep calling for more partitions until MORE
!    is returned .FALSE.
!
  implicit none
!
  integer n
!
  integer iarray(n)
  integer iff
  integer is
  integer isum
  logical more
  integer mult(n)
  integer, save :: nlast = 0
  integer npart
!
!  On the first call, set NLAST to 0.
!
  if ( .not. more ) then
    nlast = 0
  end if

  if ( n /= nlast .or. ( .not. more ) ) then
    nlast = n
    npart = 1
    iarray(npart) = n
    mult(npart) = 1
    more = mult(npart) /= n
    return
  end if

  isum = 1

  if ( iarray(npart) <= 1 ) then
    isum = mult(npart) + 1
    npart = npart - 1
  end if

  iff = iarray(npart) - 1

  if ( mult(npart) /= 1 ) then
    mult(npart) = mult(npart) - 1
    npart = npart + 1
  end if

  iarray(npart) = iff
  mult(npart) = 1 + ( isum / iff )
  is = mod ( isum, iff )

  if ( is > 0 ) then
    npart = npart + 1
    iarray(npart) = is
    mult(npart) = 1
  end if
!
!  There are more partitions, as long as we haven't just computed
!  the last one, which is N copies of 1.
!
  more = mult(npart) /= n

  return
end
subroutine intp_random ( n, iarray, mult, npart )
!
!*******************************************************************************
!
!! INTP_RANDOM selects a random partition of the integer N.
!
!
!  Discussion:
!
!    Note that some elements of the partition may be 0.  The partition is
!    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
!
!      N = sum ( 1 <= I <= N ) MULT(I) * I.
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
!    29 March 2001
!
!  Parameters:
!
!    Input, integer N, the integer to be partitioned.
!
!    Output, integer IARRAY(N), returns the number of partitions of each
!    integer from 1 to N.
!
!    Output, integer MULT(N).  MULT(I) is the multiplicity of I
!    in the chosen partition.
!
!    Output, integer NPART, the number of parts in the partition chosen,
!    that is, the number of integers I with nonzero multiplicity MULT(I).
!
  implicit none
!
  integer n
!
  integer i
  integer i1
  integer iarray(n)
  integer id
  integer is
  integer isum
  integer j
  integer m
  integer mult(n)
  integer, save :: nlast = 0
  integer npart
  real z
  real, parameter :: zhi = 1.0E+00
  real, parameter :: zlo = 0.0E+00
!
  if ( n > nlast ) then

    iarray(1) = 1
    m = nlast + 1
    nlast = n

    if ( n > 1 ) then

      do i = m, n

        isum = 0

        do id = 1, i

          is = 0
          i1 = i

          do

            i1 = i1 - id

            if ( i1 > 0 ) then
              is = is + iarray(i1)
            else
              if ( i1 == 0 ) then
                is = is + 1
              end if
              exit
            end if

          end do

          isum = isum + is * id

        end do

        iarray(i) = isum / i

      end do

    end if

  end if

  m = n
  npart = 0
  mult(1:n) = 0

  do while ( m > 0 )

    call r_random ( zlo, zhi, z )
    z = m * iarray(m) * z
    id = 0

30  continue

    id = id + 1
    i1 = m
    j = 0

    do

      j = j + 1
      i1 = i1 - id

      if ( i1 < 0 ) then
        go to 30
      end if

      if ( i1 == 0 ) then
        z = z - id
        if ( z > 0.0E+00 ) then
          go to 30
        else
          exit
        end if
      end if

      if ( i1 > 0 ) then
        z = z - id * iarray(i1)
        if ( z <= 0.0E+00 ) then
          exit
        end if
      end if

    end do

    mult(id) = mult(id) + j
    npart = npart + j
    m = i1

  end do

  return
end
subroutine involute_enum ( n, s )
!
!*******************************************************************************
!
!! INVOLUTE_ENUM enumerates the involutions of N objects.
!
!
!  Definition:
!
!    An involution is a permutation consisting only of fixed points and
!    pairwise transpositions.
!
!  Comments:
!
!    An involution is its own inverse permutation.
!
!  Recursion:
!
!    S(0) = 1
!    S(1) = 1
!    S(N) = S(N-1) + (N-1) * S(N-2)
!
!  First values:
!
!     N         S(N)
!     0           1
!     1           1
!     2           2
!     3           4
!     4          10
!     5          26
!     6          76
!     7         232
!     8         764
!     9        2620
!    10        9496
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer S(0:N), the number of involutions of 0, 1, 2, ... N
!    objects.
!
  implicit none
!
  integer n
!
  integer i
  integer s(0:n)
!
  if ( n < 0 ) then
    return
  end if

  s(0) = 1

  if ( n <= 0 ) then
    return
  end if

  s(1) = 1

  do i = 2, n
    s(i) = s(i-1) + real ( i - 1 ) * s(i-2)
  end do

  return
end
subroutine ipoly ( n, a, x0, iopt, val )
!
!*******************************************************************************
!
!! IPOLY performs operations on integer polynomials in power or factorial form.
!
!
!  Discussion:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*X**2 + ... + (AN+1)*X**N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1 + A2*(X-C) + A3*(X-C)**2 + ... + (AN+1)*(X-C)**N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2)+...
!        + (AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input/output, integer A(N), the coefficients of the polynomial.  Depending
!    on the option chosen, these coefficients may be overwritten by those
!    of a different form of the polynomial.
!
!    Input, integer X0, for IOPT = -1, 0, or positive, the value of the
!    argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Input, integer IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of the polynomial in
!    factorial form, output them in power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum form, output them
!    in factorial form.
!
!    -1: Evaluate a polynomial which has been input in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Output, integer VAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point X0.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer eps
  integer i
  integer iopt
  integer m
  integer n1
  integer val
  integer w
  integer x0
  integer z
!
  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )

  if ( iopt < -1 ) then
    n1 = n
  end if

  eps = mod ( max ( -iopt, 0 ), 2 )

  w = - n * eps

  if ( iopt > -2 ) then
    w = w + x0
  end if

  do m = 1, n1

    val = 0
    z = w

    do i = m, n
      z = z + eps
      val = a(n+m-i) + z * val
      if ( iopt /= 0 .and. iopt /= -1 ) then
        a(n+m-i) = val
      end if
    end do

    if ( iopt < 0 ) then
      w = w + 1
    end if

  end do

  return
end
subroutine ipoly_cyclo ( n, phi )
!
!*******************************************************************************
!
!! IPOLY_CYCLO computes a cyclotomic polynomial.
!
!
!  Discussion:
!
!    For N >= 1, let
!
!      I = SQRT ( - 1 )
!      L = EXP ( 2 * PI * I / N )
!
!    Then the N-th cyclotomic polynomial is defined by
!
!      PHI(N;X) = Product ( 1 <= K <= N and GCD(K,N) = 1 ) ( X - L**K )
!
!    We can use the Moebius MU function to write
!
!      PHI(N;X) = Product ( mod ( D, N ) = 0 ) ( X**D - 1 )**MU(N/D)
!
!    There is a sort of inversion formula:
!
!      X**N - 1 = Product ( mod ( D, N ) = 0 ) PHI(D;X)
!
!  Example:
!
!     N  PHI
!
!     0  1
!     1  X - 1
!     2  X + 1
!     3  X**2 + X + 1
!     4  X**2 + 1
!     5  X**4 + X**3 + X**2 + X + 1
!     6  X**2 - X + 1
!     7  X**6 + X**5 + X**4 + X**3 + X**2 + X + 1
!     8  X**4 + 1
!     9  X**6 + X**3 + 1
!    10  X**4 - X**3 + X**2 - X + 1
!
!  Reference:
!
!    Raymond Seroul,
!    Programming for Mathematicians,
!    Springer Verlag, 2000, page 269.
!
!  Modified:
!
!    08 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the cyclotomic polynomial desired.
!
!    Output, integer PHI(0:N), the N-th cyclotomic polynomial.
!
  implicit none
!
  integer n
  integer, parameter :: max_poly = 100
!
  integer d
  integer den(0:max_poly)
  integer den_n
  integer factor(0:n)
  integer mu
  integer nq
  integer nr
  integer num(0:max_poly)
  integer num_n
  integer phi(0:n)
  integer rem(0:n)
!
  num(0) = 1
  num(1:max_poly) = 0
  num_n = 0

  den(0) = 1
  den(1:max_poly) = 0
  den_n = 0

  phi(0:n) = 0

  do d = 1, n
!
!  For each divisor D of N, ...
!
    if ( mod ( n, d ) == 0 ) then

      call i_moebius ( n / d, mu )
!
!  ...multiply the numerator or denominator by (X^D-1).
!
      factor(0) = -1
      factor(1:d-1) = 0
      factor(d) = 1

      if ( mu == + 1 ) then

        if ( num_n + d > max_poly ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IPOLY_CYCLO - Fatal error!'
          write ( *, '(a)' ) '  Numerator polynomial degree too high.'
          return
        end if

        call ipoly_mul ( num_n, num, d, factor, num )

        num_n = num_n + d

      else if ( mu == -1 ) then

        if ( num_n + d > max_poly ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IPOLY_CYCLO - Fatal error!'
          write ( *, '(a)' ) '  Denominator polynomial degree too high.'
          return
        end if

        call ipoly_mul ( den_n, den, d, factor, den )

        den_n = den_n + d

      end if

    end if

  end do
!
!  PHI = NUM / DEN
!
  call ipoly_div ( num_n, num, den_n, den, nq, phi, nr, rem )

  return
end
subroutine ipoly_degree ( na, a, degree )
!
!*******************************************************************************
!
!! IPOLY_DEGREE returns the degree of a polynomial.
!
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, integer A(0:NA), the coefficients of the polynomials.
!
!    Output, integer DEGREE, the degree of A.
!
  implicit none
!
  integer na
!
  integer a(0:na)
  integer degree
!
  degree = na

  do while ( degree > 0 )

    if ( a(degree) /= 0 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine ipoly_div ( na, a, nb, b, nq, q, nr, r )
!
!*******************************************************************************
!
!! IPOLY_DIV computes the quotient and remainder of two polynomials.
!
!
!  Discussion:
!
!    Normally, the quotient and remainder would have rational coefficients.
!    This routine assumes that the special case applies that the quotient
!    and remainder are known beforehand to be integral.
!
!    The polynomials are assumed to be stored in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, integer A(0:NA), the coefficients of the polynomial to be divided.
!
!    Input, integer NB, the dimension of B.
!
!    Input, integer B(0:NB), the coefficients of the divisor polynomial.
!
!    Output, integer NQ, the degree of Q.
!    If the divisor polynomial is zero, NQ is returned as -1.
!
!    Output, integer Q(0:NA-NB), contains the quotient of A/B.
!    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
!    In any case, Q(0:NA) should be enough.
!
!    Output, integer NR, the degree of R.
!    If the divisor polynomial is zero, NR is returned as -1.
!
!    Output, integer R(0:NB-1), contains the remainder of A/B.
!    If B has full degree, R should be dimensioned R(0:NB-1).
!    Otherwise, R will actually require less space.
!
  implicit none
!
  integer na
  integer nb
!
  integer a(0:na)
  integer a2(0:na)
  integer b(0:nb)
  integer i
  integer na2
  integer nb2
  integer nq
  integer nr
  integer q(0:*)
  integer r(0:*)
!
  call ipoly_degree ( na, a, na2 )

  call ipoly_degree ( nb, b, nb2 )

  if ( b(nb2) == 0 ) then
    nq = -1
    nr = -1
    return
  end if

  a2(0:na2) = a(0:na2)

  nq = na2 - nb2
  nr = nb2 - 1

  do i = nq, 0, -1
    q(i) = a2(i+nb2) / b(nb2)
    a2(i+nb2) = 0
    a2(i:i+nb2-1) = a2(i:i+nb2-1) - q(i) * b(0:nb2-1)
  end do

  r(0:nr) = a2(0:nr)

  return
end
subroutine ipoly_mul ( na, a, nb, b, c )
!
!*******************************************************************************
!
!! IPOLY_MUL computes the product of two integer polynomials A and B.
!
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    29 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, integer A(0:NA), the coefficients of the first polynomial factor.
!
!    Input, integer NB, the dimension of B.
!
!    Input, integer B(0:NB), the coefficients of the second polynomial factor.
!
!    Output, integer C(0:NA+NB), the coefficients of A * B.
!
  implicit none
!
  integer na
  integer nb
!
  integer a(0:na)
  integer b(0:nb)
  integer c(0:na+nb)
  integer d(0:na+nb)
  integer i
  integer j
!
  d(0:na+nb) = 0

  do i = 0, na
    d(i:i+nb) = d(i:i+nb) + a(i) * b(0:nb)
  end do

  c(0:na+nb) = d(0:na+nb)

  return
end
subroutine ipoly_print ( n, a, title )
!
!*******************************************************************************
!
!! IPOLY_PRINT prints out a polynomial.
!
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, integer A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none
!
  integer n
!
  integer a(0:n)
  integer i
  integer mag
  integer n2
  character plus_minus
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call ipoly_degree ( n, a, n2 )

  if ( a(n2) < 0 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( n2 >= 2 ) then
    write ( *, '( '' p(x) = '', a1, i6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( '' p(x) = '', a1, i6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( '' p(x) = '', a1, i6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0E+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0 ) then

      if ( i >= 2 ) then
        write ( *, ' ( ''        '', a1, i6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''        '', a1, i6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''        '', a1, i6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine ivec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )
!
!*******************************************************************************
!
!! IVEC_BACKTRACK supervises a backtrack search for an integer vector.
!
!
!  Discussion:
!
!    The routine tries to construct an integer vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
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
!    24 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of positions to be filled in the vector.
!
!    Input/output, integer X(N), the partial or complete candidate vector.
!
!    Input/output, integer INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer K, if INDX=2, the current vector index being considered.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, integer STACK(MAXSTACK), a list of all current candidates for
!    all positions 1 through K.
!
!    Input, integer MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer NCAN(N), lists the current number of candidates for
!    positions 1 through K.
!
  implicit none
!
  integer n
  integer maxstack
!
  integer indx
  integer k
  integer ncan(n)
  integer nstack
  integer stack(maxstack)
  integer x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( ncan(k) > 0 ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine ivec_binary_add ( n, ivec, jvec, kvec )
!
!*******************************************************************************
!
!! IVEC_BINARY_ADD adds two binary vectors.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the length of the vectors.
!
!    Input, integer IVEC(N), JVEC(N), the vectors to be added.
!
!    Output, integer KVEC(N), the sum of the two input vectors.
!
  implicit none
!
  integer n
!
  integer i
  integer ivec(n)
  integer jvec(n)
  integer kvec(n)
  logical overflow
!
  overflow = .false.
!
  kvec(1:n) = ivec(1:n) + jvec(1:n)

  do i = n, 1, -1
    if ( kvec(i) == 2 ) then
      kvec(i) = 0
      if ( i > 1 ) then
        kvec(i-1) = kvec(i-1) + 1
      else
        overflow = .true.
      end if
    end if
  end do

  return
end
subroutine ivec_binary_complement2 ( n, ivec )
!
!*******************************************************************************
!
!! IVEC_BINARY_COMPLEMENT2 computes the two's complement of a binary vector.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the length of the vectors.
!
!    Input/output, integer IVEC(N), the vector to be two's complemented.
!
  implicit none
!
  integer n
!
  integer ivec(n)
  integer jvec(n)
  integer kvec(n)
!
  ivec(1:n) = 1 - ivec(1:n)

  jvec(1:n-1) = 0
  jvec(n) = 1

  call ivec_binary_add ( n, ivec, jvec, ivec )

  ivec(1:n) = kvec(1:n)

  return
end
subroutine ivec_binary_print ( n, a, title )
!
!*******************************************************************************
!
!! IVEC_BINARY_PRINT prints a binary integer vector, with an optional title.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  integer ihi
  integer ilo
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
  end if

  do ilo = 1, n, 80
    ihi = min ( ilo + 79, n )
    write ( *, '(80i1)' ) a(ilo:ihi)
  end do

  return
end
subroutine ivec_binary_to_i ( n, ivec, i )
!
!*******************************************************************************
!
!! IVEC_BINARY_TO_I makes an integer from a vector binary representation.
!
!
!  Example:
!
!         IVEC        I
!    --------------  --
!    0,0,...0,0,0,1   1
!    0,0,...0,0,1,0   2
!    0,0,...0,0,1,1   3
!    0,0,...0,1,0,0   4
!    0,0,...1,0,0,1   9
!    1,1,...0,1,1,1  -9
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, an integer to be represented.
!
!    Input, integer N, the dimension of the vector, which should be 32.
!
!    Output, integer IVEC(N), the binary representation.
!
  implicit none
!
  integer n
!
  integer i
  integer i_sign
  integer ivec(n)
  integer j
  integer jvec(n)
!
  jvec(1:n) = ivec(1:n)

  i_sign = 1

  if ( jvec(1) == 1 ) then
    i_sign = -1
    call ivec_binary_complement2 ( n, jvec )
  end if

  i = 0
  do j = 2, n
    i = 2 * i + jvec(j)
  end do

  i = i_sign * i

  return
end
subroutine ivec_binary_xor ( n, ivec, jvec, kvec )
!
!*******************************************************************************
!
!! IVEC_BINARY_XOR computes the exclusive OR of two binary vectors.
!
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the length of the vectors.
!
!    Input, integer IVEC(N), JVEC(N), the binary vectors to be XOR'ed.
!
!    Input, integer KVEC(N), the exclusive OR of the two vectors.
!
  implicit none
!
  integer n
!
  integer ivec(n)
  integer jvec(n)
  integer kvec(n)
!
  kvec(1:n) = mod ( ivec(1:n) + jvec(1:n), 2 )

  return
end
function ivec_descends ( n, x )
!
!*******************************************************************************
!
!! IVEC_DESCENDS determines if an integer vector is decreasing.
!
!
!  Example:
!
!    X = ( 9, 7, 7, 3, 2, 1, -8 )
!
!    IVEC_DESCEND = TRUE
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the array.
!
!    Input, integer X(N), the array to be examined.
!
!    Output, logical IVEC_DESCEND, is TRUE if the entries of X descend.
!
  implicit none
!
  integer n
!
  integer i
  logical ivec_descends
  integer x(n)
!
  ivec_descends = .false.

  do i = 1, n-1
    if ( x(i) < x(i+1) ) then
      return
    end if
  end do

  ivec_descends = .true.

  return
end
subroutine ivec_frac ( n, a, k, iafrac )
!
!*******************************************************************************
!
!! IVEC_FRAC searches for the K-th smallest element in an N-vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, integer A(N), array to search.  On output,
!    the elements of A have been somewhat rearranged.
!
!    Input, integer K, the fractile to be sought.  If K = 1, the
!    minimum entry is sought.  If K = N, the maximum is sought.
!    Other values of K search for the entry which is K-th in size.
!    K must be at least 1, and no greater than N.
!
!    Output, integer IAFRAC, the value of the K-th fractile of A.
!
  implicit none
!
  integer n
!
  integer i
  integer a(n)
  integer iafrac
  integer iryt
  integer iw
  integer ix
  integer j
  integer k
  integer left
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IVEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IVEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IVEC_FRAC  - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal K > N, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( left >= iryt ) then
      iafrac = a(k)
      exit
    end if

    ix = a(k)
    i = left
    j = iryt

    do

      if ( i > j ) then

        if ( j < k ) then
          left = i
        end if

        if ( k < i ) then
          iryt = j
        end if

        exit

      end if
!
!  Find I so that IX <= A(I).
!
      do while ( a(i) < ix )
        i = i + 1
      end do
!
!  Find J so that A(J) <= IX.
!
      do while ( ix < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call i_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine ivec_identity ( n, a )
!
!*******************************************************************************
!
!! IVEC_IDENTITY sets an integer vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n
    a(i) = i
  end do

  return
end
function ivec_imax_last ( n, x )
!
!*******************************************************************************
!
!! IVEC_MAX_LAST returns the index of the last maximal integer vector element.
!
!
!  Example:
!
!    X = ( 5, 1, 2, 5, 0, 5, 3 )
!
!    IVEC_IMAX_LAST = 6
!
!  Modified:
!
!    17 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the array.
!
!    Input, integer X(N), the array to be examined.
!
!    Output, integer IVEC_IMAX_LAST, the index of the last element of
!    X of maximal value.
!
  implicit none
!
  integer n
!
  integer i
  integer ivec_imax_last
  integer max_last
  integer x(n)
!
  ivec_imax_last = 0

  do i = 1, n
    if ( i == 1 ) then
      ivec_imax_last = 1
      max_last = x(1)
    else if ( x(i) >= max_last ) then
      ivec_imax_last = i
      max_last = x(i)
    end if
  end do

  return
end
function ivec_index ( n, a, aval )
!
!*******************************************************************************
!
!! IVEC_INDEX returns the location of the first occurrence of a given value.
!
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, integer A(N), the vector to be searched.
!
!    Input, integer AVAL, the value to be indexed.
!
!    Output, integer IVEC_INDEX, the first location in A which has the
!    value AVAL, or 0 if no such index exists.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer aval
  integer i
  integer ivec_index
!
  do i = 1, n
    if ( a(i) == aval ) then
      ivec_index = i
      return
    end if
  end do

  ivec_index = 0

  return
end
function ivec_pairwise_prime ( n, a )
!
!*******************************************************************************
!
!! IVEC_PAIRWISE_PRIME checks whether a vector of integers is pairwise prime.
!
!
!  Discussion:
!
!    Two positive integers I and J are pairwise prime if they have no common
!    factor greater than 1.
!
!  Modified:
!
!    05 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values to check.
!
!    Input, integer A(N), the vector of integers.
!
!    Output, logical IVEC_PAIRWISE_PRIME, is TRUE if the vector of integers
!    is pairwise prime.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  integer i_gcd
  logical ivec_pairwise_prime
  integer j
!
  ivec_pairwise_prime = .false.

  do i = 1, n
    do j = i+1, n
      if ( i_gcd ( a(i), a(j) ) /= 1 ) then
        return
      end if
    end do
  end do

  ivec_pairwise_prime = .true.

  return
end
subroutine ivec_print ( n, a, title )
!
!*******************************************************************************
!
!! IVEC_PRINT prints an integer vector, with an optional title.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,i10)' ) i, a(i)
  end do

  return
end
subroutine ivec_random ( n, a, alo, ahi )
!
!*******************************************************************************
!
!! IVEC_RANDOM returns a random integer vector in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, integer A(N), the vector of randomly chosen integers.
!
!    Input, integer ALO, AHI, the range allowed for the entries.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer ahi
  integer alo
  integer i
!
  do i = 1, n

    call i_random ( alo, ahi, a(i) )

  end do

  return
end
subroutine ivec_reverse ( n, a )
!
!*******************************************************************************
!
!! IVEC_REVERSE reverses the elements of an integer vector.
!
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
!
  do i = 1, n/2
    call i_swap ( a(i), a(n+1-i) )
  end do

  return
end
subroutine ivec_sort_heap_index_d ( n, a, indx )
!
!*******************************************************************************
!
!! IVEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an integer vector.
!
!
!  Discussion:
!
!    The sorting is not actually carried out.  Rather an index array is
!    created which defines the sorting.  This array may be used to sort
!    or index the array, or to sort or index related arrays keyed on the
!    original array.
!
!    Once the index array is computed, the sorting can be carried out
!    "implicitly:
!
!      A(INDX(I)), I = 1 to N is sorted,
!
!    after which A(I), I = 1 to N is sorted.
!
!  Modified:
!
!    07 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input, integer A(N), an array to be index-sorted.
!
!    Output, integer INDX(N), contains the sort index.  The
!    I-th element of the sorted array is A(INDX(I)).
!
  implicit none
!
  integer n
!
  integer a(n)
  integer aval
  integer i
  integer indx(n)
  integer indxt
  integer ir
  integer j
  integer l
!
  call ivec_identity ( n, indx )

  l = n / 2 + 1
  ir = n

  do

    if ( l > 1 ) then

      l = l - 1
      indxt = indx(l)
      aval = a(indxt)

    else

      indxt = indx(ir)
      aval = a(indxt)
      indx(ir) = indx(1)
      ir = ir - 1

      if ( ir == 1 ) then
        indx(1) = indxt
        return
      end if

    end if

    i = l
    j = l + l

    do while ( j <= ir )

      if ( j < ir ) then
        if ( a(indx(j)) > a(indx(j+1)) ) then
          j = j + 1
        end if
      end if

      if ( aval > a(indx(j)) ) then
        indx(i) = indx(j)
        i = j
        j = j + j
      else
        j = ir + 1
      end if

    end do

    indx(i) = indxt

  end do

  return
end
subroutine jfrac_to_rfrac ( m, r, s, p, q )
!
!*******************************************************************************
!
!! JFRAC_TO_RFRAC converts a J-fraction into a rational polynomial fraction.
!
!
!  Discussion:
!
!    The routine accepts a J-fraction:
!
!        R(1) / ( X + S(1)
!      + R(2) / ( X + S(2)
!      + R(3) / ...
!      + R(M) / ( X + S(M) )... ))
!
!    and returns the equivalent rational polynomial fraction:
!
!      P(1) + P(2) * X + ... + P(M) * X**(M-1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    17 April 2000
!
!  Parameters:
!
!    Input, integer M, defines the number of P, R, and S
!    coefficients, and is one less than the number of Q
!    coefficients.
!
!    Input, real R(M), S(M), the coefficients defining the J-fraction.
!
!    Output, real P(M), Q(M+1), the coefficients defining the rational
!    polynomial fraction.  The algorithm used normalizes the coefficients
!    so that Q(M+1) = 1.0.
!
  implicit none
!
  integer m
!
  real a(m,m)
  real b(m,m)
  integer i
  integer k
  real p(m)
  real q(m+1)
  real r(m)
  real s(m)
!
  a(1,1) = r(1)
  b(1,1) = s(1)

  if ( m > 1 ) then

    do k = 2, m
      a(k,k) = r(1)
      b(k,k) = b(k-1,k-1) + s(k)
    end do

    a(1,2) = r(1) * s(2)
    b(1,2) = r(2) + s(1) * s(2)

    do k = 3, m
      a(1,k) = s(k) * a(1,k-1) + r(k) * a(1,k-2)
      a(k-1,k) = a(k-2,k-1) + s(k) * r(1)
      b(1,k) = s(k) * b(1,k-1) + r(k) * b(1,k-2)
      b(k-1,k) = b(k-2,k-1) + s(k) * b(k-1,k-1) + r(k)
    end do

    do k = 4, m
      do i = 2, k-2
        a(i,k) = a(i-1,k-1) + s(k)*a(i,k-1) + r(k) * a(i,k-2)
        b(i,k) = b(i-1,k-1) + s(k)*b(i,k-1) + r(k) * b(i,k-2)
      end do
    end do

  end if

  p(1:m) = a(1:m,m)

  q(1:m) = b(1:m,m)
  q(m+1) = 1.0E+00

  return
end
subroutine josephus ( n, m, k, x )
!
!*******************************************************************************
!
!! JOSEPHUS returns the position X of the K-th man to be executed.
!
!
!  Discussion:
!
!    The classic Josephus problem concerns a circle of 41 men.
!    Every third man is killed and removed from the circle.  Counting
!    and executing continues until all are dead.  Where was the last
!    survivor sitting?
!
!    Note that the first person killed was sitting in the third position.
!    Moreover, when we get down to 2 people, and we need to count the
!    "third" one, we just do the obvious thing, which is to keep counting
!    around the circle until our count is completed.
!
!    The process may be regarded as generating a permutation of
!    the integers from 1 to N.  The permutation would be the execution
!    list, that is, the list of the executed men, by position number.
!
!  Reference:
!
!    W W Rouse Ball,
!    Mathematical Recreations and Essays,
!    Macmillan, 1962, pages 32-36.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, pages 158-159.
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 3, Sorting and Searching,
!    Addison Wesley, 1968, pages 18-19.
!
!  Modified:
!
!    03 April 2000
!
!  Parameters:
!
!    Input, integer N, the number of men.
!    N must be positive.
!
!    Input, integer M, the counting index.
!    M must not be zero.  Ordinarily, M is positive, and no greater than N.
!
!    Input, integer K, the index of the executed man of interest.
!    K must be between 1 and N.
!
!    Output, integer X, the position of the K-th man.
!    X will be between 1 and N.
!
  implicit none
!
  integer i_modp
  integer k
  integer m
  integer m2
  integer n
  integer x
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  N <= 0.'
    stop
  end if

  if ( m == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  M = 0.'
    stop
  end if

  if ( k <= 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
    write ( *, '(a)' ) '  J <= 0 or K > N.'
    stop
  end if
!
!  In case M is bigger than N, or negative, get the
!  equivalent positive value between 1 and N.
!  You can skip this operation if 1 <= M <= N.
!
  m2 = i_modp ( m, n )

  x = k * m2

  do while ( x > n )
    x = ( m2 * ( x - n ) - 1 ) / ( m2 - 1 )
  end do

  return
end
subroutine ksub_next ( n, k, iarray, more )
!
!*******************************************************************************
!
!! KSUB_NEXT generates the K subsets of an N set, one at a time.
!
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
!    Input, integer N, the size of the set from which subsets are drawn.
!
!    Input, integer K, the desired size of the subsets.  K must
!    be between 0 and N.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th element of the
!    subset.  Thus IARRAY(I) will be an integer between 1 and N.
!    Note that the routine will return the values in IARRAY
!    in sorted order: 1 <= IARRAY(1) < IARRAY(2) < ... < IARRAY(K) <= N
!
!    Input/output, logical MORE.  Set MORE = .FALSE. before first call
!    for a new sequence of subsets.  It then is set and remains
!    .TRUE. as long as the subset computed on this call is not the
!    final one.  When the final subset is computed, MORE is set to
!    .FALSE. as a signal that the computation is done.
!
  implicit none
!
  integer k
!
  integer iarray(k)
  integer j
  integer, save :: m = 0
  integer, save :: m2 = 0
  logical more
  integer n
!
  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT - Fatal error!'
    write ( *, '(a,i6)' ) 'N = ', n
    write ( *, '(a,i6)' ) 'K = ', k
    write ( *, '(a)' ) 'but 0 <= K <= N is required!'
    stop
  end if

  if ( .not. more ) then
    m2 = 0
    m = k
  else
    if ( m2 < n-m ) then
      m = 0
    end if
    m = m + 1
    m2 = iarray(k+1-m)
  end if

  do j = 1, m
    iarray(k+j-m) = m2 + j
  end do

  more = iarray(1) /= (n-k+1)

  return
end
subroutine ksub_next2 ( n, k, iarray, in, iout )
!
!*******************************************************************************
!
!! KSUB_NEXT2 computes the next K subset of an N set.
!
!
!  Discussion:
!
!    This routine uses the revolving door method.  It has no "memory".
!    It simply calculates the successor of the input set,
!    and will start from the beginning after the last set.
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
!    29 March 2001
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!    N must be positive.
!
!    Input, integer K, the size of the desired subset.  K must be
!    between 0 and N.
!
!    Input/output, integer IARRAY(K).  On input, the user must
!    supply a subset of size K in IARRAY.  That is, IARRAY must
!    contain K unique numbers, in order, between 1 and N.  On
!    output, IARRAY(I) is the I-th element of the output subset.
!    The output array is also in sorted order.
!
!    Output, integer IN, the element of the output subset which
!    was not in the input set.  Each new subset differs from the
!    last one by adding one element and deleting another.
!
!    Output, integer IOUT, the element of the input subset which
!    is not in the output subset.
!
  implicit none
!
  integer k
!
  integer iarray(k)
  integer in
  integer iout
  integer j
  integer m
  integer n
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a)' ) '  but 0 < N is required!'
    stop
  end if

  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  j = 0

  do

    if ( j > 0 .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( j > k ) then
        iarray(k) = k
        in = k
        iout = n
        return
      end if

      if ( iarray(j) /= j ) then

        iout = iarray(j)
        in = iout - 1
        iarray(j) = in

        if ( j /= 1 ) then
          in = j - 1
          iarray(j-1) = in
        end if

        return

      end if

    end if

    j = j + 1
    m = n

    if ( j < k ) then
      m = iarray(j+1) - 1
    end if

    if ( m /= iarray(j) ) then
      exit
    end if

  end do

  in = iarray(j) + 1
  iarray(j) = in
  iout = in - 1

  if ( j /= 1 ) then
    iarray(j-1) = iout
    iout = j - 1
  end if

  return
end
subroutine ksub_next3 ( n, k, iarray, more, in, iout )
!
!*******************************************************************************
!
!! KSUB_NEXT3 computes the K subsets of an N set.
!
!
!  Discussion:
!
!    The routine uses the revolving door method.
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
!    30 March 2001
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!    N must be positive.
!
!    Input, integer K, the size of the desired subsets.  K must be
!    between 0 and N.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th element of the
!    output subset.  The elements of IARRAY are sorted.
!
!    Input/output, logical MORE.  On first call, set MORE = .FALSE.
!    to signal the beginning.  MORE will be set to .TRUE., and on
!    each call, the routine will return another K-subset.
!    Finally, when the last subset has been returned,
!    MORE will be set .FALSE. and you may stop calling.
!
!    Output, integer IN, the element of the output subset which
!    was not in the input set.  Each new subset differs from the
!    last one by adding one element and deleting another.  IN is not
!    defined the first time that the routine returns, and is
!    set to zero.
!
!    Output, integer IOUT, the element of the input subset which is
!    not in the output subset.  IOUT is not defined the first time
!    the routine returns, and is set to zero.
!
  implicit none
!
  integer k
!
  integer i
  integer iarray(k)
  integer in
  integer iout
  integer j
  integer m
  logical more
  integer n
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a)' ) '  but 0 < N is required!'
    stop
  end if

  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  if ( .not. more ) then
    in = 0
    iout = 0
    call ivec_identity ( k, iarray )
    more = ( k /= n )
    return
  end if

  j = 0

  do

    if ( j > 0 .or. mod ( k, 2 ) == 0 ) then

      j = j + 1

      if ( iarray(j) /= j ) then

        iout = iarray(j)
        in = iout - 1
        iarray(j) = in

        if ( j /= 1 ) then
          in = j - 1
          iarray(j-1) = in
        end if

        if ( k /= 1 ) then
          more = ( iarray(k-1) == k-1 )
        end if

        more = ( .not. more ) .or. ( iarray(k) /= n )

        return

      end if

    end if

    j = j + 1
    m = n

    if ( j < k ) then
      m = iarray(j+1) - 1
    end if

    if ( m /= iarray(j) ) then
      exit
    end if

  end do

  in = iarray(j) + 1
  iarray(j) = in
  iout = in - 1

  if ( j /= 1 ) then
    iarray(j-1) = iout
    iout = j - 1
  end if

  if ( k /= 1 ) then
    more = ( iarray(k-1) == k-1 )
  end if

  more = ( .not. more ) .or. ( iarray(k) /= n )

  return
end
subroutine ksub_next4 ( done, iarray, k, n )
!
!*******************************************************************************
!
!! KSUB_NEXT4 generates the subsets of size K of a set of N elements.
!
!
!  Discussion:
!
!    The subsets are generated one at a time.
!
!    The routine should be used by setting DONE to .TRUE., and then calling
!    repeatedly.  Each call returns with DONE equal to .FALSE., the array
!    IARRAY contains information defining a new subset.  When DONE returns
!    equal to .TRUE., there are no more subsets.
!
!    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
!
!  Modified:
!
!    31 July 2000
!
!  Parameters:
!
!    Input/output, logical DONE.
!
!    On the first call, DONE is an input quantity with a value
!    of .TRUE. which tells the program to initialize data and
!    return the first subset.
!
!    On return, DONE is an output quantity that is .TRUE. as long as
!    the routine is returning another subset, and .FALSE. when
!    there are no more.
!
!    Input/output, integer IARRAY(K), contains information about
!    the subsets.  On the first call with DONE = .TRUE., the input contents
!    of IARRAY don't matter.  Thereafter, the input value of IARRAY
!    should be the same as the output value of the previous call.
!    In other words, leave the array alone!
!
!    On output, as long as DONE is returned .FALSE., IARRAY contains
!    information defining a subset of K elements of a set of N elements.
!    In other words, IARRAY will contain K distinct numbers (in order)
!    between 1 and N.
!
!    Input, integer K, the size of the desired subset.  K must be
!    between 0 and N.
!
!    Input, integer N, the size of the entire set.
!
  implicit none
!
  integer k
!
  logical done
  integer i
  integer iarray(k)
  integer j
  integer jsave
  integer n
!
  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_NEXT4 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if
!
!  First call:
!
  if ( done ) then
    call ivec_identity ( k, iarray )

    if ( n > 0 ) then
      done = .false.
    else
      done = .true.
    end if
!
!  Next call.
!
  else

    if ( iarray(1) < n-k+1 ) then

      done = .false.

      jsave = k

      do j = 1, k-1

        if ( iarray(j+1) > iarray(j)+1 ) then
          jsave = j
          exit
        end if

      end do

      call ivec_identity ( jsave-1, iarray )
      iarray(jsave) = iarray(jsave) + 1

    else

      done = .true.

    end if

  end if

  return
end
subroutine ksub_random ( n, k, iarray )
!
!*******************************************************************************
!
!! KSUB_RANDOM selects a random subset of size K from a set of size N.
!
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
!    01 December 2000
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!
!    Input, integer K, number of elements in desired subsets.  K must
!    be between 0 and N.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th element of the
!    output set.  The elements of IARRAY are in order.
!
  implicit none
!
  integer k
!
  integer i
  integer iarray(k)
  integer ids
  integer ihi
  integer ip
  integer ir
  integer is
  integer ix
  integer l
  integer ll
  integer m
  integer m0
  integer n
  real r
!
  if ( k < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K is required!'
    stop
  else if ( k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  do i = 1, k
    iarray(i) = ( ( i - 1 ) * n ) / k
  end do

  do i = 1, k

    do

      call i_random ( 1, n, ix )

      l = 1 + ( ix * k - 1 ) / n

      if ( ix > iarray(l) ) then
        exit
      end if

    end do

    iarray(l) = iarray(l) + 1

  end do

  ip = 0
  is = k

  do i = 1, k

    m = iarray(i)
    iarray(i) = 0

    if ( m /= ( ( i - 1 ) * n ) / k ) then
      ip = ip + 1
      iarray(ip) = m
    end if

  end do

  ihi = ip

  do i = 1, ihi
    ip = ihi + 1 - i
    l = 1 + ( iarray(ip) * k - 1 ) / n
    ids = iarray(ip) - ( ( l - 1 ) * n ) / k
    iarray(ip) = 0
    iarray(is) = l
    is = is - ids
  end do

  do ll = 1, k

    l = k + 1 - ll

    if ( iarray(l) /= 0 ) then
      ir = l
      m0 = 1 + ( ( iarray(l) - 1 ) * n ) / k
      m = ( iarray(l) * n ) / k - m0 + 1
    end if

    call i_random ( m0, m0 + m - 1, ix )

    i = l + 1

    do while ( i <= ir )

      if ( ix < iarray(i) ) then
        exit
      end if

      ix = ix + 1
      iarray(i-1) = iarray(i)
      i = i + 1

    end do

    iarray(i-1) = ix
    m = m - 1

  end do

  return
end
subroutine ksub_random2 ( n, k, iarray )
!
!*******************************************************************************
!
!! KSUB_RANDOM2 selects a random subset of size K from a set of size N.
!
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
!    25 July 2000
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!
!    Input, integer K, number of elements in desired subsets.  K must
!    be between 0 and N.
!
!    Output, integer IARRAY(K).  IARRAY(I) is the I-th element of the
!    output set.  The elements of IARRAY are in order.
!
  implicit none
!
  integer k
!
  integer i
  integer iarray(k)
  integer ic1
  integer ic2
  integer k0
  integer n
  real r
  real, parameter :: rhi = 1.0E+00
  real, parameter :: rlo = 0.0E+00
!
  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  if ( k == 0 ) then
    return
  end if

  ic1 = k
  ic2 = n
  k0 = 0
  i = 0

  do

    i = i + 1

    call r_random ( rlo, rhi, r )

    if ( real ( ic2 ) * r <= real ( ic1 ) ) then

      ic1 = ic1 - 1
      k0 = k0 + 1
      iarray(k0) = i

      if ( ic1 <= 0 ) then
        exit
      end if

    end if

    ic2 = ic2 - 1

  end do

  return
end
subroutine ksub_random3 ( n, k, iarray )
!
!*******************************************************************************
!
!! KSUB_RANDOM3 selects a random subset of size K from a set of size N.
!
!
!  Discussion:
!
!    This routine uses Floyd's algorithm.
!
!  Modified:
!
!    01 December 2000
!
!  Parameters:
!
!    Input, integer N, the size of the set from which subsets are drawn.
!
!    Input, integer K, number of elements in desired subsets.  K must
!    be between 0 and N.
!
!    Output, integer IARRAY(N).  I is an element of the subset
!    if IARRAY(I) = 1, and I is not an element if IARRAY(I)=0.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer j
  integer k
  real r
!
  if ( k < 0 .or. k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KSUB_RANDOM3 - Fatal error!'
    write ( *, '(a,i6)' ) '  N = ', n
    write ( *, '(a,i6)' ) '  K = ', k
    write ( *, '(a)' ) '  but 0 <= K <= N is required!'
    stop
  end if

  iarray(1:n) = 0

  if ( k == 0 ) then
    return
  end if

  do i = n-k+1, n

    call i_random ( 1, i, j )

    if ( iarray(j) == 0 ) then
      iarray(j) = 1
    else
      iarray(i) = 1
    end if

  end do

  return
end
subroutine matrix_product_opt ( n, rank, cost, order )
!
!*******************************************************************************
!
!! MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
!
!
!  Discussion:
!
!    The cost of multiplying an LxM matrix by an M by N matrix is
!    assessed as L*M*N.
!
!    Any particular order of multiplying a set of N matrices is equivalent
!    to parenthesizing an expression of N objects.
!
!    The actual number of ways of parenthesizing an expression
!    of N objects is C(N), the N-th Catalan number.
!
!  Reference:
!
!    Robert Sedgewick,
!    Algorithms,
!    Addison-Wesley, 1984, pages 486-489.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of matrices to be multiplied.
!
!    Input, integer RANK(N+1), the rank information for the matrices.
!    Matrix I has RANK(I) rows and RANK(I+1) columns.
!
!    Output, integer COST, the cost of the multiplication if the optimal
!    order is used.
!
!    Output, integer ORDER(N-1), indicates the order in which the N-1
!    multiplications are to be carried out.  ORDER(1) is the first
!    multiplication to do, and so on.
!
  implicit none
!
  integer, parameter :: stack_max = 100
  integer n
!
  integer best(n,n)
  integer cost
  integer cost2(n,n)
  integer cost3
  integer i
  integer i_inf
  integer i1
  integer i2
  integer i3
  integer j
  integer k
  integer order(n-1)
  integer rank(n+1)
  integer stack(stack_max)
  integer stack_num
  integer step
!
!  Initialize the cost matrix.
!
  do i = 1, n

    cost2(i,1:i) = 0
    cost2(i,i+1:n) = huge ( 1 )

  end do
!
!  Initialize the BEST matrix.
!
  best(1:n,1:n) = 0
!
!  Compute the cost and best matrices.
!
  do j = 1, n-1
    do i = 1, n-j
      do k = i+1, i+j
        cost3 = cost2(i,k-1) + cost2(k,i+j) + rank(i) * rank(k) * rank(i+j+1)
        if ( cost3 < cost2(i,i+j) ) then
          cost2(i,i+j) = cost3
          best(i,i+j) = k
        end if
      end do
    end do
  end do
!
!  Pick off the optimal cost.
!
  cost = cost2(1,n)
!
!  Backtrack to determine the optimal order.
!
  stack_num = 0

  i1 = 1
  i2 = n

  if ( i1+1 < i2 ) then
    stack_num = stack_num + 1
    stack(stack_num) = i1
    stack_num = stack_num + 1
    stack(stack_num) = i2
  end if

  step = n - 1
!
!  Take an item off the stack.
!
  do while ( stack_num > 0 )

    i3 = stack(stack_num)
    stack_num = stack_num - 1
    i1 = stack(stack_num)
    stack_num = stack_num - 1

    i2 = best(i1,i3)

    order(step) = i2 - 1
    step = step - 1
!
!  The left chunk is matrices (I1...I2-1)
!
    if ( i1 == i2-1 ) then

    else if ( i1+1 == i2-1 ) then
      order(step) = i2 - 2
      step = step - 1
    else
      stack_num = stack_num + 1
      stack(stack_num) = i1
      stack_num = stack_num + 1
      stack(stack_num) = i2 - 1
    end if
!
!  The right chunk is matrices (I2...I3)
!
    if ( i2 == i3 ) then

    else if ( i2+1 == i3 ) then
      order(step) = i2
      step = step - 1
    else
      stack_num = stack_num + 1
      stack(stack_num) = i2
      stack_num = stack_num + 1
      stack(stack_num) = i3
    end if

  end do

  return
end
subroutine moebius_matrix ( n, a, mu )
!
!*******************************************************************************
!
!! MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
!
!
!  Discussion:
!
!    This routine can be called with MATRIX and MU being the same matrix.
!    The routine will correctly compute the Moebius matrix, which
!    will, in this case, overwrite the input matrix.
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
!    Input, integer N, number of elements in the partially ordered set.
!
!    Input, integer A(N,N).  A(I,J) = 1 if I is covered by J,
!    0 otherwise.
!
!    Output, integer MU(N,N), the Moebius matrix as computed by the routine.
!
  implicit none
!
  integer n
!
  integer a(n,n)
  integer i
  integer iwork(n,2)
  integer j
  integer mu(n,n)
!
!  Compute a reordering of the elements of the partially ordered matrix.
!
  call triang ( n, a, iwork(1,1) )
!
!  Copy the matrix.
!
  mu(1:n,1:n) = a(1:n,1:n)
!
!  Apply the reordering to MU.
!
  call imat_perm2 ( n, n, mu, iwork(1,1), iwork(1,1) )
!
!  Negate the (strict) upper triangular elements of MU.
!
  do i = 1, n-1
    mu(i,i+1:n) = - mu(i,i+1:n)
  end do
!
!  Compute the inverse of MU.
!
  call imat_inverse ( n, mu, mu )
!
!  All nonzero elements are reset to 1.
!
  do i = 1, n
    do j = i, n
      if ( mu(i,j) /= 0 ) then
        mu(i,j) = 1
      end if
    end do
  end do
!
!  Invert the matrix again.
!
  call imat_inverse ( n, mu, mu )
!
!  Compute the inverse permutation.
!
  do i = 1, n
    iwork(iwork(i,1),2) = i
  end do
!
!  Unpermute the rows and columns of MU.
!
  call imat_perm2 ( n, n, mu, iwork(1,2), iwork(1,2) )

  return
end
subroutine morse_thue ( i, s )
!
!*******************************************************************************
!
!! MORSE_THUE generates a Morse_Thue number.
!
!
!  Definition:
!
!    The Morse_Thue sequence can be defined in a number of ways.
!
!    A) Start with the string containing the single letter '0'; then
!       repeatedly apply the replacement rules '0' -> '01' and
!       '1' -> '10' to the letters of the string.  The Morse_Thue sequence
!       is the resulting letter sequence.
!
!    B) Starting with the string containing the single letter '0',
!       repeatedly append the binary complement of the string to itself.
!       Thus, '0' becomes '0' + '1' or '01', then '01' becomes
!       '01' + '10', which becomes '0110' + '1001', and so on.
!
!    C) Starting with I = 0, the I-th Morse-Thue number is determined
!       by taking the binary representation of I, adding the digits,
!       and computing the remainder modulo 2.
!
!  Example:
!
!     I  binary   S
!    --  ------  --
!     0       0   0
!     1       1   1
!     2      10   1
!     3      11   0
!     4     100   1
!     5     101   0
!     6     110   0
!     7     111   1
!     8    1000   1
!
!  Modified:
!
!    17 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the index of the Morse-Thue number.
!    Normally, I is 0 or greater, but any value is allowed.
!
!    Output, integer S, the Morse-Thue number of index I.
!
  implicit none
!
  integer, parameter :: nbits = 32
!
  integer b(nbits)
  integer i
  integer i_copy
  integer s
!
  i_copy = abs ( i )
!
!  Expand I into binary form.
!
  call i_to_ivec_binary ( i_copy, nbits, b )
!
!  Sum the 1's in the binary representation.
!
  s = sum ( b(1:nbits) )
!
!  Take the value modulo 2.
!
  s = mod ( s, 2 )

  return
end
subroutine multinomial_coef1 ( nfactor, factor, ncomb )
!
!*******************************************************************************
!
!! MULTINOMIAL_COEF1 computes a multinomial coefficient.
!
!
!  Definition:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!  Formula:
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!  Method:
!
!    The log of the gamma function is used, to avoid overflow.
!
!  Modified:
!
!    16 June 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NFACTOR, the number of factors.
!
!    Input, integer FACTOR(NFACTOR), contains the factors.
!    0 <= FACTOR(I)
!
!    Output, integer NCOMB, the value of the multinomial coefficient.
!
  implicit none
!
  integer nfactor
!
  real arg
  real fack
  real facn
  integer factor(nfactor)
  real gamma_log
  integer i
  integer n
  integer ncomb
!
!  Each factor must be nonnegative.
!
  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTINOMIAL_COEF1 - Fatal error'
      write ( *, '(a,i6,a,i6)' ) '  Factor ', I, ' = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      stop
    end if

  end do
!
!  The factors sum to N.
!
  n = sum ( factor(1:nfactor) )

  arg = real ( n + 1 )
  facn = gamma_log ( arg )

  do i = 1, nfactor

    arg = real ( factor(i) + 1 )
    fack = gamma_log ( arg )
    facn = facn - fack

  end do

  ncomb = nint ( exp ( facn ) )

  return
end
subroutine multinomial_coef2 ( nfactor, factor, ncomb )
!
!*******************************************************************************
!
!! MULTINOMIAL_COEF2 computes a multinomial coefficient.
!
!
!  Definition:
!
!    The multinomial coefficient is a generalization of the binomial
!    coefficient.  It may be interpreted as the number of combinations of
!    N objects, where FACTOR(1) objects are indistinguishable of type 1,
!    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
!    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
!
!  Formula:
!
!    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
!
!  Method:
!
!    A direct method is used, which should be exact.  However, there
!    is a possibility of intermediate overflow of the result.
!
!  Modified:
!
!    01 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NFACTOR, the number of factors.
!
!    Input, integer FACTOR(NFACTOR), contains the factors.
!    0 <= FACTOR(I)
!
!    Output, integer NCOMB, the value of the multinomial coefficient.
!
  implicit none
!
  integer nfactor
!
  integer factor(nfactor)
  integer i
  integer j
  integer k
  integer ncomb
!
!  Each factor must be nonnegative.
!
  do i = 1, nfactor

    if ( factor(i) < 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MULTINOMIAL_COEF2 - Fatal error'
      write ( *, '(a,i6,a,i6)' ) '  Factor ', I, ' = ', factor(i)
      write ( *, '(a)' ) '  But this value must be nonnegative.'
      stop
    end if

  end do

  ncomb = 1
  k = 0

  do i = 1, nfactor

    do j = 1, factor(i)
      k = k + 1
      ncomb = ( ncomb * k ) / j
    end do

  end do

  return
end
subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, isorce, isink, &
  icut, node_flow )
!
!*******************************************************************************
!
!! NETWORK_FLOW_MAX finds the maximal flow and a minimal cut in a network.
!
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
!    31 July 2000
!
!  Parameters:
!
!    Input, integer NNODE, the number of nodes.
!
!    Input, integer NEDGE, the number of edges.
!
!    Input/output, integer IENDPT(2,NEDGE), the edges of the network,
!    defined as pairs of nodes.  Each edge should be listed TWICE,
!    the second time in reverse order.  On output, the edges have
!    been reordered, and so the columns of IENDPT have been rearranged.
!
!    Input/output, integer ICPFLO(2,NEDGE).  Capacities and flows.
!    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
!    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
!    been rearranged to match the reordering of IENDPT.
!
!    Input, integer ISORCE, the designated source node.
!
!    Input, integer ISINK, the designated sink node.
!
!    Output, integer ICUT(NNODE).  ICUT(I) = 1 if node I is in the
!    minimal cut set, otherwise 0.
!
!    Output, integer NODE_FLOW(NNODE).  NODE_FLOW(I) is the value of the flow
!    through node I.
!
  implicit none
!
  integer nedge
  integer nnode
!
  integer i
  integer iarray(nnode)
  integer icpflo(2,nedge)
  integer icut(nnode)
  integer idel
  integer ien1
  integer ien2
  integer iendpt(2,nedge)
  integer indx
  integer ip
  integer iparm
  integer iq
  integer ir
  integer iread
  integer irite
  integer is
  integer isink
  integer isorce
  integer isort
  integer it
  integer iwork(nnode,2)
  integer j
  integer kz
  integer lst
  integer m
  integer node_flow(nnode)
!
  iarray(1:nnode) = 0
  idel = 0

  do i = 1, nedge

    icpflo(2,i) = 0
    ip = iendpt(1,i)

    if ( ip == isorce ) then
      idel = idel + icpflo(1,i)
    end if

    iarray(ip) = iarray(ip) + 1

  end do

  node_flow(isorce) = idel
  is = 1

  do i = 1, nnode
    it = iarray(i)
    iarray(i) = is
    iwork(i,1) = is
    is = is + it
  end do

  isort = 0
  ien1 = 0
  ien2 = 0

10 continue

  indx = 0

50 continue

  do

    call sort_heap_external ( nedge, indx, ien1, ien2, is )

    if ( indx < 0 ) then

      is = iendpt(1,ien1) - iendpt(1,ien2)

      if ( is == 0 ) then
        is = iendpt(2,ien1) - iendpt(2,ien2)
      end if

    else if ( indx > 0 ) then

      do ir = 1, 2
        call i_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
        call i_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
      end do

    else

      if ( isort > 0 ) then
        return
      end if

      do i = 1, nedge
        iq = iendpt(2,i)
        iendpt(1,i) = iwork(iq,1)
        iwork(iq,1) = iwork(iq,1) + 1
      end do

      go to 100

    end if

  end do

80 continue

  iendpt(1,iendpt(1,ien1)) = ien2
  iendpt(1,iendpt(1,ien2)) = ien1

  do ir = 1, 2
    call i_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
    call i_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
  end do

  if ( indx < 0 ) then
    go to 270
  end if

  if ( indx == 0 ) then
    go to 170
  end if

  go to 50

100   continue

  indx = 0

  do i = 1, nnode

    if ( i /= isorce ) then
      node_flow(i) = 0
    end if

    iwork(i,2) = nedge + 1

    if ( i < nnode ) then
      iwork(i,2) = iarray(i+1)
    end if

    icut(i) = 0

  end do

  iread = 0
  irite = 1
  iwork(1,1) = isorce
  icut(isorce) = - 1

120   continue

  iread = iread + 1

  if ( iread <= irite ) then

    ip = iwork(iread,1)
    lst = iwork(ip,2) - 1
    i = iarray(ip) - 1

    do

      i = i + 1

      if ( i > lst ) then
        go to 120
      end if

      iq = iendpt(2,i)
      idel = icpflo(1,i) - icpflo(2,i)

      if ( icut(iq) == 0 .and. idel /= 0 ) then

        if ( iq /= isink ) then
          irite = irite + 1
          iwork(irite,1) = iq
        end if

        icut(iq) = - 1

      end if

    end do

  end if

  if ( icut(isink) == 0 ) then

    icut(1:nnode) = - icut(1:nnode)

    do i = 1, nedge
      ip = iendpt(2,iendpt(1,i))
      if ( icpflo(2,i) < 0 ) then
        node_flow(ip) = node_flow(ip) - icpflo(2,i)
      end if
      iendpt(1,i) = ip
    end do

    node_flow(isorce) = node_flow(isink)
    isort = 1
    go to 10

  end if

  icut(isink) = 1

160   continue

  iread = iread - 1

  if ( iread == 0 ) then
    go to 180
  end if

  ip = iwork(iread,1)
  ien1 = iarray(ip) - 1
  ien2 = iwork(ip,2) - 1

170   continue

  if ( ien1 /= ien2 ) then

    iq = iendpt(2,ien2)

    if ( icut(iq) <= 0 .or. icpflo(1,ien2) == icpflo(2,ien2) ) then
      ien2 = ien2 - 1
      go to 170
    end if

    iendpt(2,ien2) = - iq
    icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
    icpflo(2,ien2) = 0
    ien1 = ien1 + 1

    if ( ien1 < ien2 ) then
      go to 80
    end if

  end if

  if ( ien1 >= iarray(ip) ) then
    icut(ip) = ien1
  end if

  go to 160

180   continue

  kz = 0

  do ir = 1, irite
    if ( icut(iwork(ir,1)) > 0 ) then
      kz = kz + 1
      iwork(kz,1) = iwork(ir,1)
    end if
  end do

  indx = - 1
  m = 1

200   continue

  ip = iwork(m,1)

  if ( node_flow(ip) > 0 ) then
    go to 250
  end if

210   continue

  m = m + 1

  if ( m <= kz ) then
    go to 200
  end if

  iparm = 0

220   continue

  m = m - 1

  if ( m == 1 ) then

    do i = 1, nedge

      iq = - iendpt(2,i)

      if ( iq >= 0 ) then

        iendpt(2,i) = iq
        j = iendpt(1,i)
        icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

        idel = icpflo(2,i) - icpflo(2,j)
        icpflo(2,i) = idel
        icpflo(2,j) = - idel

      end if

    end do

    go to 100

  end if

  ip = iwork(m,1)

  if ( node_flow(ip) < 0 ) then
    go to 220
  end if

  if ( node_flow(ip) == 0 ) then

    lst = nedge + 1

    if ( ip < nnode ) then
      lst = iarray(ip+1)
    end if

    i = iwork(ip,2)
    iwork(ip,2) = lst

    do

      if ( i == lst ) then
        go to 220
      end if

      j = iendpt(1,i)
      idel = icpflo(2,j)
      icpflo(2,j) = 0
      icpflo(1,j) = icpflo(1,j) - idel
      icpflo(2,i) = icpflo(2,i) - idel
      i = i + 1

    end do

  end if

  if ( iarray(ip) > icut(ip) ) then
    go to 300
  end if

250   continue

  i = icut(ip) + 1

260  continue

  do

    i = i - 1

    if ( i < iarray(ip) ) then
      go to 290
    end if

    iq = - iendpt(2,i)

    if ( node_flow(iq) >= 0 ) then
      exit
    end if

  end do

  idel = icpflo(1,i) - icpflo(2,i)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,i) = icpflo(2,i) + idel
  node_flow(ip) = node_flow(ip) - idel
  node_flow(iq) = node_flow(iq) + idel
  iparm = 1
  ien1 = iendpt(1,i)
  ien2 = iwork(iq,2) - 1

  if ( ien1 < ien2 ) then
    go to 80
  end if

  if ( ien1 /= ien2 ) then
    go to 280
  end if

270   continue

  iwork(iq,2) = ien2

280   continue

  if ( node_flow(ip) > 0 ) then
    go to 260
  end if

  if ( icpflo(1,i) == icpflo(2,i) ) then
    i = i - 1
  end if

290   continue

  icut(ip) = i

  if ( iparm /= 0 ) then
    go to 210
  end if

300   continue

  i = iwork(ip,2)

310   continue

  j = iendpt(1,i)
  idel = icpflo(2,j)

  if ( node_flow(ip) < idel ) then
    idel = node_flow(ip)
  end if

  icpflo(2,j) = icpflo(2,j) - idel
  node_flow(ip) = node_flow(ip) - idel
  iq = iendpt(2,i)
  node_flow(iq) = node_flow(iq) + idel
  i = i + 1

  if ( node_flow(ip) > 0 ) then
    go to 310
  end if

  node_flow(ip) = - 1
  go to 220

end
subroutine nim_sum ( i, j, k )
!
!*******************************************************************************
!
!! NIM_SUM computes the Nim sum of two integers.
!
!
!  Discussion:
!
!    If K is the Nim sum of I and J, then each bit of K is the exclusive
!    OR of the corresponding bits of I and J.
!
!  Example:
!
!     I     J     K       I_2        J_2         K_2
!   ----  ----  ----  ----------  ----------  ----------
!      0     0     0           0           0           0
!      1     0     1           1           0           1
!      1     1     0           1           1           0
!      2     7     5          10         111         101
!     11    28    23        1011       11100       10111
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, J, the integers to be Nim-summed.
!
!    Input, integer K, the Nim sum of I and J.
!
  implicit none
!
  integer, parameter :: nbits = 32
!
  integer i
  integer ivec(nbits)
  integer j
  integer jvec(nbits)
  integer k
  integer kvec(nbits)
!
  call i_to_ivec_binary ( i, nbits, ivec )
  call i_to_ivec_binary ( j, nbits, jvec )
  call ivec_binary_xor ( nbits, ivec, jvec, kvec )
  call ivec_binary_to_i ( nbits, kvec, k )

  return
end
subroutine pell_basic ( d, x0, y0 )
!
!*******************************************************************************
!
!! PELL_BASIC returns the fundamental solution for Pell's basic equation.
!
!
!  Discussion:
!
!    Pell's equation has the form:
!
!      X**2 - D * Y**2 = 1
!
!    where D is a given non-square integer, and X and Y may be assumed
!    to be positive integers.
!
!  Example:
!
!     D   X0   Y0
!
!     2    3    2
!     3    2    1
!     5    9    4
!     6    5    2
!     7    8    3
!     8    3    1
!    10   19    6
!    11   10    3
!    12    7    2
!    13  649  180
!    14   15    4
!    15    4    1
!    17   33    8
!    18   17    4
!    19  170   39
!    20    9    2
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999, pages 294-307
!
!  Modified:
!
!    19 May 2000
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer D, the coefficient in Pell's equation.  D should be
!    positive, and not a perfect square.
!
!    Output, integer X0, Y0, the fundamental or 0'th solution.
!    If X0 = Y0 = 0, then the calculation was canceled because of an error.
!    Both X0 and Y0 will be nonnegative.
!
  implicit none
!
  integer, parameter :: max_term = 100
!
  integer b(0:max_term)
  integer d
  integer i
  integer n_term
  integer p
  integer pm1
  integer pm2
  integer q
  integer qm1
  integer qm2
  integer r
  integer x0
  integer y0
!
!  If these values are returned, an error has occurred.
!
  x0 = 0
  y0 = 0
!
!  Check D.
!
  if ( d <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
    write ( *, '(a)' ) '  Pell coefficient D <= 0.'
    return
  end if

  call i_sqrt ( d, q, r )

  if ( r == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
    write ( *, '(a)' ) '  Pell coefficient is a perfect square.'
    return
  end if
!
!  Find the continued fraction representation of sqrt ( D ).
!
  call i_sqrt_cf ( d, max_term, n_term, b )
!
!  If necessary, go for two periods.
!
  if ( mod ( n_term, 2 ) == 1 ) then

    do i = n_term + 1, 2*n_term
      b(i) = b(i-n_term)
    end do

    n_term = 2 * n_term

  end if
!
!  Evaluate the continued fraction using the forward recursion algorithm.
!
  pm2 = 0
  pm1 = 1
  qm2 = 1
  qm1 = 0

  do i = 0, n_term-1
    p = b(i) * pm1 + pm2
    q = b(i) * qm1 + qm2
    pm2 = pm1
    pm1 = p
    qm2 = qm1
    qm1 = q
  end do
!
!  Get the fundamental solution.
!
  x0 = p
  y0 = q

  return
end
subroutine pell_next ( d, x0, y0, xn, yn, xnp1, ynp1 )
!
!*******************************************************************************
!
!! PELL_NEXT returns the next solution of Pell's equation.
!
!
!  Discussion:
!
!    Pell's equation has the form:
!
!      X**2 - D * Y**2 = 1
!
!    where D is a given non-square integer, and X and Y may be assumed
!    to be positive integers.
!
!    To compute X0, Y0, call PELL_BASIC.
!    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
!    To compute further solutions, call again with X0, Y0 and the previous
!    solution.
!
!  Example:
!
!    ------INPUT--------  --OUTPUT--
!
!    D  X0  Y0   XN   YN  XNP1  YNP1
!
!    2   3   2    3    2    17    12
!    2   3   2   17   12    99    70
!    2   3   2   99   70   577   408
!    2   3   2  577  408  3363  2378
!
!  Reference:
!
!    Mark Herkommer,
!    Number Theory, A Programmer's Guide,
!    McGraw Hill, 1999, pages 294-307
!
!  Modified:
!
!    16 May 2000
!
!  Author:
!
!   John Burkardt
!
!  Parameters:
!
!    Input, integer D, the coefficient in Pell's equation.
!
!    Input, integer X0, Y0, the fundamental or 0'th solution.
!
!    Input, integer XN, YN, the N-th solution.
!
!    Output, integer XNP1, YNP1, the N+1-th solution.
!
  implicit none
!
  integer d
  integer x0
  integer xn
  integer xnp1
  integer y0
  integer yn
  integer ynp1
!
  xnp1 = x0 * xn + d * y0 * yn
  ynp1 = x0 * yn +     y0 * xn

  return
end
subroutine pent_enum ( n, p )
!
!*******************************************************************************
!
!! PENT_ENUM computes the N-th pentagonal number.
!
!
!  Definition:
!
!    The pentagonal number P(N) counts the number of dots in a figure of
!    N nested pentagons.  The pentagonal numbers are defined for both
!    positive and negative N.
!
!  First values:
!
!    N   P
!
!   -5  40
!   -4  26
!   -3  15
!   -2   7
!   -1   2
!    0   0
!    1   1
!    2   5
!    3  12
!    4  22
!    5  35
!
!  Formula:
!
!    P(N) = ( N * ( 3 * N - 1 ) ) / 2
!
!  Modified:
!
!    22 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the index of the pentagonal number desired.
!
!    Output, integer P, the value of the N-th pentagonal number.
!
  implicit none
!
  integer n
  integer p
!
  p = ( n * ( 3 * n - 1 ) ) / 2

  return
end
subroutine perm_ascend ( n, a, length, sub )
!
!*******************************************************************************
!
!! PERM_ASCEND computes the longest ascending subsequence of permutation.
!
!
!  Discussion:
!
!    Although this routine is intended to be applied to a permutation,
!    it will work just as well for an arbitrary vector.
!
!  Modified:
!
!    07 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the permutation.
!
!    Input, integer A(N), the permutation to be examined.
!
!    Output, integer LENGTH, the length of the longest increasing subsequence.
!
!    Output, integer SUB(N), contains in entries 1 through LENGTH
!    a longest increasing subsequence of A.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer a1
  integer a2
  integer i
  integer j
  integer k
  integer length
  integer sub(n)
  integer top(n)
  integer top_prev(n)
!
  top(1:n) = 0
  top_prev(1:n) = 0
  sub(1:n) = 0

  if ( n <= 0 ) then
    length = 0
    return
  end if

  length = 0

  do i = 1, n

    k = 0

    do j = 1, length
      if ( a(i) <= a(top(j)) ) then
        k = j
        exit
      end if
    end do

    if ( k == 0 ) then
      length = length + 1
      k = length
    end if

    top(k) = i

    if ( k > 1 ) then
      top_prev(i) = top(k-1)
    else
      top_prev(i) = 0
    end if

  end do

  j = top(length)
  sub(length) = a(j)

  do i = length-1, 1, -1
    j = top_prev(j)
    sub(i) = a(j)
  end do

  return
end
subroutine perm_break_count ( n, p, break_count )
!
!*******************************************************************************
!
!! PERM_BREAK_COUNT counts the number of "breaks" in a permutation.
!
!
!  Discussion:
!
!    We begin with a permutation of order N.  We prepend an element
!    labeled "0" and append an element labeled "N+1".  There are now
!    N+1 pairs of neighbors.  A "break" is a pair of neighbors whose
!    value differs by more than 1.  
!
!    The identity permutation has a break count of 0.  The maximum
!    break count is N+1.
!
!  Modified:
!
!    25 October 2001
!
!  Parameters:
!
!    Input, integer N, the order of the permutation.
!
!    Input, integer P(N), a permutation, in standard index form.
!
!    Output, integer BREAK_COUNT, the number of breaks in the permutation.
!
  implicit none
!
  integer n
!
  integer break_count
  integer i
  integer ierror
  integer p(n)
!
  break_count = 0
!
!  Make sure the permutation is a legal one.
!  (This is not an efficient way to do so!)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_BREAK_COUNT - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  if ( p(1) /= 1 ) then
    break_count = break_count + 1
  end if

  do i = 1, n-1
    if ( abs ( p(i+1) - p(i) ) /= 1 ) then
      break_count = break_count + 1
    end if
  end do

  if ( p(n) /= n ) then
    break_count = break_count + 1
  end if

  return
end
subroutine perm_canon_to_cycle ( n, p1, p2 )
!
!*******************************************************************************
!
!! PERM_CANON_TO_CYCLE converts a permutation from canonical to cycle form.
!
!
!  Example:
!
!    Input:
!
!      4 5 2 1 6 3
!
!    Output:
!
!      -4 5 -2 -1 6 3,
!      indicating the cycle structure
!      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, page 176.
!
!  Modified:
!
!    30 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer P1(N), the permutation, in canonical form.
!
!    Output, integer P2(N), the permutation, in cycle form.
!
  implicit none
!
  integer n
!
  integer i
  integer p1(n)
  integer p2(n)
  integer pmin
!
  p2(1:n) = p1(1:n)

  pmin = p2(1) + 1

  do i = 1, n

    if ( p2(i) < pmin ) then
      pmin = p2(i)
      p2(i) = - p2(i)
    end if

  end do

  return
end
subroutine perm_check ( n, p, ierror )
!
!*******************************************************************************
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the permutation, in standard index form.
!
!    Output, integer IERROR, error flag.
!    0, the array does represent a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none
!
  integer n
!
  integer ierror
  integer ifind
  integer iseek
  integer p(n)
!
  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end
subroutine perm_cycle ( n, p, isgn, ncycle, iopt )
!
!*******************************************************************************
!
!! PERM_CYCLE analyzes a permutation.
!
!
!  Discussion:
!
!    The routine will count cycles, find the sign of a permutation,
!    and tag a permutation.
!
!  Example:
!
!    Input:
!
!      N = 9
!      IOPT = 1
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NCYCLE = 3
!      ISGN = +1
!      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
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
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N).  On input, P describes a
!    permutation, in the sense that entry I is to be moved to P(I).
!    If IOPT = 0, then P will not be changed by this routine.
!    If IOPT = 1, then on output, P will be "tagged".  That is,
!    one element of every cycle in P will be negated.  In this way,
!    a user can traverse a cycle by starting at any entry I1 of P
!    which is negative, moving to I2 = ABS(P(I1)), then to
!    P(I2), and so on, until returning to I1.
!
!    Output, integer ISGN, the "sign" of the permutation, which is
!    +1 if the permutation is even, -1 if odd.  Every permutation
!    may be produced by a certain number of pairwise switches.
!    If the number of switches is even, the permutation itself is
!    called even.
!
!    Output, integer NCYCLE, the number of cycles in the permutation.
!
!    Input, integer IOPT, requests tagging.
!    0, the permutation will not be tagged.
!    1, the permutation will be tagged.
!
  implicit none
!
  integer n
!
  integer i
  integer i1
  integer i2
  integer ierror
  integer iopt
  integer is
  integer isgn
  integer ncycle
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_CYCLE - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  is = 1
  ncycle = n

  do i = 1, n

    i1 = p(i)

    do while ( i1 > i )
      ncycle = ncycle - 1
      i2 = p(i1)
      p(i1) = - i2
      i1 = i2
    end do

    if ( iopt /= 0 ) then
      is = - sign ( 1, p(i) )
    end if

    p(i) = sign ( p(i), is )

  end do

  isgn = 1 - 2 * mod ( n - ncycle, 2 )

  return
end
subroutine perm_cycle_to_canon ( n, p1, p2 )
!
!*******************************************************************************
!
!! PERM_CYCLE_TO_CANON converts a permutation from cycle to canonical form.
!
!
!  Example:
!
!    Input:
!
!      -6 3 1 -5, 4 -2,
!      indicating the cycle structure
!      ( 6, 3, 1 ) ( 5, 4 ) ( 2 )
!
!    Output:
!
!      4 5 2 1 6 3
!
!  Discussion:
!
!    The procedure is to "rotate" the elements of each cycle so that
!    the smallest element is first:
!
!      ( 1, 6, 3 ) ( 4, 5 ) ( 2 )
!
!    and then to sort the cycles in decreasing order of their first
!    (and lowest) element:
!
!      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
!
!    and then to drop the parentheses:
!
!      4 5 2 1 6 3
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, pages 176.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer P1(N), the permutation, in cycle form.
!
!    Output, integer P2(N), the permutation, in canonical form.
!
  implicit none
!
  integer n
!
  integer hi(n)
  integer i
  integer indx(n)
  integer j
  integer k
  integer lo(n)
  integer ncycle
  integer next
  integer nhi
  integer nlo
  integer nmin
  integer p1(n)
  integer p2(n)
  integer pmin(n)
  integer ptemp(n)
!
  p2(1:n) = p1(1:n)
!
!  Work on the next cycle.
!
  nlo = 1
  ncycle = 0

  do while ( nlo <= n )
!
!  Identify NHI, the last index in this cycle.
!
    ncycle = ncycle + 1

    nhi = nlo

    do while ( nhi < n )
      if ( p2(nhi+1) < 0 ) then
        exit
      end if
      nhi = nhi + 1
    end do
!
!  Identify the smallest value in this cycle.
!
    p2(nlo) = - p2(nlo)
    pmin(ncycle) = p2(nlo)
    nmin = nlo

    do i = nlo+1, nhi
      if ( p2(i) < pmin(ncycle) ) then
        pmin(ncycle) = p2(i)
        nmin = i
      end if
    end do
!
!  Rotate the cycle so A_MIN occurs first.
!
    ptemp(nlo+nhi+1-nmin:nhi) = p2(nlo:nmin-1)
    ptemp(nlo:nlo+nhi-nmin) = p2(nmin:nhi)

    lo(ncycle) = nlo
    hi(ncycle) = nhi
!
!  Prepare to operate on the next cycle.
!
    nlo = nhi + 1

  end do
!
!  Compute a sorting index for the cycle minima.
!
  call ivec_sort_heap_index_d ( ncycle, pmin, indx )
!
!  Copy the cycles out of the temporary array in sorted order.
!
  j = 0
  do i = 1, ncycle
    next = indx(i)
    nlo = lo(next)
    nhi = hi(next)
    do k = nlo, nhi
      j = j + 1
      p2(j) = ptemp(k)
    end do
  end do

  return
end
subroutine perm_cycle_to_index ( n, p1, p2 )
!
!*******************************************************************************
!
!! PERM_CYCLE_TO_INDEX converts a permutation from cycle to standard index form.
!
!
!  Example:
!
!    Input:
!
!      N = 9
!      P1 = -1, 2, 3, 9, -4, 6, 8, -5, 7
!
!    Output:
!
!      P2 = 2, 3, 9, 6, 7, 8, 5, 4, 1
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
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input, integer P1(N), the permutation, in cycle form.
!
!    Output, integer P2(N), the permutation, in standard index form.
!
  implicit none
!
  integer n
!
  integer j
  integer k1
  integer k2
  integer k3
  integer p1(n)
  integer p2(n)
!
  do j = 1, n

    k1 = p1(j)

    if ( k1 < 0 ) then
      k1 = - k1
      k3 = k1
    end if

    if ( j + 1 <= n ) then
      k2 = p1(j+1)
      if ( k2 < 0 ) then
        k2 = k3
      end if
    else
      k2 = k3
    end if

    p2(k1) = k2

  end do

  return
end
subroutine perm_distance ( n, a, b, k )
!
!*******************************************************************************
!
!! PERM_DISTANCE computes the Ulam metric distance of two permutations.
!
!
!  Discussion:
!
!    If we let N be the order of the permutations A and B, and L(P) be
!    the length of the longest ascending subsequence of a permutation P,
!    then the Ulam metric distance between A and B is
!
!      N - L ( A * inverse ( B ) ).
!
!  Modified:
!
!    08 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the permutation.
!
!    Input, integer A(N), B(N), the permutations to be examined.
!
!    Output, integer K, the Ulam metric distance between A and B.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer b(n)
  integer binv(n)
  integer c(n)
  integer k
  integer length
  integer sub(n)
!
  binv(1:n) = b(1:n)

  call perm_inv ( n, binv )

  call perm_mul ( n, a, binv, c )

  call perm_ascend ( n, c, length, sub )

  k = n - length

  return
end
subroutine perm_fixed_enum ( n, m, fnm )
!
!*******************************************************************************
!
!! PERM_FIXED_ENUM enumerates the permutations of N objects with M fixed.
!
!
!  Definition:
!
!    A permutation of N objects with M fixed is a permutation in which
!    exactly M of the objects retain their original positions.  If
!    M = 0, the permutation is a "derangement".  If M = N, the
!    permutation is the identity.
!
!  Formula:
!
!    F(N,M) = ( N! / M! ) * ( 1 - 1/1! + 1/2! - 1/3! ... 1/(N-M)! )
!           = COMB(N,M) * D(N-M)
!    where D(N-M) is the number of derangements of N-M objects.
!
!  Modified:
!
!    31 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!    N should be at least 1.
!
!    Input, integer M, the number of objects that retain their
!    position.  M should be between 0 and N.
!
!    Output, integer FNM, the number of derangements of N objects
!    in which M objects retain their positions.
!
  implicit none
!
  integer derange_enum
  integer fnm
  integer icnm
  integer m
  integer n
!
  if ( n <= 0 ) then

    fnm = 1

  else if ( m < 0 ) then

    fnm = 0

  else if ( m > n ) then

    fnm = 0

  else if ( m == n ) then

    fnm = 1

  else if ( n == 1 ) then

    if ( m == 1 ) then
      fnm = 1
    else
      fnm = 0
    end if

  else

    call combin2 ( n, m, icnm )

    fnm = icnm * derange_enum ( n - m )

  end if

  return
end
subroutine perm_free ( ipart, npart, ifree, nfree )
!
!*******************************************************************************
!
!! PERM_FREE reports the number of unused items in a partial permutation.
!
!
!  Discussion:
!
!    It is assumed that the N objects being permuted are the integers
!    from 1 to N, and that IPART contains a "partial" permutation, that
!    is, the NPART entries of IPART represent the beginning of a
!    permutation of all N items.
!
!    The routine returns in IFREE the items that have not been used yet.
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IPART(NPART), the partial permutation, which should
!    contain, at most once, some of the integers between 1 and
!    NPART+NFREE.
!
!    Input, integer NPART, the number of entries in IPART.  NPART may be 0.
!
!    Output, integer IFREE(NFREE), the integers between 1 and NPART+NFREE
!    that were not used in IPART.
!
!    Input, integer NFREE, the number of integers that have not been
!    used in IPART.  This is simply N - NPART.  NFREE may be zero.
!
  implicit none
!
  integer nfree
  integer npart
!
  integer i
  integer ifree(nfree)
  integer ipart(npart)
  integer j
  integer k
  integer match
  integer n
!
  n = npart + nfree

  if ( npart < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NPART < 0.'
    stop

  else if ( npart == 0 ) then

    call ivec_identity ( n, ifree )

  else if ( nfree < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
    write ( *, '(a)' ) '  NFREE < 0.'
    stop

  else if ( nfree == 0 ) then

    return

  else

    k = 0

    do i = 1, n

      match = 0

      do j = 1, npart
        if ( ipart(j) == i ) then
          match = j
          exit
        end if
      end do

      if ( match == 0 ) then

        k = k + 1

        if ( k > nfree ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
          write ( *, '(a)' ) '  The partial permutation is illegal.'
          write ( *, '(a)' ) '  It should contain, at most once, some of'
          write ( *, '(a,i6)' ) '  the integers between 1 and ', n
          stop
        end if

        ifree(k) = i

      end if

    end do

  end if

  return
end
subroutine perm_index_to_cycle ( n, p1, p2 )
!
!*******************************************************************************
!
!! PERM_INDEX_TO_CYCLE converts a permutation from standard index to cycle form.
!
!
!  Example:
!
!    Input:
!
!      N = 9
!      P1 = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      P2 = -1, 2, 3, 9, -4, 6, 8, -5, 7
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
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input, integer P1(N), the permutation, in standard index form.
!
!    Output, integer P2(N), the permutation, in cycle form.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  integer k
  integer p1(n)
  integer p2(n)
!
  i = 0
  j = 1

  do while ( j <= n )

    if ( p1(j) < 0 ) then

      j = j + 1

    else

      k = j

      i = i + 1
      p2(i) = - k

      do while ( p1(k) /= j )
        i = i + 1
        p2(i) = p1(k)
        p1(k) = - p1(k)
        k = abs ( p1(k) )
      end do

      p1(k) = - p1(k)

    end if

  end do

  p1(1:n) = abs ( p1(1:n) )

  return
end
subroutine perm_ins ( n, p, ins )
!
!*******************************************************************************
!
!! PERM_INS computes the inversion sequence of a permutation.
!
!
!  Definition:
!
!    For a given permutation P acting on objects 1 through N, the inversion
!    sequence INS is defined as:
!
!      INS(1) = 0
!      INS(I) = number of values J < I for which P(J) > P(I).
!
!  Example:
!
!    Input:
!
!      ( 3, 5, 1, 4, 2 )
!
!    Output:
!
!      ( 0, 0, 2, 1, 3 )
!
!  Note:
!
!    The original permutation can be recovered from the inversion sequence.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input, integer P(N), the permutation, in standard index form.
!    The I-th item has been mapped to P(I).
!
!    Output, integer INS(N), the inversion sequence of the permutation.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer ins(n)
  integer j
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INS - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  ins(1:n) = 0

  do i = 1, n
    do j = 1, i-1
      if ( p(j) > p(i) ) then
        ins(i) = ins(i) + 1
      end if
    end do
  end do

  return
end
subroutine perm_inv ( n, p )
!
!*******************************************************************************
!
!! PERM_INV inverts a permutation "in place".
!
!
!  Modified:
!
!    25 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!    On output, P describes the inverse permutation
!
  implicit none
!
  integer n
!
  integer i
  integer i0
  integer i1
  integer i2
  integer ierror
  integer is
  integer p(n)
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of N = ', n
    stop
  end if

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  is = 1

  do i = 1, n

    i1 = p(i)

    do while ( i1 > i )
      i2 = p(i1)
      p(i1) = - i2
      i1 = i2
    end do

    is = - sign ( 1, p(i) )
    p(i) = sign ( p(i), is )

  end do

  do i = 1, n

    i1 = - p(i)

    if ( i1 >= 0 ) then

      i0 = i

      do

        i2 = p(i1)
        p(i1) = i0

        if ( i2 < 0 ) then
          exit
        end if

        i0 = i1
        i1 = i2

      end do

    end if

  end do

  return
end
subroutine perm_inv2 ( n, p )
!
!*****************************************************************************
!
!! PERM_INV2 inverts a permutation "in place".
!
!
!  Discussion:
!
!    The routine needs no extra vector storage in order to compute the
!    inverse of a permutation.
!
!    This feature might be useful if the permutation is large.
!
!  Modified:
!
!    25 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects in the permutation.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!    On output, the inverse permutation.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer ii
  integer j
  integer k
  integer m
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_INV2 - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if
!
  do ii = 1, n

    m = n + 1 - ii
    i = p(m)

    if ( i < 0 ) then

      p(m) = - i

    else if ( i /= m ) then

      k = m

      do

        j = p(i)
        p(i) = - k

        if ( j == m ) then
          p(m) = i
          exit
        end if

        k = i
        i = j

      end do

    end if

  end do

  return
end
subroutine perm_lex ( n, p, more )
!
!*******************************************************************************
!
!! PERM_LEX generates permutations in lexical order, one at a time.
!
!
!  Example:
!
!    N = 3
!
!    1   1 2 3
!    2   1 3 2
!    3   2 1 3
!    4   2 3 1
!    5   3 1 2
!    6   3 2 1
!
!  Reference:
!
!    Mok-Kong Shen,
!    Algorithm 202: Generation of Permutations in Lexicographical Order,
!    Communications of the ACM,
!    Volume 6, September 1963, page 517.
!
!  Modified:
!
!    24 September 1998
!
!  Parameters:
!
!    Input, integer N, the number of elements being permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!
!    The user should not alter the elements of Pbetween successive
!    calls.  The routine uses the input value of P to determine
!    the output value.
!
!    Input/output, logical MORE.
!
!    On the first call, the user should set MORE = FALSE, which signals
!    the routine to do initialization.
!
!    On return, if MORE is TRUE, then another permutation has been
!    computed and returned, while if MORE is FALSE, there are no more
!    permutations.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  integer k
  logical more
  integer p(n)
  integer u
  integer w
!
!  Initialization.
!
  if ( .not. more ) then

    call ivec_identity ( n, p )
    more = .true.

  else

    if ( n <= 1 ) then
      more = .false.
      return
    end if

    w = n

    do while ( p(w) < p(w-1) )

      if ( w == 2 ) then
        more = .false.
        return
      end if

      w = w - 1
    end do

    u = p(w-1)

    do j = n, w, -1

      if ( p(j) > u ) then

        p(w-1) = p(j)
        p(j) = u

        do k = 0, ( n - w - 1 ) / 2
          call i_swap ( p(n-k), p(w+k) )
        end do

        return

      end if

    end do

  end if

  return
end
subroutine perm_mul ( n, p1, p2, p3 )
!
!*******************************************************************************
!
!! PERM_MUL "multiplies" two permutations.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the permutations.
!
!    Input, integer P1(N), P2(N), the permutations, in standard index form.
!
!    Output, integer P3(N), the product permutation.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer p1(n)
  integer p2(n)
  integer p3(n)
!
  call perm_check ( n, p1, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_MUL - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  call perm_check ( n, p2, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_MUL - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if

  p3(1:n) = p2(p1(1:n))

  return
end
subroutine perm_next ( n, p, more, even )
!
!*******************************************************************************
!
!! PERM_NEXT computes all of the permutations of N objects, one at a time.
!
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
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
!    30 March 2001
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!
!    On the first call, the input value is unimportant.
!    On subsequent calls, the input value should be the same
!    as the output value from the previous call.  In other words, the
!    user should just leave P alone.
!
!    On output, contains the "next" permutation.
!
!    Input/output, logical MORE.
!
!    Set MORE = .FALSE. before the first call.
!
!    MORE will be reset to .TRUE. and a permutation will be returned.
!
!    Each new call produces a new permutation until
!    MORE is returned .FALSE.
!
!    Input/output, logical EVEN.
!
!    The input value of EVEN should simply be its output value from the
!    previous call; (the input value on the first call doesn't matter.)
!
!    On output, EVEN is .TRUE. if the output permutation is even, that is,
!    involves an even number of transpositions.
!
  implicit none
!
  integer n
!
  logical even
  integer i
  integer i1
  integer ia
  integer id
  integer is
  integer j
  integer l
  integer m
  logical more
  integer p(n)
!
  if ( .not. more ) then

    call ivec_identity ( n, p )
    more = .true.
    even = .true.

    if ( n == 1 ) then
      more = .false.
      return
    end if

    if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
      return
    end if

    do i = 1, n-3
      if ( p(i+1) /= p(i)+1 ) then
        return
      end if
    end do

    more = .false.

  else

    if ( n == 1 ) then
      p(1) = 0
      more = .false.
      return
    end if

    if ( even ) then

      ia = p(1)
      p(1) = p(2)
      p(2) = ia
      even = .false.

      if ( p(n) /= 1 .or. p(1) /= 2 + mod ( n, 2 ) ) then
        return
      end if

      do i = 1, n-3
        if ( p(i+1) /= p(i)+1 ) then
          return
        end if
      end do

      more = .false.
      return

    else

      more = .false.

      is = 0

      do i1 = 2, n

        ia = p(i1)
        i = i1 - 1
        id = 0

        do j = 1, i
          if ( p(j) > ia ) then
            id = id + 1
          end if
        end do

        is = id + is
        if ( id /= i * mod ( is, 2 ) ) then
          more = .true.
          exit
        end if

      end do

      if ( .not. more ) then
        p(1) = 0
        return
      end if

    end if

    m = mod ( is+1, 2 ) * ( n + 1 )

    do j = 1, i

      if ( sign ( 1, p(j)-ia ) /= sign ( 1, p(j)-m ) ) then
        m = p(j)
        l = j
      end if

    end do

    p(l) = ia
    p(i1) = m
    even = .true.

  end if

  return
end
subroutine perm_next2 ( n, p, done, iactiv, idir, invers )
!
!*******************************************************************************
!
!! PERM_NEXT2 generates all the permutations of N objects.
!
!
!  Discussion:
!
!    The routine generates the permutations one at a time.  It uses a
!    particular ordering of permutations, generating them from the first
!    (which is the identity permutation) to the N!-th.  The same ordering
!    is used by the routines PERM_RANK and PERM_UNRANK.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set to be permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!
!    Input/output, logical DONE.  The user should set the input value of
!    DONE only once, before the first call to compute the permutations.
!    The user should set DONE to .TRUE., which signals the routine
!    that it is to initialize itself.
!
!    Thereafter, the routine will set DONE to .FALSE. and will
!    compute a new permutation on each call.
!
!    However, when there are no more permutations to compute, the
!    routine will not return a new permutation, but instead will
!    return DONE with the value .TRUE..  At this point, all the
!    permutations have been computed.
!
!    Workspace, integer IACTIV(N).
!
!    On each call, as long as DONE is not returned as
!    .TRUE., P will contain a new permutation.  Since this
!    permutation is computed in part from the previous one, the user
!    should not overwrite or alter the contents of P during
!    the computation.
!
!    Workspace, integer IDIR(N).
!
!    Output, integer INVERS(N).  This array, which is computed as a
!    byproduct of the main computation, contains the inverse permutation
!    to that contained in P.  In particular, INVERS(P(I)) = I.
!
  implicit none
!
  integer n
!
  logical done
  integer i
  integer iactiv(n)
  integer idir(n)
  integer invers(n)
  integer j
  integer nactiv
  integer p(n)
!
!  An input value of .FALSE. for DONE is assumed to mean a new
!  computation is beginning.
!
  if ( done ) then

    call ivec_identity ( n, p )
    invers(1:n) = p(1:n)
    idir(1:n) = -1

    iactiv(1) = 0
    iactiv(2:n) = 1
!
!  Set the DONE flag to .FALSE., signifying there are more permutations
!  to come.  Except, of course, that we must take care of the trivial case!
!
    if ( n > 1 ) then
      done = .false.
    else
      done = .true.
    end if
!
!  Otherwise, assume we are in a continuing computation
!
  else

    nactiv = 0

    do i = 1, n
      if ( iactiv(i) /= 0 ) then
        nactiv = i
      end if
    end do

    if ( nactiv <= 0 ) then

      done = .true.

    else

      j = invers(nactiv)

      p(j) = p(j+idir(nactiv))
      p(j+idir(nactiv)) = nactiv

      invers(nactiv) = invers(nactiv) + idir(nactiv)
      invers(p(j)) = j

      if ( j+2*idir(nactiv) < 1.or.j+2*idir(nactiv) > n ) then
        idir(nactiv) = - idir(nactiv)
        iactiv(nactiv) = 0
      else if ( nactiv < p(j+2*idir(nactiv)) ) then
        idir(nactiv) = - idir(nactiv)
        iactiv(nactiv) = 0
      end if

      iactiv(nactiv+1:n) = 1

    end if

  end if

  return
end
subroutine perm_next3 ( n, p, more )
!
!*******************************************************************************
!
!! PERM_NEXT3 computes all of the permutations of N objects, one at a time.
!
!
!  Discussion:
!
!    The routine is initialized by calling with MORE = TRUE, in which case
!    it returns the identity permutation.
!
!    If the routine is called with MORE = FALSE, then the successor of the
!    input permutation is computed.
!
!    Trotter's algorithm is used.
!
!  Reference:
!
!    H Trotter,
!    PERM, Algorithm 115,
!    Communications of the Association for Computing Machinery,
!    Volume 5, 1962, pages 434-435.
!
!  Modified:
!
!    25 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input/output, integer P(N), the permutation, in standard index form.
!    If MORE is .TRUE., then P is assumed to contain the
!    "previous" permutation, and on P(I) is the value
!    of the I-th object under the next permutation.
!    Otherwise, P will be set to the "first" permutation.
!
!    Input/output, logical MORE.
!    Set MORE = .FALSE. before first calling this routine.
!    MORE will be reset to .TRUE. and a permutation will be returned.
!    Each new call produces a new permutation until MORE is returned .FALSE.
!
  implicit none
!
  integer n
!
  integer i
  integer, save :: irank = 0
  integer m2
  logical more
  integer n2
  integer, save :: nfact = 0
  integer p(n)
  integer q
  integer s
  integer t
!
  if ( .not. more ) then

    call ivec_identity ( n, p )
    more = .true.
    irank = 1

    nfact = 1
    do i = 1, n
      nfact = nfact * i
    end do

  else

    n2 = n
    m2 = irank
    s = n

    do

      q = mod ( m2, n2 )
      t = mod ( m2, 2 * n2 )

      if ( q /= 0 ) then
        exit
      end if

      if ( t == 0 ) then
        s = s - 1
      end if

      m2 = m2 / n2
      n2 = n2 - 1

    end do

    if ( q == t ) then
      s = s - q
    else
      s = s + q - n2
    end if

    call i_swap ( p(s), p(s+1) )

    irank = irank + 1

    if ( irank == nfact ) then
      more = .false.
    end if

  end if

  return
end
subroutine perm_print ( n, p, title )
!
!*******************************************************************************
!
!! PERM_PRINT prints a permutation.
!
!
!  Example:
!
!    Input:
!
!      P = 7 2 4 1 5 3 6
!
!    Printed output:
!
!      "This is the permutation:"
!
!      1 2 3 4 5 6 7
!      7 2 4 1 5 3 6
!
!  Modified:
!
!    25 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer P(N), the permutation, in standard index form.
!
!    Input, character ( len = * ) TITLE, an optional title.
!    If no title is supplied, then only the permutation is printed.
!
  implicit none
!
  integer n
!
  integer i
  integer ihi
  integer ilo
  integer, parameter :: inc = 20
  integer p(n)
  character ( len = * ) title
!
  if ( len_trim ( title ) /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(a)' ) ' '
      write ( *, '(20i4)' ) ( i, i = ilo, ihi )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  else

    do ilo = 1, n, inc
      ihi = min ( n, ilo + inc - 1 )
      write ( *, '(20i4)' ) p(ilo:ihi)
    end do

  end if

  return
end
subroutine perm_random ( n, p )
!
!*******************************************************************************
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!
!  Discussion:
!
!    The routine assumes the objects are labeled 1, 2, ... N.
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
!    12 May 2002
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Output, integer P(N), a permutation of ( 1, 2, ..., N ), in standard 
!    index form.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  integer p(n)
!
  call ivec_identity ( n, p )

  do i = 1, n
    call i_random ( i, n, j )
    call i_swap ( p(i), p(j) )
  end do

  return
end
subroutine perm_random2 ( n, p )
!
!*******************************************************************************
!
!! PERM_RANDOM2 selects a random permutation of N objects.
!
!
!  Discussion:
!
!    The input values of P are used as labels; that is, the I-th object 
!    is labeled P(I).
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
!    12 May 2002
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Input/output, integer P(N), on input, a list of labels.  On output,
!    the list has been permuted randomly.
!
  implicit none
!
  integer n
!
  integer i
  integer j
  integer p(n)
!
  do i = 1, n
    call i_random ( i, n, j )
    call i_swap ( p(i), p(j) )
  end do

  return
end
subroutine perm_random3 ( n, p )
!
!*******************************************************************************
!
!! PERM_RANDOM3 selects a random permutation of N elements.
!
!
!  Reference:
!
!    K L Hoffman and D R Shier,
!    Algorithm 564,
!    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
!    ACM Transactions on Mathematical Software,
!    Volume 6, Number 4, December 1980, pages 615-617.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    James Filliben
!    National Bureau of Standards.
!
!  Parameters:
!
!    Input, integer N, the number of elements of the array.
!
!    Output, integer P(N), a permutation, in standard index form.
!
  implicit none
!
  integer n
!
  integer i
  integer iadd
  integer j
  integer p(n)
!
  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_RANDOM3 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal input value of N  = ', n
    write ( *, '(a)' ) '  N must be at least 1!'
    stop
  end if

  if ( n == 1 ) then
    p(1) = 1
    return
  end if

  call ivec_identity ( n, p )

  do i = 1, n

    call i_random ( 1, n, iadd )

    j = i + iadd

    if ( j > n ) then
      j = j - n
    end if

    if ( i /= j ) then
      call i_swap ( p(j), p(i) )
    end if

  end do

  return
end
subroutine perm_rank ( n, p, ierror, invers, irank )
!
!*******************************************************************************
!
!! PERM_RANK computes the rank of a given permutation.
!
!
!  Discussion:
!
!    This is the same as asking for the step at which PERM_NEXT2
!    would compute the permutation.  The value of the rank will be
!    between 1 and N!.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set that
!    is permuted by P.
!
!    Input, integer P(N), a permutation, in standard index form.
!
!    Output, integer IERROR, error flag.
!    0, no error.
!    nonzero, IARRAY does not represent a permutation of the integers
!    from 1 to N.  In particular, the value IERROR is missing.
!
!    Output, integer INVERS(N), the inverse permutation of P.
!    It is computed as part of the algorithm, and may be of use
!    to the user.  INVERS(P(I)) = I for each entry I.
!
!    Output, integer IRANK, the rank of the permutation.  This
!    gives the order of the given permutation in the set of all
!    the permutations on N elements.
!
  implicit none
!
  integer n
!
  integer i
  integer icount
  integer ierror
  integer invers(n)
  integer irank
  integer irem
  integer j
  integer p(n)
!
!  Make sure the permutation is a legal one.
!  (This is not an efficient way to do so!)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_RANK - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Compute the inverse permutation
!
  invers(1:n) = p(1:n)

  call perm_inv2 ( n, invers )

  irank = 0

  do i = 1, n

    icount = 0

    do j = 1, invers(i)
      if ( p(j) < i ) then
        icount = icount + 1
      end if
    end do

    if ( mod ( irank, 2 ) == 1 ) then
      irem = icount
    else
      irem = i - 1 - icount
    end if

    irank = i * irank + irem

  end do

  irank = irank + 1

  return
end
subroutine perm_sign ( n, p, p_sign )
!
!*******************************************************************************
!
!! PERM_SIGN returns the sign of a permutation.
!
!
!  Discussion:
!
!    A permutation can always be replaced by a sequence of pairwise
!    transpositions.  A given permutation can be represented by
!    many different such transposition sequences, but the number of
!    such transpositions will always be odd or always be even.
!    If the number of transpositions is even or odd, the permutation is
!    said to be even or odd.
!
!  Example:
!
!    Input:
!
!      N = 9
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      P_SIGN = +1
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
!    01 February 2000
!
!  Parameters:
!
!    Input, integer N, the number of objects permuted.
!
!    Input, integer P(N), a permutation, in standard index form.
!
!    Output, integer P_SIGN, the "sign" of the permutation.
!    +1, the permutation is even,
!    -1, the permutation is odd.
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer ivec_index
  integer j
  integer p(n)
  integer p_sign
  integer q(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_SIGN - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Make a temporary copy of the permutation.
!
  q(1:n) = p(1:n)
!
!  Start with P_SIGN indicating an even permutation.
!  Restore each element of the permutation to its correct position,
!  updating P_SIGN as you go.
!
  p_sign = 1

  do i = 1, n-1

    j = ivec_index ( n, q, i )

    if ( j /= i ) then
      call i_swap ( q(i), q(j) )
      p_sign = - p_sign
    end if

  end do

  return
end
subroutine perm_to_equiv ( n, p, npart, jarray, iarray )
!
!*******************************************************************************
!
!! PERM_TO_EQUIV computes the partition induced by a permutation.
!
!
!  Example:
!
!    Input:
!
!      N = 9
!      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
!
!    Output:
!
!      NPART = 3
!      JARRAY = 4, 3, 2
!      IARRAY = 1, 1, 1, 2  3  2  3  2, 1
!
!  Modified:
!
!    26 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects being permuted.
!
!    Input, integer P(N), a permutation, in standard index form.
!
!    Output, integer NPART, number of subsets in the partition.
!
!    Output, integer JARRAY(N).  JARRAY(I) is the number of elements
!    in the I-th subset of the partition.
!
!    Output, integer IARRAY(N).  IARRAY(I) is the class to which
!    element I belongs.
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer ierror
  integer j
  integer jarray(n)
  integer k
  integer npart
  integer p(n)
!
  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM__TO_EQUIV - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Initialize.
!
  iarray(1:n) = 0
  jarray(1:n) = 0

  npart = 0
!
!  Search for the next item J which has not been assigned a subset/orbit.
!
  do j = 1, n

    if ( iarray(j) /= 0 ) then
      cycle
    end if
!
!  Begin a new subset/orbit.
!
    npart = npart + 1
    k = j
!
!  Add the item to the subset/orbit.
!
    do

      jarray(npart) = jarray(npart) + 1
      iarray(k) = npart
!
!  Apply the permutation.  If the permuted object isn't already in the
!  subset/orbit, add it.
!
      k = p(k)

      if ( iarray(k) /= 0 ) then
        exit
      end if

    end do

  end do

  return
end
subroutine perm_to_ytb ( n, p, lambda, iarray )
!
!*******************************************************************************
!
!! PERM_TO_YTB converts a permutation to a Young tableau.
!
!
!  Discussion:
!
!    The mapping is not invertible.  In most cases, several permutations
!    correspond to the same tableau.
!
!  Example:
!
!    N = 7
!    P = 7 2 4 1 5 3 6
!
!    YTB =
!      1 2 3 6
!      4 5
!      7
!
!    LAMBDA = 4 2 1 0 0 0 0
!
!    IARRAY = 1 1 1 2 2 1 3
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
!    11 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer to be partitioned.
!
!    Input, integer P(N), a permutation, in standard index form.
!
!    Output, integer LAMBDA(N).  LAMBDA(I) is the length of the I-th row.
!
!    Output, integer IARRAY(N).  IARRAY(I) is the row containing I.
!
  implicit none
!
  integer n
!
  logical another
  integer compare
  integer i
  integer iarray(n)
  integer lambda(n)
  integer p(n)
  integer put_index
  integer put_row
  integer put_value
!
!  Initialize.
!
  lambda(1:n) = 0
  iarray(1:n) = 0
!
!  Now insert each item of the permutation.
!
  do put_index = 1, n

    put_value = p(put_index)
    put_row = 1

    do

      another = .false.

      do compare = put_value+1, n

        if ( iarray(compare) == put_row ) then
          another = .true.
          iarray(put_value) = put_row
          iarray(compare) = 0
          put_value = compare
          put_row = put_row + 1
          exit
        end if

      end do

      if ( .not. another ) then
        exit
      end if

    end do

    iarray(put_value) = put_row
    lambda(put_row) = lambda(put_row) + 1

  end do

  return
end
subroutine perm_unrank ( n, p, irank )
!
!*******************************************************************************
!
!! PERM_UNRANK "unranks" a permutation.
!
!
!  Discussion:
!
!    That is, given a rank, it computes the corresponding permutation.
!    This is the same as asking for the permutation which PERM_NEXT2
!    would compute at the IRANK-th step.
!
!    The value of the rank should be between 1 and N!.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    25 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set.
!
!    Output, integer P(N), the permutation, in standard index form.
!
!    Input, integer IRANK, the desired rank of the permutation.  This
!    gives the order of the given permutation in the set of all
!    the permutations on N elements, using the order of PERM_NEXT2.
!
  implicit none
!
  integer n
!
  integer i
  integer icount
  integer iprev
  integer irank
  integer irem
  integer j
  integer jdir
  integer jrank
  integer nfact
  integer p(n)
!
  p(1:n) = 0

  nfact = 1

  do i = 1, n
    nfact = nfact * i
  end do

  if ( irank < 1 .or. irank > nfact ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PERM_UNRANK - Fatal error!'
    write ( *, '(a)' ) '  Illegal input value for IRANK.'
    write ( *, '(a,i6)' ) '  IRANK must be between 1 and ', nfact
    write ( *, '(a,i6)' ) '  but the input value is ', irank
    stop
  end if

  jrank = irank - 1

  do i = 1, n

    iprev = n + 1 - i
    irem = mod ( jrank, iprev )
    jrank = jrank / iprev

    if ( mod ( jrank, 2 ) == 1 ) then
      j = 0
      jdir = 1
    else
      j = n + 1
      jdir = -1
    end if

    icount = 0

    do

      j = j + jdir

      if ( p(j) == 0 ) then
        icount = icount + 1
      end if

      if ( icount > irem ) then
        exit
      end if

    end do

    p(j) = iprev

  end do

  return
end
subroutine pord_check ( n, a, ierror )
!
!*******************************************************************************
!
!! PORD_CHECK checks a matrix representing a partial ordering.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements in the set.  The number
!    of rows and columns used in IZETA.
!
!    Input, integer A(N,N).  Contains the description of the
!    partial ordering.  A(I,J) = 1 if I is less than J
!    in the partial ordering, 0 otherwise.
!
!    Output, integer IERROR, error flag.
!    0, no errors detected.
!    1, N <= 0.
!    2, A(I,J) > 0 and A(J,I) > 0 for some I and J.
!
  implicit none
!
  integer n
!
  integer a(n,n)
  integer i
  integer ierror
  integer j
!
  ierror = 0

  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
    write ( *, '(a,i6)' ) '  N must be positive, but N = ', n
    ierror = 1
    return
  end if

  do i = 1, n
    do j = i+1, n

      if ( a(i,j) > 0 ) then
        if ( a(j,i) > 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
          write ( *, '(a,i6)' ) '  For indices I = ', i
          write ( *, '(a,i6)' ) '  and J = ', j
          write ( *, '(a,i6)' ) '  A(I,J) = ', a(i,j)
          write ( *, '(a,i6)' ) '  A(J,I) = ', a(j,i)
          ierror = 2
          return
        end if
      end if

    end do
  end do

  return
end
subroutine power_mod ( a, n, m, x )
!
!*******************************************************************************
!
!! POWER_MOD computes mod ( A**N, M ).
!
!
!  Discussion:
!
!    Some programming tricks are used to speed up the computation, and to
!    allow computations in which A**N is much too large to store in a
!    real word.
!
!    First, for efficiency, the power A**N is computed by determining
!    the binary expansion of N, then computing A, A**2, A**4, and so on
!    by repeated squaring, and multiplying only those factors that
!    contribute to A**N.
!
!    Secondly, the intermediate products are immediately "mod'ed", which
!    keeps them small.
!
!    For instance, to compute mod ( A**13, 11 ), we essentially compute
!
!       13 = 1 + 4 + 8
!
!       A**13 = A * A**4 * A**8
!
!       mod ( A**13, 11 ) = mod ( A, 11 ) * mod ( A**4, 11 ) * mod ( A**8, 11 ).
!
!    Fermat's little theorem says that if P is prime, and A is not divisible
!    by P, then ( A**(P-1) - 1 ) is divisible by P.
!
!  Modified:
!
!    29 July 1999
!
!  Parameters:
!
!    Input, integer A, the base of the expression to be tested.
!    A should be nonnegative.
!
!    Input, integer N, the power to which the base is raised.
!    N should be nonnegative.
!
!    Input, integer M, the divisor against which the expression is tested.
!    M should be positive.
!
!    Output, integer X, the remainder when A**N is divided by M.
!
  implicit none
!
  integer a
  integer a_square
  integer d
  integer m
  integer n
  integer ncopy
  integer x
!
  x = - 1

  if ( a < 0 ) then
    return
  end if

  if ( m < 0 ) then
    return
  end if

  if ( n < 0 ) then
    return
  end if
!
!  A_SQUARE contains the successive squares of A.
!
  a_square = a
  ncopy = n
  x = 1

  do while ( ncopy > 0 )

    d = mod ( ncopy, 2 )

    if ( d == 1 ) then
      x = mod ( x * a_square, m )
    end if

    a_square = mod ( a_square * a_square, m )
    ncopy = ( ncopy - d ) / 2

  end do

  return
end
subroutine power_series1 ( a, c, n, alpha )
!
!*******************************************************************************
!
!! POWER_SERIES1 computes a power series for a function F(Z) = (1+H(Z))**ALPHA.
!
!
!  Discussion:
!
!    The power series for H(Z) is given.
!
!    The form of the power series are:
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
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
!    Output, real A(N), the power series coefficients for F(Z).
!
!    Input, real C(N), the power series coefficients for H(Z).
!
!    Input, integer N, the number of terms in the power series.
!
!    Input, real ALPHA, the exponent of 1+H(Z) in the definition of F(Z).
!
  implicit none
!
  integer n
!
  real a(n)
  real alpha
  real c(n)
  integer i
  integer j
  real v
!
  do j = 1, n

    v = 0.0E+00

    do i = 1, j-1
      v = v + a(i) * c(j-i) * ( alpha * ( j - i ) - i )
    end do

    a(j) = ( alpha * c(j) + v / real(j) )

  end do

  return
end
subroutine power_series2 ( a, c, n )
!
!*******************************************************************************
!
!! POWER_SERIES2 computes the power series for a function F(Z) = EXP(H(Z)) - 1.
!
!
!  Discussion:
!
!    The power series for H(Z) is given.
!
!    The power series have the form:
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
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
!    Output, real A(N), the power series coefficients for F(Z).
!
!    Input, real C(N), the power series coefficients for H(Z).
!
!    Input, integer N, the number of terms in the power series.
!
  implicit none
!
  integer n
!
  real a(n)
  real c(n)
  integer i
  integer j
  real v
!
  do j = 1, n

    v = 0.0E+00

    do i = 1, j-1
      v = v + a(i) * c(j-i) * real ( j - i )
    end do

    a(j) = c(j) + v / real ( j )

  end do

  return
end
subroutine power_series3 ( a, b, c, n )
!
!*******************************************************************************
!
!! POWER_SERIES3 computes the power series for a function F(Z) = G(H(Z)).
!
!
!  Discussion:
!
!    The power series for G and H are given.
!
!    We assume that
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
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
!    02 November 2001
!
!  Parameters:
!
!    Output, real A(N), the power series for F.
!
!    Input, real B(N), the power series for G.
!
!    Input, real C(N), the power series for H.
!
!    Input, integer N, the number of terms in the power series.
!
  implicit none
!
  integer n
!
  real a(n)
  real b(n)
  real c(n)
  integer i
  integer iq
  integer j
  integer m
  real r
  real v
  real work(n)
!
  work(1:n) = b(1) * c(1:n)
!
!  Search for IQ, the index of the first nonzero entry in C.
!
  iq = 0

  do i = 1, n

    if ( c(i) /= 0.0E+00 ) then
      iq = i
      exit
    end if

  end do

  if ( iq /= 0 ) then

    m = 1

    do

      m = m + 1

      if ( m * iq > n ) then
        exit
      end if

      if ( b(m) == 0.0E+00 ) then
        cycle
      end if

      r = b(m) * c(iq)**m
      work(m*iq) = work(m*iq) + r

      do j = 1, n-m*iq

        v = 0.0E+00
        do i = 1, j-1
          v = v + a(i) * c(j-i+iq) * real ( m * ( j - i ) - i )
        end do

        a(j) = ( real ( m ) * c(j) + v / real ( j ) ) / c(iq)

      end do

      do i = 1, n-m*iq
        work(i+m*iq) = work(i+m*iq) + a(i) * r
      end do

    end do

  end if

  a(1:n) = work(1:n)

  return
end
subroutine power_series4 ( a, b, c, n )
!
!*******************************************************************************
!
!! POWER_SERIES4 computes the power series for a function G(Z) = F ( 1/H(Z) ).
!
!
!  Discussion:
!
!    POWER_SERIES4 is given the power series for the functions F and H.
!
!    We assume that
!
!      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
!      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
!      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
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
!    Input, real A(N).  Power series for F.
!
!    Output, real B(N).  Power series for G.
!
!    Input, real C(N).  Power series for H.  For this problem, C(1)
!    may not be 0.0.
!
!    Input, integer N, the number of terms in the power series.
!
  implicit none
!
  integer n
!
  real a(n)
  real b(n)
  real c(n)
  integer i
  integer l
  integer m
  real s
  real t
  real work(n)
!
  if ( c(1) == 0.0E+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER_SERIES4 - Fatal error!'
    write ( *, '(a)' ) '  C(1) is zero.'
    stop
  end if

  t = 1.0E+00

  do i = 1, n
    t = t / c(1)
    b(i) = a(i) * t
    work(i) = c(i) * t
  end do

  do m = 2, n
    s = - work(m)
    do i = m, n
      do l = i, n
        b(l) = b(l) + s * b(l+1-m)
        work(l) = work(l) + s * work(l+1-m)
      end do
    end do
  end do

  return
end
function prime ( n )
!
!*******************************************************************************
!
!! PRIME returns any of the first MAXPRIME prime numbers.
!
!
!  Note:
!
!    MAXPRIME is 1500, and the largest prime stored is 12553.
!
!  Modified:
!
!    21 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964, pages 870-873.
!
!    Daniel Zwillinger,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
!    Input, integer N, the index of the desired prime number.
!    N = -1 returns MAXPRIME, the index of the largest prime available.
!    N = 0 is legal, returning PRIME = 1.
!    It should generally be true that 0 <= N <= MAXPRIME.
!
!    Output, integer PRIME, the N-th prime.  If N is out of range, PRIME
!    is returned as 0.
!
  implicit none
!
  integer, parameter :: maxprime = 1500
!
  integer i
  integer n
  integer, save, dimension ( maxprime ) :: npvec
  integer prime
!
  data ( npvec(i), i = 1, 100 ) / &
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /

  data ( npvec(i), i = 101, 200 ) / &
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /

  data ( npvec(i), i = 201, 300 ) / &
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /

  data ( npvec(i), i = 301, 400 ) / &
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /

  data ( npvec(i), i = 401, 500 ) / &
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /

  data ( npvec(i), i = 501, 600 ) / &
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /

  data ( npvec(i), i = 601, 700 ) / &
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /

  data ( npvec(i), i = 701, 800 ) / &
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /

  data ( npvec(i), i = 801, 900 ) / &
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /

  data ( npvec(i), i = 901, 1000 ) / &
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /

  data ( npvec(i), i = 1001, 1100 ) / &
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /

  data ( npvec(i), i = 1101, 1200 ) / &
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /

  data ( npvec(i), i = 1201, 1300 ) / &
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /

  data ( npvec(i), i = 1301, 1400 ) / &
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
    10861,10867,10883,10889,10891,10903,10909,19037,10939,10949, &
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /

  data ( npvec(i),i = 1401, 1500 ) / &
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /
!
  if ( n == -1 ) then
    prime = maxprime
  else if ( n == 0 ) then
    prime = 1
  else if ( n <= maxprime ) then
    prime = npvec(n)
  else
    prime = 0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIMES - Fatal error!'
    write ( *, '(a)' ) '  Illegal prime index N = ', n
    write ( *, '(a,i6)' ) '  but N must be between 0 and MAXPRIME =',  maxprime
    stop
  end if

  return
end
subroutine pythag_triple_next ( i, j, a, b, c )
!
!*******************************************************************************
!
!! PYTHAG_TRIPLE_NEXT computes the next Pythagorean triple.
!
!
!  Modified:
!
!    08 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J, the generators.
!    On first call, set I = J = 0.  On repeated calls, leave I and J
!    at their output values from the previous call.
!
!    Output, integer A, B, C, the next Pythagorean triple.
!    A, B, and C are positive integers which have no common factors,
!    and A**2 + B**2 = C**2.
!
  implicit none
!
  integer a
  integer b
  integer c
  integer i
  integer j
!
!  I starts at 2 and increases;
!
!  J starts out at 2 if I is odd, or 1 if I is even, increases by 2,
!    but is always less than I.
!
  if ( i == 0 .and. j == 0 ) then
    i = 2
    j = 1
  else if ( j + 2 < i ) then
    j = j + 2
  else
    i = i + 1
    j = mod ( i, 2 ) + 1
  end if

  a = i**2 - j**2
  b = 2 * i * j
  c = i**2 + j**2

  return
end
function r_agm ( a, b )
!
!*******************************************************************************
!
!! R_AGM finds the arithmetic-geometric mean of two numbers.
!
!
!  Discussion:
!
!    The AGM of (A,B) is produced by the following iteration:
!
!      (A,B) -> ( (A+B)/2, SQRT(A*B) ).
!
!    The sequence of successive values of (A,B) quickly converge to the AGM.
!
!  Modified:
!
!    02 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, B, the numbers whose AGM is desired.  A and B should
!    both be non-negative.
!
!    Output, real R_AGM, the AGM of the two numbers.
!
  implicit none
!
  real a
  real a1
  real a2
  real b
  real b1
  real b2
  real tol
  real r_agm
!
  if ( a < 0.0E+00 ) then
    r_agm = -1.0E+00
    return
  end if

  if ( b < 0.0E+00 ) then
    r_agm = -1.0E+00
    return
  end if

  if ( a == 0.0E+00 .or. b == 0.0E+00 ) then
    r_agm = 0.0E+00
    return
  end if

  if ( a == b ) then
    r_agm = a
    return
  end if

  tol = epsilon ( tol ) * ( a + b + 1.0E+00 )

  a1 = a
  b1 = b

  do

    a2 = ( a1 + b1 ) / 2.0E+00
    b2 = sqrt ( a1 * b1 )

    if ( abs ( a2 - b2 ) < tol ) then
      r_agm = ( a2 + b2 ) / 2.0E+00
      exit
    end if

    a1 = a2
    b1 = b2

  end do

  return
end
function r_pi ( )
!
!*******************************************************************************
!
!! R_PI returns the value of pi.
!
!
!  Modified:
!
!    04 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real R_PI, the value of pi.
!
  implicit none
!
  real r_pi
!
  r_pi = 3.14159265358979323846264338327950288419716939937510E+00

  return
end
subroutine r_random ( rlo, rhi, r )
!
!*******************************************************************************
!
!! R_RANDOM returns a random real in a given range.
!
!
!  Modified:
!
!    26 October 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RLO, RHI, the minimum and maximum values.
!
!    Output, real R, the randomly chosen value.
!
  implicit none
!
  real r
  real rhi
  real rlo
  real t
!
!  Pick a random number in (0,1).
!
  call random_number ( harvest = t )
!
!  Set R.
!
  r = ( 1.0E+00 - t ) * rlo + t * rhi

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine r_to_cfrac ( r, n, a, p, q )
!
!*******************************************************************************
!
!! R_TO_CFRAC converts a real value to a continued fraction.
!
!
!  Discussion:
!
!    The routine is given a real number R.  It computes a sequence of
!    continued fraction approximations to R, returning the results as
!    simple fractions of the form P(I) / Q(I).
!
!  Example:
!
!    X = 2 * PI
!    N = 7
!
!    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
!    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
!    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
!
!    (This ignores roundoff error, which will cause later terms to differ).
!
!  Reference:
!
!    Norman Richert,
!    Strang's Strange Figures,
!    American Mathematical Monthly,
!    Volume 99, Number 2, February 1992, pages 101-107.
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the real value.
!
!    Input, integer N, the number of convergents to compute.
!
!    Output, integer A(0:N), the partial quotients.
!
!    Output, integer P(-1:N), Q(-1:N), the numerators and denominators
!    of the continued fraction approximations.
!
  implicit none
!
  integer n
!
  integer a(0:n)
  integer i
  integer p(-1:n)
  integer q(-1:n)
  real r
  real r_copy
  real x(0:n)
!
  if ( r == 0.0E+00 ) then
    a(0:n) = 0
    p(-1:n) = 0
    q(-1:n) = 1
    return
  end if

  r_copy = abs ( r )

  p(-1) = 1
  q(-1) = 0

  p(0) = int ( r_copy )
  q(0) = 1
  x(0) = r_copy
  a(0) = int ( x(0) )

  do i = 1, n
    x(i) = 1.0E+00 / ( x(i-1) - real ( a(i-1) ) )
    a(i) = int ( x(i) )
    p(i) = a(i) * p(i-1) + p(i-2)
    q(i) = a(i) * q(i-1) + q(i-2)
  end do

  if ( r < 0.0E+00 ) then
    p(-1:n) = - p(-1:n)
  end if

  return
end
subroutine r_to_dec ( rval, itop, ibot )
!
!*******************************************************************************
!
!! R_TO_DEC converts a real value to a decimal fraction form.
!
!
!  Discussion:
!
!    The routine computes ITOP and IBOT, so that approximately:
!
!      RVAL = ITOP * 10 ** IBOT
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RVAL, the real number whose decimal
!    representation is desired.
!
!    Output, integer ITOP, IBOT, form the decimal
!    representation of RVAL, approximately.
!
!    ITOP is an integer, strictly between -10**NDIG and 10**NDIG.
!    IBOT is an integer exponent of 10.
!
  implicit none
!
  integer dec_digit
  integer ibot
  integer itop
  real rtop
  real rval
  real ten1
  real ten2
!
  dec_digit = 0
  call i_data ( 'GET', 'DEC_DIGIT', dec_digit )
!
!  Special cases.
!
  if ( rval == 0.0E+00 ) then
    itop = 0
    ibot = 0
    return
  end if
!
!  Factor RVAL = RTOP * 10**IBOT
!
  rtop = rval
  ibot = 0
!
!  Now normalize so that 10**(DEC_DIGIT-1) <= ABS(RTOP) < 10**(DEC_DIGIT)
!
  ten1 = 10.0E+00**( dec_digit - 1 )
  ten2 = 10.0E+00**dec_digit

  do while ( abs ( rtop ) < ten1 )
    rtop = rtop * 10.0E+00
    ibot = ibot - 1
  end do

  do while ( abs ( rtop ) >= ten2 )
    rtop = rtop / 10.0E+00
    ibot = ibot + 1
  end do
!
!  ITOP is the integer part of RTOP, rounded.
!
  itop = nint ( rtop )
!
!  Now divide out any factors of ten from ITOP.
!
  if ( itop /= 0 ) then
    do while ( 10 * ( itop / 10 ) == itop )
      itop = itop / 10
      ibot = ibot + 1
    end do
  end if

  return
end
subroutine r_to_rat ( a, iatop, iabot, ndig )
!
!*******************************************************************************
!
!! R_TO_RAT converts a real value to a rational value.
!
!
!  Discussion:
!
!    The rational value (IATOP/IABOT) is essentially computed by truncating
!    the decimal representation of the real value after a given number of
!    decimal digits.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A, the real value to be converted.
!
!    Output, integer IATOP, IABOT, the numerator and denominator
!    of the rational value that approximates A.
!
!    Input, integer NDIG, the number of decimal digits used.
!
  implicit none
!
  real a
  real factor
  integer i_gcd
  integer iabot
  integer iatop
  integer ibot
  integer ifac
  integer itemp
  integer itop
  integer jfac
  integer ndig
!
  factor = 10.0E+00**ndig

  if ( ndig > 0 ) then
    ifac = 10**ndig
    jfac = 1
  else
    ifac = 1
    jfac = 10**(-ndig)
  end if

  itop = nint ( a * factor ) * jfac
  ibot = ifac
!
!  Factor out the greatest common factor.
!
  itemp = i_gcd ( itop, ibot )

  iatop = itop / itemp
  iabot = ibot / itemp

  return
end
subroutine r_to_s_left ( rval, s )
!
!*******************************************************************************
!
!! R_TO_S_LEFT represents a real using 14 left_justified characters.
!
!
!  Modified:
!
!    22 January 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real RVAL, a real number.
!
!    Output, character ( len = * ) S, a left-justified character variable
!    containing the representation of RVAL.
!
  implicit none
!
  character ( len = 14 ) chrtmp
  integer i
  real rval
  character ( len = * ) s
!
!  We can't seem to write directly into the string because of compiler
!  quibbles.
!
  if ( real ( int ( rval ) ) == rval .and. abs ( rval ) < 1.0E+13 ) then

    write ( chrtmp, '(i14)' ) int ( rval )

  else

    write ( chrtmp, '(g14.6)' ) rval

  end if

  do i = 1, len ( chrtmp )
    if ( chrtmp(i:i) /= ' ' ) then
      s = chrtmp(i:)
      return
    end if
  end do

  s = ' '

  return
end
subroutine random_initialize ( seed )
!
!*******************************************************************************
!
!! RANDOM_INITIALIZE initializes the FORTRAN 90 random number seed.
!
!
!  Discussion:
!
!    If you don't initialize the random number generator, its behavior
!    is not specified.  If you initialize it simply by:
!
!      call random_seed
!
!    its behavior is not specified.  On the DEC ALPHA, if that's all you
!    do, the same random number sequence is returned.  In order to actually
!    try to scramble up the random number generator a bit, this routine
!    goes through the tedious process of getting the size of the random
!    number seed, making up values based on the current time, and setting
!    the random number seed.
!
!  Modified:
!
!    03 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED.
!    If SEED is zero on input, then you're asking this routine to come up
!    with a seed value, which is returned as output.
!    If SEED is nonzero on input, then you're asking this routine to
!    use the input value of SEED to initialize the random number generator.
!
  implicit none
!
  integer count
  integer count_max
  integer count_rate
  integer i
  integer seed
  integer, allocatable :: seed_vector(:)
  integer seed_size
  real t
!
!  Initialize the random number seed.
!
  call random_seed ( )
!
!  Determine the size of the random number seed.
!
  call random_seed ( size = seed_size )
!
!  Allocate a seed of the right size.
!
  allocate ( seed_vector(seed_size) )

  if ( seed /= 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) '  Initialize RANDOM_NUMBER with user SEED = ', seed

  else

    call system_clock ( count, count_rate, count_max )

    seed = count

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RANDOM_INITIALIZE'
    write ( *, '(a,i12)' ) &
      '  Initialize RANDOM_NUMBER with arbitrary SEED = ', seed

  end if
!
!  Now set the seed.
!
  seed_vector(1:seed_size) = seed

  call random_seed ( put = seed_vector(1:seed_size) )
!
!  Free up the seed space.
!
  deallocate ( seed_vector )
!
!  Call the random number routine a bunch of times.
!
  do i = 1, 100
    call random_number ( harvest = t )
  end do

  return
end
subroutine rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
!
!*******************************************************************************
!
!! RAT_ADD adds two rational values.
!
!
!  Discussion:
!
!    The routine computes
!
!      ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
!
!    while trying to avoid integer overflow.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the first value to add.
!
!    Input, integer ITOP2, IBOT2, the second value to add.
!
!    Output, integer ITOP, IBOT, the sum.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.  The addition of the two values
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
!    Output, integer ITOP, the numerator of the result.
!
!    Input, integer ITOP1, ITOP2, the numerators of the
!    two rational values to be added.
!
  implicit none
!
  integer i_gcd
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer i_max
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer jbot1
  integer jbot2
  integer jbot3
  integer jtop1
  integer jtop2
!
  i_max = huge ( i_max )

  ierror = 0

  if ( itop1 == 0 ) then
    itop = itop2
    ibot = ibot2
    return
  else if ( itop2 == 0 ) then
    itop = itop1
    ibot = ibot1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Compute the greatest common factor of the two denominators,
!  and factor it out.
!
  jbot3 = i_gcd ( jbot1, jbot2 )
  jbot1 = jbot1 / jbot3
  jbot2 = jbot2 / jbot3
!
!  The fraction may now be formally written as:
!
!    (jtop1*jbot2 + jtop2*jbot1) / (jbot1*jbot2*jbot3)
!
!  Check the tops for overflow.
!
  if ( abs ( real ( jtop1 ) * real ( jbot2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  jtop1 = jtop1 * jbot2

  if ( abs ( real ( jtop2 ) * real ( jbot1 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  jtop2 = jtop2 * jbot1

  if ( abs ( real ( jtop1 ) + real ( jtop2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational sum.'
    itop = 0
    stop
  end if

  itop = jtop1 + jtop2
!
!  Check the bottom for overflow.
!
  if ( abs ( real ( jbot1 ) * real ( jbot2 ) * real ( jbot3 ) ) > &
    real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational sum.'
    ibot = 1
    stop
  end if

  ibot = jbot1 * jbot2 * jbot3
!
!  Put the fraction in lowest terms.
!
  itemp = i_gcd ( itop, ibot )
  itop = itop / itemp
  ibot = ibot / itemp
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if

  return
end
subroutine rat_div ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
!
!*******************************************************************************
!
!! RAT_DIV divides one rational value by another.
!
!
!  Discussion:
!
!    The routine computes
!
!      ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
!
!    while avoiding integer overflow.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the numerator.
!
!    Input, integer ITOP2, IBOT2, the denominator.
!
!    Output, integer ITOP, IBOT, the result.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.  One of the quantities IBOT1, IBOT2,
!    or ITOP2 is zero, or the result of the division
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
  implicit none
!
  integer i_gcd
  integer i_max
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer jbot1
  integer jbot2
  integer jtop1
  integer jtop2
!
  ierror = 0

  i_max = huge ( i_max )

  if ( ibot1 == 0 .or. itop2 == 0 .or. ibot2 == 0 ) then
    ierror = 1
    return
  end if

  if ( itop1 == 0 ) then
    itop = 0
    ibot = 1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!  Implicitly invert the divisor fraction here.  The rest of
!  the code will be a multiply operation.
!
  jbot1 = ibot1
  jbot2 = itop2
  jtop1 = itop1
  jtop2 = ibot2
!
!  Get rid of all common factors in top and bottom.
!
  itemp = i_gcd ( jtop1, jbot1 )
  jtop1 = jtop1 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop1, jbot2 )
  jtop1 = jtop1 / itemp
  jbot2 = jbot2 / itemp
  itemp = i_gcd ( jtop2, jbot1 )
  jtop2 = jtop2 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop2, jbot2 )
  jtop2 = jtop2 / itemp
  jbot2 = jbot2 / itemp
!
!  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
!
!  Check the top for overflow.
!
  if ( abs ( real ( jtop1 ) * abs ( jtop2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational fraction.'
    itop = 0
    return
  end if

  itop = jtop1 * jtop2
!
!  Check the bottom IBOT1*ITOP2 for overflow.
!
  if ( abs ( real ( jbot1 ) * real ( jbot2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational fraction.'
    ibot = 1
    return
  end if

  ibot = jbot1 * jbot2
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy.
!
  return
end
subroutine rat_farey ( n, max_frac, num_frac, a, b )
!
!*******************************************************************************
!
!! RAT_FAREY computes the N-th row of the Farey fraction table.
!
!
!  Example:
!
!    N = 5
!
!    NUM_FRAC = 11
!    A =  0  1  1  1  2  1  3  2  3  4  1
!    B =  1  5  4  3  5  2  5  3  4  5  1
!
!  Discussion:
!
!    In this form of the Farey fraction table, fractions in row N lie between
!    0 and 1, are in lowest terms, and have a denominator that is no greater
!    than N.  Row N is computed directly, and does not require the computation
!    of previous rows.
!
!    The data satisfy the relationship:
!
!      A(K+1) * B(K) - A(K) * B(K+1) = 1
!
!    The number of items in the N-th row is roughly N**2 / PI**2.
!
!  Reference:
!
!    Donald Knuth,
!    The Art of Computer Programming,
!    Volume 1, Fundamental Algorithms,
!    Addison Wesley, 1968, page 157.
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the desired row number.  N must be positive.
!
!    Input, integer MAX_FRAC, the maximum number of fractions to compute.
!
!    Output, integer NUM_FRAC, the number of fractions computed.
!
!    Output, integer A(MAX_FRAC), B(MAX_FRAC), contains the NUM_FRAC
!    numerators and denominators of the N-th row of the Farey fraction table.
!
  implicit none
!
  integer max_frac
!
  integer a(max_frac)
  integer b(max_frac)
  integer c
  integer k
  integer n
  integer num_frac
!
  if ( n <= 0 ) then
    num_frac = 0
    return
  end if

  if ( max_frac <= 0 ) then
    num_frac = k
    return
  end if

  k = 1
  a(k) = 0
  b(k) = 1

  if ( max_frac <= 1 ) then
    num_frac = k
    return
  end if

  k = 2
  a(k) = 1
  b(k) = n

  do while ( k < max_frac )

    if ( a(k) == 1 .and. b(k) == 1 ) then
      exit
    end if

    k = k + 1
    c = int ( ( b(k-2) + n ) / b(k-1) )
    a(k) = c * a(k-1) - a(k-2)
    b(k) = c * b(k-1) - b(k-2)

  end do

  num_frac = k

  return
end
subroutine rat_farey2 ( n, a, b )
!
!*******************************************************************************
!
!! RAT_FAREY2 computes the next row of the Farey fraction table.
!
!
!  Example:
!
!    Input:
!
!      N = 3
!      A =  0  1  1  2  1
!      B =  1  3  2  3  1
!
!    Output:
!
!      A =  0  1  1  2  1  3  2  3  1
!      B =  1  4  3  5  2  5  3  4  1
!
!  Discussion:
!
!    In this form of the Farey fraction table, fractions in row N lie between
!    0 and 1, and are in lowest terms.  For every adjacent pair of input
!    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
!    and inserted between them.
!
!    The number of items in the N-th row is 1+2**(N-1).
!
!  Modified:
!
!    31 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the input row number.  N must be nonnegative.
!    If N is zero, then the input is ignored, and the entries of
!    row 1 are computed directly.
!
!    Input/output, integer A(1+2**N), B(1+2**N).
!    On input, entries 1 through 1+2**(N-1) contain the entries of row N.
!    On output, entries 1 through 1+2**N contain the entries of row N+1.
!
  implicit none
!
  integer n
!
  integer a(1+2**n)
  integer b(1+2**n)
  integer i

  if ( n == 0 ) then
    a(1) = 0
    b(1) = 1
    a(2) = 1
    b(2) = 1
    return
  end if
!
!  Shift the current data.
!
  do i = 1+2**(n-1), 1, -1
    a(2*i-1) = a(i)
    b(2*i-1) = b(i)
  end do
!
!  Compute the mediants.
!
  do i = 2, 2**n, 2
    a(i) = a(i-1) + a(i+1)
    b(i) = b(i-1) + b(i+1)
  end do

  return
end
subroutine rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )
!
!*******************************************************************************
!
!! RAT_MUL multiplies two fractions.
!
!
!  Discussion:
!
!    The routine computes
!
!      ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
!
!    while avoiding integer overflow.
!
!  Modified:
!
!    01 February 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ITOP1, IBOT1, the first factor.
!
!    Input, integer ITOP2, IBOT2, the second factor.
!
!    Output, integer ITOP, IBOT, the product.
!
!    Output, integer IERROR.
!    0, no error occurred.
!    1, an error occurred.  The multiplication of the two values
!    requires a numerator or denominator larger than the
!    maximum legal integer.
!
  implicit none
!
  integer i_gcd
  integer i_max
  integer ibot
  integer ibot1
  integer ibot2
  integer ierror
  integer itemp
  integer itop
  integer itop1
  integer itop2
  integer jbot1
  integer jbot2
  integer jtop1
  integer jtop2
!
  ierror = 0

  i_max = huge ( i_max )

  if ( itop1 == 0 .or. itop2 == 0 ) then
    itop = 0
    ibot = 1
    return
  end if
!
!  Make copies of the input arguments, since we will change them.
!
  jbot1 = ibot1
  jbot2 = ibot2
  jtop1 = itop1
  jtop2 = itop2
!
!  Get rid of all common factors in top and bottom.
!
  itemp = i_gcd ( jtop1, jbot1 )
  jtop1 = jtop1 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop1, jbot2 )
  jtop1 = jtop1 / itemp
  jbot2 = jbot2 / itemp
  itemp = i_gcd ( jtop2, jbot1 )
  jtop2 = jtop2 / itemp
  jbot1 = jbot1 / itemp
  itemp = i_gcd ( jtop2, jbot2 )
  jtop2 = jtop2 / itemp
  jbot2 = jbot2 / itemp
!
!  The fraction (ITOP1*ITOP2)/(IBOT1*IBOT2) is in lowest terms.
!
!  Check the top ITOP1*ITOP2 for overflow.
!
  if ( abs ( real ( jtop1 ) * real ( jtop2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
    write ( *, '(a)' ) '  Overflow of top of rational product.'
    itop = 0
    return
  end if

  itop = jtop1 * jtop2
!
!  Check the bottom IBOT1*IBOT2 for overflow.
!
  if ( abs ( real ( jbot1 ) * real ( jbot2 ) ) > real ( i_max ) ) then
    ierror = 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
    write ( *, '(a)' ) '  Overflow of bottom of rational product.'
    ibot = 1
    return
  end if

  ibot = jbot1 * jbot2
!
!  The bottom should be positive.
!
  if ( ibot < 0 ) then
    ibot = - ibot
    itop = - itop
  end if
!
!  The fraction is ITOP/IBOT with no loss of accuracy.
!
  return
end
subroutine rat_sum_formula ( n, a, b )
!
!*******************************************************************************
!
!! RAT_SUM_FORMULA computes the formulas for sums of powers of integers.
!
!
!  Example:
!
!    N = 6
!
!        1    2    3    4    5    6    7
!    -----------------------------------
!    0 | 1    0    0    0    0    0    0
!      |
!    1 | 1    1    0    0    0    0    0
!      | 2    2
!      |
!    2 | 1    1    1    0    0    0    0
!      | 3    2    6
!      |
!    3 | 1    1    1    0    0    0    0
!      | 4    2    4
!      | 
!    4 | 1    1    1    0   -1    0    0
!      | 5    2    3        30
!      |
!    5 | 1    1    5    0   -1    0    0
!      | 6    2   12        12
!      |
!    6 | 1    1    1    0   -1    0    1
!      | 7    2    2         6        42
!
!    The interpretation of row 2, for instance, is:
!
!      sum ( 1 <= I <= N ) I**2 = 1/3 N**3 + 1/2 N**2 + 1/6 N
!
!    This suggests that a more sensible way to display the table
!    is to reverse the order of the entries in the row, so that
!    the entry in column J is the coeeficient of N**J, which is
!    not the case now.
!
!  Reference:
!
!    Robert Owens,
!    Sums of Powers of Integers,
!    Mathematics Magazine,
!    Volume 65, Number 1, February 1992, pages 38-40.
!
!  Modified:
!
!    11 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows of coefficients to compute.
!
!    Output, integer A(0:N,N+1), B(0:N,N+1), the numerator and denominator
!    of the coefficients.
!
  implicit none
!
  integer n
!
  integer a(0:n,1:n+1)
  integer asum
  integer b(0:n,1:n+1)
  integer bsum
  integer i
  integer ierror
  integer j
!
  a(0,1) = 1
  b(0,1) = 1
  a(0,2:n+1) = 0
  b(0,2:n+1) = 1

  do i = 1, n

    asum = 0
    bsum = 0
!
!  Subdiagonal entries are multiples of entries above them.
!
    do j = 1, i

      call rat_mul ( a(i-1,j), b(i-1,j), i, i+2-j, a(i,j), b(i,j), ierror )

      call rat_add ( asum, bsum, a(i,j), b(i,j), asum, bsum, ierror )

    end do
!
!  Diagonal entry is 1 - sum of previous entries in row.
!
    asum = - asum
    call rat_add ( 1, 1, asum, bsum, a(i,i+1), b(i,i+1), ierror )
!
!  Superdiagonal entries are zero.
!
    a(i,i+2:n+1) = 0
    b(i,i+2:n+1) = 1

  end do

  return
end
subroutine rat_to_cfrac ( a, ierror, m, n, ip, iq )
!
!*******************************************************************************
!
!! RAT_TO_CFRAC converts a rational value to a continued fraction.
!
!
!  Discussion:
!
!    The routine is given a rational number represented by IP/IQ, and
!    computes the monic or "simple" continued fraction representation
!    with integer coefficients of the number:
!
!      A(1) + 1/ (A(2) + 1/ (A(3) + ... + 1/A(N) ...))
!
!    The user must dimension A to a value M which is "large enough".
!    The actual number of terms needed in the continued fraction
!    representation cannot be known beforehand.
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Output, integer A(M), contains the continued fraction
!    representation of the number.
!
!    Output, integer IERROR, error indicator.  0 if no error,
!    1 if there was an error, namely, M is not large enough.
!
!    Input, integer M, the dimension of A.  If M is not great
!    enough, the algorithm may run out of space.
!
!    Output, integer N, the actual number of entries used in A.
!
!    Input, integer IP, IQ, the numerator and denominator of the
!    rational value whose continued fraction representation is
!    desired.
!
  implicit none
!
  integer m
!
  integer a(m)
  integer ierror
  integer ip
  integer iq
  integer jp
  integer jq
  integer n
!
  jp = ip
  jq = iq

  n = 0

  do

    n = n + 1

    if ( n > m ) then
      ierror = 1
      return
    end if

    a(n) = jp / jq
    jp = mod ( jp, jq )

    if ( jp == 0 ) then
      return
    end if

    n = n + 1

    if ( n > m ) then
      ierror = 1
      return
    end if

    a(n) = jq / jp
    jq = mod ( jq, jp )

    if ( jq == 0 ) then
      exit
    end if

  end do

  return
end
subroutine rat_to_r ( a, iatop, iabot )
!
!*******************************************************************************
!
!! RAT_TO_R converts rational values to real values.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real A, the value of the rational quantity.
!
!    Input, integer IATOP, IABOT, the rational quantity
!    (IATOP/IABOT) that is to be converted.
!
  implicit none
!
  real a
  integer iabot
  integer iatop
!
  if ( iabot == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RAT_TO_R - Warning!'
    write ( *, '(a)' ) '  The input fraction to be converted had a'
    write ( *, '(a)' ) '  zero denominator.'
    a = 0.0E+00
  else
    a = real ( iatop ) / real ( iabot )
  end if

  return
end
subroutine rat_to_s_left ( ival, jval, s )
!
!*******************************************************************************
!
!! RAT_TO_S_LEFT returns a left-justified representation of IVAL/JVAL.
!
!
!  Discussion:
!
!    If the ratio is negative, a minus sign precedes IVAL.
!    A slash separates IVAL and JVAL.
!
!  Modified:
!
!    13 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, JVAL, the two integers whose
!    ratio IVAL/JVAL is to be represented.
!
!    Note that if IVAL is nonzero and JVAL is 0, STRING will
!    be returned as "Inf" or "-Inf" (Infinity), and if both
!    IVAL and JVAL are zero, STRING will be returned as "NaN"
!    (Not-a-Number).
!
!    Output, character ( len = * ) S, a left-justified string
!    containing the representation of IVAL/JVAL.
!
  implicit none
!
  integer ival
  integer ival2
  integer jval
  integer jval2
  character ( len = * ) s
  character ( len = 22 ) s2
!
!  Take care of simple cases right away.
!
  if ( ival == 0 ) then

    if ( jval /= 0 ) then
      s2 = '0'
    else
      s2= 'NaN'
    end if

  else if ( jval == 0 ) then

    if ( ival > 0 ) then
      s2 = 'Inf'
    else
      s2 = '-Inf'
    end if
!
!  Make copies of IVAL and JVAL.
!
  else

    ival2 = ival
    jval2 = jval

    if ( jval2 == 1 ) then
      write ( s2, '(i11)' ) ival2
    else
      write ( s2, '(i11, ''/'', i10)' ) ival2, jval2
    end if

    call s_blank_delete ( s2 )

  end if

  s = s2

  return
end
subroutine ratmat_det ( n, iatop, iabot, idtop, idbot, ierror )
!
!*******************************************************************************
!
!! RATMAT_DET finds the determinant of an N by N matrix of rational entries.
!
!
!  Discussion:
!
!    The brute force method is used.
!
!    This routine should only be used for small matrices, since this
!    calculation requires the summation of N! products of N numbers.
!
!  Modified:
!
!    17 April 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of A.
!
!    Input, integer IATOP(N,N), IABOT(N,N), the numerators
!    and denominators of the entries of the matrix.
!
!    Output, integer IDTOP, IDBOT, the determinant of the matrix,
!    expressed as IDTOP/IDBOT.
!
!    Output, integer IERROR.
!    0, the determinant was computed.
!    1, an overflow error occurred, and the determinant was not
!    computed.
!
  implicit none
!
  integer n
!
  logical even
  integer i
  integer iabot(n,n)
  integer iatop(n,n)
  integer iarray(n)
  integer ibot
  integer ibot1
  integer ibot2
  integer idbot
  integer idtop
  integer ierror
  integer itop
  integer itop1
  integer itop2
  logical more
!
  ierror = 0

  more = .false.
  idtop = 0
  idbot = 1

  do

    call perm_next ( n, iarray, more, even )

    if ( even ) then
      itop = 1
    else
      itop = -1
    end if

    ibot = 1

    do i = 1, n

      itop1 = itop
      ibot1 = ibot
      itop2 = iatop(i,iarray(i))
      ibot2 = iabot(i,iarray(i))

      call rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

      if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
        write ( *, '(a)' ) '  An overflow occurred.'
        write ( *, '(a)' ) '  The determinant calculation cannot be done'
        write ( *, '(a)' ) '  for this matrix.'
        idtop = 0
        idbot = 1
        return
      end if

    end do

    itop1 = itop
    ibot1 = ibot

    itop2 = idtop
    ibot2 = idbot

    call rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

    if ( ierror == 0 ) then
      idtop = itop
      idbot = ibot
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
      write ( *, '(a)' ) '  An overflow occurred.'
      write ( *, '(a)' ) '  The determinant calculation cannot be done'
      write ( *, '(a)' ) '  for this matrix.'
      idtop = 0
      idbot = 1
      return
    end if

    if ( .not. more ) then
      exit
    end if

  end do
!
!  The bottom should be positive.
!
  if ( idbot < 0 ) then
    idbot = - idbot
    idtop = - idtop
  end if

  return
end
subroutine ratmat_print ( m, n, a, b, title )
!
!*******************************************************************************
!
!! RATMAT_PRINT prints out rational vectors or matrices.
!
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in the matrix.
!
!    Input, integer A(M,N), B(M,N), the current rational or decimal
!    matrix.
!
!    Input, character ( len = * ) TITLE, a label for the object being printed.
!
  implicit none
!
  integer m
  integer n
!
  integer a(m,n)
  integer b(m,n)
  character ( len = 10 ) chrtmp1
  character ( len = 10 ) chrtmp2
  character ( len = 10 ) chrtmp3
  character ( len = 40 ) format1
  character ( len = 40 ) format2
  integer i
  integer ichi
  integer iclo
  integer imax
  integer imin
  integer ione
  integer itemp
  integer izhi
  integer izlo
  integer j
  integer jmax
  integer jmin
  integer kmax
  integer, parameter :: ncolum = 80
  integer none
  integer npline
  character ( len = 100 ) output
  character ( len = * ) title
!
!  Figure out how many rationals we can get in NCOLUM columns.
!
  kmax = 3

  do i = 1, m
    do j = 1, n

      itemp = abs ( a(i,j) )

      do while ( itemp >= 10**(kmax-2) )
        kmax = kmax + 1
      end do

      itemp = abs ( b(i,j) )

      do while ( itemp > 10**(kmax-2) )
        kmax = kmax + 1
      end do

    end do
  end do

  kmax = kmax + 1
  npline = ncolum / kmax
!
!  Create the formats.
!
  call i_to_s_left ( npline, chrtmp2 )
  call i_to_s_left ( kmax, chrtmp3 )

  format1 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format1 )

  format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

  call s_blank_delete ( format2 )

  do jmin = 1, n, npline

    jmax = min ( jmin+npline-1, n )

    write ( *, '(a)' ) ' '

    if ( jmin == 1 ) then
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
    end if

    if ( jmin > 1 .or. jmax < n ) then
      write ( output, * ) 'Columns ', jmin, ' to ', jmax
      call s_blanks_delete ( output )
      write ( *, '(a)' ) trim ( output )
      write ( *, '(a)' ) ' '
    end if

    do i = 1, m

      write ( *, format1 ) a(i,jmin:jmax)
      write ( output, format1 ) b(i,jmin:jmax)
!
!  Delete each denominator that is 1.  If all are 1, don't
!  even print out the line.
!
      none = 0

      do j = jmin, jmax

        if ( b(i,j) == 1 ) then
          ione = (j-jmin+1) * kmax
          output(ione:ione) = ' '
        else
          none = 1
        end if

      end do

      write ( *, '(a)' ) trim ( output )

      if ( jmax == n .and. i == m ) then
      else
        write ( *, '(a)' ) ' '
      end if

    end do

  end do

  return
end
subroutine regro_next ( done, n, v, vmax )
!
!*******************************************************************************
!
!! REGRO_NEXT computes restricted growth functions one at a time.
!
!
!  Definition:
!
!    A restricted growth function on N is a vector (V(1), ..., V(N) )
!    of values V(I) between 1 and N, satisfying the requirements:
!      V(1) = 1;
!      V(I) <= 1 + max ( V(1), V(2), ..., V(I-1) ).
!
!  Comments:
!
!    The number of restricted growth functions on N is equal to
!    the Bell number B(N).
!
!    There is a bijection between restricted growth functions on N
!    and set partitions of N.
!
!  Examples:
!
!    The 15 restricted growth functions for N = 4 are:
!
!    (1111), (1112), (1121), (1122), (1123),
!    (1211), (1212), (1213), (1221), (1222),
!    (1223), (1231), (1232), (1233), (1234).
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986, page 19.
!
!  Modified:
!
!    15 July 1998
!
!  Parameters:
!
!    Input/output, logical DONE.
!    On first call, set DONE to .TRUE., and then do not alter it.
!
!    On output, DONE will be .FALSE. if the routine has computed another
!    restricted growth function, or .TRUE. if all the restricted
!    growth functions have been returned.
!
!    Input, integer N, the number of components in the restricted growth
!  function.
!
!    Input/output, integer V(N).  The user need not set this quantity
!    before the initial call, and should not alter it between successive
!    calls.  On each return from the routine, with DONE = .FALSE.,
!    V will contain the componentwise values of the next restricted
!    growth function.
!
!    Input/output, integer VMAX(N).  The user need not set this quantity
!    before the initial call, and should not alter it between calls.
!    VMAX(I) records the largest value that component V(I) could take,
!    given the values of components 1 through I-1.
!
  implicit none
!
  integer n
!
  logical done
  integer i
  integer j
  integer v(n)
  integer vmax(n)
!
!  First call:
!
  if ( done ) then

    v(1:n) = 1

    vmax(1) = 1
    vmax(2:n) = 2

    done = .false.
!
!  Later calls.
!
  else

    j = n

    do

      if ( j == 1 ) then
        done = .true.
        return
      end if

      if ( v(j) /= vmax(j) ) then
        exit
      end if

      j = j - 1

    end do

    v(j) = v(j) + 1

    do i = j+1, n

      v(i) = 1

      if ( v(j) == vmax(j) ) then
        vmax(i) = vmax(j) + 1
      else
        vmax(i) = vmax(j)
      end if

    end do

  end if

  return
end
subroutine rfrac_to_cfrac ( m, p, q, t, ierror )
!
!*******************************************************************************
!
!! RFRAC_TO_CFRAC converts a rational polynomial fraction to a continued fraction.
!
!
!  Discussion:
!
!    That is, it accepts
!
!      P(1) + P(2) * X + ... + P(M) * X**(M-1)
!      -------------------------------------------------------
!      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!    and returns the equivalent continued fraction:
!
!      1 / (T(1) + X/(T(2) + X/(...T(2*M-1)+X/(T(2*M) ... )))
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    17 April 2000
!
!  Parameters:
!
!    Input, integer M, defines the number of P coefficients,
!    and is one less than the number of Q coefficients, and one
!    half the number of T coefficients.
!
!    Input, real P(M), Q(M+1), the coefficients defining the rational
!    polynomial fraction.
!
!    Output, real T(2*M), the coefficients defining the continued fraction.
!
!    Output, integer IERROR, error flag.
!    0, no error;
!    nonzero, the algorithm broke down at some point with a zero divisor.
!
  implicit none
!
  integer m
!
  real a(m+1,2*m+1)
  integer i
  integer ierror
  integer k
  real p(m)
  real q(m+1)
  real t(2*m)
  real ta
!
  ierror = 0

  a(1:m+1,1) = q(1:m+1)
  a(1:m,  2) = p(1:m)

  t(1) = a(1,1) / a(1,2)
  ta = a(m+1,1)

  do i = 1, m
    a(m-i+1,2*i+1) = ta
  end do

  do k = 1, 2*m-2

    do i = 1, (2*m-k)/2
      a(i,k+2) = a(i+1,k) - t(k) * a(i+1,k+1)
    end do

    if ( a(1,k+2) == 0.0E+00 ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RFRAC_TO_CFRAC - Fatal error!'
      write ( *, '(a,i6)' ) '  A(1,K+2) is zero for K = ', k
      stop
    end if

    t(k+1) = a(1,k+1) / a(1,k+2)

  end do

  t(2*m) = a(1,2*m) / a(1,2*m+1)

  return
end
subroutine rfrac_to_jfrac ( m, p, q, r, s, ierror )
!
!*******************************************************************************
!
!! RFRAC_TO_JFRAC converts a rational polynomial fraction to a J fraction.
!
!
!  Discussion:
!
!    The routine accepts
!
!    P(1) + P(2) * X + ... + P(M) * X**(M-1)
!    -------------------------------------------------------
!    Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
!
!    and returns the equivalent J-fraction:
!
!    R(1)/ (X + S(1) + R(2)/(X + S(2) + R(3)/ ... + R(M)/(X+S(M))... ))
!
!  Reference:
!
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Modified:
!
!    17 April 2000
!
!  Parameters:
!
!    Input, integer M, defines the number of P, R, and S coefficients,
!    and is one less than the number of Q coefficients.
!
!    Input, real P(M), Q(M+1), the coefficients defining the rational
!    polynomial fraction.
!
!    Output, real R(M), S(M), the coefficients defining the
!    J-fraction.
!
!    Output, integer IERROR, error flag.
!    IERROR is zero on return, unless the algorithm broke down at
!    some point with a zero divisor, in which case IERROR = 1.
!
  implicit none
!
  integer m
!
  real a(m+1,m+1)
  integer i
  integer ierror
  integer k
  real p(m)
  real q(m+1)
  real r(m)
  real s(m)
!
  ierror = 0

  a(1:m+1,1) = q(1:m+1)
  a(1:m,  2) = p(1:m)

  if ( m > 1 ) then

    r(1) = a(m,2) / a(m+1,1)
    s(1) = ( r(1) * a(m,1) - a(m-1,2) ) / a(m,2)

    do k = 1, m-2

      a(1,k+2) = r(k) * a(1,k) - s(k) * a(1,k+1)

      do i = 2, m-k

        a(i,k+2) = r(k) * a(i,k) - a(i-1,k+1) - s(k) * a(i,k+1)

        if ( a(m-k,k+2) == 0.0E+00 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RFRAC_TO_JFRAC - Fatal error!'
          write ( *, '(a,i6)' ) '  A(M-K,K+2) = 0 for K=', k
          stop
        end if

        r(k+1) = a(m-k,k+2) / a(m-k+1,k+1)
        s(k+1) = ( r(k+1) * a(m-k,k+1) - a(m-k-1,k+2) ) / a(m-k,k+2)

      end do

    end do

    a(1,m+1) = r(m-1) * a(1,m-1) - s(m-1) * a(1,m)

  end if

  r(m) = a(1,m+1) / a(2,m)
  s(m) = a(1,m) / a(2,m)

  return
end
function rise ( x, n )
!
!*******************************************************************************
!
!! RISE computes the rising factorial function [X]^N.
!
!
!  Discussion:
!
!    [X}^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
!
!    Note that the number of ways of arranging N objects in M ordered
!    boxes is [M}^N.  (Here, the ordering in each box matters).  Thus,
!    2 objects in 2 boxes have the following 6 possible arrangements:
!
!      -/12, 1/2, 12/-, -/21, 2/1, 21/-.
!
!    Moreover, the number of non-decreasing maps from a set of
!    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
!    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
!
!      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
!      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the argument of the rising factorial function.
!
!    Input, integer N, the order of the rising factorial function.
!    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
!    negative, a "falling" factorial will be computed.
!
!    Output, real RISE, the value of the rising factorial function.
!
  implicit none
!
  real arg
  integer i
  integer n
  real rise
  real x
!
  rise = 1.0E+00

  arg = x

  if ( n > 0 ) then

    do i = 1, n
      rise = rise * arg
      arg = arg + 1.0E+00
    end do

  else if ( n < 0 ) then

    do i = -1, n
      rise = rise * arg
      arg = arg - 1.0E+00
    end do

  end if

  return
end
subroutine rmat_det ( n, a, det )
!
!*******************************************************************************
!
!! RMAT_DET finds the determinant of a real N by N matrix.
!
!
!  Discussion:
!
!    A brute force calculation is made.
!
!    This routine should only be used for small matrices, since this
!    calculation requires the summation of N! products of N numbers.
!
!  Modified:
!
!    10 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of rows and columns of A.
!
!    Input, real A(N,N), the matrix whose determinant is desired.
!
!    Output, real DET, the determinant of the matrix.
!
  implicit none
!
  integer n
!
  real a(n,n)
  real det
  logical even
  integer i
  integer iarray(n)
  logical more
  real term
!
  more = .false.
  det = 0.0E+00

  do

    call perm_next ( n, iarray, more, even )

    if ( even ) then
      term = 1.0E+00
    else
      term = -1.0E+00
    end if

    do i = 1, n
      term = term * a(i,iarray(i))
    end do

    det = det + term

    if ( .not. more ) then
      exit
    end if

  end do

  return
end
subroutine rmat_perm ( n, a, p )
!
!*******************************************************************************
!
!! RMAT_PERM permutes the rows and columns of a square real matrix.
!
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
!    10 June 2002
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(N,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(N), a permutation to be applied to the rows
!    and columns.  P(I) is the new number of row and column I.
!
  implicit none
!
  integer n
!
  real a(n,n)
  integer i
  integer i1
  integer is
  real it
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(n)
!
  call perm_cycle ( n, p, is, nc, 1 )

  do i = 1, n

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( p(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            it = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( p(j1) )

              call r_swap ( a(i1,j1), it )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( p(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do
!
!  Restore the positive signs of the data.
!
  p(1:n) = abs ( p(1:n) )

  return
end
subroutine rmat_perm2 ( m, n, a, p, q )
!
!*******************************************************************************
!
!! RMAT_PERM2 permutes rows and columns of a real rectangular matrix, in place.
!
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
!    28 October 1999
!
!  Parameters:
!
!    Input, integer M, number of rows in the matrix.
!
!    Input, integer N, number of columns in the matrix.
!
!    Input/output, real A(M,N).
!    On input, the matrix to be permuted.
!    On output, the permuted matrix.
!
!    Input, integer P(M), the row permutation.  P(I) is the new number of row I.
!
!    Input, integer Q(N), the column permutation.  Q(I) is the new number of
!    column I.  Note that the routine allows you to pass a single array as both
!    P and Q.
!
  implicit none
!
  integer m
  integer n
!
  real a(m,n)
  integer i
  integer i1
  integer is
  integer j
  integer j1
  integer j2
  integer k
  integer lc
  integer nc
  integer p(m)
  integer q(n)
  real t
!
  call perm_cycle ( m, p, is, nc, 1 )

  if ( q(1) > 0 ) then
    call perm_cycle ( n, q, is, nc, 1 )
  end if

  do i = 1, m

    i1 = - p(i)

    if ( i1 > 0 ) then

      lc = 0

      do

        i1 = p(i1)
        lc = lc + 1

        if ( i1 <= 0 ) then
          exit
        end if

      end do

      i1 = i

      do j = 1, n

        if ( q(j) <= 0 ) then

          j2 = j
          k = lc

          do

            j1 = j2
            t = a(i1,j1)

            do

              i1 = abs ( p(i1) )
              j1 = abs ( q(j1) )

              call r_swap ( a(i1,j1), t )

              if ( j1 /= j2 ) then
                cycle
              end if

              k = k - 1

              if ( i1 == i ) then
                exit
              end if

            end do

            j2 = abs ( q(j2) )

            if ( k == 0 ) then
              exit
            end if

          end do

        end if

      end do

    end if

  end do

  p(1:m) = abs ( p(1:m) )

  if ( q(1) <= 0 ) then

    q(1:n) = abs ( q(1:n) )

  end if

  return
end
subroutine rmat_permanent ( n, a, perm )
!
!*******************************************************************************
!
!! RMAT_PERMANENT computes the permanent of a real matrix.
!
!
!  Discussion:
!
!    The permanent function is similar to the determinant.  Recall that
!    the determinant of a matrix may be defined as the sum of all the
!    products:
!
!      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
!
!    where I is any permutation of the columns of the matrix, and S is the
!    sign of the permutation.  By contrast, the permanent function is
!    the (unsigned) sum of all the products
!
!      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
!
!    where I is any permutation of the columns of the matrix.  The only
!    difference is that there is no permutation sign multiplying each summand.
!
!    Symbolically, then, the determinant of a 2 by 2 matrix
!
!      a b
!      c d
!
!    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
!
!
!    The permanent is invariant under row and column permutations.
!    If a row or column of the matrix is multiplied by S, then the
!      permanent is likewise multiplied by S.
!    If the matrix is square, then the permanent is unaffected by
!      transposing the matrix.
!    Unlike the determinant, however, the permanent does change if
!      one row is added to another, and it is not true that the
!      permanent of the product is the product of the permanents.
!
!
!    Note that if A is a matrix of all 1's and 0's, then the permanent
!    of A counts exactly which permutations hit exactly 1's in the matrix.
!    This fact can be exploited for various combinatorial purposes.
!
!    For instance, setting the diagonal of A to 0, and the offdiagonals
!    to 1, the permanent of A counts the number of derangements of N
!    objects.
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
!    Input, integer N, number of rows and columns in matrix.
!
!    Input, real A(N,N).  The matrix whose permanent is desired.
!
!    Output, real PERM, the value of the permanent function of A.
!
  implicit none
!
  integer n
!
  real a(n,n)
  integer i
  integer iadd
  integer iwork(n)
  integer j
  logical more
  integer ncard
  real p
  real perm
  real prod
  real sgn
  real work(n)
  real z
!
  more = .false.

  do i = 1, n
    work(i) = a(i,n)
    do j = 1, n
      work(i) = work(i) - 0.5E+00 * a(i,j)
    end do
  end do

  p = 0.0E+00
  sgn = - 1.0E+00

  do

    sgn = - sgn
    call sub_next ( n-1, iwork, more, ncard, iadd )

    if ( ncard /= 0 ) then
      z = real ( 2 * iwork(iadd) - 1 )
      work(1:n) = work(1:n) + z * a(1:n,iadd)
    end if

    p = p + sgn * product ( work )

    if ( .not. more ) then
      exit
    end if

  end do

  perm = p * real ( 4 * mod ( n, 2 ) - 2 )

  return
end
subroutine rmat_print ( m, n, a, title )
!
!*******************************************************************************
!
!! RMAT_PRINT prints a real matrix.
!
!
!  Modified:
!
!    24 March 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows in A.
!
!    Input, integer N, the number of columns in A.
!
!    Input, real A(M,N), the matrix to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer m
  integer n
!
  real a(m,n)
  integer i
  integer j
  integer jhi
  integer jlo
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) title
  end if

  do jlo = 1, n, 5
    jhi = min ( jlo + 4, n )
    write ( *, '(a)' ) ' '
    write ( *, '(6x,5(i7,7x))' ) ( j, j = jlo, jhi )
    write ( *, '(a)' ) ' '
    do i = 1, m
      write ( *, '(i6,5g14.6)' ) i, a(i,jlo:jhi)
    end do
  end do

  return
end
subroutine rpoly ( n, a, x0, iopt, val )
!
!*******************************************************************************
!
!! RPOLY performs operations on real polynomials in power or factorial form.
!
!
!  Discussion:
!
!    The power sum form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*X**2 + ... + (AN+1)*X**N
!
!    The Taylor expansion at C has the form
!
!      P(X) = A1 + A2*(X-C) + A3*(X-C)**2+... + (AN+1)*(X-C)**N
!
!    The factorial form of a polynomial is
!
!      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2) + ...
!        + (AN+1)*(X)*(X-1)*...*(X-N+1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of coefficients in the polynomial
!    (in other words, the polynomial degree + 1)
!
!    Input/output, real A(N), the coefficients of the polynomial.  Depending
!    on the option chosen, these coefficients may be overwritten by those
!    of a different form of the polynomial.
!
!    Input, real X0, for IOPT = -1, 0, or positive, the value of the
!    argument at which the polynomial is to be evaluated, or the
!    Taylor expansion is to be carried out.
!
!    Input, integer IOPT, a flag describing which algorithm is to
!    be carried out:
!
!    -3: Reverse Stirling.  Input the coefficients of the polynomial in
!    factorial form, output them in power sum form.
!
!    -2: Stirling.  Input the coefficients in power sum
!    form, output them in factorial form.
!
!    -1: Evaluate a polynomial which has been input
!    in factorial form.
!
!    0:  Evaluate a polynomial input in power sum form.
!
!    1 or more:  Given the coefficients of a polynomial in
!    power sum form, compute the first IOPT coefficients of
!    the polynomial in Taylor expansion form.
!
!    Output, real VAL, for IOPT = -1 or 0, the value of the
!    polynomial at the point X0.
!
  implicit none
!
  integer n
!
  real a(n)
  real eps
  integer i
  integer iopt
  integer m
  integer n1
  real val
  real w
  real x0
  real z
!
  n1 = min ( n, iopt )
  n1 = max ( 1, n1 )

  if ( iopt < -1 ) then
    n1 = n
  end if

  eps = real ( mod ( max ( -iopt, 0 ), 2 ) )

  w = - real ( n ) * eps

  if ( iopt > -2 ) then
    w = w + x0
  end if

  do m = 1, n1

    val = 0.0E+00
    z = w

    do i = m, n
      z = z + eps
      val = a(n+m-i) + z * val
      if ( iopt /= 0 .and. iopt /= -1 ) then
        a(n+m-i) = val
      end if
    end do

    if ( iopt < 0 ) then
      w = w + 1.0E+00
    end if

  end do

  return
end
subroutine rpoly_degree ( na, a, degree )
!
!*******************************************************************************
!
!! RPOLY_DEGREE returns the degree of a polynomial.
!
!
!  Discussion:
!
!    The degree of a polynomial is the index of the highest power
!    of X with a nonzero coefficient.
!
!    The degree of a constant polynomial is 0.  The degree of the
!    zero polynomial is debatable, but this routine returns the
!    degree as 0.
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real A(0:NA), the coefficients of the polynomials.
!
!    Output, integer DEGREE, the degree of A.
!
  implicit none
!
  integer na
!
  real a(0:na)
  integer degree
!
  degree = na

  do while ( degree > 0 )

    if ( a(degree) /= 0.0E+00 ) then
      return
    end if

    degree = degree - 1

  end do

  return
end
subroutine rpoly_div ( na, a, nb, b, nq, q, nr, r )
!
!*******************************************************************************
!
!! RPOLY_DIV computes the quotient and remainder of two polynomials.
!
!
!  Discussion:
!
!    The polynomials are assumed to be stored in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    23 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real A(0:NA), the coefficients of the polynomial to be divided.
!
!    Input, integer NB, the dimension of B.
!
!    Input, real B(0:NB), the coefficients of the divisor polynomial.
!
!    Output, integer NQ, the degree of Q.
!    If the divisor polynomial is zero, NQ is returned as -1.
!
!    Output, real Q(0:NA-NB), contains the quotient of A/B.
!    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
!    In any case, Q(0:NA) should be enough.
!
!    Output, integer NR, the degree of R.
!    If the divisor polynomial is zero, NR is returned as -1.
!
!    Output, real R(0:NB-1), contains the remainder of A/B.
!    If B has full degree, R should be dimensioned R(0:NB-1).
!    Otherwise, R will actually require less space.
!
  implicit none
!
  integer na
  integer nb
!
  real a(0:na)
  real a2(0:na)
  real b(0:nb)
  integer i
  integer na2
  integer nb2
  integer nq
  integer nr
  real q(0:*)
  real r(0:*)
!
  call rpoly_degree ( na, a, na2 )
  call rpoly_degree ( nb, b, nb2 )

  if ( b(nb2) == 0.0E+00 ) then
    nq = -1
    nr = -1
    return
  end if

  a2(0:na) = a(0:na)

  nq = na2 - nb2
  nr = nb2 - 1

  do i = nq, 0, -1
    q(i) = a2(i+nb2) / b(nb2)
    a2(i+nb2) = 0.0E+00
    a2(i:i+nb2-1) = a2(i:i+nb2-1) - q(i) * b(0:nb2-1)
  end do

  r(0:nr) = a2(0:nr)

  return
end
subroutine rpoly_f2p ( n, a )
!
!*******************************************************************************
!
!! RPOLY_F2P converts a real polynomial from factorial form to power sum form.
!
!
!  Discussion:
!
!    The factorial form is
!
!      p(x) =   a(1)
!             + a(2) * x
!             + a(3) * x*(x-1)
!             ...
!             + a(n) * x*(x-1)*...*(x-(n-2))
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, real N, the dimension of A.
!
!    Input/output, real A(N), on input, the polynomial
!    coefficients in factorial form.  On output, the polynomial
!    coefficients in power sum form.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer m
  real val
  real w
  real z
!
  w = - n

  do m = 1, n

    val = 0.0E+00
    z = w

    do i = m, n
      z = z + 1.0E+00
      val = a(n+m-i) + z * val
      a(n+m-i) = val
    end do

    w = w + 1.0E+00

  end do

  return
end
subroutine rpoly_fval ( n, a, x, val )
!
!*******************************************************************************
!
!! RPOLY_FVAL evaluates a real polynomial in factorial form.
!
!
!  Discussion:
!
!    The factorial form of a polynomial is:
!
!      p(x) = a(1)
!           + a(2)  *x
!           + a(3)  *x*(x-1)
!           +...
!           + a(n-1)*x*(x-1)*(x-2)...*(x-(n-3))
!           + a(n)  *x*(x-1)*(x-2)...*(x-(n-3))*(x-(n-2))
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Output, real VAL, the value of the polynomial at X.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  real val
  real x
!
  val = 0.0E+00
  do i = 1, n
    val = a(n+1-i) + ( x - n + i ) * val
  end do

  return
end
subroutine rpoly_mul ( na, a, nb, b, c )
!
!*******************************************************************************
!
!! RPOLY_MUL computes the product of two real polynomials A and B.
!
!
!  Discussion:
!
!    The polynomials are in power sum form.
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real A(0:NA), the coefficients of the first polynomial factor.
!
!    Input, integer NB, the dimension of B.
!
!    Input, real B(0:NB), the coefficients of the second polynomial factor.
!
!    Output, real C(0:NA+NB), the coefficients of A * B.
!
  implicit none
!
  integer na
  integer nb
!
  real a(0:na)
  real b(0:nb)
  real c(0:na+nb)
  real d(0:na+nb)
  integer i
  integer j
!
  d(0:na+nb) = 0.0E+00

  do i = 0, na
    d(i:i+nb) = d(i:i+nb) + a(i) * b(0:nb)
  end do

  c(0:na+nb) = d(0:na+nb)

  return
end
subroutine rpoly_n2p ( n, a, xarray )
!
!*******************************************************************************
!
!! RPOLY_N2P converts a real polynomial from Newton form to power sum form.
!
!
!  Discussion:
!
!    This is done by shifting all the Newton abscissas to zero.
!
!    Actually, what happens is that the abscissas of the Newton form
!    are all shifted to zero, which means that A is the power sum
!    polynomial description and A, XARRAY is the Newton polynomial
!    description.  It is only because all the abscissas are shifted to
!    zero that A can be used as both a power sum and Newton polynomial
!    coefficient array.
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(1) + a(2)*x + ... + a(n-1)*x**(n-2) + a(n)*x**(n-1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input/output, real A(N).  On input, the coefficients
!    of the polynomial in Newton form, and on output, the coefficients
!    in power sum form.
!
!    Input/output, real XARRAY(N).  On input, the abscissas of
!    the Newton form of the polynomial.  On output, these values
!    have all been set to zero.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  real xarray(n)
  real, parameter :: zero = 0.0E+00
!
  do i = 1, n
    call rpoly_nx ( n, a, xarray, zero )
  end do

  return
end
subroutine rpoly_nval ( n, a, xarray, x, val )
!
!*******************************************************************************
!
!! RPOLY_NVAL evaluates a real polynomial in Newton form.
!
!
!  Definition:
!
!    The Newton form of a polynomial is;
!
!      p(x) = a(1)
!           + a(2)  *(x-x1)
!           + a(3)  *(x-x1)*(x-x2)
!           +...
!           + a(n-1)*(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))
!           + a(n)  *(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))*(x-x(n-1))
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real A(N), the coefficients of the polynomial.
!    A(1) is the constant term.
!
!    Input, real XARRAY(N-1), the N-1 points X which are part
!    of the definition of the polynomial.
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Output, real VAL, the value of the polynomial at X.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  real val
  real x
  real xarray(n-1)
!
  val = a(n)
  do i = n-1, 1, -1
    val = a(i) + ( x - xarray(i) ) * val
  end do

  return
end
subroutine rpoly_nx ( n, a, xarray, x )
!
!*******************************************************************************
!
!! RPOLY_NX replaces one of the base points in a polynomial in Newton form.
!
!
!  Discussion:
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input/output, real A(N), the polynomial coefficients of the Newton form.
!
!    Input/output, real XARRAY(N), the set of abscissas that
!    are part of the Newton form of the polynomial.  On output,
!    the abscissas have been shifted up one index, so that
!    the first location now holds X, and the original contents
!    of XARRAY(N) have been discarded.
!
!    Input, real X, the new point to be shifted into XARRAY.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  real x
  real xarray(n)
!
  do i = n-1, 1, -1
    a(i) = a(i) + ( x - xarray(i) ) * a(i+1)
  end do

  do i = n, 2, -1
    xarray(i) = xarray(i-1)
  end do

  xarray(1) = x

  return
end
subroutine rpoly_p2f ( n, a )
!
!*******************************************************************************
!
!! RPOLY_P2F converts a real polynomial from power sum form to factorial form.
!
!
!  Discussion:
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!    The factorial form is
!
!      p(x) =   a(1)
!             + a(2) * x
!             + a(3) * x*(x-1)
!             ...
!             + a(n) * x*(x-1)*...*(x-(n-2))
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, real N, the dimension of A.
!
!    Input/output, real A(N), on input, the polynomial
!    coefficients in the power sum form, on output, the polynomial
!    coefficients in factorial form.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer m
  real val
!
  do m = 1, n
    val = 0.0E+00
    do i = m, n
      val = a(n+m-i) + real ( m - 1 ) * val
      a(n+m-i) = val
    end do
  end do

  return
end
subroutine rpoly_p2n ( n, a, xarray )
!
!*******************************************************************************
!
!! RPOLY_P2N converts a real polynomial from power sum form to Newton form.
!
!
!  Discussion:
!
!    This is done by shifting all the Newton abscissas from zero.
!
!    The power sum form of a polynomial is:
!
!      p(x) = a(1) + a(2)*x + ... + a(n-1)*x**(n-2) + a(n)*x**(n-1)
!
!    The Newton form of a polynomial is described by an array of N coefficients
!    A and N abscissas X:
!
!      p(x) =   a(1)
!             + a(2) * (x-x(1))
!             + a(3) * (x-x(1)) * (x-x(2))
!             ...
!             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
!
!    X(N) does not occur explicitly in the formula for the evaluation of p(x),
!    although it is used in deriving the coefficients A.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input/output, real A(N).  On input, the coefficients
!    of the polynomial in power sum form, and on output, the
!    coefficients in Newton form.
!
!    Input/output, real XARRAY(N).  On input, the desired abscissas of
!    the Newton form of the polynomial.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  real xarray(n)
  real work(n)
!
  work(1:n) = 0.0E+00

  do i = n, 1, -1
    call rpoly_nx ( n, a, work, xarray(i) )
  end do

  return
end
subroutine rpoly_p2t ( n, a, x )
!
!*******************************************************************************
!
!! RPOLY_P2T converts a real polynomial from power sum form to Taylor form.
!
!
!  Discussion:
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!    The Taylor form is
!
!      p(x) =   a(1)
!             + a(2) * (x-x0)
!             + a(3) * (x-x0)**2
!             ...
!             + a(n) * (x-x0)**(n-1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, real N, the dimension of A.
!
!    Input/output, real A(N), on input, the coefficients in
!    power sum form, and on output, the coefficients in Taylor form.
!
!    Input, real X, the point at which the Taylor form of the
!    polynomial is to be based.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer m
  real val
  real x
!
  do m = 1, n
    val = 0.0E+00
    do i = m, n
      val = a(n+m-i) + x * val
      a(n+m-i) = val
    end do
  end do

  return
end
subroutine rpoly_power ( na, a, p, b )
!
!*******************************************************************************
!
!! RPOLY_POWER computes a positive integer power of a polynomial.
!
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    20 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NA, the dimension of A.
!
!    Input, real A(0:NA), the polynomial to be raised to the power.
!
!    Input, integer P, the nonnegative power to which A is raised.
!
!    Output, real B(0:P*NA), the power of the polynomial.
!
  implicit none
!
  integer na
  integer p
!
  real a(0:na)
  real b(0:p*na)
  integer i
  integer j
  integer nonzer
!
!  Zero out B.
!
  b(0:p*na) = 0.0E+00
!
!  Search for the first nonzero element in A.
!
  nonzer = 0

  do i = 0, na
    if ( a(i) /= 0.0E+00 ) then
      nonzer = i
      exit
    end if
  end do

  if ( nonzer == 0 ) then
    return
  end if

  b(0) = a(nonzer)**p

  do i = 1, p*(na-nonzer)

    if ( i + nonzer <= na ) then
      b(i) = real ( i * p ) * b(0) * a(i+nonzer)
    else
      b(i) = 0.0E+00
    end if

    do j = 1, i-1

      if ( j+nonzer <= na ) then
        b(i) = b(i) - real ( i - j ) * a(j+nonzer) * b(i-j)
      end if

      if ( i-j+nonzer <= na ) then
        b(i) = b(i) + real ( i - j ) * real ( p ) * b(j) * a(i-j+nonzer)
      end if

    end do

    b(i) = b(i) / ( real ( i ) * a(nonzer) )

  end do
!
!  Shift B up.
!
  do i = p*nonzer, p*na
    b(i) = b(i-p*nonzer)
  end do

  do i = 0, p * nonzer-1
    b(i) = 0.0E+00
  end do

  return
end
subroutine rpoly_print ( n, a, title )
!
!*******************************************************************************
!
!! RPOLY_PRINT prints out a polynomial.
!
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    06 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real A(0:N), the polynomial coefficients.
!    A(0) is the constant term and
!    A(N) is the coefficient of X**N.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none
!
  integer n
!
  real a(0:n)
  integer i
  real mag
  integer n2
  character plus_minus
  character ( len = * ) title
!
  if ( len_trim ( title ) > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  call rpoly_degree ( n, a, n2 )

  if ( a(n2) < 0.0E+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n2) )

  if ( n2 >= 2 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n2
  else if ( n2 == 1 ) then
    write ( *, '( '' p(x) = '', a1, g14.6, '' * x'' )' ) plus_minus, mag
  else if ( n2 == 0 ) then
    write ( *, '( '' p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n2-1, 0, -1

    if ( a(i) < 0.0E+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0E+00 ) then

      if ( i >= 2 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i
      else if ( i == 1 ) then
        write ( *, ' ( ''        '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 0 ) then
        write ( *, ' ( ''        '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine rpoly_pval ( n, a, x, val )
!
!*******************************************************************************
!
!! RPOLY_PVAL evaluates a real polynomial in power sum form.
!
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
!
!  Modified:
!
!    08 December 1999
!
!  Parameters:
!
!    Input, integer N, the dimension of A.
!
!    Input, real A(0:N), the coefficients of the polynomial.
!    A(0) is the constant term.
!
!    Input, real X, the point at which the polynomial is to be evaluated.
!
!    Output, real VAL, the value of the polynomial at X.
!
  implicit none
!
  integer n
!
  real a(0:n)
  integer i
  real val
  real x
!
  val = 0.0E+00
  do i = n, 0, -1
    val = val * x + a(i)
  end do

  return
end
subroutine rpoly_t2p ( n, a, x )
!
!*******************************************************************************
!
!! RPOLY_T2P converts a real polynomial from Taylor form to power sum form.
!
!
!  Discussion:
!
!    The Taylor form is
!
!      p(x) =   a(1)
!             + a(2) * (x-x0)
!             + a(3) * (x-x0)**2
!             ...
!             + a(n) * (x-x0)**(n-1)
!
!    The power sum form is
!
!      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, real N, the dimension of A.
!
!    Input/output, real A(N).  On input, the coefficients in Taylor form,
!    and on output, the coefficients in power sum form.
!
!    Input, real X, the point at which the Taylor form polynomial is based.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  integer j
  real x
!
  do i = n, 1, -1
    do j = i, n-1
      a(j) = a(j) - a(j+1) * x
    end do
  end do

  return
end
subroutine rvec_backtrack ( n, x, indx, k, nstack, stack, maxstack, ncan )
!
!*******************************************************************************
!
!! RVEC_BACKTRACK supervises a backtrack search for a real vector.
!
!
!  Discussion:
!
!    The routine tries to construct a real vector one index at a time,
!    using possible candidates as supplied by the user.
!
!    At any time, the partially constructed vector may be discovered to be
!    unsatisfactory, but the routine records information about where the
!    last arbitrary choice was made, so that the search can be
!    carried out efficiently, rather than starting out all over again.
!
!    First, call the routine with INDX = 0 so it can initialize itself.
!
!    Now, on each return from the routine, if INDX is:
!      1, you've just been handed a complete candidate vector;
!         Admire it, analyze it, do what you like.
!      2, please determine suitable candidates for position X(K).
!         Return the number of candidates in NCAN(K), adding each
!         candidate to the end of STACK, and increasing NSTACK.
!      3, you're done.  Stop calling the routine;
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
!    24 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of positions to be filled in the vector.
!
!    Input/output, real X(N), the partial or complete candidate vector.
!
!    Input/output, integer INDX, a communication flag.
!    On input,
!      0 to start a search.
!    On output:
!      1, a complete output vector has been determined and returned in X(1:N);
!      2, candidates are needed for position X(K);
!      3, no more possible vectors exist.
!
!    Output, integer K, if INDX=2, the current vector index being considered.
!
!    Input/output, integer NSTACK, the current length of the stack.
!
!    Input, real STACK(MAXSTACK), a list of all current candidates for
!    all positions 1 through K.
!
!    Input, integer MAXSTACK, the maximum length of the stack.
!
!    Input/output, integer NCAN(N), lists the current number of candidates for
!    positions 1 through K.
!
  implicit none
!
  integer n
  integer maxstack
!
  integer indx
  integer k
  integer ncan(n)
  integer nstack
  real stack(maxstack)
  real x(n)
!
!  If this is the first call, request a candidate for position 1.
!
  if ( indx == 0 ) then
    k = 1
    nstack = 0
    indx = 2
    return
  end if
!
!  Examine the stack.
!
  do
!
!  If there are candidates for position K, take the first available
!  one off the stack, and increment K.
!
!  This may cause K to reach the desired value of N, in which case
!  we need to signal the user that a complete set of candidates
!  is being returned.
!
    if ( ncan(k) > 0 ) then

      x(k) = stack(nstack)
      nstack = nstack - 1

      ncan(k) = ncan(k) - 1

      if ( k /= n ) then
        k = k + 1
        indx = 2
      else
        indx = 1
      end if

      exit
!
!  If there are no candidates for position K, then decrement K.
!  If K is still positive, repeat the examination of the stack.
!
    else

      k = k - 1

      if ( k <= 0 ) then
        indx = 3
        exit
      end if

    end if

  end do

  return
end
subroutine rvec_frac ( n, a, k, afrac )
!
!*******************************************************************************
!
!! RVEC_FRAC searches for the K-th smallest entry in an N-vector.
!
!
!  Discussion:
!
!    Hoare's algorithm is used.
!
!  Modified:
!
!    17 July 2000
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Input/output, real A(N).
!
!    On input, A is the array to search.
!
!    On output, the elements of A have been somewhat rearranged.
!
!    Input, integer K, the fractile to be sought.  If K = 1, the minimum
!    entry is sought.  If K = N, the maximum is sought.  Other values
!    of K search for the entry which is K-th in size.  K must be at
!    least 1, and no greater than N.
!
!    Output, real AFRAC, the value of the K-th fractile of A.
!
  implicit none
!
  integer n
!
  real a(n)
  real afrac
  integer i
  integer iryt
  integer j
  integer k
  integer left
  real w
  real x
!
  if ( n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RVEC_FRAC - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of N = ', n
    stop
  end if

  if ( k <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RVEC_FRAC - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal nonpositive value of K = ', k
    stop
  end if

  if ( k > n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RVEC_FRAC - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal K > N, K = ', k
    stop
  end if

  left = 1
  iryt = n

  do

    if ( left >= iryt ) then
      afrac = a(k)
      exit
    end if

    x = a(k)
    i = left
    j = iryt

    do

      if ( i > j ) then
        if ( j < k ) then
          left = i
        end if
        if ( k < i ) then
          iryt = j
        end if
        exit
      end if
!
!  Find I so that X< = A(I)
!
      do while ( a(i) < x )
        i = i + 1
      end do
!
!  Find J so that A(J) < =  X
!
      do while ( x < a(j) )
        j = j - 1
      end do

      if ( i <= j ) then
        call r_swap ( a(i), a(j) )
        i = i + 1
        j = j - 1
      end if

    end do

  end do

  return
end
subroutine rvec_identity ( n, a )
!
!*******************************************************************************
!
!! RVEC_IDENTITY sets a real vector to the identity vector A(I)=I.
!
!
!  Modified:
!
!    09 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, real A(N), the array to be initialized.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
!
  do i = 1, n
    a(i) = real ( i )
  end do

  return
end
subroutine rvec_print ( n, a, title )
!
!*******************************************************************************
!
!! RVEC_PRINT prints a real vector.
!
!
!  Modified:
!
!    16 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none
!
  integer n
!
  real a(n)
  integer i
  character ( len = * ) title
!
  if ( title /= ' ' ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  return
end
subroutine rvec_random ( n, a, alo, ahi )
!
!*******************************************************************************
!
!! RVEC_RANDOM returns a random real vector in a given range.
!
!
!  Modified:
!
!    01 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Output, real A(N), the vector of randomly chosen values.
!
!    Input, real ALO, AHI, the range allowed for the entries.
!
  implicit none
!
  integer n
!
  real a(n)
  real ahi
  real alo
  integer i
!
  do i = 1, n

    call r_random ( alo, ahi, a(i) )

  end do

  return
end
subroutine s_blank_delete ( s )
!
!*******************************************************************************
!
!! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
!
!
!  Comment:
!
!    All TAB characters are also removed.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
!
  character c
  integer iget
  integer iput
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  iput = 0

  do iget = 1, len ( s )

    c = s(iget:iget)

    if ( c /= ' ' .and. c /= TAB ) then
      iput = iput + 1
      s(iput:iput) = c
    end if

  end do

  s(iput+1:) = ' '

  return
end
subroutine s_blanks_delete ( s )
!
!*******************************************************************************
!
!! S_BLANKS_DELETE replaces consecutive blanks by one blank.
!
!
!  Discussion:
!
!    The remaining characters are left justified and right padded with blanks.
!    TAB characters are converted to spaces.
!
!  Modified:
!
!    26 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) S, the string to be transformed.
!
  implicit none
!
  integer i
  integer j
  character newchr
  character oldchr
  character ( len = * ) s
  character, parameter :: TAB = char ( 9 )
!
  j = 0
  newchr = ' '

  do i = 1, len ( s )

    oldchr = newchr
    newchr = s(i:i)

    if ( newchr == TAB ) then
      newchr = ' '
    end if

    s(i:i) = ' '

    if ( oldchr /= ' ' .or. newchr /= ' ' ) then
      j = j + 1
      s(j:j) = newchr
    end if

  end do

  return
end
function s_eqi ( s1, s2 )
!
!*******************************************************************************
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none
!
  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2
!
  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine schroeder ( n, s )
!
!*******************************************************************************
!
!! SCHROEDER generates the Schroeder numbers.
!
!
!  Definition:
!
!    The Schroeder number S(N) counts the number of ways to insert
!    parentheses into an expression of N items, with two or more items within
!    a parenthesis.
!
!    Note that the Catalan number C(N) counts the number of ways
!    to legally arrange a set of N left and N right parentheses.
!
!  Example:
!
!    N = 4
!
!    1234
!    12(34)
!    1(234)
!    1(2(34))
!    1(23)4
!    1((23)4)
!    (123)4
!    (12)34
!    (12)(34)
!    (1(23))4
!    ((12)3)4
!
!  First Values:
!
!           1
!           1
!           3
!          11
!          45
!         197
!         903
!        4279
!       20793
!      103049
!      518859
!     2646723
!    13648869
!    71039373
!
!  Formula:
!
!    S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
!    where P(N)(X) is the N-th Legendre polynomial.
!
!  Recursion:
!
!    S(1) = 1
!    S(2) = 1
!    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
!
!  Reference:
!
!    R P Stanley,
!    Hipparchus, Plutarch, Schroeder, and Hough,
!    American Mathematical Monthly,
!    Volume 104, Number 4, 1997, pages 344-350.
!
!    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
!    On the Second Number of Plutarch,
!    American Mathematical Monthly, May 1998, page 446.
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of Schroeder numbers desired.
!
!    Output, integer S(N), the Schroeder numbers.
!
  implicit none
!
  integer n
!
  integer i
  integer s(n)
!
  if ( n <= 0 ) then
    return
  end if

  s(1) = 1

  if ( n <= 1 ) then
    return
  end if

  s(2) = 1

  if ( n <= 2 ) then
    return
  end if

  do i = 3, n
    s(i) = ( ( 6 * i - 9 ) * s(i-1) - ( i - 3 ) * s(i-2) ) / i
  end do

  return
end
subroutine sort_heap_external ( n, indx, i, j, isgn )
!
!*******************************************************************************
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into linear order.
!
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    19 May 1999
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if I > J;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    ISGN => 0 means I is greater than or equal to J.
!
  implicit none
!
  integer i
  integer indx
  integer isgn
  integer j
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    n1 = n
    k = n / 2
    k1 = k
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i = i + 1
      end if

      j = k1
      k1 = i
      indx = - 1
      return

    end if

    if ( isgn > 0 ) then
      indx = 2
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        indx = 0
      else
        i = n1
        n1 = n1 - 1
        j = 1
        indx = 1
      end if

      return

    end if

    k = k - 1
    k1 = k
!
!  INDX > 0, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i = 2 * k1

    if ( i == n1 ) then
      j = k1
      k1 = i
      indx = - 1
      return
    else if ( i <= n1 ) then
      j = i + 1
      indx = - 2
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    indx = 0
  else
    i = n1
    n1 = n1 - 1
    j = 1
    indx = 1
  end if

  return
end
subroutine sub_by_size_next ( n, iarray, size, more )
!
!*******************************************************************************
!
!! SUB_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
!
!
!  Example:
!
!    N = 4:
!
!    1 2 3 4
!    1 2 3
!    1 2 4
!    1 3 4
!    1 3
!    1 4
!    2 3
!    1
!    2
!    3
!    (the empty set)
!
!  Discussion:
!
!    The subsets are returned in decreasing order of size, with the
!    empty set last.
!
!    For a given size K, the K subsets are returned in lexicographic order.
!
!  Modified:
!
!    02 December 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the set.
!
!    Output, integer IARRAY(N).  The entries IARRAY(1:SIZE) contain
!    the elements of the subset.  The elements are given in ascending
!    order.
!
!    Output, integer SIZE, the number of elements in the subset.
!
!    Input/output, logical MORE.  Set MORE = .FALSE. before first call
!    for a new sequence of subsets.  It then is set and remains
!    .TRUE. as long as the subset computed on this call is not the
!    final one.  When the final subset is computed, MORE is set to
!    .FALSE. as a signal that the computation is done.
!
  implicit none
!
  integer n
!
  integer iarray(n)
  logical more
  logical, save :: more2 = .false.
  integer size
!
  if ( .not. more ) then
    more = .true.
    more2 = .false.
    size = n
  else if ( .not. more2 ) then
    size = size - 1
  end if
!
!  Compute the next subset of size SIZE.
!
  if ( size > 0 ) then
    call ksub_next ( n, size, iarray, more2 )
  else if ( size == 0 ) then
    more = .false.
  end if

  return
end
subroutine sub_lex ( n, k, a, jmp, ndim )
!
!*******************************************************************************
!
!! SUB_LEX generates the subsets of a set of N elements, one at a time.
!
!
!  Discussion:
!
!    The subsets are generated in lexical order.  The routine can also be
!    forced to generate only those subsets whose size is no greater than
!    some user-specified maximum.
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
!    Input, integer N, the order of the main set from which subsets
!    are chosen.
!
!    Input/output, integer K.  On first call, the user must set K = 0 as
!    a startup signal to the program.  Thereafter, the routine returns
!    the size of the computed subset in K.  On the last return,
!    the empty set is returned and K is 0, which is a signal to
!    the user that the computation is complete.
!
!    Output, integer A(NDIM).  A(I) is the I-th element of the
!    subset, listed in increasing order, with 0's in entries
!    beyond entry K.
!
!    Input, logical JMP.  In the simplest case, set JMP = .FALSE. for
!    a normal computation.  But to jump over supersets of the input set,
!    set JMP = .TRUE..  Setting JMP=K.EQ.3 before every new call
!    will, for example, force all the subsets returned
!    to have cardinality 3 or less.
!
!    Input, integer NDIM, the allowed storage for A.  If NDIM < N,
!    JMP must be used to avoid creation of a subset too large to store in A.
!
  implicit none
!
  integer ndim
!
  integer a(ndim)
  integer is
  logical jmp
  integer k
  integer n
!
  if ( k == 0 ) then

    if ( jmp ) then
      return
    end if

    is = 0
    k = k + 1
    a(1) = 1

  else if ( a(k) /= n ) then

    is = a(k)

    if ( .not. jmp ) then
      k = k + 1
    end if

    a(k) = is + 1

  else

    k = k - 1

    if ( k /= 0 ) then
      a(k) = a(k) + 1
    end if

  end if

  return
end
subroutine sub_next ( n, a, more, ncard, iadd )
!
!*******************************************************************************
!
!! SUB_NEXT generates all subsets of a set of order N, one at a time.
!
!
!  Discussion:
!
!    It generates the subsets one at a time, by adding or subtracting
!    exactly one element on each step.
!
!    The user should set MORE = .FALSE. and the value of N before
!    the first call.  On return, the user may examine A which contains
!    the definition of the new subset, and must check .MORE., because
!    as soon as it is .FALSE. on return, all the subsets have been
!    generated and the user probably should cease calling.
!
!    The first set returned is the empty set.
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
!    Input, integer N, the order of the total set from which
!    subsets will be drawn.
!
!    Output, integer A(N).  On each return, the Gray code for the newly
!    generated subset.  A(I) = 0 if element I is in the subset, 1 otherwise.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but once
!    all the subsets have been generated, MORE will be
!    reset .FALSE. on return and you should stop calling the program.
!
!    Output, integer NCARD, the cardinality of the set returned,
!    which may be any value between 0 (the empty set) and N (the
!    whole set).
!
!    Output, integer IADD, the element which was added or removed to the
!    previous subset to generate the current one.  Exception:
!    the empty set is returned on the first call, and IADD is set to 0.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  integer iadd
  logical more
  integer ncard
!
!  First set returned is the empty set.
!
  if ( .not. more ) then

    a(1:n) = 0

    iadd = 0
    ncard = 0
    more = .true.

  else

    iadd = 1

    if ( mod ( ncard, 2 ) /= 0 ) then

      do

        iadd = iadd + 1
        if ( a(iadd-1) /= 0 ) then
          exit
        end if

      end do

    end if

    a(iadd) = 1 - a(iadd)
    ncard = ncard + 2 * a(iadd) - 1
!
!  Last set returned is the singleton A(N).
!
    if ( ncard == a(n) ) then
      more = .false.
    end if

  end if

  return
end
subroutine sub_random ( n, a )
!
!*******************************************************************************
!
!! SUB_RANDOM selects a random subset of an N-set.
!
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
!    01 December 2000
!
!  Parameters:
!
!    Input, integer N, the size of the full set.
!
!    Output, integer A(N).  A vector to hold the information about
!    the set chosen.  On return, if A(I) = 1, then
!    I is in the random subset, otherwise, A(I) = 0
!    and I is not in the random subset.
!
  implicit none
!
  integer n
!
  integer a(n)
  integer i
  real r
!
  do i = 1, n
    call i_random ( 0, 1, a(i) )
  end do

  return
end
subroutine sub_rank ( a, irank, k )
!
!*******************************************************************************
!
!! SUB_RANK computes the rank of a subset of an N set.
!
!
!  Discussion:
!
!    The routine accepts an array representing a subset of size K from a set
!    of size N, and returns the rank (or order) of that subset.  This
!    is the same order in which routine SUB_NEXT2 would produce that subset.
!    Note the value of N is not input, and is not, in fact,
!    needed.
!
!  Parameters:
!
!    Input, integer A(K), contains K distinct numbers between
!    1 and N, in order.
!
!    Output, integer IRANK, the rank of this subset.
!
!    Input, integer K, the number of elements in the subset.
!
  implicit none
!
  integer k
!
  integer a(k)
  integer i
  integer iprod
  integer irank
  integer j
!
  irank = 0

  do i = 1, k

    iprod = 1

    do j = i+1, a(i)-1
      iprod = iprod * j
    end do

    do j = 1, a(i)-i-1
      iprod = iprod / j
    end do

    if ( a(i) == 1 ) then
      iprod = 0
    end if

    irank = irank + iprod

  end do

  irank = irank + 1

  return
end
subroutine sub_unrank ( a, irank, k )
!
!*******************************************************************************
!
!! SUB_UNRANK returns the subset of a given rank.
!
!
!  Discussion:
!
!    SUB_UNRANK is given a rank and returns the corresponding subset of K
!    elements of a set of N elements.  It uses the same ranking that
!    SUB_NEXT2 uses to generate all the subsets one at a time.  Note that
!    the value of N itself is not input, nor is it needed.
!
!  Modified:
!
!    24 July 2000
!
!  Parameters:
!
!    Output, integer A(K), K distinct integers in order between
!    1 and N, which define the subset.
!
!    Input, integer IRANK, the rank of the desired subset.
!    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
!    subsets, so IRANK must be between 1 and that value.
!
!    Input, integer K, the number of elements in the subset.
!
  implicit none
!
  integer k
!
  integer a(k)
  integer i
  integer ii
  integer ip
  integer iprod
  integer irank
  integer jrank
!
  jrank = irank - 1

  do ii = 1, k

    i = k + 1 - ii
    ip = i - 1
    iprod = 1

    do

      ip = ip + 1

      if ( ip /= i ) then
        iprod = ( ip * iprod ) / ( ip - i )
      end if

      if ( iprod > jrank ) then
        exit
      end if

    end do

    if ( ip /= i ) then
      iprod = ( ( ip - i ) * iprod ) / ip
    end if

    jrank = jrank - iprod
    a(i) = ip

  end do

  return
end
subroutine thue_binary_next ( n, thue )
!
!*******************************************************************************
!
!! THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
!
!
!  Discussion:
!
!    Thue demonstrated that arbitrarily long sequences of 0's and
!    1's could be generated which had the "cubefree" property.  In
!    other words, for a given string S, there was no substring W
!    such that S contained "WWW".  In fact, a stronger result holds:
!    if "a" is the first letter of W, it is never the case that S
!    contains the substring "WWa". 
!
!    In this example, the digits allowed are binary, that is, just
!    "0" and "1".  The replacement rules are:
!
!    "0" -> "01"
!    "1" -> "10"
!
!    This routine produces the next binary Thue sequence in a given series.
!    However, the input sequence must be a Thue sequence in order for
!    us to guarantee that the output sequence will also have the
!    cubic nonrepetition property.
!
!    Also, enough space must be set aside in THUE to hold the
!    output sequence.  This will always be twice the input
!    value of N.
!
!  Modified:
!
!    05 November 2001
!
!  Parameters:
!
!    Input/output, integer N.  On input, the length of the input sequence.
!    On output, the length of the output sequence.
!
!    Input, integer THUE(N).  On input, the initial Thue sequence, and on
!    output, the result of applying the substitution rules once.
!
  implicit none
!
  integer n
!
  integer i
  integer n_out
  integer thue(*)
  integer thue_out(2*n)
!
  n_out = 0

  do i = 1, n

    if ( thue(i) == 0 ) then
      n_out = n_out + 1
      thue_out(n_out) = 0
      n_out = n_out + 1
      thue_out(n_out) = 1
    else if ( thue(i) == 1 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 0
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'THUE_BINARY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The input sequence contains a non-binary digit'
      write ( *, '(a,i6,a,i6)' ) '  THUE(', i, ') = ', thue(i)
      stop
    end if

  end do

  n = n_out
  thue(1:n) = thue_out(1:n)

  return
end
subroutine thue_ternary_next ( n, thue )
!
!*******************************************************************************
!
!! THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
!
!
!  Discussion:
!
!    Thue was interested in showing that there were arbitrarily long
!    sequences of digits which never displayed a pair of contiguous
!    repetitions of any length.  That is, there was no occurrence of 
!    "00" or "1010" or "121121", anywhere in the string.  This makes
!    the string "squarefree".
!
!    To do this, he demonstrated a way to start with a single digit,
!    and to repeatedly apply a series of transformation rules to each 
!    digit of the sequence, deriving nonrepeating sequences of ever 
!    greater length.
!
!    In this example, the digits allowed are ternary, that is, just
!    "0", "1" and "2".  The replacement rules are:
!
!    "0" -> "12"
!    "1" -> "102"
!    "2" -> "0"
!
!    This routine produces the next Thue sequence in a given series.
!    However, the input sequence must be a Thue sequence in order for
!    us to guarantee that the output sequence will also have the 
!    nonrepetition property.
!
!    Also, enough space must be set aside in THUE to hold the
!    output sequence.  This will never be more than 3 times the input
!    value of N.
!
!  Reference:
!
!    Brian Hayes,
!    Third Base,
!    American Scientist, 
!    Volume 89, Number 6, pages 490-494, November-December 2001.
!
!  Modified:
!
!    28 October 2001
!
!  Parameters:
!
!    Input/output, integer N.  On input, the length of the input sequence.
!    On output, the length of the output sequence.
!
!    Input, integer THUE(N).  On input, the initial Thue sequence, and on
!    output, the result of applying the substitution rules once.
!
  implicit none
!
  integer n
!
  integer i
  integer n_out
  integer thue(*)
  integer thue_out(3*n)
!
  n_out = 0

  do i = 1, n

    if ( thue(i) == 0 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 2
    else if ( thue(i) == 1 ) then
      n_out = n_out + 1
      thue_out(n_out) = 1
      n_out = n_out + 1
      thue_out(n_out) = 0
      n_out = n_out + 1
      thue_out(n_out) = 2
    else if ( thue(i) == 2 ) then
      n_out = n_out + 1
      thue_out(n_out) = 0
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'THUE_TERNARY_NEXT - Fatal error!'
      write ( *, '(a)' ) '  The input sequence contains a non-ternary digit'
      write ( *, '(a,i6,a,i6)' ) '  THUE(', i, ') = ', thue(i)
      stop
    end if

  end do

  n = n_out
  thue(1:n) = thue_out(1:n)

  return
end
!!$subroutine timestamp ( )
!!$!
!!$!*******************************************************************************
!!$!
!!$!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!!$!
!!$!
!!$!  Example:
!!$!
!!$!    May 31 2001   9:45:54.872 AM
!!$!
!!$!  Modified:
!!$!
!!$!    31 May 2001
!!$!
!!$!  Author:
!!$!
!!$!    John Burkardt
!!$!
!!$!  Parameters:
!!$!
!!$!    None
!!$!
!!$  implicit none
!!$!
!!$  character ( len = 8 ) ampm
!!$  integer d
!!$  character ( len = 8 ) date
!!$  integer h
!!$  integer m
!!$  integer mm
!!$  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
!!$    'January  ', 'February ', 'March    ', 'April    ', &
!!$    'May      ', 'June     ', 'July     ', 'August   ', &
!!$    'September', 'October  ', 'November ', 'December ' /)
!!$  integer n
!!$  integer s
!!$  character ( len = 10 )  time
!!$  integer values(8)
!!$  integer y
!!$  character ( len = 5 ) zone
!!$!
!!$  call date_and_time ( date, time, zone, values )
!!$
!!$  y = values(1)
!!$  m = values(2)
!!$  d = values(3)
!!$  h = values(5)
!!$  n = values(6)
!!$  s = values(7)
!!$  mm = values(8)
!!$
!!$  if ( h < 12 ) then
!!$    ampm = 'AM'
!!$  else if ( h == 12 ) then
!!$    if ( n == 0 .and. s == 0 ) then
!!$      ampm = 'Noon'
!!$    else
!!$      ampm = 'PM'
!!$    end if
!!$  else
!!$    h = h - 12
!!$    if ( h < 12 ) then
!!$      ampm = 'PM'
!!$    else if ( h == 12 ) then
!!$      if ( n == 0 .and. s == 0 ) then
!!$        ampm = 'Midnight'
!!$      else
!!$        ampm = 'AM'
!!$      end if
!!$    end if
!!$  end if
!!$
!!$  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!!$    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )
!!$
!!$  return
!!$end
subroutine triang ( n, izeta, p )
!
!*******************************************************************************
!
!! TRIANG renumbers elements in accordance with a partial ordering.
!
!
!  Discussion:
!
!    TRIANG is given a partially ordered set, where the partial ordering
!    is defined by a matrix IZETA, where element I is partially less than
!    or equal to element J if and only if IZETA(I,J) = 1.
!
!    TRIANG renumbers the elements with a permutation P so that if
!    element I is partially less than element J in the partial ordering,
!    then P(I) < P(J) in the usual, numerical ordering.
!
!    In other words, the elements are relabeled so that their labels
!    reflect their ordering.  This is equivalent to relabeling the
!    matrix so that, on unscrambling it, the matrix would be upper
!    triangular.
!
!    Calling IMAT_PERM or RMAT_PERM with P used for both the row
!    and column permutations applied to matrix IZETA will result in
!    an upper triangular matrix.
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
!    Input, integer N, the number of elements in the set.  The number
!    of rows and columns used in IZETA.
!
!    Input, integer IZETA(N,N), contains the description of the
!    partial ordering.  To avoid confusion with the true "less than"
!    relationship, symbolized by "<" or ".LT.", let us use
!    "I .PLT. J" to mean that I is less than J in the partial ordering.
!
!    The diagonal entries of IZETA are meaningless.
!
!    Consider any two distinct elements, I and J.  There are three
!    possibilities:
!
!    * I and J are not related in the partial ordering.  That is, neither
!      I .PLT. J nor J .PLT. I are true.  In that case,
!      set IZETA(I,J) = IZETA(J,I) = 0.
!
!    * I and J are related, and I .PLT. J.  In that case, set
!      IZETA(I,J) = 1, IZETA(J,I) = 0.
!
!    * I and J are related, and J .PLT. I.  In that case, set
!      IZETA(I,J) = 0, IZETA(J,I) = 1
!
!    Note that it should never happen that IZETA(I,J) and IZETA(J,I)
!    are both nonzero.  This would represent incorrect data.
!
!    Thus, the values of IZETA(I,J) are:
!      0, for diagonal elements, or for unrelated elements, or if J .PLT. I.
!      1, (or any nonzero value) if I .PLT. J.
!
!    Output, integer P(N), a permutation of the elements that reflects
!    their partial ordering.  P(I) is the new label of element I, with
!    the property that if IZETA(I,J) = 1, that is, I .PLT. J,
!    then P(I) < P(J) (in the usual ordering).
!
  implicit none
!
  integer n
!
  integer i
  integer ierror
  integer iq
  integer ir

  integer it
  integer izeta(n,n)
  integer l
  integer m
  integer p(n)
!
!  Make sure IZETA represents a partially ordered set.  In other words,
!  if IZETA(I,J) = 1, then IZETA(J,I) must NOT be 1.
!
  call pord_check ( n, izeta, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANG - Fatal error!'
    write ( *, '(a)' ) '  The matrix IZETA does not represent a'
    write ( *, '(a)' ) '  partial ordering.'
    stop
  end if

  m = 0
  l = 0
  p(1:n) = 0
!
!  Find the next value of M for which P(M) is 0.
!
10 continue

  do

    m = m + 1

    if ( p(m) == 0 ) then
      exit
    end if

    if ( m == n ) then
      return
    end if

  end do

  it = m + 1
  ir = m + 1

  do

    if ( ir <= n ) then

      if ( p(ir) == 0 .and. izeta(ir,m) /= 0 ) then
        p(ir) = m
        m = ir
        ir = it
      else
        ir = ir + 1
      end if

    else

      l = l + 1
      iq = p(m)
      p(m) = l

      if ( iq == 0 ) then
        if ( m == n ) then
          exit
        end if
        go to 10
      end if

      ir = m + 1
      m = iq

    end if

  end do

  return
end
subroutine tuple_next ( m, n, k, x )
!
!*******************************************************************************
!
!! TUPLE_NEXT computes the next element of a tuple space.
!
!
!  Discussion:
!
!    The elements are N vectors.  Each entry is constrained to lie
!    between 1 and M.  The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!  Examples:
!
!    N = 2, M = 3
!
!    K   X
!    1   1 1
!    2   1 2
!    3   1 3
!    4   2 1
!    5   2 2
!    6   2 3
!    7   3 1
!    8   3 2
!    9   3 3
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, integer M, the maximum entry.
!
!    Input, integer N, the number of components.
!
!    Input/output, integer K, counts the elements.
!    On first call, set K to 0.  Thereafter, K will indicate the
!    order of the element returned.  When there are no more elements,
!    K will be returned as 0.
!
!    Input/output, integer X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none
!
  integer n
!
  integer i
  integer k
  integer m
  integer x(n)
!
  if ( m < 1 ) then
    return
  end if

  if ( k <= 0 ) then

    x(1:n) = 1
    k = 1

  else

    k = k + 1
    i = n

    do

      if ( x(i) < m ) then
        x(i) = x(i) + 1
        exit
      end if

      x(i) = 1

      if ( i == 1 ) then
        k = 0
        exit
      end if

      i = i - 1

    end do

  end if

  return
end
subroutine tuple_next_ge ( m, n, k, x )
!
!*******************************************************************************
!
!! TUPLE_NEXT_GE computes the next "nondecreasing" element of a tuple space.
!
!
!  Discussion:
!
!    The elements are N vectors.  Each element is constrained to lie
!    between 1 and M, and to have components that are nondecreasing.
!    That is, for an element X, and any positive K,
!      X(I) <= X(I+K)
!
!    The elements are produced one at a time.
!    The first element is
!      (1,1,...,1)
!    and the last element is
!      (M,M,...,M)
!    Intermediate elements are produced in lexicographic order.
!
!  Examples:
!
!    N = 3, M = 3
!
!    K   X
!    1   1 1 1
!    2   1 1 2
!    3   1 1 3
!    4   1 2 2
!    5   1 2 3
!    6   1 3 3
!    7   2 2 2
!    8   2 2 3
!    9   2 3 3
!   10   3 3 3
!
!  Modified:
!
!    14 August 2001
!
!  Parameters:
!
!    Input, integer M, the maximum entry.
!
!    Input, integer N, the number of components.
!
!    Input/output, integer K, counts the elements.
!    On first call, set K to 0.  Thereafter, K will indicate the
!    order of the element returned.  When there are no more elements,
!    K will be returned as 0.
!
!    Input/output, integer X(N), on input the previous tuple.
!    On output, the next tuple.
!
  implicit none
!
  integer n
!
  integer i
  integer k
  integer m
  integer x(n)
!
  if ( m < 1 ) then
    return
  end if

  if ( k <= 0 ) then
    x(1:n) = 1
    k = 1
    return
  end if

  do i = n, 1, -1

    if ( x(i) < m ) then
      x(i) = x(i) + 1
      x(i+1:n) = x(i)
      k = k + 1
      return
    end if

  end do

  k = 0
  x(1:n) = 0

  return
end
subroutine tuple_next2 ( n, xmin, xmax, x, rank )
!
!*******************************************************************************
!
!! TUPLE_NEXT2 computes the next element of an integer tuple space.
!
!
!  Discussion:
!
!    The elements X are N vectors.
!
!    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
!
!    The elements are produced one at a time.
!
!    The first element is
!      (XMIN(1), XMIN(2), ..., XMIN(N)),
!    the second is (probably)
!      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
!    and the last element is
!      (XMAX(1), XMAX(2), ..., XMAX(N))
!
!    Intermediate elements are produced in a lexicographic order, with
!    the first index more important than the last, and the ordering of
!    values at a fixed index implicitly defined by the sign of
!    XMAX(I) - XMIN(I).
!
!  Examples:
!
!    N = 2,
!    XMIN = (/ 1, 10 /)
!    XMAX = (/ 3,  8 /)
!
!    RANK    X
!    ----  -----
!      1   1 10
!      2   1  9
!      3   1  8
!      4   2 10
!      5   2  9
!      6   2  8
!      7   3 10
!      8   3  9
!      9   3  8
!
!  Modified:
!
!    21 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components.
!
!    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
!    These values are minimum and maximum only in the sense of the lexicographic
!    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
!    than XMAX(I).
!
!    Input/output, integer X(N), on input the previous tuple.
!    On output, the next tuple.
!
!    Input/output, integer RANK, the rank of the item.  On first call,
!    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
!    there are no more items in the sequence.
!
  implicit none
!
  integer n
!
  integer i
  integer rank
  integer x(n)
  integer xmin(n)
  integer xmax(n)
!
  if ( rank < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank > product ( 1 + abs ( xmax(1:n) - xmin(1:n) ) ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
    write ( *, '(a,i6)' ) '  Illegal value of RANK = ', rank
    stop
  end if

  if ( rank == 0 ) then
    x(1:n) = xmin(1:n)
    rank = 1
    return
  end if

  rank = rank + 1
  i = n

  do

    if ( x(i) /= xmax(i) ) then
      x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
      exit
    end if

    x(i) = xmin(i)

    if ( i == 1 ) then
      rank = 0
      exit
    end if

    i = i - 1

  end do

  return
end
subroutine vec_next ( n, iarray, more, base )
!
!*******************************************************************************
!
!! VEC_NEXT generates all N-vectors of integers modulo a given base.
!
!
!  Examples:
!
!    N = 2, BASE = 3
!
!    0   0
!    0   1
!    0   2
!    1   0
!    1   1
!    1   2
!    2   0
!    2   1
!    2   2
!
!  Comment:
!
!    The vectors are produced in lexical order, starting with
!    (0,0,...,0), (0,0,...,1), ... through (BASE-1,BASE-1,...,BASE-1).
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input, integer N, the size of the vectors to be used.
!
!    Output, integer IARRAY(N).  On each return, IARRAY
!    will contain entries in the range 0 to BASE-1.
!
!    Input/output, logical MORE.  Set this variable .FALSE. before
!    the first call.  Normally, MORE will be returned .TRUE. but
!    once all the vectors have been generated, MORE will be
!    reset .FALSE. and you should stop calling the program.
!
!    Input, integer BASE, the base to be used.  BASE = 2 will
!    give vectors of 0's and 1's, for instance.
!
  implicit none
!
  integer n
!
  integer base
  integer i
  integer iarray(n)
  integer, save :: kount = 0
  integer, save :: last = 0
  logical more
  integer nn
!
  if ( .not. more ) then

    kount = 1
    last = base**n
    more = .true.
    iarray(1:n) = 0

  else

    kount = kount + 1

    if ( kount == last ) then
      more = .false.
    end if

    iarray(n) = iarray(n) + 1

    do i = 1, n

      nn = n - i

      if ( iarray(nn+1) < base ) then
        return
      end if

      iarray(nn+1) = 0

      if ( nn /= 0 ) then
        iarray(nn) = iarray(nn) + 1
      end if

    end do

  end if

  return
end
subroutine vec_next2 ( done, iactiv, iarray, idir, imax, n )
!
!*******************************************************************************
!
!! VEC_NEXT2 computes the elements of a product space.
!
!
!  Discussion:
!
!    The elements are produced one at a time.
!
!    This routine handles the case where the number of degrees of freedom may
!    differ from one component to the next.
!
!  Examples:
!
!    N = 2, IMAX = ( 2, 3 )
!
!    0 0
!    0 1
!    0 2
!    1 2
!    1 1
!    1 0
!
!  Comments:
!
!    A method similar to the Gray code is used, so that successive
!    elements returned by this routine differ by only a single element.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    15 April 1999
!
!  Parameters:
!
!    Input/output, logical DONE.  On the first call, the user must
!    set DONE to .TRUE.  This signals the program to initialize data.
!    On every return, if DONE is .FALSE., the program has computed
!    another entry, which is contained in IARRAY.  If DONE is .TRUE.,
!    then there are no more entries, and the program should not be
!    called for any more.
!
!    Workspace, integer IACTIV(N).
!
!    Input/output, integer IARRAY(N).  On the first call, the input value
!    of IARRAY doesn't matter.  Thereafter, it should be the same as
!    its output value from the previous call.  On output, if DONE
!    is .FALSE., then IARRAY contains the next element of the space.
!
!    Workspace, integer IDIR(N).
!
!    Input, integer IMAX(N), contains the number of degrees of
!    freedom of each component.  The output values of IARRAY will
!    satisfy 0 <= IARRAY(I) < IMAX(I).
!
!    Input, integer N, the number of components.
!
  implicit none
!
  integer n
!
  logical done
  integer i
  integer iactiv(n)
  integer iarray(n)
  integer idir(n)
  integer imax(n)
  integer ip
  integer maxact
!
!  The user is calling for the first time.
!
  if ( done ) then

    done = .false.

    do i = 1, n

      iarray(i) = 0
      idir(i) = 1

      if ( imax(i) < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VEC_NEXT2 - Warning!'
        write ( *, '(a,i6)' ) '  For index I = ',i
        write ( *, '(a,i6)' ) '  the nonpositive value of IMAX(I) = ',  imax(i)
        write ( *, '(a)' ) '  which was reset to 1!'
        imax(i) = 1
        iactiv(i) = 0
      else if ( imax(i) == 1 ) then
        iactiv(i) = 0
      else
        iactiv(i) = 1
      end if

    end do
!
!  The user is asking for the next vector.
!
  else
!
!  Find the maximum active index.
!
    maxact = 0

    do i = 1, n
      if ( iactiv(i) /= 0 ) then
        maxact = i
      end if
    end do

    if ( maxact /= 0 ) then

      done = .false.
      ip = maxact
      iarray(ip) = iarray(ip) + idir(ip)

      if ( iarray(ip) == imax(ip)-1 .or. iarray(ip) == 0 ) then

        idir(ip) = - idir(ip)
        iactiv(ip) = 0

      end if

      do i = ip+1, n
        if ( imax(i) > 1 ) then
          iactiv(i) = i
        end if
      end do

    else
      done = .true.
    end if

  end if

  return
end
subroutine vec_random ( n, iarray, base )
!
!*****************************************************************************
!
!! VEC_RANDOM selects a random N-vector of integers modulo a given base.
!
!
!  Modified:
!
!    15 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the size of the vector to be generated.
!
!    Output, integer IARRAY(N), a list of N random values between
!    0 and BASE-1.
!
!    Input, integer BASE, the base to be used.
!
  implicit none
!
  integer n
!
  integer base
  integer i
  integer iarray(n)
  integer ival
!
  do i = 1, n
    call i_random ( 0, base-1, ival )
    iarray(i) = ival
  end do

  return
end
subroutine vec_rank ( iarray, imax, n, rank )
!
!*******************************************************************************
!
!! VEC_RANK computes the rank of a product space element.
!
!
!  Discussion:
!
!    The rank applies only to the elements as produced by the routine
!    VEC_NEXT2.
!
!  Examples:
!
!    N = 2, IMAX = ( 2, 3 ), IARRAY = ( 1, 2 ),
!
!    RANK = 4.
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    27 July 1998
!
!  Parameters:
!
!    Input, integer IARRAY(N), the product space element, with the
!    property that 0 <= IARRAY(I) < IMAX(I) for each entry I.
!
!    Input, integer IMAX(N), contains the number of degrees of
!    freedom of each component.  The output values of IARRAY will
!    satisfy 0 <= IARRAY(I) < IMAX(I).
!
!    Input, integer N, the number of components.
!
!    Output, integer RANK, the rank, or order, of the element in
!    the list of all elements.  The rank count begins at 1.
!
  implicit none
!
  integer n
!
  integer c
  integer i
  integer iarray(n)
  integer imax(n)
  integer rank
!
  rank = 0

  do i = 1, n

    if ( mod ( rank, 2 ) == 1 ) then
      c = imax(i) - iarray(i) - 1
    else
      c = iarray(i)
    end if

    rank = imax(i) * rank + c

  end do

  rank = rank + 1

  return
end
subroutine vec_unrank ( iarray, imax, n, rank )
!
!*******************************************************************************
!
!! VEC_UNRANK computes the product space element of a given rank.
!
!
!  Discussion:
!
!    The rank applies only to the elements as produced by the routine
!    VEC_NEXT2.
!
!  Examples:
!
!    N = 2, IMAX = ( 2, 3 ), RANK = 4.
!
!    IARRAY = ( 1, 2 ).
!
!  Reference:
!
!    Dennis Stanton and Dennis White,
!    Constructive Combinatorics,
!    Springer Verlag, New York, 1986.
!
!  Modified:
!
!    27 July 1998
!
!  Parameters:
!
!    Output, integer IARRAY(N), the product space element of the given rank.
!
!    Input, integer IMAX(N), contains the number of degrees of
!    freedom of each component.  The output values of IARRAY will
!    satisfy 0 <= IARRAY(I) < IMAX(I).
!
!    Input, integer N, the number of components.
!
!    Input, integer RANK, the desired rank, or order, of the element in
!    the list of all elements.  The rank count begins at 1 and extends
!    to MAXRANK = Product ( I = 1 to N ) IMAX(I).
!
  implicit none
!
  integer n
!
  integer i
  integer iarray(n)
  integer imax(n)
  integer rank
  integer s
!
  s = rank - 1

  do i = n, 1, -1

    iarray(i) = mod ( s, imax(i) )
    s = s / imax(i)

    if ( mod ( s, 2 ) == 1 ) then
      iarray(i) = imax(i) - iarray(i) - 1
    end if

  end do

  return
end
subroutine ytb_enum ( n, ytb_num )
!
!*******************************************************************************
!
!! YTB_ENUM enumerates the Young tableau of size N.
!
!
!  Discussion:
!
!    If A(N) is the number of Young tableau of size N, then A(1) = 1,
!    A(2) = 2, and
!
!    A(N) = A(N-1) + (N-1) * A(N-2).
!
!  Modified:
!
!    05 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer which is to be partitioned.
!
!    Output, integer YTB_NUM, the number of Young tableau of N.
!
  implicit none
!
  integer a1
  integer a2
  integer a3
  integer i
  integer n
  integer ytb_num
!
  if ( n <= 0 ) then
    ytb_num = 0
  else if ( n == 1 ) then
    ytb_num = 1
  else if ( n == 2 ) then
    ytb_num = 2
  else
    a2 = 1
    a3 = 2
    do i = 3, n
      a1 = a2
      a2 = a3
      a3 = a2 + ( i - 1 ) * a1
    end do
    ytb_num = a3
  end if

  return
end
subroutine ytb_next ( n, lambda, iarray, more )
!
!*******************************************************************************
!
!! YTB_NEXT computes the next Young tableau for a given shape.
!
!
!  Discussion:
!
!    When the routine is called with MORE = .FALSE. (the first time), and
!    if LAMBDA on this call has M parts, with M<N, then the user
!    must also make sure that LAMBDA(M+1) = 0.
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
!    21 March 2001
!
!  Parameters:
!
!    Input, integer N, the integer which is to be partitioned.
!
!    Output, integer LAMBDA(N).  LAMBDA(I) is the I-th part of the partition.
!
!    Output, integer IARRAY(N).  IARRAY(I) is the row containing I
!    in the output tableau.
!
!    Input/output, logical MORE.  Set MORE .FALSE. before first call.
!    It is reset to .TRUE. as the program returns a new tableau
!    on each call, until the last tableau is computed, when
!    the program also sets MORE = .FALSE.
!
  implicit none
!
  integer n
!
  integer i
  integer ii
  integer ir
  integer it
  integer j
  integer l
  logical more
  integer lambda(n)
  integer iarray(n)
  integer isave
!
  it = n

  if ( more ) then

    lambda(1) = 1
    lambda(2:n) = 0

    isave = 0

    do i = 2, n

      lambda(iarray(i)) = lambda(iarray(i)) + 1

      if ( iarray(i) < iarray(i-1) ) then
        isave = i
        exit
      end if

    end do

    if ( isave == 0 ) then
      more = .false.
      return
    end if

    it = lambda(1+iarray(isave))

    do ii = 1, n

      i = n + 1 - ii

      if ( lambda(i) == it ) then
        iarray(isave) = i
        lambda(i) = lambda(i) - 1
        it = isave-1
        exit
      end if

    end do

  end if

  l = 1
  ir = 1

  do

    if ( ir > n ) then
      exit
    end if

    if ( lambda(ir) /= 0 ) then
      iarray(l) = ir
      lambda(ir) = lambda(ir) - 1
      l = l + 1
      ir = ir + 1
      cycle
    end if

    if ( l > it ) then
      exit
    end if

    ir = 1

  end do

  if ( n == 1 ) then
    more = .false.
    return
  end if

  do j = 2, n
    if ( iarray(j) < iarray(j-1) ) then
      more = .true.
      return
    end if
  end do

  more = .false.

  return
end
subroutine ytb_print ( n, iarray )
!
!*******************************************************************************
!
!! YTB_PRINT prints a Young tableau.
!
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the integer that is partitioned.
!
!    Input, integer IARRAY(N), describes the Young tableau.
!    IARRAY(I) is the row of the tableau on which I occurs.
!
  implicit none
!
  integer n
!
  integer iarray(n)
  integer j
  integer jarray(n)
  integer row_i
  integer row_length
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Young tableau:'
  write ( *, '(a)' ) ' '

  row_i = 0

  do

    row_i = row_i + 1

    row_length = 0

    do j = 1, n

      if ( iarray(j) == row_i ) then
        row_length = row_length + 1
        jarray(row_length) = j
      end if

    end do

    if ( row_length <= 0 ) then
      exit
    end if

    write ( *, '(20i4)' ) jarray(1:row_length)

  end do

  return
end
subroutine ytb_random ( n, lambda, iarray )
!
!*******************************************************************************
!
!! YTB_RANDOM selects a random Young tableau of a given shape.
!
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
!    01 December 2000
!
!  Parameters:
!
!    Input, integer N, the integer which has been partitioned.
!
!    Input, integer LAMBDA(N).  N = LAMBDA(1)+LAMBDA(2)+...is the
!    partition of N.
!
!    Output, integer IARRAY(N).  The vector describing the Young tableau.
!
  implicit none
!
  integer n
  real r
!
  integer i
  integer iarray(n)
  integer ih
  integer j
  integer l
  integer lambda(n)
  integer m
!
  iarray(1:n) = 0

  i = 0
  l = 0

  do

    i = i + 1
    do j = 1, lambda(i)
      iarray(j) = iarray(j) + 1
      l = l + 1
    end do

    if ( l >= n ) then
      exit
    end if

  end do

  do m = 1, n

    do

      call i_random ( 1, iarray(1), i )
      call i_random ( 1, lambda(1), j )

      if ( i <= iarray(j) .and. j <= lambda(i) ) then
        exit
      end if

    end do

    do

      ih = iarray(j) + lambda(i) - i - j

      if ( ih == 0 ) then
        exit
      end if

      call i_random ( 1, ih, l )

      if ( l <= lambda(i)-j ) then
        j = j + l
      else
        i = l - lambda(i) + i + j
      end if

    end do

    lambda(i) = lambda(i) - 1
    iarray(j) = iarray(j) - 1
    iarray(n+1-m) = i

  end do

  do i = 1, n
    lambda(iarray(i)) = lambda(iarray(i)) + 1
  end do

  return
end
