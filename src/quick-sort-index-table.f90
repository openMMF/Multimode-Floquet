!!$PROGRAM SORTING
!!$
!!$  INTEGER, PARAMETER :: N = 4
!!$  DOUBLE PRECISION   :: V(N)
!!$  INTEGER            :: INDEX_T(N)
!!$
!!$  DO I = 1,N
!!$     INDEX_T(I) = I
!!$     V(I) = 1.0*ran()
!!$     write(*,*) V(i)
!!$  END DO
!!$
!!$  
!!$  CALL QUICK_SORT_I_T(v,INDEX_T,N)
!!$
!!$  write(*,*) 
!!$  DO I = 1,N
!!$     write(*,*) V(i),INDEX_T(i),V(INDEX_T(i))
!!$  END DO
!!$END PROGRAM

SUBROUTINE QUICK_SORT_I_T(v,index_t,N)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  
  !INTEGER, DIMENSION(N),INTENT(INOUT) :: v
  DOUBLE PRECISION, DIMENSION(N),INTENT(INOUT) :: v
  INTEGER, DIMENSION(N),INTENT(INOUT) :: index_t

  INTEGER, PARAMETER :: NN=2500, NSTACK=500
  real :: a, cpu
  integer :: trial
  INTEGER :: k,i,j,jstack,l,r,istack(NSTACK),indice
  
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then
        do j=l+1,r
           indice = index_t(j)
           a=v(j)
           do i=j-1,l,-1
              if (v(index_t(i)) <= a) exit
              index_t(i+1) = index_t(i)
           end do
           index_t(i+1) = indice
        end do
        if (jstack == 0) exit
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+r)/2
        call swap(index_t(k),index_t(l+1))
        if (v(index_t(l)) > v(index_t(r))) call swap(index_t(l),index_t(r))
        if (v(index_t(l+1))>v(index_t(r))) call swap(index_t(l+1),index_t(r))        
        if (v(index_t(l))>v(index_t(l+1))) call swap(index_t(l),index_t(l+1))
        
        i=l+1
        j=r
        indice = index_t(l+1)
        a=v(indice)
        do
           do
              i=i+1
              if (v(index_t(i)) >= a) exit
           end do
           do
              j=j-1
              if (v(index_t(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index_t(i),index_t(j))
        end do
        index_t(l+1)=index_t(j)
        index_t(j)=indice
        jstack=jstack+2
        if (jstack > NSTACK) then
           print *, "NSTACK too small", jstack
        end if
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
  
END SUBROUTINE QUICK_SORT_I_T

SUBROUTINE QUICK_SORT_INTEGERS(v,index_t,N)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, DIMENSION(N),INTENT(INOUT) :: v
  INTEGER, DIMENSION(N),INTENT(INOUT) :: index_t

  INTEGER, PARAMETER :: NN=10000, NSTACK=8000
  real :: a, cpu
  integer :: trial
  INTEGER :: k,i,j,jstack,l,r,istack(NSTACK),indice

  DO j=1,N
     index_t(j) = j
  END DO
  
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then
        do j=l+1,r
           indice = index_t(j)
           a=v(j)
           do i=j-1,l,-1
              if (v(index_t(i)) <= a) exit
              index_t(i+1) = index_t(i)
           end do
           index_t(i+1) = indice
        end do
        if (jstack == 0) exit
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+r)/2
        call swap(index_t(k),index_t(l+1))
        if (v(index_t(l)) > v(index_t(r))) call swap(index_t(l),index_t(r))
        if (v(index_t(l+1))>v(index_t(r))) call swap(index_t(l+1),index_t(r))        
        if (v(index_t(l))>v(index_t(l+1))) call swap(index_t(l),index_t(l+1))
        
        i=l+1
        j=r
        indice = index_t(l+1)
        a=v(indice)
        do
           do
              i=i+1
              if (v(index_t(i)) >= a) exit
           end do
           do
              j=j-1
              if (v(index_t(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index_t(i),index_t(j))
        end do
        index_t(l+1)=index_t(j)
        index_t(j)=indice
        jstack=jstack+2
        if (jstack > NSTACK) then
           print *, "NSTACK too small", jstack
        end if
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
  
END SUBROUTINE QUICK_SORT_INTEGERS

subroutine swap(a,b)
  !real, intent(inout) :: a, b
  !real :: c
  integer, intent(inout) :: a,b
  integer :: c
  c = a; a = b; b = c
end subroutine swap

