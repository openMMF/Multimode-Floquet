DOUBLE PRECISION FUNCTION  D_H(N_SITES,N_PARTICLES,stats) 

  !Calculates the number of states of N_particles bosos in N_sites lattice sites
  ! Wikipedia:
  ! D_B  = (N_particles + N_sites - 1)!)/(N_particles! * (L-1)!)))
  ! D_F  = N_sites!/(N_particles! (N_sites-N_particles)!)))
  IMPLICIT NONE

  INTEGER, INTENT(IN):: N_SITES,N_PARTICLES
  !DOUBLE PRECISION D_H

  INTEGER          :: N_,K_,i
  DOUBLE PRECISION :: Num,Den


  IF(PRESENT(stats)) THEN
     SELECT CASE (stats)
     CASE("B")

        N_ = N_SITES+N_PARTICLES-1;
        K_ = N_PARTICLES;

        Num = 1.0;
        do while (n_.gt.0) then
           Num = Num*n_;
           n_ = n_-1;
        end do
        Den = 1.0;
        do while (k_.gt.0) then !i = 1,k_-1
           Den = Den*k_;
           k_ = k_-1;
        end do

     CASE("F")
        N_ = N_SITES
        Num = 1.0;
        do while (n_.gt.0) then
           Num = Num*n_
           n_ = n_-1;
        end do
        k_ = N_SITES - N_PARTICLES
        Den = 1.0;
        do while (k_.gt.0)
           Den = Den*n_;
           k_ = k_-1;
        end do
        k_ = N_PARTICLES
        do while (k_.gt.0) then
           Den = Den*k_;
           k_ = k_-1;
        end do

     CASE DEFAULT
        D_H = N_SITEX
     END SELECT

     D_H = Num/Den;

  ELSE

     D_H = N_SITES

  END IF

END FUNCTION D_H
