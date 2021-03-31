
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      INTEGER I
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION D
      DOUBLE PRECISION PI

      D = PAR(1)

      PI = 4 * atan(1.0_16);

      DO I = 2,NDIM-1
            F(I) = D*(U(I+1)-2*U(I)+U(I-1)) - U(I)*(1 - U(I)*U(I))
      END DO

      ! Dirichlet BCs
      F(1)    = D*(U(2) - 2*U(1) + -1) - U(1)*(1 - U(1)*U(1))
      F(NDIM) = D*(-1 - 2*U(NDIM) + U(NDIM-1)) - U(NDIM)*(1 - U(NDIM)*U(NDIM))

END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values

      IMPLICIT NONE
      INTEGER NDIM
      INTEGER OFFSET
      INTEGER I
      DOUBLE PRECISION U(NDIM), PAR(*), T
      DOUBLE PRECISION W
      DOUBLE PRECISION PI

      PI = 4 * atan(1.0_16);

! Initialize the equation parameters
      PAR(1)=0

! Initialize the solution, kink-antikink

      DO I = 1,NDIM
            U(I) = -1
      END DO

      OFFSET = 1
      DO I = NDIM/2 - OFFSET, NDIM/2+OFFSET
            U(I) = 1
      END DO

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
