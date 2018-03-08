    ! ***************************************************************
    ! * Given an N x N matrix A, this routine replaces it by the LU *
    ! * decomposition of a rowwise permutation of itself. A and N *
    ! * are input. INDX is an output vector which records the row *
    ! * permutation effected by the partial pivoting; D is output *
    ! * as -1 or 1, depending on whether the number of row inter- *
    ! * changes was even or odd, respectively. This routine is used *
    ! * in combination with LUBKSB to solve linear equations or to *
    ! * invert a matrix. Return code is 1, if matrix is singular. *
    ! ***************************************************************
     Subroutine DLUDCMP(A,N,INDX,D,CODE)
     IMPLICIT NONE
     integer, parameter :: nmax = 10000
     real*8, parameter :: tiny = 1.d-30

     real*8, intent(inout), dimension(N,N) :: A
     integer, intent(in) :: N
     integer, intent(out) :: D, CODE
     integer, intent(out), dimension(N) :: INDX
     !f2py depend(N) A, indx

     real*8 :: AMAX, DUM, SUMM, VV(NMAX)
     INTEGER :: i, j, k, imax

     D=1; CODE=0

     DO I=1,N
       AMAX=0.d0
       DO J=1,N
         IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
       END DO ! j loop
       IF(AMAX.EQ.0.d0) THEN
         CODE = 1
         write(*,*) i, amax
         RETURN

       END IF
       VV(I) = 1.d0 / AMAX
     END DO ! i loop

     DO J=1,N
       DO I=1,J-1
         SUMM = A(I,J)
         DO K=1,I-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
         SUMM = A(I,J)
         DO K=1,J-1
           SUMM = SUMM - A(I,K)*A(K,J)
         END DO ! k loop
         A(I,J) = SUMM
         DUM = VV(I)*ABS(SUMM)
         IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
         END IF
       END DO ! i loop
       
       IF(J.NE.IMAX) THEN
         DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         END DO ! k loop
         D = -D
         VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(A(J,J).EQ.0.d0) THEN
         write(*,*) 'Warning: Newton direction is not precise', j
         A(J,J) = TINY
       END IF

       IF(J.NE.N) THEN
         DUM = 1.d0 / A(J,J)
         DO I=J+1,N
           A(I,J) = A(I,J)*DUM
         END DO ! i loop
       END IF
     END DO ! j loop

     RETURN
END subroutine DLUDCMP


! ******************************************************************
! * Solves the set of N linear equations A . X = B. Here A is *
! * input, not as the matrix A but rather as its LU decomposition, *
! * determined by the routine LUDCMP. INDX is input as the permuta-*
! * tion vector returned by LUDCMP. B is input as the right-hand *
! * side vector B, and returns with the solution vector X. A, N and*
! * INDX are not modified by this routine and can be used for suc- *
! * cessive calls with different right-hand sides. This routine is *
! * also efficient for plain matrix inversion. *
! ******************************************************************
 Subroutine DLUBKSB(A, N, INDX, B)
 integer, intent(in) :: N
 real*8, intent(in), dimension(N,N) :: A
 integer, intent(in), dimension(N) :: INDX
 real*8, intent(inout), dimension(N) :: B
 !f2py depend(N) A, INDX, B

 real*8 SUMM

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUMM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUMM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUMM
 END DO ! i loop

 DO I=N,1,-1
   SUMM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUMM / A(I,I)
 END DO ! i loop

 RETURN
END subroutine DLUBKSB
