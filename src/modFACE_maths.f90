 module modFACE_maths
    use modFACE_precision
    use modFACE_header
     use modFACE_error
    implicit none
    contains
    subroutine solve_linear_system(a,b,n)
    integer,intent(in):: n
    real(DP):: a(n,n)
    real(DP),intent(inout):: b(n)
    integer :: indx(n), d, code
    !     Externals from the LAPACK library
      external dgesv
!       DGESV computes the solution to a real system of linear equations
!    A * X = B,
! where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!
! The LU decomposition with partial pivoting and row interchanges is
! used to factor A as
!    A = P * L * U,
! where P is a permutation matrix, L is unit lower triangular, and U is
! upper triangular.  The factored form of A is then used to solve the
! system of equations A * X = B.
!Parameters
![in]    N
!          N is INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
![in]    NRHS
!          NRHS is INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
![in,out]    A
!          A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
![in]    LDA
!          LDA is INTEGER
 !         The leading dimension of the array A.  LDA >= max(1,N).
![out]   IPIV
!          IPIV is INTEGER array, dimension (N)
 !         The pivot indices that define the permutation matrix P;
 !         row i of the matrix was interchanged with row IPIV(i).
![in,out]    B
 !         B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
![in]    LDB
!          LDB is INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
![out]   INFO
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!                has been completed, but the factor U is exactly
!                singular, so the solution could not be computed.
    !call DLUDCMP(a,n,indx,d,code)
    call dgesv(n, 1, a, n, indx, b, n, code)
        if (code .ne. 0) then
            call face_error("Singular Jacobian! code=",code)
        endif
    end subroutine solve_linear_system


subroutine which_eq(j_eq)

        integer,intent(in):: j_eq
        integer j_idx,k_idx

        if (solve_heat_eq .eq. "yes") then

            if (j_eq.ge.neq-ngrd+1.and.j_eq.le.neq) then

                j_idx=ngrd+1-(neq-j_eq)
                write(iout,*) "eq: heat eq at j=",j_idx

            elseif (j_eq.lt.neq-ngrd+1) then


                k_idx=floor(real(j_eq)/real(ngrd+3))
                j_idx=j_eq-k_idx*(ngrd+3)+1
                if  (j_idx.eq.1) then
                    write(iout,*) "eq: spc k = ",k_idx," at left surface "
                elseif  (j_idx.eq.ngrd+3) then
                    write(iout,*) "eq: spc k = ",k_idx," at right surface"
                elseif (j_idx.lt.ngrd+3.and.j_idx.gt.1) then
                    write(iout,*) "eq: volume spc k = ",k_idx," at j = ",j_idx-1
                else
                    call face_error("cannot find the equation ",j_eq)
                endif
            endif

            elseif (solve_heat_eq .eq. "no") then

                k_idx=floor(real(j_eq)/real(ngrd+3))
                j_idx=j_eq-k_idx*(ngrd+3)+1
                if  (j_idx.eq.1) then
                    write(iout,*) "eq: spc k = ",k_idx," at left surface "
                elseif  (j_idx.eq.ngrd+3) then
                    write(iout,*) "eq: spc k = ",k_idx," at right surface"
                elseif (j_idx.lt.ngrd+3.and.j_idx.gt.1) then
                    write(iout,*) "eq: volume spc k = ",k_idx," at j = ",j_idx-1
                else
                    call face_error("cannot find the equation ;", j_eq,"/",neq,"; jdix=",j_idx)
                endif
            else

                call face_error("j_eq too large j_eq=",j_eq,"; neq= ",neq)
            endif
    end subroutine which_eq
   end module modFACE_maths
