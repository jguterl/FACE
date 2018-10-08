
module modFACE_solverf90
    use modFACE_precision
    use modFACE_header
    use modFACE_functions
    use modFACE_output
    use modFACE_error
    use modFACE_compute
    use modFACE_maths
    implicit none
contains
    subroutine jac(u,f,fdot,norm)
        integer i, j
        real(DP):: u(neq), um(neq), up(neq)
        real(DP):: f(neq), fm(neq), fp(neq)
        real(DP):: fdot(neq,neq)
        real(DP):: eps, ddu, uud, amax, norm
        real(DP):: ftran
        parameter (eps=1.d-3)
        parameter (ftran=1.d-0)

        do i=1,neq
            um(i)=u(i)
            up(i)=u(i)
        enddo

        do i=1,neq
            ddu=u(i)*eps
            um(i)=u(i)-ddu
            up(i)=u(i)+ddu
            call compute_f(um,fm)
            call compute_f(up,fp)
            uud=1.d0/(ddu+ddu)
            do j=1,neq
                fdot(j,i)=(fp(j)-fm(j))*uud
            enddo
            um(i)=u(i)
            up(i)=u(i)
        enddo
        !     --- pseudo-transient continuation ---
        if (steady_state .eq. "yes") then
            do i=1,neq
                fdot(i,i)=fdot(i,i)-ftran*norm
            enddo
        endif
        !     --- check for singular Jacobian ---
        do i=1,neq
            amax=0.d0
            do j=1,neq
                if (abs(fdot(i,j)) .gt. amax) amax=abs(fdot(i,j))
            enddo
            if (amax .eq. 0.d0) then
                write (iout,*) '***warning: underflow in Jacobian raw ', i
                call which_eq(i)
                fdot(i,i)=-1.d+99
            endif
        enddo

    end subroutine jac

    subroutine build_vector(u,du)

        real(DP):: u(:), du(:)
        integer i,k,j
        i=0

        do k=1,nspc
            i=i+1
            u (i)=dsrfl(ndt-1,k)
            du(i)=0.d0
            do j=0,ngrd
                i=i+1
                u (i)=dens(ndt-1,j,k)
                du(i)=0.d0
            enddo
            i=i+1
            u (i)=dsrfr(ndt-1,k)
            du(i)=0.d0
        enddo
        if (solve_heat_eq .eq. "yes") then
            do j=0,ngrd
                i=i+1
                u (i)=temp(ndt-1,j)
                du(i)=0.d0
            enddo
        endif
        if (i.ne.neq) then
            call face_error("mismatch in vector size: neq=",neq, 'size(vec)=',i)

        endif
    end subroutine build_vector

    subroutine dsolve(du,u,f,fdot)

        real(DP)::  du(neq), u(neq), f(neq), fdot(neq,neq)
        real(DP):: a(neq,neq), b(neq)
        integer indx(neq)
        integer n, d, code
        integer i, j

        n=neq
        do i=1,neq
            b(i)=-f(i)
            do j=1,neq
                a(i,j)=fdot(i,j)
            enddo
        enddo
        ! solve A.X=B where A (N,N), B(N) X(N). Result is returned as b
        call solve_linear_system(a,b,n)

        do i=1,neq
            if (u(i) .gt. 1.d-10) then
                du(i)=b(i)
            else
                if (b(i) .gt. 0.d0) then
                    du(i)=b(i)
                else
                    du(i)=0.d0
                endif
            endif
        enddo

    !      do i=1,neq
    !!!       c(i)=0.d0
    !!!       do j=1,neq
    !!!        c(i)=c(i)+fdot(i,j)*du(j)
    !!!       enddo
    !!!       c(i)=c(i)+f(i)
    !!!       if (abs(c(i)) .gt. 1.d-6) then
    !!!        write (*,*) '***warning: loss of direction precision ',
    !!!     +              i, c(i)
    !!!       endif
    !!!      enddo

    end subroutine dsolve
end module modFACE_solverf90
