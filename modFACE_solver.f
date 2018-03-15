
      module modFACE_solver
       use modFACE_precision
      use modFACE_header
      use modFACE_functions
      use modFACE_output
      use modFACE_error
       use modFACE_compute
       use modFACE_solverf90
       implicit none

       contains


            subroutine newton_solver
      integer i,ii
      integer iter_solver_max, cntl, idx

      real(DP) eps, udspl, fdspl, gdspl, fstp
      real(DP):: u(neq), du(neq), f(neq), fdot(neq,neq)
      real(DP)  unew(neq), unorm(neq)
      real(DP) norm, normnew

      parameter (iter_solver_max=100)
      parameter (eps=3.d-3, fstp=1.d-1)
      parameter (udspl=9.d-1, fdspl=9.d0, gdspl=1.d-3)

       iter_solver=0

      call build_vector(u,du)

      norm=compute_fnorm(u)

      call compute_f(u,f)

      if (verbose_debug) then
      do ii=1,neq
      write(iout,*) 'i=',ii,'u=',u(ii),'f(i)=',f(ii)
      enddo
      endif

5199  call jac(u,f,fdot,norm)
      call dsolve(du,u,f,fdot)
c
c     --- newton step reduction ---
      norm=compute_fnorm(u)
      cntl=0
5198  idx=0
      do i=1,neq
       unew (i)=u(i)+du(i)
       unorm(i)=du(i)/u(i)
c       if (abs(unorm(i)) .lt. 1.d-30) then
c        write (*,*) '***warning: displacement vector is too small at '
c     +              , i
c       endif
       if ((-unorm(i) .gt. udspl) .or.
     +     ( unorm(i) .gt. 1.d0/(1.d0-udspl))) then
        idx=i
       endif
      enddo
      if (idx .ne. 0) then
       do i=1,neq
        du(i)=fstp*du(i)
       enddo
       cntl=cntl+1
       if (verbose_step)write (*,*) cntl, ' var step reduction',
     +idx
       goto 5198
      endif
      normnew=compute_fnorm(unew)
      if ((normnew/norm-1.d0) .gt. fdspl) then
       do i=1,neq
        du(i)=fstp*du(i)
       enddo
       cntl=cntl+1
c       write (*,*) cntl, ' norm step reduction', normnew
       goto 5198
      else
       do i=1,neq
        u(i)=unew(i)
       enddo
       norm=normnew
      endif
c
      call compute_f(u,f)
c
c     --- check convergence ---
      iter_solver=iter_solver+1
      if (iter_solver .lt. iter_solver_max) then
      if (verbose_step) then
      write (iout,*)'-- Newton iter# ',iter_solver,
     +' norm ', norm
      endif

       if (norm .gt. eps) then
       goto 5199
      endif
       ! if norm<eps then do nothing and exit if at line 107
      else
       if (verbose_step) then
       write (iout,*) 'iter# ', iter_solver,
     * ' norm ', norm
       call face_warning('implicit precision was not reached!')
      endif

      endif

      call get_density_values(u)
      normf=norm
      end subroutine newton_solver

      end module modFACE_solver
