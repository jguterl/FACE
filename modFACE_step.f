      module modFACE_step
      use modFACE_solver
      use modFACE_header
      use modFACE_functions
      use modFACE_output

      implicit none
      contains


      subroutine step()
c     ******************************************************************
c     * advance problem in time                                        *
c     *                                                                *
c     * Author: Roman D. Smirnov                                       *
c     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
c     *                                                                *
c     ******************************************************************

c
      integer i, j, k,ii
      integer cntmax, cntl, idx
      real(DP) told

      real(DP) eps, udspl, fdspl, gdspl, fstp
      real(DP):: u(neq), du(neq), f(neq), fdot(neq,neq)
      real(DP)  unew(neq), unorm(neq)
      real(DP) norm, normnew

      parameter (cntmax=100)
      parameter (eps=3.d-3, fstp=1.d-1)
      parameter (udspl=9.d-1, fdspl=9.d0, gdspl=1.d-3)

      call shift_array()
c     ------------------------------------------------------------------
c      implicit step
c     ------------------------------------------------------------------
      told=time
      time=told+dt_face
c
      delta=1.d0/(1.d0+2.d0*pi*nucut*dt_face)

      call source_modification
      call coeff_modification


c
      cnt=0
c
      call build_vector(u,du)

      norm=fnorm(u)
      call func(u,f)
      if (verbose_debug) then
      do ii=1,neq
      write(iout,*) 'i=',ii,'u=',u(ii),'f(i)=',f(ii)
      enddo
      endif

5199  call jac(u,f,fdot,norm)
      call dsolve(du,u,f,fdot)
c
c     --- newton step reduction ---
      norm=fnorm(u)
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
      normnew=fnorm(unew)
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
      call func(u,f)
c
c     --- check convergence ---
      cnt=cnt+1
      if (cnt .lt. cntmax) then
       if (verbose_step) write (iout,*) 'newton iter# ', cnt,
     +' norm ', norm
       if (norm .gt. eps) goto 5199
      else
       if (verbose_step) write (iout,*) 'iter# ', cnt, ' norm ', norm
       if (verbose_step) write (iout,*)
     +'Warning: implicit precision was not reached!'
      endif
c
      i=0
      do k=1,nspc
       i=i+1
       dsrfl(ndt,k)=u(i)
       do j=0,ngrd
        i=i+1
        dens(ndt,j,k)=u(i)
       enddo
       i=i+1
       dsrfr(ndt,k)=u(i)
      enddo
      if (solve_heat_eq .eq. "yes") then
       do j=0,ngrd
        i=i+1
        temp(ndt,j)=u(i)
       enddo
      endif

      call final_check()
      normf=norm
      end subroutine step

      subroutine final_check
      integer i,j,k
      ! 1) check that density are below their max values and above zero
      !  2) check that density are not nan
      ! 3) check that temperature >0 and not nan

      !1) to be completed
!      2)
       if (finalcheck) then
       do i=1,ndt
         do j=0,ngrd
          do k=1,nspc
         if  (dens(i,j,k).lt.0) then
         write(iout,*) 'dens<0',i,j,k
         stop
         elseif (isnan(dens(i,j,k))) then
         write(iout,*) 'dens is nan',i,j,k
         stop
         endif
         enddo
         enddo
         enddo

         do i=1,ndt
          do k=1,nspc
         if  (dsrfl(i,k).lt.0) then
         write(iout,*) 'dsrfl<0',i,k
         stop
         elseif (isnan(dsrfl(i,k))) then
         write(iout,*) 'dsrfl is nan',i,k
         stop
         endif
         enddo
         enddo

          do i=1,ndt
          do k=1,nspc
         if  (dsrfr(i,k).lt.0) then
         write(iout,*) 'dsrfr<0',i,k
         stop
         elseif (isnan(dsrfr(i,k))) then
         write(iout,*) 'dsrfr is nan',i,k
         stop
         endif
         enddo
         enddo
         endif

      !2)to be completed
      end subroutine

      end module
