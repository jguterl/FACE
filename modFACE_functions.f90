!     ******************************************************************
!     * process terms                                                  *
!     *                                                                *
!     * Author: Roman D. Smirnov                                       *
!     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
!     *                                                                *
!     ******************************************************************
!
      module modFACE_functions
      use modFACE_header
      implicit none

      contains

      real(DP) function source(j,k)
      integer j, k, n

      source=0.d0

      if (k .eq. 1) then
       source=srate(1)*exp(-0.5d0*abs((x(j)-rdpth(1))/sigma(1))**2.d0)
       return
      endif
      ! create empty traps
      do n=2,nspc-1,2
       if (k .eq. n) then
        source=srate(n)*(1.d0-erf((x(j)-rdpth(n))/(sqrt2*sigma(n))))
       return
       endif
      enddo
      return
      end

!     Rate of induced-detrapping: detrappping rate from trap n+1 into 1 + n
     !due to collisions of filled traps n with incoming ion
      real(DP) function srcbin(j,k,l)
!
      integer j, k, l, n
      real(DP):: s, rs, crsc
      ! why?
!      data crsc /3.d-20/
      data crsc /0.d0/


      srcbin=0
      do n=3,nspc,2
       if (l .eq. n) then
        s=crsc*inflx(1)
        rs=s*(1.d0-erf((x(j)-rdpth(n))/(sqrt2*sigma(n))))
       if ((k .eq. 1) .or. (k .eq. n-1)) then
        srcbin=+rs ! detrap from n -> create 1 free H
        srcbin=+rs ! detrap from n -> create n-1 empty trap
        return
        elseif (k .eq. n) then
        srcbin=-rs ! detrap from trap n -> filled trap n disapears
        return
        endif
       endif
      enddo

      return
      end function srcbin
!
!
!    linear cap function to model saturation of total H
      real(DP) function cspcs(i,j,k)
      integer i, j, k
      if (dens(i,j,k) .lt. densm(k)) then
       cspcs=1.d0-dens(i,j,k)/densm(k)
      else
       cspcs=0.d0
      endif
      return
      end

      ! linear cap function to model saturation of material with traps or H
      ! i= index time step
      ! j = index grid
      ! k= species
      real(DP) function csours(i,j,k)
      integer i, j, k, n
      csours=1.0

      if (k .eq. 1) then
       csours=cspcs(i,j,1)
      endif

      do n=2,nspc-1,2
       if (k .eq. n) then
        csours=(1.d0-(dens(i,j,n)+dens(i,j,n+1))/(densm(n)+densm(n+1)))
       endif
      enddo

      return

      end

        ! linear cap function to model saturation of material with traps or H
      ! i= index time step
      ! j = index grid
      ! k,l= species
      real(DP) function csrbin(i,j,k,l)
      integer i, j, k, l, n
      csrbin=1.0
      do n=3,nspc,2
       if (l.eq.n) then
       if ((k.eq.1).or.(k.eq.n-1).or.(k.eq.n)) then
        csrbin=cspcs(i,j,1)
       endif
       endif
      enddo
      return
      end


      ! implantation rate
      real(DP) function srate(k)
      integer k, n
      real(DP) r, sig, crsc
!      data crsc /3.d-20/
       crsc=0.d0

       srate=0.0

        if (k .eq. 1) then
         r  =rdpth(1)
         sig=sqrt2*sigma(1)
         srate=inflx(1)/(0.5d0*sqrt(pi)*(1.d0+erf(r/sig))*sig)
         return
        endif
        ! create empty traps
        do n=2,nspc-1,2
         if (k .eq. n) then
          srate=crsc*inflx(1)/lambda3c
          return
          endif
        enddo

      return
      end


     ! implantation depth
     ! TODO use gsigma and gaussian def
      real(DP) function rdpth(k)
      integer k
      rdpth=1.d-8
      return
      end


     ! sigma of implantation profile
      real(DP) function sigma(k)
      integer k
      sigma=3.d-9
      return
      end

      ! pre-expoential factor of trapping process
      ! H + n->n+1
      real(DP) function kbinar(k,l,m)
      integer k, l, m, n
       kbinar=0.0
       do n=2,nspc-1,2
        if(m.eq.1) then
       if (l.eq.n) then
       if ((k.eq.1).or.(k.eq.n)) then
        kbinar=-nu(n)*lambda3c
        kbinar=-nu(n)*lambda3c
        return
        elseif (k.eq.n+1) then
        kbinar=+nu(n)*lambda3c
        return
        endif
! this if for multiple H in traps
!  tmp(1  ,n+1,1)=+rk(n,2)
!        tmp(n  ,n+1,1)=+rk(n,2)
!       tmp(n+1,n+1,1)=-rk(n,2)
            endif
              endif
       enddo

      return
      end


     !activation energy of trapping process
 !    (1) + (n) -> (n+1)
      real(DP) function ebinar(k,l,m)
      integer k, l, m, n

      ebinar=0.0
       do n=2,nspc-1,2
       if (m.eq.1) then
       if (l.eq.n) then
       if ((k.eq.1).or.(k.eq.n).or.(k.eq.n+1)) then
        ebinar=etr(n)
        return
        endif
        endif
        endif
!        tmp(1  ,n+1,1)=edtr(n+1)
!        tmp(n  ,n+1,1)=edtr(n+1)
!        tmp(n+1,n+1,1)=edtr(n+1)
       enddo
      return
      end


!     linear cap for trapping
      real(DP) function cbinar(i,j,k,l,m)

      integer i, j, k, l, m
      cbinar=1.d0
!      do n=2,nspc-1,2
!
!       if ((l .eq. n+1) .and. (m .eq. 1)) then
!        if ((k .eq. 1) .or. (k .eq. n) .or. (k .eq. n+1)) then
!         c0=cspcs(i,j,1)
!         tmp(1  ,n+1,1)=c0
!         tmp(n  ,n+1,1)=c0
!         tmp(n+1,n+1,1)=c0
!        endif
!       endif
!
!      enddo
!
      return
      end

      !prex-exponential factor for thermal (activated) detrapping
      !(n+1)-> (1) + (n)
      real(DP) function ktherm(k,l)
      implicit none
      integer k, l, n
       ktherm=0.d0
       do n=3,nspc,2
       if (l .eq. n) then
       if ((k .eq. 1) .or. (k .eq. n-1)) then
        ktherm=+nu(n)
        ktherm=+nu(n)
        return
        elseif (k .eq. n) then
        ktherm=-nu(n)
        return
        endif
        endif
       enddo

      return

      end function ktherm


      !> Activation energy for thermal (activated) detrapping
  !! @param aggr information about the aggregates
  !! @todo Handle special case
  !! (n) -> (1) +(n-1)
      real(DP) function etherm(k,l)
      implicit none

      integer ::k,l !< index of species
      integer::n
      etherm=0
       do n=3,nspc,2
           if (l .eq. n) then
               if ((k .eq. 1) .or. (k .eq. n-1) .or. (k .eq. n)) then
        etherm=edtr(n)
        return
        endif
        endif
       enddo
      return
      end


!>     linear cap for thermal detrappinng
!!limited by amount of free H in material)
      real(DP) function ctherm(i,j,k,l)
      integer i, j, k, l, n

      ctherm=1.d0

      do n=3,nspc,2
       if (l .eq. n) then
        ctherm=cspcs(i,j,1)
       endif
      enddo
      return
      end

       subroutine print_source(ifile)
      integer j,k,ifile
      write(ifile,*) "#source rate (function source(j,k))"
      write(ifile,*)"xgrid ",(namespc(k)//" ",k=1,nspc)
      do j=1,ngrd
      write(ifile,*) x(j),(source(j,k),k=1,nspc)
      enddo
      end subroutine print_source

      end module modFACE_functions
