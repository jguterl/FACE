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
      use modFACE_error
      implicit none

      contains

          real(DP) function source(j,k)
              ! Source rate of incoming species
              integer j, k
              real(DP) :: source_rate
              source=0.d0
              source_rate=srate(k)
              if (source_rate.gt.0d0) then

                  if (implantation_model(k).eq.'G'.and.implantation_width(k).ne.0d0) then

                      source=source_rate*exp(-0.5d0*abs((x(j)-implantation_depth(k))/implantation_width(k))**2.d0)
                      return

                  elseif  (implantation_model(k).eq.'S') then

                      if ((j_implantation_depth(k).gt.0).and.j.le.j_implantation_depth(k)-1) then
                          source=source_rate
                      else
                          source=0d0
                      endif
                      return

                  elseif  (implantation_model(k).eq.'E'.and.implantation_width(k).ne.0d0) then

                      source=source_rate*(1.d0-erf((x(j)-implantation_depth(k))/(sqrt2*implantation_width(k))))
                      return

                  else

                      if (implantation_width(k).eq.0d0) then
                          call face_error("implantation_width(k)=0 with gaussian-like implantation model: k=",k)
                      else
                          call face_error("Unknown implantation model : ",implantation_model(k),"; k=",k)
                      endif

                  endif

              elseif (source_rate.lt.0d0) then

                  call face_error("source rate < 0 for k=",k," ; srate(k)=",srate(k))
              endif

              return

          end function source

               ! implantation rate
          real(DP) function srate(k)
              integer k
              real(DP) r, sig
              !      data crsc /3.d-20/

              srate=0.0
              if (inflx(k).gt.0d0) then

                  if (implantation_model(k).eq.'G'.and.implantation_width(k).ne.0d0) then

                      r  =implantation_depth(k)
                      sig=sqrt2*implantation_width(k)
                      srate=inflx(k)/(0.5d0*sqrt(pi)*(1.d0+erf(r/sig))*sig)
                      return

                  elseif (implantation_model(k).eq.'E'.and.implantation_width(k).ne.0d0) then

                      call face_error("implentation model E required proper normalization in FACE source code")
                      srate=inflx(k)
                      return

                  elseif  (implantation_model(k).eq.'S') then

                      if (j_implantation_depth(k).gt.0) then
                          srate=inflx(k)/(x(j_implantation_depth(k))-x(0))
                      else
                      srate=0d0
                      endif
                      return

                  else
                      call face_error("Unknown implantation model : ",implantation_model(k),"; k=",k)
                  endif

              elseif (inflx(k).lt.0d0) then

                  call face_error("inflx(k) < 0 : k=",k," ;inflx(k)=",inflx(k))

              endif


              return
          end function srate



      real(DP) function srcbin(j,k,l)
      ! Rate of induced-detrapping: detrappping rate from trap n+1 into 1 + n
      ! due to collisions of filled traps n with incoming ion
      integer j, k, l, n
      real(DP):: s, rs, crsc
      ! why?
!      data crsc /3.d-20/
      data crsc /0.d0/


      srcbin=0.d0
      do n=3,nspc,2
       if (l .eq. n) then
        s=crsc*inflx(n)
        rs=s*(1.d0-erf((x(j)-implantation_depth(n))/(sqrt2*implantation_width(n))))
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




      ! pre-expoential factor of trapping process
      ! H + n->n+1
      real(DP) function kbinar(k,l,m)
      integer k, l, m, n
       kbinar=0.0
       do n=2,nspc-1,2
        if(m.eq.1) then
       if (l.eq.n) then
       if ((k.eq.1).or.(k.eq.n)) then
        kbinar=-nu(n)*lambda**3*cvlm
        kbinar=-nu(n)*lambda**3*cvlm
        return
        elseif (k.eq.n+1) then
        kbinar=+nu(n)*lambda**3*cvlm
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




       subroutine print_source(ifile)
      integer j,k,ifile
      write(ifile,*) "#source rate (function source(j,k))"
      write(ifile,*)"xgrid ",(namespc(k)//" ",k=1,nspc)
      do j=1,ngrd
      write(ifile,*) x(j),(source(j,k),k=1,nspc)
      enddo
      end subroutine print_source

      end module modFACE_functions
