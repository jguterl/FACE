module modFACE_cap
use modFACE_header
implicit none
contains
    !    linear cap function to model saturation of total H
    real(DP) function cspcs(i,j,k)
        integer i, j, k
        if (active_cap) then
            if (dens(i,j,k) .lt. densm(k)) then
                cspcs=1.d0-dens(i,j,k)/densm(k)
            else
                cspcs=0.d0
            endif
        else
            cspcs=1.d0
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
      if (active_cap) then
      if (k .eq. 1) then
       csours=cspcs(i,j,1)
      endif

      do n=2,nspc-1,2
       if (k .eq. n) then
        csours=(1.d0-(dens(i,j,n)+dens(i,j,n+1))/(densm(n)+densm(n+1)))
       endif
      enddo
      endif
      return
      end

        ! linear cap function to model saturation of material with traps or H
      ! i= index time step
      ! j = index grid
      ! k,l= species
      real(DP) function csrbin(i,j,k,l)
      integer i, j, k, l, n
      csrbin=1.0
      if (active_cap) then
      do n=3,nspc,2
       if (l.eq.n) then
       if ((k.eq.1).or.(k.eq.n-1).or.(k.eq.n)) then
        csrbin=cspcs(i,j,1)
       endif
       endif
      enddo
      endif
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



      !>     linear cap for thermal detrappinng
!!limited by amount of free H in material)
      real(DP) function ctherm(i,j,k,l)
      integer i, j, k, l, n

      ctherm=1.d0
      if (active_cap) then
          do n=3,nspc,2
              if (l .eq. n) then
                  ctherm=cspcs(i,j,1)
              endif
          enddo
      endif

      return
      end


      subroutine cap_srf_flx(k,i,cabs_l, cdes_l, cb_l, cads_l,cabs_r, cdes_r, cb_r, cads_r)
            real(DP) :: cabs_l, cdes_l, cb_l, cads_l
      real(DP) :: cabs_r, cdes_r, cb_r, cads_r
      integer  :: k,i

        ! computing cap factor to mimic saturation by hydrogen
        if (dsrfl(i,k) .lt. dsrfm(k)) then
            cabs_l=1.d0-dsrfl(i,k)/dsrfm(k)
        else
            cabs_l=0.d0
        endif

        if (dsrfr(i,k) .lt. dsrfm(k)) then
            cabs_r=1.d0-dsrfr(i,k)/dsrfm(k)
        else
            cabs_r=0.d0
        endif

        if (dsrfl(i,k) .gt. 0.d0) then
            cdes_l=1.d0
        else
            cdes_l=0.d0
        endif

        if (dsrfr(i,k) .gt. 0.d0) then
            cdes_r=1.d0
        else
            cdes_r=0.d0
        endif

        if ((dsrfl(i,k) .gt. 0.d0) .and. (dens(i,0,k) .lt. densm(k))) then
            cb_l=1.d0-dens(i,0,k)/densm(k)
        else
            cb_l=0.d0
        endif

        if ((dsrfr(i,k) .gt. 0.d0) .and.(dens(i,ngrd,k) .lt. densm(k))) then
            cb_r=1.d0-dens(i,ngrd,k)/densm(k)
        else
            cb_r=0.d0
        endif

        cads_l=1.d0
        cads_r=1.d0
        end subroutine cap_srf_flx
      end module modFACE_cap
