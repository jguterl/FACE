
 module modFACE_solver
    use modFACE_precision
    use modFACE_header
    use modFACE_functions
    use modFACE_output
    use modFACE_error
    implicit none

contains

    subroutine func(u,f)

        integer i, j, k
        real(DP):: u(neq), f(neq)


        call update(u)

        i=0
        do k=1,nspc

            i=i+1
            if (steady_state .eq. "no") then
                !     --- 1st order BDF ---
                if (order_solver.eq.1) then
                    f(i)=u(i)-a11*dsrfl(ndt-1,k)
                    f(i)=f(i)-a12*Gsrf_l (ndt  ,k)*dt_face
                !     --- 2nd order BDF ---
                elseif (order_solver.eq.2) then
                    f(i)=u(i)-a21*dsrfl(ndt-1,k) &
                        -a22*dsrfl(ndt-2,k)
                    f(i)=f(i)-a23*Gsrf_l (ndt  ,k)*dt_face
                !     --- 5th order BDF ---
                elseif (order_solver.eq.5) then
                    f(i)=u(i)-a51*dsrfl(ndt-1,k)&
                        -a52*dsrfl(ndt-2,k)&
                        -a53*dsrfl(ndt-3,k)&
                        -a54*dsrfl(ndt-4,k)
                    f(i)=f(i)-a55*Gsrf_l (ndt  ,k)*dt_face
                else
                    write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
                endif
            else
                if (dsrfl(ndt,k) .gt. 1.d0) then
                    f(i)=Gsrf_l(ndt,k)
                else
                    f(i)=Gsrf_l(ndt,k)*dsrfl(ndt,k)
                endif
            endif

            do j=0,ngrd
                !     --- step ---
                i=i+1
                if (steady_state .eq. "no") then
                    !     --- 1st order BDF ---
                    if (order_solver.eq.1) then
                        f(i)=u(i)-a11*dens(ndt-1,j,k)
                        f(i)=f(i)-a12*rtd (ndt  ,j,k)*dt_face
                    !     --- 2nd order BDF ---
                    elseif (order_solver.eq.2) then
                        f(i)=u(i)-a21*dens(ndt-1,j,k)&
                            -a22*dens(ndt-2,j,k)
                        f(i)=f(i)-a23*rtd (ndt  ,j,k)*dt_face
                    !     --- 5th order BDF ---
                    elseif (order_solver.eq.5) then
                        f(i)=u(i)-a51*dens(ndt-1,j,k)&
                            -a52*dens(ndt-2,j,k)&
                            -a53*dens(ndt-3,j,k)&
                            -a54*dens(ndt-4,j,k)
                        f(i)=f(i)-a55*rtd (ndt  ,j,k)*dt_face
                    else
                        write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
                    endif
                else
                    if (dens(ndt,j,k) .gt. 1.d0) then
                        f(i)=rtd(ndt,j,k)
                    else
                        f(i)=rtd(ndt,j,k)*dens(ndt,j,k)
                    endif
                endif
            enddo

            i=i+1
            if (steady_state .eq. "no") then
            !     --- 1st order BDF ---
            if (order_solver.eq.1) then
                f(i)=u(i)-a11*dsrfr(ndt-1,k)
                f(i)=f(i)-a12*Gsrf_r (ndt  ,k)*dt_face
            !     --- 2nd order BDF ---
            elseif (order_solver.eq.2) then
                f(i)=u(i)-a21*dsrfr(ndt-1,k)&
                    -a22*dsrfr(ndt-2,k)
                f(i)=f(i)-a23*Gsrf_r (ndt  ,k)*dt_face
            !     --- 5th order BDF ---
            elseif (order_solver.eq.5) then
                f(i)=u(i)-a51*dsrfr(ndt-1,k)&
                    -a52*dsrfr(ndt-2,k)&
                    -a53*dsrfr(ndt-3,k)&
                    -a54*dsrfr(ndt-4,k)
                f(i)=f(i)-a55*Gsrf_r (ndt  ,k)*dt_face
            else
                write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
            endif
        else
        if (dsrfr(ndt,k) .gt. 1.d0) then
            f(i)=Gsrf_r(ndt,k)
        else
            f(i)=Gsrf_r(ndt,k)*dsrfr(ndt,k)
        endif
    endif

enddo

if (solve_heat_eq .eq. "yes") then
    do j=0,ngrd
        i=i+1
        if (steady_state .eq. "no") then

            !     --- 1st order BDF ---
        if (order_solver.eq.1) then
            f(i)=u(i)-a11*temp(ndt-1,j)
            f(i)=f(i)-a12*rtt (ndt  ,j)*dt_face
        !     --- 2nd order BDF ---
        elseif (order_solver.eq.2) then
            f(i)=u(i)-a21*temp(ndt-1,j)&
                -a22*temp(ndt-2,j)
            f(i)=f(i)-a23*rtt (ndt  ,j)*dt_face

        !     --- 5th order BDF ---
        elseif (order_solver.eq.5) then
            f(i)=u(i)-a51*temp(ndt-1,j) &
                -a52*temp(ndt-2,j)  &
                -a53*temp(ndt-3,j)  &
                -a54*temp(ndt-4,j)
            f(i)=f(i)-a55*rtt (ndt  ,j)*dt_face
        else
            write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
        endif
    else
        if (temp(ndt,j) .gt. 1.d0) then
            f(i)=rtt(ndt,j)
        else
            f(i)=rtt(ndt,j)*temp(ndt,j)
        endif
    endif
enddo

endif
!
if (i .ne. neq) then
    write (*,*) '***Error: wrong number of functions!', i
    stop
endif
end subroutine func
    !
    !
    !
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine jac(u,f,fdot,norm)
        implicit none


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
            call func(um,fm)
            call func(up,fp)
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
                fdot(i,i)=-1.d+99
            endif
        enddo

    end
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
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

        call DLUDCMP(a,n,indx,d,code)
        if (code .ne. 1) then
            call DLUBKSB(a,n,indx,b)
        else
            write (*,*) '***error: singular Jacobian!'
            stop
        endif

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

    end



    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    real(DP) function fnorm(u)

        integer i, idx
        real(DP):: u(neq), f(neq)
        real(DP):: norm, mxel, tmp

        call func(u,f)
        norm=0.d0
        mxel=0.d0
        do i=1,neq
            if (u(i) .gt. 1.d-10) then
                tmp=abs(f(i)/u(i))
            else
                tmp=0.d0
            endif
            norm=norm+tmp*tmp
            if (tmp .gt. mxel) then
                mxel=tmp
                idx=i
            endif
        enddo
        !      write (*,*) 'Max norm element ', mxel, ' at eq', idx
        fnorm=sqrt(norm)

        return
    end



    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine flx_update()
        implicit none

        integer k
        real(DP):: trel, tmp
        real(DP):: cenr
        parameter (cenr=1.d0)

        ! pulsed plasma flow
        if (pulsed_flux.eq."yes") then
            ! check that the pulse period is not zero
            if (tpulse.le.0) then
            write(iout,*) "ERROR: Pulsed incoming plasma flux activated but pulse_period =0"
            stop 'Exiting FACE'
            endif

            trel=time-tpulse*int(time/tpulse)
            if (trel .le. t1) then ! phase 1 of pulse
                tmp=trel/t1
                do k=1,nspc
                    inflx(k)=inflx_in(k)+(inflx_in_max(k)-inflx_in(k))*tmp
                enddo
                rad =rad_min +(rad_max -rad_min )*tmp
                cero=cero_min+(cero_max-cero_min)*tmp
            elseif (trel .le. t2) then ! phase 2 of pusle
                do k=1,nspc
                    inflx(k)=inflx_in_max(k)
                enddo
                rad =rad_max
                cero=cero_max
            elseif (trel .le. t3) then ! phase 3 of pulse
                tmp=(trel-t2)/(t3-t2)
                do k=1,nspc
                    inflx(k)=inflx_in_max(k)+(inflx_in(k)-inflx_in_max(k))*tmp
                enddo
                rad =rad_max +(rad_min -rad_max )*tmp
                cero=cero_max+(cero_min-cero_max)*tmp
            endif

        elseif (pulsed_flux.eq."no") then ! no pulse
            do k=1,nspc
                inflx(k)=inflx_in(k)
            enddo
            rad =rad_min
            cero=cero_min
        else
        call face_error("Unknown mode for pulsed_flux")
        endif
        cero=cero+gamero*inflx(1)*lambda3c ! sputtering
        ! TODO modif Q flux here
        qflx=0.d0
        do k=1,nspc
            qflx=qflx+cenr*ee*enrg(k)*inflx(k)
        enddo
        qflx=qflx+rad
        qflx=qflx-cero*ee*qform/lambda3c
        qflx=qflx-emiss*sigma_sb*(temp(ndt,0)**4.d0-temp(ndt,ngrd)**4.d0)

    end subroutine flx_update

    subroutine shift_array
    ! shifting time array down in time (current time: ndt, past time: ndt-1,ndt-2,...)
        integer i,j,k,l
        ! volume
        do i=1,ndt-1
            do j=0,ngrd
                temp(i,j)=temp(i+1,j)
                rtt (i,j)=rtt (i+1,j)
                flxt(i,j)=flxt(i+1,j)
                do k=1,nspc
                    dens(i,j,k)=dens(i+1,j,k)
                    flx (i,j,k)=flx (i+1,j,k)
                    ero_flx (i,j,k)=ero_flx (i+1,j,k)
                    src (i,j,k)=src (i+1,j,k)
                    srs (i,j,k)=srs (i+1,j,k)
                    cdif(i,j,k)=cdif(i+1,j,k)
                    rct (i,j,k)=rct (i+1,j,k)
                    rtd (i,j,k)=rtd (i+1,j,k)
                    do l=1,nspc
                        srb(i,j,k,l)=srb(i+1,j,k,l)
                    enddo
                enddo
            enddo
        enddo
        ! surface
        do i=1,ndt-1
            do k=1,nspc
                dsrfl(i,k)=dsrfl(i+1,k)
                Gsrf_l (i,k)=Gsrf_l (i+1,k)
                Gabs_l  (i,k)=Gabs_l  (i+1,k)
                Gdes_l  (i,k)=Gdes_l  (i+1,k)
                Gb_l  (i,k)=Gb_l  (i+1,k)
                Gads_l  (i,k)=Gads_l  (i+1,k)
                dsrfr(i,k)=dsrfr(i+1,k)
                Gsrf_r (i,k)=Gsrf_r (i+1,k)
                Gabs_r  (i,k)=Gabs_r  (i+1,k)
                Gdes_r  (i,k)=Gdes_r  (i+1,k)
                Gb_r  (i,k)=Gb_r  (i+1,k)
                Gads_r  (i,k)=Gads_r  (i+1,k)
                jout (i,k)=jout (i+1,k)
            enddo
        enddo

    end subroutine shift_array

    subroutine update_reaction(k)
        integer j,k,l,m
        !     --- reactions ---

        do j=0,ngrd
            rct(ndt,j,k)=0.d0
            if (k .eq. 1) then
                do l=1,nspc
                    m=1
                    if (kbin(k,l,m) .ne. 0.d0) then
                        rct(ndt,j,k)=rct(ndt,j,k) +kbin(k,l,m)*cbinar(ndt,j,k,l,m)*dens(ndt,j,l)*dens(ndt,j,m)
                    endif
                    if (nuth(k,l) .ne. 0.d0) then
                        rct (ndt,j,k)=rct (ndt,j,k)+nuth(k,l)*ctherm(ndt,j,k,l)*dens(ndt,j,l)
                    endif
                enddo
            elseif (k .eq. nspc) then
                do l=k-1,k
                    m=1
                    if (kbin(k,l,m) .ne. 0.d0) then
                        rct(ndt,j,k)=rct(ndt,j,k)+kbin(k,l,m)*cbinar(ndt,j,k,l,m)*dens(ndt,j,l)*dens(ndt,j,m)
                    endif
                    if (nuth(k,l) .ne. 0.d0) then
                        rct (ndt,j,k)=rct (ndt,j,k)+nuth(k,l)*ctherm(ndt,j,k,l)*dens(ndt,j,l)
                    endif
                enddo
            else
                do l=k-1,k+1
                    m=1
                    if (kbin(k,l,m) .ne. 0.d0) then
                        rct(ndt,j,k)=rct(ndt,j,k)+kbin(k,l,m)*cbinar(ndt,j,k,l,m)*dens(ndt,j,l)*dens(ndt,j,m)
                    endif
                    if (nuth(k,l) .ne. 0.d0) then
                        rct (ndt,j,k)=rct (ndt,j,k)+nuth(k,l)*ctherm(ndt,j,k,l)*dens(ndt,j,l)
                    endif
                enddo
            endif
        enddo
    end subroutine update_reaction

    subroutine update_surface(k)
        integer k
        real(DP) tmp
        !     --- surface ---
        ! reducing qch when close to surface saturation
        tmp=1.d2*(dsrfl(ndt,k)/dsrfm(k)-1.d0)
        qchtl(k)=qchl(k)*0.5d0*(1.d0-erf(tmp))
        tmp=1.d2*(dsrfr(ndt,k)/dsrfm(k)-1.d0)
        qchtr(k)=qchr(k)*0.5d0*(1.d0-erf(tmp))

        ! calculate rates of surface processes
        ! - left surface
        Kabs_l(k)=j0(k)*K0abs_l(k)*exp(-eekb* echl(k)/temp(ndt,   0))
        Kdes_l(k)=2.d0 *K0des_l(k)*exp(-2.d0*eekb*(echl(k)+qchtl(k))/temp(ndt,   0))
        Kb_l(k)=        K0b_l(k)  *exp(-     eekb*(ebl (k)+qchtl(k))/temp(ndt,   0))
        Kads_l(k)=      K0ads_l(k)*exp(-     eekb*(ebl (k)-esl  (k))/temp(ndt,   0))

        ! - right surface
        Kabs_r(k)=j0(k)*K0abs_r(k)*exp(-     eekb* echr(k)          /temp(ndt,ngrd))
        Kdes_r(k)=2.d0 *K0des_r(k)*exp(-2.d0*eekb*(echr(k)+qchtr(k))/temp(ndt,ngrd))
        Kb_r(k)=        K0b_r(k)  *exp(-     eekb*(ebr (k)+qchtr(k))/temp(ndt,ngrd))
        Kads_r(k)=      K0ads_r(k)*exp(-     eekb*(ebr (k)-esr  (k))/temp(ndt,ngrd))

        ! calculate surface fluxes using rates and density
        ! - left surface
        Gabs_l (ndt,k)=Kabs_l(k)
        Gdes_l (ndt,k)=Kdes_l(k)*dsrfl(ndt,k)*dsrfl(ndt,k)
        Gb_l (ndt,k)  =Kb_l(k)  *dsrfl(ndt,k)
        Gads_l (ndt,k)=Kads_l(k)*dens(ndt,0   ,k)
        ! - right surface
        Gabs_r (ndt,k)=Kabs_r(k)                           ! Gabsorp=K(gas)
        Gdes_r (ndt,k)=Kdes_r(k)*dsrfr(ndt,k)*dsrfr(ndt,k) ! Gdesorp=K*ns^2
        Gb_r (ndt,k)  =Kb_r(k)  *dsrfr(ndt,k)              ! Gbulk  =K*ns
        Gads_r (ndt,k)=Kads_r(k)*dens(ndt,ngrd,k)          ! Gadsorb=K*nb

        ! apply cap factor to mimic effects of saturation
        call set_cap_factor_surface(k,ndt)

        ! calculate effective desorptiopn and heat fluxes
        if (solve_heat_eq .eq. "yes") then
            jout(ndt,k)=jout(ndt,k)+Gdes_l(ndt,k)
            qflx=qflx+jout(ndt,k)*(ee*esl(k)-2.d0*kb*temp(ndt,0))
        endif

        ! - net flux onto surface
        Gsrf_l(ndt,k)=Gabs_l(ndt,k)-Gdes_l(ndt,k)-Gb_l(ndt,k)+Gads_l(ndt,k)
        Gsrf_r(ndt,k)=Gabs_r(ndt,k)-Gdes_r(ndt,k)-Gb_r(ndt,k)+Gads_r(ndt,k)
    end subroutine update_surface

    subroutine update_source(k)
        integer j,k,l
        real(DP) csrs, csrb,tmp
        !     --- sources ---
        if (solve_heat_eq .eq. "no") then
            do j=0,ngrd
                src(ndt,j,k)=0.d0
                if (srs(ndt,j,k) .ne. 0.d0) then
                    src (ndt,j,k)=src(ndt,j,k)+srs(ndt,j,k)*csours(ndt,j,k)
                endif
                do l=1,nspc
                    if (srb(ndt,j,k,l) .ne. 0.d0) then
                        src (ndt,j,k  )=src(ndt,j,k  )+srb(ndt,j,k,l)*dens(ndt,j,l)*csrbin(ndt,j,k,l)
                    endif
                enddo
            enddo
        elseif (solve_heat_eq .eq. "yes") then
            jout(ndt,k)=0.d0
            do j=0,ngrd
                src(ndt,j,k)=0.d0
                if (srs(ndt,j,k) .ne. 0.d0) then
                    csrs         =csours(ndt,j,k)
                    src (ndt,j,k)=src   (ndt,j,k)+srs   (ndt,j,k)*csrs
                    jout(ndt,  k)=jout  (ndt,  k)+srs(ndt,j,k)*(1.d0-csrs)*dx(j)
                endif
                do l=1,nspc
                    if (srb(ndt,j,k,l) .ne. 0.d0) then
                        csrb           =csrbin(ndt,j,k,l)
                        src (ndt,j,k  )=src   (ndt,j,k  )+srb   (ndt,j,k,l)*dens(ndt,j,l)*csrb
                        jout(ndt,  k  )=jout  (ndt,  k  )+srb   (ndt,j,k,l)*dens(ndt,j,l)*(1.d0-csrb)*dx(j)
                    endif
                enddo
            enddo
        else
            call face_error("Unknown option for solve_heat_eq:", solve_heat_eq)
        endif
    end subroutine update_source

    subroutine update_gradient_dens(k)

        integer k,j
        do j=0,ngrd-1
            flx (ndt,j,k)=(dens(ndt,j,k)-dens(ndt,j+1,k))/dx(j)
        enddo
        flx (ndt,ngrd,k)=flx (ndt,ngrd-1,k)

    end subroutine update_gradient_dens

    subroutine check_value_inputs
        integer ::j,k,l,m
        real(DP)::tmp
        character(500)::myfmt
        !     --- check for negative densities ---
        do k=1,nspc
            if (dsrfl(ndt,k) .lt. 0.d0) then
                write(myfmt,*) "(' *** error: negative left surface density of species ',i2,"//&
                    "'  the density at previous step was ', 1pe19.9e4, ' m^-2')"
                write (iout, myfmt) k,  dsrfl(ndt-1,k)
                stop
            endif


            if (dsrfr(ndt,k) .lt. 0.d0) then
                write(myfmt,*) "(' *** error: negative right surface density of species ',i2,"//&
                    "'  the density at previous step was ', 1pe19.9e4, ' m^-2')"
                write (iout, myfmt) k,  dsrfr(ndt-1,k)
                stop
            endif

            do j=0,ngrd
                if (dens(ndt,j,k) .lt. 0.d0) then
                    write(myfmt,*) "(' *** error: negative density of species ',i2, ' in cell number ', i6,",&
                        "'  the density at previous step was ', 1pe19.9e4, ' m^-3')"
                    write (iout, myfmt) k, j, dens(ndt-1,j,k)
                    stop
                endif

            enddo
        enddo

        !     --- check for density exceeding maximum ---
         do j=0,ngrd
        if (dens(ndt,j,1) .gt. densm(1)) then
            write (iout,*) '***warning: hydrogen density exceeded the limit',j, dens(ndt,j,1)
            dens(ndt,j,1)=densm(1)
        endif
        !
        do k=2,(nspc-1),2
            if ((dens(ndt,j,k)+dens(ndt,j,k+1)) .gt.(densm(k)+densm(k+1))) then
                write (iout,*) '***warning: trap density exceeded the limit',k/2, j, dens(ndt,j,k)+dens(ndt,j,k+1)
                tmp=(densm(k)+densm(k+1))/(dens(ndt,j,k)+dens(ndt,j,k+1))
                dens(ndt,j,k  )=tmp*dens(ndt,j,k  )
                dens(ndt,j,k+1)=tmp*dens(ndt,j,k+1)
            endif
        enddo
        enddo

        !    --- check for negative temperatures and update coefficients ---
        if (solve_heat_eq .eq. "yes") then
            do j=0,ngrd
                if (temp(ndt,j) .lt. 0.d0) then
                    write(myfmt,*)"' *** error: negative temperature in cell ',i4,"//&
                        "'  the temperature at previous step was ', 1pe19.9e4, ' K'"
                    write (iout, myfmt) j,  temp(ndt-1,j)
                    stop 'Exiting FACE'
                endif

            enddo
            do k=1,nspc
                do j=0,ngrd
                    cdif(ndt,j,k)=cdif0(k)*exp(-eekb*edif(k)/temp(ndt,j))
                    do l=1,nspc
                        m=1
                        kbin(k,l,m)=kbin0(k,l,m)*exp(-eekb*ebin(k,l,m)/temp(ndt,j))
                        nuth(k,l  )=nuth0(k,l  )*exp(-eekb*eth (k,l  )/temp(ndt,j))
                    enddo
                enddo
            enddo
        endif
    end subroutine check_value_inputs

    subroutine source_modification
        integer::j,k,l
        !     --- source modification ---
        if (solve_heat_eq .eq. "no") then
            do k=1,nspc
                do j=0,ngrd
                    srs(ndt,j,k)=srs(ndt-1,j,k)
                    do l=1,nspc
                        srb(ndt,j,k,l)=srb(ndt-1,j,k,l)
                    enddo
                enddo
            enddo
        elseif (solve_heat_eq .eq. "yes") then
            call flx_update()
            do k=1,nspc
                do j=0,ngrd
                    srs(ndt,j,k)=source(j,k)
                    do l=1,nspc
                        srb(ndt,j,k,l)=srcbin(j,k,l)
                    enddo
                enddo
            enddo
        else
        write(iout,*) "ERROR: Unknown option for solve_heat_eq:", solve_heat_eq
        stop
        endif

    end subroutine source_modification


    subroutine  coeff_modification
        integer j,k,l,m,n
        !     --- temperature and coefficients ---
        if (solve_heat_eq .eq. "no") then

            if (framp .eq. "none") then

                if (time .lt. tramp0) then
                    do j=0,ngrd
                        temp(ndt,j)=temp0
                    enddo
                elseif (time .le. tramp1) then
                    do j=0,ngrd
                        temp(ndt,j)=temp0+dtemp*(time-tramp0)
                    enddo
                else
                    do j=0,ngrd
                        temp(ndt,j)=temp1
                    enddo
                endif

            else ! use temperature ramp file

                if (rtime(1) .gt. time) then
                    do j=0,ngrd
                        temp(ndt,j)=rtemp(1)
                    enddo
                elseif (rtime(nramp) .gt. time) then
                    do j=0,ngrd
                        temp(ndt,j)=rtemp(nramp)
                    enddo
                else
                    do n=2,nramp-1
                        if (rtime(n) .gt. time) then
                            do j=0,ngrd
                                temp(ndt,j)=rtemp(n-1)+(rtemp(n)-rtemp(n-1))*(time-rtime(n-1))/(rtime(n)-rtime(n-1))
                            enddo
                            exit
                        endif
                    enddo
                endif
            endif ! framp

            ! WARNING: check consistency here wqith orioginal FACE
            do k=1,nspc
                do j=0,ngrd
                    cdif(ndt,j,k)=cdif0(k)*exp(-eekb*edif(k)/temp(ndt,j))
                    do l=1,nspc
                        m=1
                        kbin(k,l,m)=kbin0(k,l,m)*exp(-eekb*ebin(k,l,m)/temp(ndt,j))
                        nuth(k,l  )=nuth0(k,l  )*exp(-eekb*eth (k,l  )/temp(ndt,j))
                    enddo
                enddo
            enddo
        endif
    end subroutine coeff_modification

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
            write(iout,*)" ERROR: mismatch in vector size: neq=",neq, 'size(vec)=',i
            stop
        endif
    end subroutine build_vector

    subroutine update(u)
   integer i, j, k
      real(DP) u(neq)

!     ------------------------------------------------------------------
!      update equation terms
!     ------------------------------------------------------------------
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

      call check_value_inputs

      do k=1,nspc
      call update_source(k)    !
      call update_gradient_dens(k)      ! set flx (ndt,j,k)=(n(ndt,j+1,k)-n(ndt,j,k)/dx(j))
      call update_reaction(k)  !
!     write(iout,*) "784: Gb_l",Gb_l(ndt,k)
      call update_surface(k)


!     --- low-pass filter ---
       Gsrf_l(ndt,k)=delta*Gsrf_l(ndt-1,k)+(1.d0-delta)*Gsrf_l(ndt,k)
       Gsrf_r(ndt,k)=delta*Gsrf_r(ndt-1,k)+(1.d0-delta)*Gsrf_r(ndt,k)

!     --- erosion ---
       ero_flx(ndt,0,k)=-cero*flx(ndt,0,k) ! until now, flx is the gradient (see line 764)
       do j=1,ngrd-1
!        ero(ndt,j,k)=-cero*(dx(j)*flx(ndt,j-1,k)+dx(j-1)*flx(ndt,j,k))/(dx(j-1)+dx(j))
        ero_flx(ndt,j,k)=-cero*flx(ndt,j,k)
       enddo
!       ero(ndt,ngrd,k)=-cero*flx(ndt,ngrd,k)
       ero_flx(ndt,ngrd,k)=0.d0 ! no erosion on the right side of the material (only left side is facing plasma)

       do j=0,ngrd
        rtd(ndt,j,k)=ero_flx(ndt,j,k)
       enddo

!     --- diffusion ---
       flx(ndt,0,k)=cdif(ndt,0,k)*flx(ndt,0,k)
       rtd(ndt,0,k)=rtd(ndt,0,k)+(Gb_l(ndt,k)-Gads_l(ndt,k)-flx(ndt,0,k))*2.d0/dx(0)

       do j=1,ngrd-1
        flx(ndt,j,k)=cdif(ndt,j,k)*flx(ndt,j,k)
        rtd(ndt,j,k)=rtd(ndt,j,k)+(flx(ndt,j-1,k)-flx(ndt,j,k))/(0.5d0*(dx(j-1)+dx(j)))
       enddo
       flx(ndt,ngrd,k)=cdif(ndt,ngrd,k)*flx(ndt,ngrd,k)
       rtd(ndt,ngrd,k)=rtd(ndt,ngrd,k)+(Gb_r(ndt,k)-Gads_r(ndt,k)+flx(ndt,ngrd,k))*2.d0/dx(ngrd-1)
       do j=0,ngrd
!     --- sources ---
        rtd(ndt,j,k)=rtd(ndt,j,k)+src(ndt,j,k)
!     --- reactions ---
        rtd(ndt,j,k)=rtd(ndt,j,k)+rct(ndt,j,k)
!     --- low-pass filter ---
        rtd(ndt,j,k)=delta*rtd(ndt-1,j,k)+(1.d0-delta)*rtd(ndt,j,k)
       enddo
      enddo
!
!     --- heat conduction ---
      if (solve_heat_eq .eq. "yes") then
       do j=0,ngrd-1
        flxt(ndt,j)=(temp(ndt,j)-temp(ndt,j+1))/dx(j)
       enddo
       flxt(ndt,ngrd)=flxt(ndt,ngrd-1)

       erot(ndt,0)=-cero*flxt(ndt,0)
       do j=1,ngrd-1
!        erot(ndt,j)=-cero*(dx(j)*flxt(ndt,j-1)+dx(j-1)*flxt(ndt,j))
!      +                  /(dx(j-1)+dx(j))
        erot(ndt,j)=-cero*flxt(ndt,j)
       enddo
!       erot(ndt,ngrd)=-cero*flxt(ndt,ngrd)
       erot(ndt,ngrd)=0.d0
!
       do j=0,ngrd
        flxt(ndt,j)=thcond*flxt(ndt,j)
       enddo
!
       rtt(ndt,0)=(qflx-flxt(ndt,0))*2.d0/dx(0)/rhocp+erot(ndt,0)
       do j=1,ngrd-1
        rtt(ndt,j)=(flxt(ndt,j-1)-flxt(ndt,j))/(0.5d0*(dx(j-1)+dx(j))*rhocp)+erot(ndt,j)
       enddo
       rtt(ndt,ngrd)=0.d0
!
       do j=0,ngrd
!     --- low-pass filter ---
        rtt(ndt,j)=delta*rtt(ndt-1,j)+(1.d0-delta)*rtt(ndt,j)
       enddo
      endif
!
      end subroutine update

      subroutine set_cap_factor_surface(k,i)
      ! setting cap factor to mimic saturation by hydrogen
      real(DP) :: c1l, c2l, c3l, c4l
      real(DP) :: c1r, c2r, c3r, c4r
      integer  :: k,i

        ! left surface
        if (dsrfl(i,k) .lt. dsrfm(k)) then
            c1l=1.d0-dsrfl(i,k)/dsrfm(k)
        else
            c1l=0.d0
        endif

        if (dsrfr(i,k) .lt. dsrfm(k)) then
            c1r=1.d0-dsrfr(i,k)/dsrfm(k)
        else
            c1r=0.d0
        endif

        if (dsrfl(i,k) .gt. 0.d0) then
            c2l=1.d0
        else
            c2l=0.d0
        endif

        if (dsrfr(i,k) .gt. 0.d0) then
            c2r=1.d0
        else
            c2r=0.d0
        endif

        if ((dsrfl(i,k) .gt. 0.d0) .and. (dens(i,0,k) .lt. densm(k))) then
            c3l=1.d0-dens(i,0,k)/densm(k)
        else
            c3l=0.d0
        endif

        if ((dsrfr(i,k) .gt. 0.d0) .and.(dens(i,ngrd,k) .lt. densm(k))) then
            c3r=1.d0-dens(i,ngrd,k)/densm(k)
        else
            c3r=0.d0
        endif

        c4l=1.d0
        c4r=1.d0
        Gabs_l (i,k)=Gabs_l (i,k) *c1l
        Gdes_l (i,k)=Gdes_l (i,k) *c2l
        Gb_l (i,k)=Gb_l (i,k)     *c3l
        Gads_l (i,k)=Gads_l (i,k) *c4l
        ! right surface
        Gabs_r (i,k)=Gabs_r (i,k) *c1r
        Gdes_r (i,k)=Gdes_r (i,k) *c2r
        Gb_r (i,k)=Gb_r (i,k)     *c3r
        Gads_r (i,k)=Gads_r (i,k) *c4r

      end subroutine set_cap_factor_surface

end module modFACE_solver
