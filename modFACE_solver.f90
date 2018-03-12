
 module modFACE_solver
    use modFACE_precision
    use modFACE_header
    use modFACE_functions
    use modFACE_output
    implicit none

contains

    subroutine func(u,f)

        integer i, j, k
        real(DP):: u(neq), f(neq)


        call update(u)

        i=0
        do k=1,nspc

            i=i+1
            if (stdst .eq. "no") then
                !     --- 1st order BDF ---
                if (order_solver.eq.1) then
                    f(i)=u(i)-a11*dsrfl(ndt-1,k)
                    f(i)=f(i)-a12*rtsl (ndt  ,k)*dt_face
                !     --- 2nd order BDF ---
                elseif (order_solver.eq.2) then
                    f(i)=u(i)-a21*dsrfl(ndt-1,k) &
                        -a22*dsrfl(ndt-2,k)
                    f(i)=f(i)-a23*rtsl (ndt  ,k)*dt_face
                !     --- 5th order BDF ---
                elseif (order_solver.eq.5) then
                    f(i)=u(i)-a51*dsrfl(ndt-1,k)&
                        -a52*dsrfl(ndt-2,k)&
                        -a53*dsrfl(ndt-3,k)&
                        -a54*dsrfl(ndt-4,k)
                    f(i)=f(i)-a55*rtsl (ndt  ,k)*dt_face
                else
                    write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
                endif
            else
                if (dsrfl(ndt,k) .gt. 1.d0) then
                    f(i)=rtsl(ndt,k)
                else
                    f(i)=rtsl(ndt,k)*dsrfl(ndt,k)
                endif
            endif

            do j=0,ngrd
                !     --- step ---
                i=i+1
                if (stdst .eq. "no") then
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
            if (stdst .eq. "no") then
            !     --- 1st order BDF ---
            if (order_solver.eq.1) then
                f(i)=u(i)-a11*dsrfr(ndt-1,k)
                f(i)=f(i)-a12*rtsr (ndt  ,k)*dt_face
            !     --- 2nd order BDF ---
            elseif (order_solver.eq.2) then
                f(i)=u(i)-a21*dsrfr(ndt-1,k)&
                    -a22*dsrfr(ndt-2,k)
                f(i)=f(i)-a23*rtsr (ndt  ,k)*dt_face
            !     --- 5th order BDF ---
            elseif (order_solver.eq.5) then
                f(i)=u(i)-a51*dsrfr(ndt-1,k)&
                    -a52*dsrfr(ndt-2,k)&
                    -a53*dsrfr(ndt-3,k)&
                    -a54*dsrfr(ndt-4,k)
                f(i)=f(i)-a55*rtsr (ndt  ,k)*dt_face
            else
                write(iout,*) "ERROR: order of solver not implanted. order=",order_solver
            endif
        else
        if (dsrfr(ndt,k) .gt. 1.d0) then
            f(i)=rtsr(ndt,k)
        else
            f(i)=rtsr(ndt,k)*dsrfr(ndt,k)
        endif
    endif

enddo

if (solve_heat_eq .eq. "yes") then
    do j=0,ngrd
        i=i+1
        if (stdst .eq. "no") then

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
        if (stdst .eq. "yes") then
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
        if (tp.gt.0) then
        trel=time-tp*int(time/tp)
        if (trel .le. t1) then ! phase 1 of pulse
            tmp=trel/t1
            do k=1,nspc
                inflx(k)=inflx_min(k)+(inflx_max(k)-inflx_min(k))*tmp
            enddo
            rad =rad_min +(rad_max -rad_min )*tmp
            cero=cero_min+(cero_max-cero_min)*tmp
        elseif (trel .le. t2) then ! phase 2 of pusle
            do k=1,nspc
                inflx(k)=inflx_max(k)
            enddo
            rad =rad_max
            cero=cero_max
        elseif (trel .le. t3) then ! phase 3 of pulse
            tmp=(trel-t2)/(t3-t2)
            do k=1,nspc
                inflx(k)=inflx_max(k)+(inflx_min(k)-inflx_max(k))*tmp
            enddo
            rad =rad_max +(rad_min -rad_max )*tmp
            cero=cero_max+(cero_min-cero_max)*tmp
        endif
        else ! no pulse
            do k=1,nspc
                inflx(k)=inflx_min(k)
            enddo
            rad =rad_min
            cero=cero_min
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

    subroutine shift_array()
        integer i,j,k,l
        !    ------------------------------------------------------------------
        !      shifting time array down
        !     ------------------------------------------------------------------
        do i=1,ndt-1
            do j=0,ngrd
                temp(i,j)=temp(i+1,j)
                rtt (i,j)=rtt (i+1,j)
                flxt(i,j)=flxt(i+1,j)
                do k=1,nspc
                    dens(i,j,k)=dens(i+1,j,k)
                    flx (i,j,k)=flx (i+1,j,k)
                    ero (i,j,k)=ero (i+1,j,k)
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
        do i=1,ndt-1
            do k=1,nspc
                dsrfl(i,k)=dsrfl(i+1,k)
                rtsl (i,k)=rtsl (i+1,k)
                j1l  (i,k)=j1l  (i+1,k)
                j2l  (i,k)=j2l  (i+1,k)
                j3l  (i,k)=j3l  (i+1,k)
                j4l  (i,k)=j4l  (i+1,k)
                dsrfr(i,k)=dsrfr(i+1,k)
                rtsr (i,k)=rtsr (i+1,k)
                j1r  (i,k)=j1r  (i+1,k)
                j2r  (i,k)=j2r  (i+1,k)
                j3r  (i,k)=j3r  (i+1,k)
                j4r  (i,k)=j4r  (i+1,k)
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
        real(DP) c1l, c2l, c3l, c4l,tmp
        real(DP) c1r, c2r, c3r, c4r
        !     --- surface ---
        tmp=1.d2*(dsrfl(ndt,k)/dsrfm(k)-1.d0)
        qchtl(k)=qchl(k)*0.5d0*(1.d0-erf(tmp))
        !write(iout,*) "443: tmp",tmp
        k1l(k)=j0(k)*r1l(k) *exp(-eekb* echl(k)/temp(ndt,   0))
        k2l(k)=2.d0 *r2l(k)*exp(-2.d0*eekb*(echl(k)+qchtl(k))/temp(ndt,   0))
        !write(iout,*) "445: r3l",r3l
        !write(iout,*) "446: temp(ndt,   0)",temp(ndt,   0)
        !write(iout,*) "447: ebl",ebl
        !write(iout,*) "447: qchtl",qchtl
        k3l(k)=      r3l(k)*exp(-     eekb*(ebl (k)+qchtl(k))/temp(ndt,   0))
        k4l(k)=      r4l(k)*exp(-     eekb*(ebl (k)-esl  (k))/temp(ndt,   0))

        tmp=1.d2*(dsrfr(ndt,k)/dsrfm(k)-1.d0)
        qchtr(k)=qchr(k)*0.5d0*(1.d0-erf(tmp))
        k1r(k)=j0(k)*r1r(k)*exp(-     eekb* echr(k)          /temp(ndt,ngrd))
        k2r(k)=2.d0 *r2r(k)*exp(-2.d0*eekb*(echr(k)+qchtr(k))/temp(ndt,ngrd))
        k3r(k)=      r3r(k)*exp(-     eekb*(ebr (k)+qchtr(k))/temp(ndt,ngrd))
        k4r(k)=      r4r(k)*exp(-     eekb*(ebr (k)-esr  (k))/temp(ndt,ngrd))

        if (dsrfl(ndt,k) .lt. dsrfm(k)) then
            c1l=1.d0-dsrfl(ndt,k)/dsrfm(k)
        else
            c1l=0.d0
        endif
        if (dsrfr(ndt,k) .lt. dsrfm(k)) then
            c1r=1.d0-dsrfr(ndt,k)/dsrfm(k)
        else
            c1r=0.d0
        endif

        if (dsrfl(ndt,k) .gt. 0.d0) then
            c2l=dsrfl(ndt,k)*dsrfl(ndt,k)
        else
            c2l=0.d0
        endif
        if (dsrfr(ndt,k) .gt. 0.d0) then
            c2r=dsrfr(ndt,k)*dsrfr(ndt,k)
        else
            c2r=0.d0
        endif

        if ((dsrfl(ndt,k) .gt. 0.d0) .and. (dens(ndt,0,k) .lt. densm(k))) then
            c3l=dsrfl(ndt,k)*(1.d0-dens(ndt,0,k)/densm(k))
        else
            c3l=0.d0
        endif
        if ((dsrfr(ndt,k) .gt. 0.d0) .and.(dens(ndt,ngrd,k) .lt. densm(k))) then
            c3r=dsrfr(ndt,k)*(1.d0-dens(ndt,ngrd,k)/densm(k))
        else
            c3r=0.d0
        endif

        c4l=dens(ndt,0   ,k)
        c4r=dens(ndt,ngrd,k)
         !write(iout,*) "490: c3l",c3l
         !write(iout,*) "490: k3l",k3l
        j1l (ndt,k)=k1l(k)*c1l
        j2l (ndt,k)=k2l(k)*c2l
        j3l (ndt,k)=k3l(k)*c3l
        j4l (ndt,k)=k4l(k)*c4l
        j1r (ndt,k)=k1r(k)*c1r
        j2r (ndt,k)=k2r(k)*c2r
        j3r (ndt,k)=k3r(k)*c3r
        j4r (ndt,k)=k4r(k)*c4r

        if (solve_heat_eq .eq. "yes") then
            jout(ndt,k)=jout(ndt,k)+j2l(ndt,k)
            qflx=qflx+jout(ndt,k)*(ee*esl(k)-2.d0*kb*temp(ndt,0))
        endif

        !    --- surface ---
        rtsl(ndt,k)=j1l(ndt,k)-j2l(ndt,k)-j3l(ndt,k)+j4l(ndt,k)
        rtsr(ndt,k)=j1r(ndt,k)-j2r(ndt,k)-j3r(ndt,k)+j4r(ndt,k)
    end subroutine update_surface

    subroutine update_source(k)
        integer j,k,l
        real(DP) csrs, csrb,tmp
        !     --- sources ---
        if (solve_heat_eq .eq. "no") then
            do j=0,ngrd
                src(ndt,j,k)=0.d0
                if (srs(ndt,j,k) .ne. 0.d0) then
                    csrs         =csours(ndt,j,k)
                    tmp          =srs   (ndt,j,k)*csrs
                    src (ndt,j,k)=src   (ndt,j,k)+tmp
                endif
                do l=1,nspc
                    if (srb(ndt,j,k,l) .ne. 0.d0) then
                        csrb           =csrbin(ndt,j,k,l)
                        tmp            =srb   (ndt,j,k,l)*dens(ndt,j,l)*csrb
                        src (ndt,j,k  )=src   (ndt,j,k  )+tmp
                    endif
                enddo
            enddo
        elseif (solve_heat_eq .eq. "yes") then
            jout(ndt,k)=0.d0
            do j=0,ngrd
                src(ndt,j,k)=0.d0
                if (srs(ndt,j,k) .ne. 0.d0) then
                    csrs         =csours(ndt,j,k)
                    tmp          =srs   (ndt,j,k)*csrs
                    src (ndt,j,k)=src   (ndt,j,k)+tmp
                    jout(ndt,  k)=jout  (ndt,  k)+srs(ndt,j,k)*(1.d0-csrs)*dx(j)
                endif
                do l=1,nspc
                    if (srb(ndt,j,k,l) .ne. 0.d0) then
                        csrb           =csrbin(ndt,j,k,l)
                        tmp            =srb   (ndt,j,k,l)*dens(ndt,j,l)*csrb
                        src (ndt,j,k  )=src   (ndt,j,k  )+tmp
                        jout(ndt,  k  )=jout  (ndt,  k  )+srb   (ndt,j,k,l)*dens(ndt,j,l)*(1.d0-csrb)*dx(j)
                    endif
                enddo
            enddo
        else
            write(iout,*) "ERROR: Unknown option for solve_heat_eq:", solve_heat_eq
            stop
        endif
    end subroutine update_source

    subroutine update_flux(k)
        integer k,j
        do j=0,ngrd-1
            flx (ndt,j,k)=(dens(ndt,j,k)-dens(ndt,j+1,k))/dx(j)
        enddo
        flx (ndt,ngrd,k)=flx (ndt,ngrd-1,k)
    end subroutine update_flux

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
      call update_source(k)
      call update_flux(k)
      call update_reaction(k)
!     write(iout,*) "784: j3l",j3l(ndt,k)
      call update_surface(k)


!     --- low-pass filter ---
       rtsl(ndt,k)=delta*rtsl(ndt-1,k)+(1.d0-delta)*rtsl(ndt,k)
       rtsr(ndt,k)=delta*rtsr(ndt-1,k)+(1.d0-delta)*rtsr(ndt,k)
!
!     --- erosion ---
       ero(ndt,0,k)=-cero*flx(ndt,0,k)
       do j=1,ngrd-1
!        ero(ndt,j,k)=-cero*(dx(j)*flx(ndt,j-1,k)+dx(j-1)*flx(ndt,j,k))
!     +                    /(dx(j-1)+dx(j))
        ero(ndt,j,k)=-cero*flx(ndt,j,k)
       enddo
!       ero(ndt,ngrd,k)=-cero*flx(ndt,ngrd,k)
       ero(ndt,ngrd,k)=0.d0
!
       do j=0,ngrd
        rtd(ndt,j,k)=ero(ndt,j,k)
  !      if (isnan(rtd(ndt,j,k))) then
 !       write(iout,*) "804: rtd is nan"
 !       stop
  !      endif
       enddo
!     --- diffusion ---
       flx(ndt,0,k)=cdif(ndt,0,k)*flx(ndt,0,k)
!       write(iout,*) "811: j3l",j3l(ndt,k)
       rtd(ndt,0,k)=rtd(ndt,0,k)+(j3l(ndt,k)-j4l(ndt,k)-flx(ndt,0,k))*2.d0/dx(0)
!       if (isnan(rtd(ndt,0,k))) then
!        write(iout,*) "813: rtd is nan",flx(ndt,0,k),dx(0),j4l(ndt,k),j3l(ndt,k)
 !       endif
       do j=1,ngrd-1
        flx(ndt,j,k)=cdif(ndt,j,k)*flx(ndt,j,k)
        rtd(ndt,j,k)=rtd(ndt,j,k)+(flx(ndt,j-1,k)-flx(ndt,j,k))/(0.5d0*(dx(j-1)+dx(j)))
       enddo
       flx(ndt,ngrd,k)=cdif(ndt,ngrd,k)*flx(ndt,ngrd,k)
       rtd(ndt,ngrd,k)=rtd(ndt,ngrd,k)+(j3r(ndt,k)-j4r(ndt,k)+flx(ndt,ngrd,k))*2.d0/dx(ngrd-1)
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

end module modFACE_solver
