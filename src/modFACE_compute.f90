module modFACE_compute
    use modFACE_header
    use modFACE_cap
    use modFACE_error
    use modFACE_functions
    implicit none
contains

    subroutine compute_f(u,f)

        integer i, j, k
        real(DP),intent(in):: u(:)
        real(DP),intent(out):: f(:)

        call compute_eq_rates(u)
        ! ** check if NaN exist in u
        call check_isNaN(u,"u compute_f")

        ! ** build f
        i=0
        ! *** loop over each species
        do k=1,nspc
            i=i+1

            ! **** left surface (dsrf)
            if (steady_state .eq. "no") then
                !     --- 1st order BDF ---
                if (order_solver.eq.1) then
                    f(i)=u(i)-a11*dsrfl(ndt-1,k)
                    f(i)=f(i)-a12*Gsrf_l (ndt  ,k)*dt_face

                !     --- 2nd order BDF ---
                elseif (order_solver.eq.2) then
                    f(i)=u(i)-a21*dsrfl(ndt-1,k) &
                        -a22*dsrfl(ndt-2,k)
                    !write(iout,'(a,i3,a,es12.3)',advance="no") "i=",i,"; f(i)=",f(i)
                    f(i)=f(i)-a23*Gsrf_l (ndt  ,k)*dt_face
                    !write(iout,'(a,es12.3)') " Gsrf_l (ndt  ,k)*dt_face =",Gsrf_l (ndt  ,k)*dt_face
                !     --- 5th order BDF ---
                elseif (order_solver.eq.5) then
                    f(i)=u(i)-a51*dsrfl(ndt-1,k)&
                        -a52*dsrfl(ndt-2,k)&
                        -a53*dsrfl(ndt-3,k)&
                        -a54*dsrfl(ndt-4,k)
                    f(i)=f(i)-a55*Gsrf_l (ndt  ,k)*dt_face
                else
                    call face_error("ERROR: order of solver not implemented. order=",order_solver)
                endif

            elseif(steady_state .eq. "yes") then
                if (dsrfl(ndt,k) .gt. 1.d0) then
                    f(i)=Gsrf_l(ndt,k)
                else
                    f(i)=Gsrf_l(ndt,k)*dsrfl(ndt,k)
                endif
            else
                call face_error("Unknown mode for steady_state:",steady_state)
            endif

            ! **** bulk (dens)
            do j=0,ngrd
                !     --- step ---
                i=i+1
                if (steady_state .eq. "no") then
                    !     --- 1st order BDF ---
                    if (order_solver.eq.1) then
                        f(i)=u(i)-a11*dens(ndt-1,j,k)
                        f(i)=f(i)-a12*rate_d (ndt  ,j,k)*dt_face
                       ! write(iout,*) 'rate_d(ndt,0,k)*dt_face',rate_d(ndt,0,k)*dt_face,'dens(ndt-1,j,k)=',dens(ndt-1,j,k)
                    !     --- 2nd order BDF ---
                    elseif (order_solver.eq.2) then
                        f(i)=u(i)-a21*dens(ndt-1,j,k)&
                                 -a22*dens(ndt-2,j,k)
                        f(i)=f(i)-a23*rate_d (ndt  ,j,k)*dt_face
                    !     --- 5th order BDF ---
                    elseif (order_solver.eq.5) then
                        f(i)=u(i)-a51*dens(ndt-1,j,k)&
                            -a52*dens(ndt-2,j,k)&
                            -a53*dens(ndt-3,j,k)&
                            -a54*dens(ndt-4,j,k)
                        f(i)=f(i)-a55*rate_d (ndt  ,j,k)*dt_face

                    else
                        call face_error("ERROR: order of solver not implemented. order=",order_solver)
                    endif
                elseif (steady_state .eq. "yes") then
                    if (dens(ndt,j,k) .gt. 1.d0) then
                        f(i)=rate_d(ndt,j,k)
                    else
                        f(i)=rate_d(ndt,j,k)*dens(ndt,j,k)
                    endif

                else
                    call face_error("Unknown mode for steady_state:",steady_state)
                end if
            enddo

            ! **** right surface (dsrfr)
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
                     !   write(iout,'(a,i3,a,es12.3)',advance="no") "i=",i,"; f(i)=",f(i)
                    f(i)=f(i)-a23*Gsrf_r (ndt  ,k)*dt_face
                    ! write(iout,'(a,es12.3)') " Gsrf_r (ndt  ,k)*dt_face =",Gsrf_r (ndt  ,k)*dt_face
                !     --- 5th order BDF ---
                elseif (order_solver.eq.5) then
                    f(i)=u(i)-a51*dsrfr(ndt-1,k)&
                        -a52*dsrfr(ndt-2,k)&
                        -a53*dsrfr(ndt-3,k)&
                        -a54*dsrfr(ndt-4,k)
                    f(i)=f(i)-a55*Gsrf_r (ndt  ,k)*dt_face

                else
                    call face_error("ERROR: order of solver not implemented. order=",order_solver)
                endif
            elseif(steady_state .eq. "yes") then
                if (dsrfr(ndt,k) .gt. 1.d0) then
                    f(i)=Gsrf_r(ndt,k)
                else
                    f(i)=Gsrf_r(ndt,k)*dsrfr(ndt,k)
                endif
            else
                call face_error("Unknown mode for steady_state:",steady_state)
            endif

        enddo
        ! ** temperature
        if (solve_heat_eq) then
            do j=0,ngrd
                i=i+1
                if (steady_state .eq. "no") then

                        !     --- 1st order BDF ---
                    if (order_solver.eq.1) then
                        f(i)=u(i)-a11*temp(ndt-1,j)
                        f(i)=f(i)-a12*rate_t (ndt  ,j)*dt_face
                    !     --- 2nd order BDF ---
                    elseif (order_solver.eq.2) then
                        f(i)=u(i)-a21*temp(ndt-1,j)&
                            -a22*temp(ndt-2,j)
                        f(i)=f(i)-a23*rate_t (ndt  ,j)*dt_face

                    !     --- 5th order BDF ---
                    elseif (order_solver.eq.5) then
                        f(i)=u(i)-a51*temp(ndt-1,j) &
                            -a52*temp(ndt-2,j)  &
                            -a53*temp(ndt-3,j)  &
                            -a54*temp(ndt-4,j)
                        f(i)=f(i)-a55*rate_t (ndt  ,j)*dt_face
                    else
                        call face_error("ERROR: order of solver not implemented. order=",order_solver)
                    endif
                elseif(steady_state .eq. "yes") then
                    if (temp(ndt,j) .gt. 1.d0) then
                        f(i)=rate_t(ndt,j)
                    else
                        f(i)=rate_t(ndt,j)*temp(ndt,j)
                    endif
                else
                    call face_error("Unknown mode for steady_state:",steady_state)
                endif
            enddo

        endif
        ! ** check if the size of f is consistant with the number of equations
        if (i .ne. neq) then
            call face_error('Wrong number of functions size(f)=', i, ' neq=',neq)
        endif
        ! ** check if NaN exist in f
        call check_isNaN(f,"f compute_f")

    end subroutine compute_f

    subroutine compute_gradient_dens(k)

        integer k,j
         flx (ndt,0,k)=(dens(ndt,0,k)-dens(ndt,1,k))/dx(0)
        do j=1,ngrd-1

            flx (ndt,j,k)=(dens(ndt,j,k)-dens(ndt,j+1,k))/dx(j)
          !  flx (ndt,j,k)=(dens(ndt,j-1,k)-dens(ndt,j+1,k))/(dx(j-1)+dx(j))
        enddo
        flx (ndt,ngrd,k)=flx (ndt,ngrd-1,k)

    end subroutine compute_gradient_dens

    subroutine compute_gradient_temp

        integer j
        do j=0,ngrd-1
            qflx(ndt,j)=(temp(ndt,j)-temp(ndt,j+1))/dx(j)
        enddo
        qflx(ndt,ngrd)=qflx(ndt,ngrd-1)

    end subroutine compute_gradient_temp

    subroutine compute_ero_flx(k)

        integer j,k

        ero_flx(ndt,0,k)=-cero*flx(ndt,0,k) ! until now, flx is the density gradient (see line 764)
        do j=1,ngrd-1
            !        ero(ndt,j,k)=-cero*(dx(j)*flx(ndt,j-1,k)+dx(j-1)*flx(ndt,j,k))/(dx(j-1)+dx(j))
            ero_flx(ndt,j,k)=-cero*flx(ndt,j,k)
        enddo
        !       ero(ndt,ngrd,k)=-cero*flx(ndt,ngrd,k)
        ero_flx(ndt,ngrd,k)=0.d0 ! no erosion on the right side of the material (only left side is facing plasma)

    end subroutine compute_ero_flx

    subroutine compute_diff_flx(k)
        integer j,k
        do j=0,ngrd
            dif_flx(ndt,j,k)=cdif(ndt,j,k)*flx(ndt,j,k)
            flx(ndt,j,k)=cdif(ndt,j,k)*flx(ndt,j,k)
        enddo

    end subroutine compute_diff_flx

    subroutine compute_ero_qflx

        integer j

        ero_qflx(ndt,0)=-cero*qflx(ndt,0)
        do j=1,ngrd-1
            !        ero_qflx(ndt,j)=-cero*(dx(j)*qflx(ndt,j-1)+dx(j-1)*qflx(ndt,j))/(dx(j-1)+dx(j))
            ero_qflx(ndt,j)=-cero*qflx(ndt,j)
        enddo
        !       ero_qflx(ndt,ngrd)=-cero*qflx(ndt,ngrd)
        ero_qflx(ndt,ngrd)=0.d0

    end subroutine compute_ero_qflx

    subroutine compute_thermal_qflx

        integer j
        do j=0,ngrd
            qflx(ndt,j)=thcond*qflx(ndt,j)
        enddo

    end subroutine compute_thermal_qflx

    subroutine compute_arrhenius_coeffs

        integer j,k,l,m

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



    end subroutine compute_arrhenius_coeffs

    subroutine compute_source
        integer::j,l,k
        !     --- source modification ---
        call compute_inflx
!        if (.not.solve_heat_eq) then
!
!            do k=1,nspc
!                do j=0,ngrd
!                    srs(ndt,j,k)=srs(ndt-1,j,k)
!                    do l=1,nspc
!                        srb(ndt,j,k,l)=srb(ndt-1,j,k,l)
!                    enddo
!                enddo
!            enddo
!
!        elseif (solve_heat_eq) then
            do k=1,nspc
                do j=0,ngrd
                    srs(ndt,j,k)=source(j,k)
                    do l=1,nspc
                        srb(ndt,j,k,l)=srcbin(j,k,l)
                    enddo
                enddo
            enddo
 !       endif

    end subroutine compute_source

    subroutine compute_source_rate(k)
        integer::j,l
        integer,intent(in) :: k
        real(DP) csrs, csrb
        !     --- sources ---
        if (.not.solve_heat_eq) then
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
        elseif (solve_heat_eq) then
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

        endif
    end subroutine compute_source_rate

    subroutine compute_surface_flx(k)
        integer::k
        real(DP):: tmp,Edes_rc,Edes_lc
        !     --- surface ---

        ! reducing qch when close to surface saturation
        if (active_cap) then
            tmp=1.d2*(dsrfl0(k)/dsrfm(k)-1.d0)
            Edes_lc=Edes_l(k)*0.5d0*(1.d0-erf(tmp))
            tmp=1.d2*(dsrfr0(k)/dsrfm(k)-1.d0)
            Edes_rc=Edes_r(k)*0.5d0*(1.d0-erf(tmp))
        else
            Edes_lc=Edes_l(k)
            Edes_rc=Edes_r(k)
        endif
        !        tmp=1.d2*(dsrfl(ndt,k)/dsrfm(k)-1.d0)
        !        Edes_lc=Edes_l(k)*0.5d0*(1.d0-erf(tmp))
        !        tmp=1.d2*(dsrfr(ndt,k)/dsrfm(k)-1.d0)
        !        Edes_rc=Edes_r(k)*0.5d0*(1.d0-erf(tmp))
        ! calculate rates of surface processes
        ! - left surface
        if (left_surface_model(k).eq."S") then

            Kabs_l(k)=j0(k)*K0abs_l(k)*exp(-  ee*Eabs_l(k) /(kb*temp(ndt,0)))
            Kdes_l(k)=2.d0 *K0des_l(k)*exp(-  ee*Edes_lc /(kb*temp(ndt,0)))
            Kb_l(k)=        K0b_l(k)  *exp(-  ee*Eb_l(k)   /(kb*temp(ndt,0)))
            Kads_l(k)=      K0ads_l(k) *exp(-  ee*Eads_l(k) /(kb*temp(ndt,0)))
        elseif (left_surface_model(k).eq."B") then
            Kabs_l(k)=min_rate_surface
            Kdes_l(k)=2.d0 *K0des_l(k)*exp(-  ee*Edes_lc /(kb*temp(ndt,0)))
            Kb_l(k)= min_rate_surface
            Kads_l(k)= min_rate_surface
        elseif (left_surface_model(k).eq."N") then
            Kabs_l(k)=0d0
            Kdes_l(k)=0d0
            Kb_l(k)=0d0
            Kads_l(k)=0d0
        else
            call face_error("unknown left surface model:",left_surface_model(k))
        endif



        ! - right surface
        if (right_surface_model(k).eq."S") then

            Kabs_r(k)=j0(k)*K0abs_r(k)*exp(-  ee*Eabs_r(k) /(kb*temp(ndt,0)))
            Kdes_r(k)=2.d0 *K0des_r(k)*exp(-  ee*Edes_rc /(kb*temp(ndt,0)))
            Kb_r(k)=        K0b_r(k)  *exp(-  ee*Eb_r(k)   /(kb*temp(ndt,0)))
            Kads_r(k)=      K0ads_r(k)*exp(-  ee*Eads_r(k) /(kb*temp(ndt,0)))
        elseif (right_surface_model(k).eq."B") then
            Kabs_r(k)=0d0
            Kdes_r(k)=2.d0 *K0des_r(k)*exp(-  ee*Edes_rc /(kb*temp(ndt,0)))
            Kb_r(k)= 0d0
            Kads_r(k)=0d0
        elseif (right_surface_model(k).eq."N") then
            Kabs_r(k)=0d0
            Kdes_r(k)=0d0
            Kb_r(k)=0d0
            Kads_r(k)=0d0
        else
            call face_error("unknown right surface model:",right_surface_model(k))
        endif

                ! calculate surface fluxes using rates and density
                ! - left surface


        if ((left_surface_model(k).eq."S")) then
            Gabs_l (ndt,k)=Kabs_l(k)
            Gdes_l (ndt,k)=Kdes_l(k)*dsrfl(ndt,k)**order_desorption_left(k)

            Gb_l (ndt,k)  =Kb_l(k)  *dsrfl(ndt,k)
            Gads_l (ndt,k)=Kads_l(k)*dens(ndt,0   ,k)

            if (verbose_surface) then
                write(iout,*) 'Kdesl=', Kdes_l(k), 'Kdesl=', K0des_l(k),';ns^order=',dsrfl(ndt,k)**order_desorption_left(k)
                write(iout,*) 'Kadsl=', Kads_l(k), ';K0ads_l(k)',K0ads_l(k),'dens(ndt,0   ,k)=',dens(ndt,0   ,k)
            endif
       elseif (left_surface_model(k).eq."N") then
            Gabs_l (ndt,k)=0d0
            Gdes_l (ndt,k)=0d0

            Gb_l (ndt,k)  =dsrfl(ndt,k)
            Gads_l (ndt,k)=0d0
        elseif (left_surface_model(k).eq."B") then

            Gabs_l (ndt,k)=0d0
            Gdes_l (ndt,k)=Kdes_l(k)*dens(ndt,0   ,k)**order_desorption_left(k)
            if (verbose_surface) then
                write(iout,*) 'Kdesl=', Kdes_l(k),'K0desl=', K0des_l(k), 'n^order=',dens(ndt,0,k)**order_desorption_left(k)
            endif
            Gb_l (ndt,k)  =dsrfl(ndt,k)
            Gads_l (ndt,k)=0d0

        endif

        ! - right surface
        if ((right_surface_model(k).eq."S"))then
            Gabs_r (ndt,k)=Kabs_r(k)                           ! Gabsorp=K(gas)
            Gdes_r (ndt,k)=Kdes_r(k)*dsrfr(ndt,k)**order_desorption_right(k)          ! Gdesorp=K*ns^2
            Gb_r (ndt,k)  =Kb_r(k)  *dsrfr(ndt,k)              ! Gbulk  =K*ns
            Gads_r (ndt,k)=Kads_r(k)*dens(ndt,ngrd,k)          ! Gadsorb=K*nb
        elseif (right_surface_model(k).eq."N") then
              Gabs_r (ndt,k)=0d0                           ! Gabsorp=K(gas)
            Gdes_r (ndt,k)=Kdes_r(k)*dsrfr(ndt,k)**order_desorption_right(k)          ! Gdesorp=K*ns^2
            Gb_r (ndt,k)  =0d0               ! Gbulk  =K*ns
            Gads_r (ndt,k)=0d0
        elseif (right_surface_model(k).eq."B") then
            Gabs_r (ndt,k)=0d0                          ! Gabsorp=K(gas)
            Gdes_r (ndt,k)=Kdes_r(k)*dens(ndt,ngrd,k)**order_desorption_right(k)          ! Gdesorp=K*ns^2
            Gb_r (ndt,k)  =dsrfr(ndt,k)             ! Gbulk  =K*ns
            Gads_r (ndt,k)=0d0      ! Gadsorb=K*n
        endif
        ! apply cap factor to mimic effects of saturation
        if (active_cap) then
            call compute_cap_factor_surface(k,ndt)
        endif

        ! calculate effective desorptiopn and heat fluxes
        if (solve_heat_eq) then
            jout(ndt,k)=jout(ndt,k)+Gdes_l(ndt,k)
            qflx_in=qflx_in+jout(ndt,k)*(ee*Eads_l(k)-2.d0*kb*temp(ndt,0))
        endif

        ! - net flux onto surface
        if ((left_surface_model(k).eq."S"))then
            Gsrf_l(ndt,k)=Gabs_l(ndt,k)-Gdes_l(ndt,k)-Gb_l(ndt,k)+Gads_l(ndt,k)
            elseif (left_surface_model(k).eq."N") then
            Gsrf_l(ndt,k)=-Gb_l(ndt,k)
        elseif (left_surface_model(k).eq."B") then
            Gsrf_l(ndt,k)=-Gb_l(ndt,k)
        endif

        if ((right_surface_model(k).eq."S")) then
        Gsrf_r(ndt,k)=Gabs_r(ndt,k)-Gdes_r(ndt,k)-Gb_r(ndt,k)+Gads_r(ndt,k)
        elseif (right_surface_model(k).eq."N") then
            Gsrf_r(ndt,k)=-Gb_r(ndt,k)
        elseif (right_surface_model(k).eq."B") then
            Gsrf_r(ndt,k)=-Gb_r(ndt,k)
        endif

        if (verbose_surface) then
            write(iout,*) 'Gsrf_l(ndt,k)=',Gsrf_l(ndt,k),'; Gsrf_r(ndt,k)=',Gsrf_r(ndt,k),&
                'dsrfl(ndt,k)=',dsrfl(ndt,k),'dsrfr(ndt,k)=',dsrfr(ndt,k)
        endif

        !     --- low-pass filter ---
        delta=0d0
        Gsrf_l(ndt,k)=delta*Gsrf_l(ndt-1,k)+(1.d0-delta)*Gsrf_l(ndt,k)
        Gsrf_r(ndt,k)=delta*Gsrf_r(ndt-1,k)+(1.d0-delta)*Gsrf_r(ndt,k)

    end subroutine compute_surface_flx

    subroutine compute_reaction_rate(k)
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
    end subroutine compute_reaction_rate

    subroutine compute_temp_rate
        integer j
        rate_t(ndt,0)=(qflx_in-qflx(ndt,0))*2.d0/dx(0)/rhocp+ero_qflx(ndt,0)
        do j=1,ngrd-1
            rate_t(ndt,j)=(qflx(ndt,j-1)-qflx(ndt,j))/(0.5d0*(dx(j-1)+dx(j))*rhocp)+ero_qflx(ndt,j)
        enddo
        rate_t(ndt,ngrd)=0.d0

        do j=0,ngrd
            !     --- low-pass filter ---
            rate_t(ndt,j)=delta*rate_t(ndt-1,j)+(1.d0-delta)*rate_t(ndt,j)
        enddo



    end subroutine compute_temp_rate

    subroutine compute_dens_rate(k)
        integer j,k
        do j=0,ngrd
            rate_d(ndt,j,k)=ero_flx(ndt,j,k)
        enddo

        if (left_surface_model(k).eq."S") then
         !rate_d(ndt,0,k)=rate_d(ndt,0,k)+(Gb_l(ndt,k)-Gads_l(ndt,k)-dif_flx(ndt,0,k))/dx(0)
         rate_d(ndt,0,k)=rate_d(ndt,0,k)+(Gb_l(ndt,k)-Gads_l(ndt,k)-dif_flx(ndt,0,k))/dx(0)
        elseif (left_surface_model(k).eq."N") then
            !rate_d(ndt,0,k)=rate_d(ndt,0,k)+(Gb_l(ndt,k)-Gads_l(ndt,k)-dif_flx(ndt,0,k))*2.d0/dx(0)
            rate_d(ndt,0,k)=rate_d(ndt,0,k)-dif_flx(ndt,0,k)/dx(0)
        elseif(left_surface_model(k).eq."B") then
            !rate_d(ndt,0,k)=rate_d(ndt,0,k)+(-Gdes_l(ndt,k)-dif_flx(ndt,0,k))*2.d0/dx(0)
             !rate_d(ndt,0,k)=rate_d(ndt,0,k)+(-Gdes_l(ndt,k)-dif_flx(ndt,0,k))/dx(0)
             rate_d(ndt,0,k)=rate_d(ndt,0,k)+(-Gdes_l(ndt,k)-dif_flx(ndt,0,k))/dx(0)
        !write(iout,*) 'rate_d(ndt,0,k)',rate_d(ndt,0,k),'Gdes_l(ndt,k)',Gdes_l(ndt,k),'dif_flx(ndt,0,k)',dif_flx(ndt,0,k)
        endif

        do j=1,ngrd-1
            !rate_d(ndt,j,k)=rate_d(ndt,j,k)+(dif_flx(ndt,j-1,k)-dif_flx(ndt,j,k))/(0.5d0*(dx(j-1)+dx(j)))
            rate_d(ndt,j,k)=rate_d(ndt,j,k)+(dif_flx(ndt,j-1,k)-dif_flx(ndt,j,k))/dx(j-1)
            !rate_d(ndt,j,k)=rate_d(ndt,j,k)+(dif_flx(ndt,j-1,k)-dif_flx(ndt,j+1,k))/((dx(j-1)+dx(j)))
        enddo

        if (right_surface_model(k).eq."S") then
        !rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+(Gb_r(ndt,k)-Gads_r(ndt,k)+dif_flx(ndt,ngrd,k))/dx(ngrd-1)
        rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+(Gb_r(ndt,k)-Gads_r(ndt,k)+dif_flx(ndt,ngrd,k))/dx(ngrd-1)
        elseif(right_surface_model(k).eq."N") then
            !rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+(Gb_r(ndt,k)-Gads_r(ndt,k)+dif_flx(ndt,ngrd,k))*2.d0/dx(ngrd-1)
            rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+dif_flx(ndt,ngrd,k)/dx(ngrd-1)
        elseif(right_surface_model(k).eq."B") then
        !rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+(-Gdes_r(ndt,k)+dif_flx(ndt,ngrd,k))/dx(ngrd-1)
            rate_d(ndt,ngrd,k)=rate_d(ndt,ngrd,k)+(-Gdes_r(ndt,k)+dif_flx(ndt,ngrd,k))/dx(ngrd-1)

        endif
        do j=0,ngrd
            !     --- sources ---
            rate_d(ndt,j,k)=rate_d(ndt,j,k)+src(ndt,j,k)
            !     --- reactions ---
            rate_d(ndt,j,k)=rate_d(ndt,j,k)+rct(ndt,j,k)
            !     --- low-pass filter ---
            delta=0d0
            rate_d(ndt,j,k)=delta*rate_d(ndt-1,j,k)+(1.d0-delta)*rate_d(ndt,j,k)

        enddo
    end subroutine compute_dens_rate


    subroutine compute_cap_factor_surface(k,i)
        integer,intent(in)::k,i
        real(DP) :: cabs_l, cdes_l, cb_l, cads_l
        real(DP) :: cabs_r, cdes_r, cb_r, cads_r
        if (active_cap) then
            call cap_srf_flx(k,i,cabs_l, cdes_l, cb_l, cads_l,cabs_r, cdes_r, cb_r, cads_r)

            ! multiply surfaces fluxes by cap factors
            ! - left surface
            Gabs_l (i,k)=Gabs_l (i,k) *cabs_l
            Gdes_l (i,k)=Gdes_l (i,k) *cdes_l
            Gb_l (i,k)  =Gb_l (i,k)   *cb_l
            Gads_l (i,k)=Gads_l (i,k) *cads_l
            ! - right surface
            Gabs_r (i,k)=Gabs_r (i,k) *cabs_r
            Gdes_r (i,k)=Gdes_r (i,k) *cdes_r
            Gb_r (i,k)  =Gb_r (i,k)   *cb_r
            Gads_r (i,k)=Gads_r (i,k) *cads_r
        endif
    end subroutine compute_cap_factor_surface

    subroutine compute_inflx()
        implicit none

        integer k
        real(DP):: trel, tmp
        real(DP):: cenr
        parameter (cenr=1.d0)


        ! pulsed plasma flow
        do k=1,nspc
            if (inflx_in_pulse(k).eq."R") then
                ! check that the pulse period is not zero
                if (inflx_in_pulse_period(k).le.0) then
                    call face_error('Pulsed incoming plasma flux activated but pulse_period =0')
                endif

                trel=time-inflx_in_pulse_period(k)*int(time/inflx_in_pulse_period(k))
                if (trel .le. t1) then ! phase 1 of pulse
                    tmp=trel/t1

                    inflx(k)=inflx_in(k)+(inflx_in_max(k)-inflx_in(k))*tmp

                    rad =rad_min +(rad_max -rad_min )*tmp
                    cero=cero_min+(cero_max-cero_min)*tmp
                elseif (trel .le. t2) then ! phase 2 of pusle

                    inflx(k)=inflx_in_max(k)

                    rad =rad_max
                    cero=cero_max
                elseif (trel .le. t3) then ! phase 3 of pulse
                    tmp=(trel-t2)/(t3-t2)

                    inflx(k)=inflx_in_max(k)+(inflx_in(k)-inflx_in_max(k))*tmp

                    rad =rad_max +(rad_min -rad_max )*tmp
                    cero=cero_max+(cero_min-cero_max)*tmp
                endif
            elseif (inflx_in_pulse(k).eq."N") then ! no pulse
                inflx(k)=inflx_in(k)

            elseif (inflx_in_pulse(k).eq."S") then
                if (time.ge.inflx_in_pulse_starttime(k)) then
                    if (inflx_in_pulse_period(k).le.0d0) then
                        call face_error(' Gamma pulse_period <=0')
                    endif
                    inflx(k)=inflx_in(k)+(inflx_in_max(k)-inflx_in(k))*&
                    sin(2d0*pi*(time-inflx_in_pulse_starttime(k))/inflx_in_pulse_period(k))

                else
                inflx(k)=inflx_in(k)
                endif

            else
                call face_error("Unknown mode for Gamma_in_pulse k=",k)
            endif
        enddo

        rad =rad_min
        cero=cero_min


        if (gamero.ne.0d0) cero=cero+gamero*inflx(1)*(lambda**3*cvlm) ! sputtering

        ! TODO modif Q flux here
        qflx_in=0.d0
        do k=1,nspc
            qflx_in=qflx_in+ee*enrg(k)*inflx(k)
        enddo
        if (rad.ne.0d0)  qflx_in=qflx_in+rad

        if (cero.ne.0d0)  qflx_in=qflx_in-cero*ee*qform/(lambda**3*cvlm)
        if (emiss.ne.0d0) qflx_in=qflx_in-emiss*sigma_sb*(temp(ndt,0)**4.d0-temp(ndt,ngrd)**4.d0)

    end subroutine compute_inflx

    subroutine  compute_temperature
        integer j,n
        !     --- compute new temperature and update arrhenius coeffs in heat eq is not solved.
        ! if heat eq is solved, then arrhenius coefficients are computed within cpomputation of f functions

        if (.not.solve_heat_eq) then
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


            call compute_arrhenius_coeffs
        endif
    end subroutine compute_temperature

    subroutine check_isNaN(u,str)
        ! update equation terms
        real(DP),intent(in):: u(:)
        character(*) :: str
        integer :: i
        do i=1,neq
            if (isnan(u(i))) then
                call face_error(' ; location:',str,'NaN found: at i=' ,i )
            endif
        enddo
    end subroutine check_isNaN

    subroutine compute_eq_rates(u)
        ! update equation terms
        real(DP),intent(in):: u(:)
        integer :: k
        !      if (verbose_debug) write (iout,*) '-- compute_eq_rates'
        ! ** get values of density and temperature from vector u

        call get_density_values(u)

        call check_isNaN(u,'after get_density')
        ! ** check if calculated new values are positive and below max values

        call check_positivity_max

        ! ** update all equations rates (rate_t (temperature), rate_d (density) and surface flux)
        if (solve_heat_eq) then
            call compute_arrhenius_coeffs ! compute trapping,detrapping and diffusion coeff which depend on temperature. Ifheaqt eq not solved, this is done upfront in routine step
        endif
        ! *** density
        do k=1,nspc
            call compute_source_rate(k)             ! compute source term for species
            call compute_gradient_dens(k)      ! set flx (ndt,j,k)=(n(ndt,j+1,k)-n(ndt,j,k)/dx(j))
            call compute_reaction_rate(k)           !
            call compute_surface_flx(k)
            call compute_ero_flx(k)
            call compute_diff_flx(k)
            call compute_dens_rate(k)
        enddo

        ! *** heat conduction
        if (solve_heat_eq) then
            call compute_gradient_temp ! compute qflx (ndt,j,k)=(n(ndt,j+1,k)-n(ndt,j,k)/dx(j))
            call compute_ero_qflx      ! compute ero_qflx=cero*qflx
            call compute_thermal_qflx  ! compute qflx=kappa*qflx
            call compute_temp_rate    !
        endif

    end subroutine compute_eq_rates

    real(DP) function compute_fnorm(u) result(fnorm)
        integer i, idx
        real(DP),intent(in):: u(:)
        real(DP):: f(neq)
        real(DP):: norm, mxel, tmp

        call compute_f(u,f)

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
          !    write (iout,*) 'Max norm element ', mxel, ' at eq', idx
        fnorm=sqrt(norm)
    end function compute_fnorm

    subroutine check_positivity_max
        integer ::j,k
        real(DP)::tmp
        !     --- check for negative densities ---
        do k=1,nspc
            if (dsrfl(ndt,k) .lt. 0.d0) then
                call face_error("left srf dens <0 for species k= ",k,&
                " dsrfl=",dsrfl(ndt,k),"previous step: dsrfl(ndt-1)=",dsrfl(ndt-1,k))
            endif

            if (dsrfl(ndt,k) .gt. dsrfm(k)) then
                call face_error("left srf dens >max for species k= ",k,&
                " dsrfl=",dsrfl(ndt,k),"dsrfm=",dsrfm(k))
            endif

            if (dsrfr(ndt,k) .gt. dsrfm(k)) then
                call face_error("right srf dens >max for species k= ",k,&
                " dsrfr=",dsrfr(ndt,k),"dsrfm=",dsrfm(k))
            endif

            if (dsrfr(ndt,k) .lt. 0.d0) then
                call face_error("right srf dens <0 for species k= ",k,&
                " dsrfr=",dsrfr(ndt,k),"previous step: dsrf(ndt-1)=",dsrfr(ndt-1,k))
            endif

            do j=0,ngrd
                if (dens(ndt,j,k) .lt. 0.d0) then
                    call face_error("dens <0 for species k= ",k," cell j=",j,&
                    " dens=",dens(ndt,j,k),"previous step: dens(ndt-1)=",dens(ndt-1,j,k))
                endif

            enddo
        enddo

        !     --- check for density exceeding maximum ---
        do j=0,ngrd
            if (dens(ndt,j,1) .gt. densm(1)) then
                call face_warning('hydrogen density exceeded the limit j=',j, " dens=",dens(ndt,j,1))
                dens(ndt,j,1)=densm(1)
            endif
            !
            do k=2,(nspc-1),2
                if ((dens(ndt,j,k)+dens(ndt,j,k+1)) .gt.(densm(k)+densm(k+1))) then
                    call face_warning('Trap density exceeded the limit k=',k,' j=' ,j,&
                        ' dens tr(k)+otr(k+1)=',dens(ndt,j,k)+dens(ndt,j,k+1))
                    tmp=(densm(k)+densm(k+1))/(dens(ndt,j,k)+dens(ndt,j,k+1))
                    dens(ndt,j,k  )=tmp*dens(ndt,j,k  )
                    dens(ndt,j,k+1)=tmp*dens(ndt,j,k+1)
                endif
            enddo
        enddo

        !    --- check for negative temperatures and update coefficients ---
        if (solve_heat_eq) then
            do j=0,ngrd
                if (temp(ndt,j) .lt. 0.d0) then
                    call face_error("negative temperature in cell",j,"T=",temp(ndt,j)," previous step:T(ndt-1)=",temp(ndt-1,j))
                endif

            enddo
        endif

    end subroutine check_positivity_max

    subroutine get_density_values(u)
        integer ::i,j,k
        real(DP),intent(in):: u(:)
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
        !      write(iout,*) "i=",i, " ;dsrfr(ndt,k)=" , dsrfr(ndt,k)

        enddo
        if (solve_heat_eq) then
            do j=0,ngrd
                i=i+1
                temp(ndt,j)=u(i)
            enddo

        endif

    end subroutine get_density_values

    subroutine compute_trace_flux
        ! sum up the outgassing flux left and right over time to estimate the average outgassing flux over the simulation
        ! end verify if <Gdes_l>\approx Gdes_l(end)
        real(DP):: int_src
        integer k,j

        do k=1,nspc
            int_src=0.d0
            do j=0,ngrd-1
                int_src=int_src+source(j,k)*dx(j)
            enddo
            trace_flux(k)%sum_inflx=trace_flux(k)%sum_inflx+int_src*dt_face
            trace_flux(k)%sum_Gdes_l=trace_flux(k)%sum_Gdes_l+Gdes_l(ndt,k)*dt_face
            trace_flux(k)%sum_Gdes_r=trace_flux(k)%sum_Gdes_r+Gdes_r(ndt,k)*dt_face
            trace_flux(k)%sig_Gdes_l=trace_flux(k)%sig_Gdes_l+Gdes_l(ndt,k)**2
            trace_flux(k)%sig_Gdes_r=trace_flux(k)%sig_Gdes_r+Gdes_r(ndt,k)**2
            trace_flux(k)%max_Gdes_l=max(trace_flux(k)%max_Gdes_l,Gdes_l(ndt,k))
            trace_flux(k)%max_Gdes_r=max(trace_flux(k)%max_Gdes_r,Gdes_r(ndt,k))
            trace_flux(k)%min_Gdes_l=min(trace_flux(k)%min_Gdes_l,Gdes_l(ndt,k))
            trace_flux(k)%min_Gdes_r=min(trace_flux(k)%min_Gdes_r,Gdes_r(ndt,k))
        enddo
    end subroutine compute_trace_flux

    subroutine compute_onthefly_inventory
        ! sum up the outgassing flux left and right over time to estimate the average outgassing flux over the simulation
        ! end verify if <Gdes_l>\approx Gdes_l(end)
        real(DP):: int_src,int_dens,int_dsrf,int_des
        integer k,j
         if (print_onthefly_inventory) then
        do k=1,nspc
            int_src=0.d0
            do j=0,ngrd-1
                int_src=int_src+source(j,k)*dx(j)*dt_face
            enddo
            int_dens=0.d0

            do j=0,ngrd-1
                int_dens=int_dens+dens(ndt,j,k)*dx(j)
            enddo
            int_dsrf=0d0
            if ((left_surface_model(k).eq."S")) then
            int_dsrf=dsrfl(ndt,k)+int_dsrf
            endif

            if ((right_surface_model(k).eq."S")) then
            int_dsrf=dsrfr(ndt,k)+int_dsrf
            endif

            int_des=Gdes_l(ndt,k)*dt_face+Gdes_r(ndt,k)*dt_face

            onthefly_inventory(k)%net_int_dens=int_dens-onthefly_inventory(k)%int_dens
            onthefly_inventory(k)%net_int_dsrf=int_dsrf-onthefly_inventory(k)%int_dsrf
            onthefly_inventory(k)%int_dens=int_dens
            onthefly_inventory(k)%int_src=int_src
            onthefly_inventory(k)%int_dsrf=int_dsrf
            onthefly_inventory(k)%int_des=int_des
        enddo
        endif

    end subroutine compute_onthefly_inventory

    subroutine compute_dt_update

        integer  ::i, j, k
        real(DP) :: tmp,dt
        real(DP),parameter:: rtfm=5.d-1, rtfp=5.d-1
        logical :: adjust_reduction_factor=.false.
        tmp=0d0

        if ((variable_timestep)) then

        if (adjust_reduction_factor) then
        if ((iter_solver.gt.2)) then
        reduction_factor_dt=reduction_factor_dt/15d0
        else
        reduction_factor_dt=reduction_factor_dt*1.2d0
        endif

        endif


        if (reduction_factor_dt.gt.1d0) then
        reduction_factor_dt=1d0
        endif


            dt=end_time-time

            ! bulk
            do k=1,nspc
                do j=0,ngrd
                    if (rate_d(ndt,j,k) .ne. 0.d0) then
                        if (rate_d(ndt,j,k) .lt. 0.d0) then
                            tmp=-rtfm*dens(ndt,j,k)/rate_d(ndt,j,k)*reduction_factor_dt
                        else
                            tmp= rtfp*dens(ndt,j,k)/rate_d(ndt,j,k)*reduction_factor_dt
                        endif
                    else
                        tmp=dt
                    endif

                    if ((dt .gt. tmp) .and. (tmp .ne. 0.d0)) then
                        dt=tmp
                    endif
                enddo

                ! left surface
                if ((left_surface_model(k).eq."S")) then
                    if (Gsrf_l (ndt  ,k) .ne. 0.d0) then
                        if (Gsrf_l (ndt  ,k) .lt. 0.d0) then
                            tmp=-rtfm*dsrfl(ndt,k)/Gsrf_l (ndt  ,k)*reduction_factor_dt
                        else
                            tmp= rtfp*dsrfl(ndt,k)/Gsrf_l (ndt  ,k)*reduction_factor_dt
                        endif
                    else
                        tmp=dt
                    endif
                    if ((dt .gt. tmp) .and. (tmp .ne. 0.d0)) then
                        dt=tmp
                    endif
                endif
                ! right surface
                if ((right_surface_model(k).eq."S") ) then
                    if (Gsrf_r (ndt  ,k) .ne. 0.d0) then
                        if (Gsrf_r (ndt  ,k) .lt. 0.d0) then
                            tmp=-rtfm*dsrfr(ndt,k)/Gsrf_r (ndt  ,k)*reduction_factor_dt
                        else
                            tmp= rtfp*dsrfr(ndt,k)/Gsrf_r (ndt  ,k)*reduction_factor_dt
                        endif
                    else
                        tmp=dt
                    endif
                    if ((dt .gt. tmp) .and. (tmp .ne. 0.d0)) then
                        dt=tmp
                    endif
                endif
            enddo



            ! heat equation
            if (solve_heat_eq) then
                do j=0,ngrd-1
                    if (rate_t(ndt,j) .ne. 0.d0) then
                        if (rate_t(ndt,j) .lt. 0.d0) then
                            tmp=-rtfm*temp(ndt,j)/rate_t(ndt,j)*reduction_factor_dt
                        else
                            tmp= rtfp*temp(ndt,j)/rate_t(ndt,j)*reduction_factor_dt
                        endif
                    else
                        tmp=dt
                    endif
                    if ((dt .gt. tmp) .and. (tmp .ne. 0.d0)) then
                        dt=tmp
                    endif
                enddo
            endif


            if (dt .lt. min_dt_face) dt=min_dt_face



        do k=1,nspc
            if (inflx_in_pulse(k).ne."N") then
                if (time.ge.inflx_in_pulse_starttime(k)) then
                    if (dt>inflx_in_pulse_period(k)/30d0) then
                        dt=inflx_in_pulse_period(k)/30d0
                    endif
                endif
            endif
        enddo


        dt_face=dt


    endif
end subroutine compute_dt_update
end module modFACE_compute
