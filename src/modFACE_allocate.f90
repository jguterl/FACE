      module modFACE_allocate
      use modFACE_header
      use modFACE_output
      implicit none
          interface init_zero
        module procedure init_zero_r
        module procedure init_zero_rarr
        module procedure init_zero_rarr2
        module procedure init_zero_rarr3
        module procedure init_zero_rarr4
        module procedure init_zero_rarr5
        module procedure init_zero_s
        module procedure init_zero_sarr
        module procedure init_zero_i
        module procedure init_zero_iarr
    end interface init_zero
      contains
         subroutine init_zero_rarr2(r)
        real(DP),allocatable::r(:,:)
!        integer:: n1,n2,i,j
         if (.not.allocated(r))  then
         call face_error('trying to set value to non-allocated variable')
         endif
!        n1=size(r,1)
!        n2=size(r,2)
!        do i=1,n1
!        do j=1,n2
        r=0.d0
!        enddo
!         enddo

    end subroutine init_zero_rarr2

      subroutine init_zero_rarr3(r)
        real(DP),allocatable::r(:,:,:)
!        integer:: n1,n2,n3,i,j,k
         if (.not.allocated(r))  then
         call face_error("trying to set value to non-allocated variable")
         endif
!        n1=size(r,1)
!        n2=size(r,2)
!        n3=size(r,3)
!        do i=1,n1
!        do j=1,n2
!        do k=1,n3
        r=0.d0
!        enddo
!         enddo
!         enddo

    end subroutine init_zero_rarr3

     subroutine init_zero_rarr4(r)
        real(DP),allocatable::r(:,:,:,:)
!        integer:: n1,n2,n3,n4,i,j,k,l
         if (.not.allocated(r))  then
         call face_error("trying to set value to non-allocated variable")
         endif
        r=0.d0

    end subroutine init_zero_rarr4

     subroutine init_zero_rarr5(r)
        real(DP),allocatable::r(:,:,:,:,:)
!        integer:: n1,n2,n3,n4,n5,i,j,k,l,m
         if (.not.allocated(r))  then
         call face_error("trying to set value to non-allocated variable")
         endif
        r=0.d0

    end subroutine init_zero_rarr5

    subroutine init_zero_rarr(r)
        real(DP),allocatable::r(:)
          if (.not.allocated(r))  then
         call face_error("trying to set value to non-allocated variable")
         endif
        r=0
    end subroutine init_zero_rarr

    subroutine init_zero_r(r)
        real(DP)::r
        r=0
    end subroutine init_zero_r


    subroutine init_zero_sarr(s)
        character(*)::s(:)
        s=''
    end subroutine init_zero_sarr

    subroutine init_zero_s(s)
        character(*)::s
        s=''
    end subroutine init_zero_s


    subroutine init_zero_i(i)
        integer::i
        i=0
    end subroutine init_zero_i

    subroutine init_zero_iarr(i)
        integer::i(:)
        i=0
    end subroutine init_zero_iarr

      subroutine alloc_variables()
          ! grid
          if (.not.allocated(x)) then
              allocate(x(0:ngrd))
              call init_zero(x)
          else
              call face_error('x already allocated')

          endif

          if (.not.allocated(dx)) then
              allocate(dx(0:ngrd))
              call init_zero(dx)
          else
              call face_error('dx already allocated')

          endif

          ! temperature
          if (.not.allocated(temp)) then
          allocate( temp(ndt,0:ngrd))
          call init_zero(temp)
          else
              call face_error('dx already allocated')
          endif
          ! surface

          allocate(inflx(nspc))

          ! bulk
          if (.not.allocated(dens)) then
              allocate(dens(ndt,0:ngrd,nspc))
              call init_zero(dens)
          else
              call face_error('dens already allocated')
          endif

          if (.not.allocated(flx)) then
              allocate(flx (ndt,0:ngrd,nspc))
              call init_zero(flx)
          else
              call face_error('flx already allocated')
          endif

          if (.not.allocated(rate_d)) then
              allocate(rate_d (ndt,0:ngrd,nspc))
              call init_zero(rate_d)
          else
              call face_error('rate_d already allocated')
          endif

          if (.not.allocated(cdif)) then
              allocate(cdif(ndt,0:ngrd,nspc))
              call init_zero(cdif)
          else
              call face_error('cdif already allocated')

          endif

          if (.not.allocated(rct)) then
              allocate(rct (ndt,0:ngrd,nspc))
              call init_zero(rct)
          else
             call face_error(' rct already allocated')
          endif
          if (.not.allocated(ero_flx)) then
              allocate(ero_flx (ndt,0:ngrd,nspc))
              call init_zero(ero_flx)
          else
              call face_error('ero_flx already allocated')

          endif

          if (.not.allocated(dif_flx)) then
              allocate(dif_flx (ndt,0:ngrd,nspc))
              call init_zero(dif_flx)
          else
              call face_error('dif_flx already allocated')

          endif

          if (.not.allocated(rate_t)) then
              allocate(rate_t(ndt,0:ngrd))
              call init_zero(rate_t)
          else
              call face_error('rate_t already allocated')

          endif
          if (.not.allocated(qflx)) then
              allocate(qflx(ndt,0:ngrd))
              call init_zero(qflx)
          else
              call face_error('qflx already allocated')

          endif
          if (.not.allocated(ero_qflx)) then
              allocate(ero_qflx(ndt,0:ngrd))
              call init_zero(ero_qflx)
          else
              call face_error('ero_qflx already allocated')

          endif

          if (.not.allocated(srs)) then
              allocate(srs(ndt,0:ngrd,nspc))
              call init_zero(srs)
          else
              call face_error('srs already allocated')

          endif

          if (.not.allocated(srb)) then
              allocate(srb (ndt,0:ngrd,nspc,nspc))
              call init_zero(srb)
          else
              call face_error('srb already allocated')

          endif
          if (.not.allocated(src)) then
              allocate(src (ndt,0:ngrd,nspc))
              call init_zero(src)
          else
              call face_error('src already allocated')

          endif

          if (.not.allocated(src_profile)) then
              allocate(src_profile (0:ngrd,nspc))
              call init_zero(src_profile)
          else
              call face_error('src_profile already allocated')

          endif

         ! reaction parameters
         if (.not.allocated(kbin0)) then
     allocate(kbin0(nspc,nspc,nspc))
         call init_zero(kbin0)
          else
              call face_error('kbin0 already allocated')

          endif
     if (.not.allocated(ebin)) then
      allocate(ebin(nspc,nspc,nspc))
          call init_zero(ebin)
          else
              call face_error('ebin already allocated')
          endif
      if (.not.allocated(kbin)) then
      allocate(kbin(ndt,0:ngrd,nspc,nspc,nspc))
          call init_zero(kbin)
          else
              call face_error('kbin already allocated')
          endif


      if (.not.allocated(nuth0)) then
       allocate(nuth0(nspc,nspc))
           call init_zero(nuth0)
          else
              call face_error('nuth0 already allocated')
          endif

       if (.not.allocated(eth)) then
      allocate(eth(nspc,nspc))
          call init_zero(eth)
          else
              call face_error('eth already allocated')
          endif
      if (.not.allocated(nuth)) then
      allocate(nuth(ndt,0:ngrd,nspc,nspc))
          call init_zero(nuth)
          else
              call face_error('nuth already allocated')

          endif

      ! surface
      if (.not.allocated(Kabs_l)) then
      allocate(Kabs_l(nspc))
                call init_zero(Kabs_l)
          else
              call face_error('Kabs_l already allocated')
          endif
          if (.not.allocated(Kdes_l)) then
      allocate(Kdes_l(nspc))
                call init_zero(Kdes_l)
          else
             call face_error('Kdes_l already allocated')
          endif
          if (.not.allocated(Kb_l)) then
      allocate(Kb_l(nspc))
                call init_zero(Kb_l)
          else
             call face_error('Kb_l already allocated')

          endif

          if (.not.allocated(Kads_l)) then
      allocate(Kads_l(nspc))
                call init_zero(Kads_l)
          else
              call face_error('Kads_l already allocated')
          endif


          if (.not.allocated(Kabs_r)) then
      allocate(Kabs_r(nspc))
                call init_zero(Kabs_r)
          else
              call face_error('Kabs_r already allocated')
          endif

          if (.not.allocated(Kdes_r)) then
      allocate(Kdes_r(nspc))
                call init_zero(Kdes_r)
          else
              call face_error('Kdes_r already allocated')
          endif

          if (.not.allocated(Kb_r)) then
      allocate(Kb_r(nspc))
                call init_zero(Kb_r)
          else
              call face_error('Kb_r already allocated')
          endif

          if (.not.allocated(Kads_r)) then
      allocate(Kads_r(nspc))
                call init_zero(Kads_r)
          else
              call face_error('Kads_r already allocated')
          endif

          if (.not.allocated(K0abs_l)) then
      allocate(K0abs_l(nspc))
                call init_zero(K0abs_l)
          else
             call face_error('K0abs_l already allocated')
          endif

          if (.not.allocated(K0des_l)) then
      allocate(K0des_l(nspc))
                call init_zero(K0des_l)
          else
              call face_error('K0des_l already allocated')
          endif

          if (.not.allocated(K0b_l)) then
      allocate(K0b_l(nspc))
                call init_zero(K0b_l)
          else
              call face_error('K0b_l already allocated')
          endif

          if (.not.allocated(K0ads_l)) then
      allocate(K0ads_l(nspc))
                call init_zero(K0ads_l)
          else
              call face_error('K0ads_l already allocated')
          endif

          if (.not.allocated(K0abs_r)) then
      allocate(K0abs_r(nspc))
                call init_zero(K0abs_r)
          else
              call face_error('K0abs_r already allocated')
          endif

          if (.not.allocated(K0des_r)) then
      allocate(K0des_r(nspc))
                call init_zero(K0des_r)
          else
              call face_error('K0des_r already allocated')
          endif

          if (.not.allocated(K0b_r)) then
      allocate(K0b_r(nspc))
                call init_zero(K0b_r)
          else
              call face_error('K0b_r already allocated')
          endif

          if (.not.allocated(K0ads_r)) then
      allocate(K0ads_r(nspc))
                call init_zero(K0ads_r)
          else
              call face_error('K0ads_r already allocated')
          endif

          if (.not.allocated(j0)) then
      allocate(j0(nspc))
                call init_zero(j0)
          else
              call face_error('j0 already allocated')
          endif

          if (.not.allocated(dsrfl)) then
      allocate(dsrfl(ndt,nspc))
                call init_zero(dsrfl)
          else
              call face_error('dsrfl already allocated')
          endif


          if (.not.allocated(Gsrf_l)) then
      allocate( Gsrf_l(ndt,nspc))
                call init_zero(Gsrf_l)
          else
              call face_error('Gsrf_l already allocated')
          endif

          if (.not.allocated(dsrfr)) then
      allocate(dsrfr(ndt,nspc))
                call init_zero(dsrfr)
          else
              call face_error('dsrfr already allocated')
          endif

          if (.not.allocated(Gsrf_r)) then
      allocate( Gsrf_r(ndt,nspc))
                call init_zero(Gsrf_r)
          else
             call face_error('Gsrf_r already allocated')
          endif

          if (.not.allocated(Gabs_l)) then
      allocate(Gabs_l(ndt,nspc))
                call init_zero(Gabs_l)
          else
              call face_error('Gabs_l already allocated')
          endif

          if (.not.allocated(Gdes_l)) then
      allocate( Gdes_l(ndt,nspc))
                call init_zero(Gdes_l)
          else
              call face_error('Gdes_l already allocated')
          endif


          if (.not.allocated(Gb_l)) then
      allocate( Gb_l(ndt,nspc))
                call init_zero(Gb_l)
          else
              call face_error('Gb_l already allocated')
          endif

          if (.not.allocated(Gads_l)) then
      allocate( Gads_l(ndt,nspc))
                call init_zero(Gads_l)
          else
              call face_error('Gads_l already allocated')
          endif

          if (.not.allocated(Gabs_r)) then
      allocate(Gabs_r(ndt,nspc))
                call init_zero(Gabs_r)
          else
              call face_error('Gabs_r already allocated')
          endif

          if (.not.allocated(Gdes_r)) then
      allocate( Gdes_r(ndt,nspc))
                call init_zero(Gdes_r)
          else
              call face_error('Gdes_r already allocated')

          endif



          if (.not.allocated(Gb_r)) then
      allocate(Gb_r(ndt,nspc))
                call init_zero(Gb_r)
          else
              call face_error('Gb_r already allocated')

          endif

          if (.not.allocated(Gads_r)) then
      allocate(Gads_r(ndt,nspc))
                call init_zero(Gads_r)
          else
              call face_error('Gads_r already allocated')

          endif

          if (.not.allocated(jout)) then
      allocate(jout(ndt,nspc))
                call init_zero(jout)
          else
              call face_error('jout already allocated')
          endif

               if (.not.allocated(init_inventory)) then
      allocate(init_inventory(nspc))

          else
              call face_error('init_inventory already allocated')
          endif



               if (.not.allocated(final_inventory)) then
      allocate(final_inventory(nspc))

          else
              call face_error('final_inventory already allocated')
          endif

               if (.not.allocated(trace_flux)) then
      allocate(trace_flux(nspc))

          else
              call face_error('trace_outgassing already allocated')
          endif
          if (.not.allocated(onthefly_inventory)) then
      allocate(onthefly_inventory(nspc))

          else
              call face_error('onthefly_inventory already allocated')
          endif



      end subroutine alloc_variables

      subroutine alloc_input_species()

          allocate(namespc(nspc))
          call init_zero(namespc)

          allocate(dens0(nspc))
          call init_zero(dens0)

          allocate(gprof(nspc))
          call init_zero(gprof)

          allocate(gxmax(nspc))
          call init_zero(gxmax)

          allocate(gsigm(nspc))
          call init_zero(gsigm)

          allocate(densm(nspc))
          call init_zero(densm)

          allocate(cdif0(nspc))
          call init_zero(cdif0)

          allocate(edif   (nspc))
          call init_zero(edif)

          allocate(etr   (nspc))
          call init_zero(etr)

          allocate(edtr   (nspc))
          call init_zero(edtr)
          allocate(left_surface_model(nspc))
          call init_zero(left_surface_model)
          allocate(right_surface_model(nspc))
          call init_zero(right_surface_model)
          allocate(left_surface_model_string(nspc))
          call init_zero(left_surface_model_string)
          allocate(right_surface_model_string(nspc))
          call init_zero(right_surface_model_string)


          allocate(order_desorption_left(nspc))
          call init_zero(order_desorption_left)

          allocate(order_desorption_right(nspc))
          call init_zero(order_desorption_right)
          allocate(order_desorption_left_sat(nspc))
          call init_zero(order_desorption_left_sat)

          allocate(order_desorption_right_sat(nspc))
          call init_zero(order_desorption_right_sat)

          allocate(dsrfm(nspc))
          call init_zero(dsrfm)

          allocate(dsrfl0(nspc))
          call init_zero(dsrfl0)

          allocate(dsrfr0(nspc))
          call init_zero(dsrfr0)

          allocate(Eabs_l(nspc))
          call init_zero(Eabs_l)

          allocate(Edes_l(nspc))
          call init_zero(Edes_l)

          allocate(Edes_lsat(nspc))
          call init_zero(Edes_lsat)

          allocate(Eb_l(nspc))
          call init_zero(Eb_l)

          allocate(Eads_l(nspc))
          call init_zero(Eads_l)

          allocate(Eabs_r(nspc))
          call init_zero(Eabs_r)

          allocate(Edes_r(nspc))
          call init_zero(Edes_r)
          allocate(Edes_rsat(nspc))
          call init_zero(Edes_rsat)

          allocate(Eb_r(nspc))
          call init_zero(Eb_r)

          allocate(Eads_r(nspc))
          call init_zero(Eads_r)

          allocate(nu(nspc))
          call init_zero(nu)
          ! implantation
          allocate(implantation_model(nspc))
          call init_zero(implantation_model)

          allocate(implantation_depth(nspc))
          call init_zero(implantation_depth)

           allocate(diagnostic_depth(nspc))
          call init_zero(diagnostic_depth)

          allocate(j_implantation_depth(nspc))
          call init_zero(j_implantation_depth)

          allocate(j_diagnostic_depth(nspc))
          call init_zero(j_diagnostic_depth)

          allocate(implantation_width(nspc))
          call init_zero(implantation_width)

          allocate(enrg(nspc))
          call init_zero(enrg)
allocate(Gamma_in_max(nspc))

          call init_zero(Gamma_in_max)
          allocate(Gamma_in_pulse(nspc))
          allocate(Gamma_in_pulse_string(nspc))
          allocate(Gamma_in_pulse_period(nspc))
          allocate(Gamma_in_pulse_duration(nspc))
          call init_zero(Gamma_in_pulse_period)
          allocate(Gamma_in_pulse_starttime(nspc))
          call init_zero(Gamma_in_pulse_starttime)
          allocate(Gamma_in_base(nspc))
          call init_zero(Gamma_in_base)


          allocate(gas_pressure  (nspc))
          call init_zero(gas_pressure)

          allocate( gas_temp(nspc))
          call init_zero(gas_temp)

          allocate( mass(nspc))
          call init_zero(mass)

        if (verbose_input) write(iout,*) 'Allocation of input parameters for species: OK'
    end subroutine alloc_input_species

subroutine deallocate_variables()
        if (allocated(namespc)) deallocate(namespc)
        if (allocated(dens0)) deallocate(dens0)
        if (allocated(dens)) deallocate(dens)
        if (allocated(flx)) deallocate(flx)
        if (allocated(qflx)) deallocate(qflx)
        if (allocated(rate_d)) deallocate(rate_d)
        if (allocated(rate_t)) deallocate(rate_t)
        if (allocated(cdif)) deallocate(cdif)
if (allocated(rct)) deallocate(rct)
if (allocated(ero_flx)) deallocate(ero_flx)
if (allocated(dif_flx)) deallocate(dif_flx)
if (allocated(ero_qflx)) deallocate(ero_qflx)
        if (allocated(gprof)) deallocate(gprof)
        if (allocated(gxmax)) deallocate(gxmax)
        if (allocated(gsigm)) deallocate(gsigm)
        if (allocated(densm)) deallocate(densm)

        if (allocated(cdif0)) deallocate(cdif0)
        if (allocated(edif)) deallocate(edif   )
        if (allocated(etr )) deallocate(etr   )
        if (allocated(edtr )) deallocate(edtr   )


        if (allocated(Eabs_l)) deallocate(Eabs_l)
        if (allocated(Edes_l)) deallocate(Edes_l)
        if (allocated(Edes_lsat)) deallocate(Edes_lsat)
         if (allocated(Eb_l)) deallocate(Eb_l)
        if (allocated(Eads_l))  deallocate(Eads_l)
         if (allocated(Eabs_r)) deallocate(Eabs_r)
        if (allocated(Edes_r))  deallocate(Edes_r)
        if (allocated(Edes_rsat))  deallocate(Edes_rsat)
        if (allocated(Eb_r))  deallocate(Eb_r)
         if (allocated(Eads_r)) deallocate(Eads_r)

         if (allocated(nu)) deallocate(nu)
      if (allocated(temp)) deallocate(temp)
         ! implantation

         if (allocated(implantation_model)) deallocate(implantation_model)
         if (allocated( implantation_depth)) deallocate( implantation_depth)
         if (allocated(j_implantation_depth)) deallocate( j_implantation_depth)
         if (allocated( diagnostic_depth)) deallocate( diagnostic_depth)
         if (allocated(j_diagnostic_depth)) deallocate( j_diagnostic_depth)
         if (allocated( implantation_width)) deallocate( implantation_width)
         if (allocated(enrg)) deallocate(enrg)
         if (allocated(inflx)) deallocate(inflx)
         if (allocated(Gamma_in_max)) deallocate(Gamma_in_max)
         if (allocated(Gamma_in_base)) deallocate(Gamma_in_base)
         if (allocated(Gamma_in_pulse_period)) deallocate(Gamma_in_pulse_period)
         if (allocated(Gamma_in_pulse_duration)) deallocate(Gamma_in_pulse_duration)
         if (allocated(Gamma_in_pulse_starttime)) deallocate(Gamma_in_pulse_starttime)
         if (allocated(Gamma_in_pulse)) deallocate(Gamma_in_pulse)
         if (allocated(Gamma_in_pulse_string)) deallocate(Gamma_in_pulse_string)
         if (allocated(gas_pressure)) deallocate(gas_pressure  )
         if (allocated(gas_temp)) deallocate(gas_temp)
         if (allocated(mass)) deallocate( mass)
            if (allocated(x))     deallocate(x)
       if (allocated(dx)) deallocate(dx)

       if (allocated(srs)) deallocate(srs)
       if (allocated(srb)) deallocate(srb)
       if (allocated(jout)) deallocate(jout)
       if (allocated(src)) deallocate(src)
       if (allocated(src_profile)) deallocate(src_profile)
       if (allocated(j0)) deallocate(j0)
       ! reaction
       if (allocated(nuth0)) deallocate(nuth0)
       if (allocated(nuth)) deallocate(nuth)
       if (allocated(kbin0)) deallocate(kbin0)
       if (allocated(kbin)) deallocate(kbin)
       if (allocated(ebin)) deallocate(ebin)
       if (allocated(eth)) deallocate(eth)
! surface
        if (allocated(dsrfm)) deallocate(dsrfm)
        if (allocated(dsrfl0)) deallocate(dsrfl0)
        if (allocated(dsrfr0)) deallocate(dsrfr0)
        if (allocated(dsrfl)) deallocate(dsrfl)
        if (allocated(dsrfr)) deallocate(dsrfr)
       if (allocated(right_surface_model)) deallocate(right_surface_model)
       if (allocated(left_surface_model)) deallocate(left_surface_model)
       if (allocated(right_surface_model_string)) deallocate(right_surface_model_string)
       if (allocated(left_surface_model_string)) deallocate(left_surface_model_string)

       if (allocated(order_desorption_left)) deallocate(order_desorption_left)
       if (allocated(order_desorption_right)) deallocate(order_desorption_right)
       if (allocated(order_desorption_left_sat)) deallocate(order_desorption_left_sat)
       if (allocated(order_desorption_right_sat)) deallocate(order_desorption_right_sat)


      if (allocated(Kabs_l)) deallocate(Kabs_l)
      if (allocated(Kdes_l)) deallocate(Kdes_l)
      if (allocated(Kb_l)) deallocate(Kb_l)
      if (allocated(Kads_l)) deallocate(Kads_l)
      if (allocated(Kabs_r)) deallocate(Kabs_r)
      if (allocated(Kdes_r)) deallocate(Kdes_r)
      if (allocated(Kb_r)) deallocate(Kb_r)
      if (allocated(Kads_r)) deallocate(Kads_r)
      if (allocated(K0abs_l)) deallocate(K0abs_l)
      if (allocated(K0des_l)) deallocate(K0des_l)
      if (allocated(K0b_l)) deallocate(K0b_l)
      if (allocated(K0ads_l)) deallocate(K0ads_l)
      if (allocated(K0abs_r)) deallocate(K0abs_r)
      if (allocated(K0des_r)) deallocate(K0des_r)
      if (allocated(K0b_r)) deallocate(K0b_r)
      if (allocated(K0ads_r)) deallocate(K0ads_r)

      if (allocated(Gads_l)) deallocate(Gads_l)
      if (allocated(Gads_r)) deallocate(Gads_r)
      if (allocated(Gdes_l)) deallocate(Gdes_l)
      if (allocated(Gdes_r)) deallocate(Gdes_r)
      if (allocated(Gabs_l)) deallocate(Gabs_l)
      if (allocated(Gabs_r)) deallocate(Gabs_r)
      if (allocated(Gb_l)) deallocate(Gb_l)
      if (allocated(Gb_r)) deallocate(Gb_r)
       if (allocated(Gsrf_l)) deallocate(Gsrf_l)
      if (allocated(Gsrf_r)) deallocate(Gsrf_r)


      if (allocated(trace_flux)) deallocate(trace_flux)
       if (allocated(onthefly_inventory)) deallocate(onthefly_inventory)
       if (allocated(init_inventory)) deallocate(init_inventory)
       if (allocated(final_inventory)) deallocate(final_inventory)
    end subroutine deallocate_variables


       end module modFACE_allocate
