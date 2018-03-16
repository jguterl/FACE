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
          if (.not.allocated(rate_t)) then
              allocate(rate_t(ndt,0:ngrd))
              call init_zero(rate_t)
          else
              call face_error('srs already allocated')

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
      allocate(kbin(nspc,nspc,nspc))
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
      allocate(nuth(nspc,nspc))
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

          if (.not.allocated(qchtr)) then
      allocate(qchtr(nspc))
                call init_zero(qchtr)
          else
              call face_error('qchtr already allocated')
          endif

          if (.not.allocated(qchtl)) then
      allocate(qchtl(nspc))
                call init_zero(qchtl)
          else
              call face_error('qchtl already allocated')
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

          allocate(dsrfm(nspc))
          call init_zero(dsrfm)

          allocate(dsrfl0(nspc))
          call init_zero(dsrfl0)

          allocate(dsrfr0(nspc))
          call init_zero(dsrfr0)

          allocate(echl(nspc))
          call init_zero(echl)

          allocate(qchl(nspc))
          call init_zero(qchl)

          allocate(ebl(nspc))
          call init_zero(ebl)

          allocate(esl(nspc))
          call init_zero(esl)

          allocate(echr(nspc))
          call init_zero(echr)

          allocate(qchr(nspc))
          call init_zero(qchr)

          allocate(ebr(nspc))
          call init_zero(ebr)

          allocate(esr(nspc))
          call init_zero(esr)

          allocate(nu(nspc))
          call init_zero(nu)

          allocate(enrg(nspc))
          call init_zero(enrg)

          allocate(inflx_in_max(nspc))
          call init_zero(inflx_in_max)

          allocate(inflx_in(nspc))
          call init_zero(inflx_in)

          allocate(gas_pressure  (nspc))
          call init_zero(gas_pressure)

          allocate( gas_temp(nspc))
          call init_zero(gas_temp)

          allocate( mass(nspc))
          call init_zero(mass)

        if (verbose_input) write(iout,*) 'Allocation of input parameters for species: OK'
    end subroutine alloc_input_species

subroutine deallocate_variables()
        deallocate(namespc)
        deallocate(dens0)
        deallocate(gprof)
        deallocate(gxmax)
        deallocate(gsigm)
        deallocate(densm)

        deallocate(cdif0)
        deallocate(edif   )
        deallocate(etr   )
        deallocate(edtr   )

        deallocate(dsrfm)
        deallocate(dsrfl0)
        deallocate(dsrfr0)
        deallocate(echl)
        deallocate(qchl)
        deallocate(ebl)
        deallocate(esl)
        deallocate(echr)
        deallocate(qchr)
        deallocate(ebr)
        deallocate(esr)

        deallocate(nu)
        deallocate(enrg)
        deallocate(inflx)
        deallocate(inflx_in_max)
        deallocate(inflx_in)
        deallocate(gas_pressure  )
        deallocate(gas_temp)
        deallocate( mass)
               deallocate(x)
      deallocate(dx)

      deallocate(srs)
      deallocate(srb)
      deallocate(jout)
      deallocate(src)


      deallocate(Kabs_l)
      deallocate(Kdes_l)
      deallocate(Kb_l)
      deallocate(Kads_l)
      deallocate(Kabs_r)
      deallocate(Kdes_r)
      deallocate(Kb_r)
      deallocate(Kads_r)
      deallocate(K0abs_l)
      deallocate(K0des_l)
      deallocate(K0b_l)
      deallocate(K0ads_l)
      deallocate(K0abs_r)
      deallocate(K0des_r)
      deallocate(K0b_r)
      deallocate(K0ads_r)
      deallocate(temp)
      deallocate(Gads_l)
      deallocate(Gads_r)
      deallocate(Gdes_l)
      deallocate(Gdes_r)
      deallocate(Gabs_l)

      deallocate(Gb_l)
      deallocate(Gb_r)
      deallocate(trace_flux)
       deallocate(init_inventory)
       deallocate(final_inventory)
    end subroutine deallocate_variables


       end module modFACE_allocate
