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
         write(iout,*) 'ERROR: trying to set value to non-allocated variable'
         stop 'Exiting FACE...'
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
         write(iout,*) 'ERROR: trying to set value to non-allocated variable'
         stop 'Exiting FACE...'
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
         write(iout,*) 'ERROR: trying to set value to non-allocated variable'
         stop 'Exiting FACE...'
         endif
        r=0.d0

    end subroutine init_zero_rarr4

    subroutine init_zero_rarr(r)
        real(DP),allocatable::r(:)
          if (.not.allocated(r))  then
         write(iout,*) 'ERROR: trying to set value to non-allocated variable'
         stop 'Exiting FACE...'
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
              write(iout,*) 'ERROR: x already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(dx)) then
              allocate(dx(0:ngrd))
              call init_zero(dx)
          else
              write(iout,*) 'ERROR: dx already allocated'
              STOP 'Exiting FACE...'
          endif

          ! temperature
          if (.not.allocated(temp)) then
          allocate( temp(ndt,0:ngrd))
          call init_zero(temp)
          else
              write(iout,*) 'ERROR: dx already allocated'
              STOP 'Exiting FACE...'
          endif
          ! surface

          allocate(inflx(nspc))

          ! bulk
          if (.not.allocated(dens)) then
              allocate(dens(ndt,0:ngrd,nspc))
              call init_zero(dens)
          else
              write(iout,*) 'ERROR: dens already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(flx)) then
              allocate(flx (ndt,0:ngrd,nspc))
              call init_zero(flx)
          else
              write(iout,*) 'ERROR: flx already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(rtd)) then
              allocate(rtd (ndt,0:ngrd,nspc))
              call init_zero(rtd)
          else
              write(iout,*) 'ERROR: rtd already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(cdif)) then
              allocate(cdif(ndt,0:ngrd,nspc))
              call init_zero(cdif)
          else
              write(iout,*) 'ERROR: cdif already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(rct)) then
              allocate(rct (ndt,0:ngrd,nspc))
              call init_zero(rct)
          else
              write(iout,*) 'ERROR: rct already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(ero)) then
              allocate(ero (ndt,0:ngrd,nspc))
              call init_zero(ero)
          else
              write(iout,*) 'ERROR: ero already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(rtt)) then
              allocate(rtt(ndt,0:ngrd))
              call init_zero(rtt)
          else
              write(iout,*) 'ERROR: srs already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(flxt)) then
              allocate(flxt(ndt,0:ngrd))
              call init_zero(flxt)
          else
              write(iout,*) 'ERROR: flxt already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(erot)) then
              allocate(erot(ndt,0:ngrd))
              call init_zero(erot)
          else
              write(iout,*) 'ERROR: srs already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(srs)) then
              allocate(srs(ndt,0:ngrd,nspc))
              call init_zero(srs)
          else
              write(iout,*) 'ERROR: srs already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(srb)) then
              allocate(srb (ndt,0:ngrd,nspc,nspc))
              call init_zero(srb)
          else
              write(iout,*) 'ERROR: srb already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(src)) then
              allocate(src (ndt,0:ngrd,nspc))
              call init_zero(src)
          else
              write(iout,*) 'ERROR: src already allocated'
              STOP 'Exiting FACE...'
          endif

         ! reaction parameters
         if (.not.allocated(kbin0)) then
     allocate(kbin0(nspc,nspc,nspc))
         call init_zero(kbin0)
          else
              write(iout,*) 'ERROR: kbin0 already allocated'
              STOP 'Exiting FACE...'
          endif
     if (.not.allocated(ebin)) then
      allocate(ebin(nspc,nspc,nspc))
          call init_zero(ebin)
          else
              write(iout,*) 'ERROR: ebin already allocated'
              STOP 'Exiting FACE...'
          endif
      if (.not.allocated(kbin)) then
      allocate(kbin(nspc,nspc,nspc))
          call init_zero(kbin)
          else
              write(iout,*) 'ERROR: kbin already allocated'
              STOP 'Exiting FACE...'
          endif


      if (.not.allocated(nuth0)) then
       allocate(nuth0(nspc,nspc))
           call init_zero(nuth0)
          else
              write(iout,*) 'ERROR: nuth0 already allocated'
              STOP 'Exiting FACE...'
          endif

       if (.not.allocated(eth)) then
      allocate(eth(nspc,nspc))
          call init_zero(eth)
          else
              write(iout,*) 'ERROR: eth already allocated'
              STOP 'Exiting FACE...'
          endif
      if (.not.allocated(nuth)) then
      allocate(nuth(nspc,nspc))
          call init_zero(nuth)
          else
              write(iout,*) 'ERROR: nuth already allocated'
              STOP 'Exiting FACE...'
          endif

      ! surface
      if (.not.allocated(k1l)) then
      allocate(k1l(nspc))
                call init_zero(k1l)
          else
              write(iout,*) 'ERROR: k1l already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(k2l)) then
      allocate(k2l(nspc))
                call init_zero(k2l)
          else
              write(iout,*) 'ERROR: k2l already allocated'
              STOP 'Exiting FACE...'
          endif
          if (.not.allocated(k3l)) then
      allocate(k3l(nspc))
                call init_zero(k3l)
          else
              write(iout,*) 'ERROR: k3l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(k4l)) then
      allocate(k4l(nspc))
                call init_zero(k4l)
          else
              write(iout,*) 'ERROR: k4l already allocated'
              STOP 'Exiting FACE...'
          endif


          if (.not.allocated(k1r)) then
      allocate(k1r(nspc))
                call init_zero(k1r)
          else
              write(iout,*) 'ERROR: k1r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(k2r)) then
      allocate(k2r(nspc))
                call init_zero(k2r)
          else
              write(iout,*) 'ERROR: k2r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(k3r)) then
      allocate(k3r(nspc))
                call init_zero(k3r)
          else
              write(iout,*) 'ERROR: k3r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(k4r)) then
      allocate(k4r(nspc))
                call init_zero(k4r)
          else
              write(iout,*) 'ERROR: k4r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r1l)) then
      allocate(r1l(nspc))
                call init_zero(r1l)
          else
              write(iout,*) 'ERROR: r1l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r2l)) then
      allocate(r2l(nspc))
                call init_zero(r2l)
          else
              write(iout,*) 'ERROR: r2l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r3l)) then
      allocate(r3l(nspc))
                call init_zero(r3l)
          else
              write(iout,*) 'ERROR: r3l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r4l)) then
      allocate(r4l(nspc))
                call init_zero(r4l)
          else
              write(iout,*) 'ERROR: r4l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r1r)) then
      allocate(r1r(nspc))
                call init_zero(r1r)
          else
              write(iout,*) 'ERROR: r1r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r2r)) then
      allocate(r2r(nspc))
                call init_zero(r2r)
          else
              write(iout,*) 'ERROR: r2r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r3r)) then
      allocate(r3r(nspc))
                call init_zero(r3r)
          else
              write(iout,*) 'ERROR: r3r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(r4r)) then
      allocate(r4r(nspc))
                call init_zero(r4r)
          else
              write(iout,*) 'ERROR: r4r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j0)) then
      allocate(j0(nspc))
                call init_zero(j0)
          else
              write(iout,*) 'ERROR: j0 already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(qchtr)) then
      allocate(qchtr(nspc))
                call init_zero(qchtr)
          else
              write(iout,*) 'ERROR: qchtr already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(qchtl)) then
      allocate(qchtl(nspc))
                call init_zero(qchtl)
          else
              write(iout,*) 'ERROR: qchtl already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(dsrfl)) then
      allocate(dsrfl(ndt,nspc))
                call init_zero(dsrfl)
          else
              write(iout,*) 'ERROR: dsrfl already allocated'
              STOP 'Exiting FACE...'
          endif


          if (.not.allocated(rtsl)) then
      allocate( rtsl(ndt,nspc))
                call init_zero(rtsl)
          else
              write(iout,*) 'ERROR: rtsl already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(dsrfr)) then
      allocate(dsrfr(ndt,nspc))
                call init_zero(dsrfr)
          else
              write(iout,*) 'ERROR: dsrfr already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(rtsr)) then
      allocate( rtsr(ndt,nspc))
                call init_zero(rtsr)
          else
              write(iout,*) 'ERROR: rtsr already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j1l)) then
      allocate(j1l(ndt,nspc))
                call init_zero(j1l)
          else
              write(iout,*) 'ERROR: j1l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j2l)) then
      allocate( j2l(ndt,nspc))
                call init_zero(j2l)
          else
              write(iout,*) 'ERROR: j2l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j3l)) then
      allocate( j3l(ndt,nspc))
                call init_zero(j3l)
          else
              write(iout,*) 'ERROR: j3l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j4l)) then
      allocate( j4l(ndt,nspc))
                call init_zero(j4l)
          else
              write(iout,*) 'ERROR: j4l already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j1r)) then
      allocate(j1r(ndt,nspc))
                call init_zero(j1r)
          else
              write(iout,*) 'ERROR: j1r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j2r)) then
      allocate( j2r(ndt,nspc))
                call init_zero(j2r)
          else
              write(iout,*) 'ERROR: j2r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j3r)) then
      allocate(j3r(ndt,nspc))
                call init_zero(j3r)
          else
              write(iout,*) 'ERROR:j3r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(j4r)) then
      allocate(j4r(ndt,nspc))
                call init_zero(j4r)
          else
              write(iout,*) 'ERROR: j4r already allocated'
              STOP 'Exiting FACE...'
          endif

          if (.not.allocated(jout)) then
      allocate(jout(ndt,nspc))
                call init_zero(jout)
          else
              write(iout,*) 'ERROR: jout already allocated'
              STOP 'Exiting FACE...'
          endif



      end subroutine alloc_variables

      subroutine alloc_input_species()

          allocate(namespc(nspc))
          call init_zero(namespc)

          allocate(dens0(nspc))
          call init_zero(namespc)

          allocate(gprof(nspc))
          call init_zero(namespc)

          allocate(gxmax(nspc))
          call init_zero(namespc)

          allocate(gsigm(nspc))
          call init_zero(namespc)

          allocate(densm(nspc))
          call init_zero(namespc)

          allocate(cdif0(nspc))
          call init_zero(namespc)

          allocate(edif   (nspc))
          call init_zero(namespc)

          allocate(etr   (nspc))
          call init_zero(namespc)

          allocate(edtr   (nspc))
          call init_zero(namespc)

          allocate(dsrfm(nspc))
          call init_zero(namespc)

          allocate(dsrfl0(nspc))
          call init_zero(namespc)

          allocate(dsrfr0(nspc))
          call init_zero(namespc)

          allocate(echl(nspc))
          call init_zero(namespc)

          allocate(qchl(nspc))
          call init_zero(namespc)

          allocate(ebl(nspc))
          call init_zero(namespc)

          allocate(esl(nspc))
          call init_zero(namespc)

          allocate(echr(nspc))
          call init_zero(namespc)

          allocate(qchr(nspc))
          call init_zero(namespc)

          allocate(ebr(nspc))
          call init_zero(namespc)

          allocate(esr(nspc))
          call init_zero(namespc)

          allocate(nu(nspc))
          call init_zero(namespc)

          allocate(enrg(nspc))
          call init_zero(namespc)

          allocate(inflx_max(nspc))
          call init_zero(namespc)

          allocate(inflx_min(nspc))
          call init_zero(namespc)

          allocate(prg  (nspc))
          call init_zero(namespc)

          allocate( tg(nspc))
          call init_zero(namespc)

          allocate( mass(nspc))
          call init_zero(namespc)

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
        deallocate(inflx_max)
        deallocate(inflx_min)
        deallocate(prg  )
        deallocate( tg)
        deallocate( mass)
               deallocate(x)
      deallocate(dx)

      deallocate(srs)
      deallocate(srb)
      deallocate(jout)
      deallocate(src)


      deallocate(k1l)
      deallocate(k2l)
      deallocate(k3l)
      deallocate(k4l)
      deallocate(k1r)
      deallocate(k2r)
      deallocate(k3r)
      deallocate(k4r)
      deallocate(r1l)
      deallocate(r2l)
      deallocate(r3l)
      deallocate(r4l)
      deallocate(r1r)
      deallocate(r2r)
      deallocate(r3r)
      deallocate(r4r)
      deallocate(temp)
    end subroutine deallocate_variables


       end module modFACE_allocate
