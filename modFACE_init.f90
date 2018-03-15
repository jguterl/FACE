      module modFACE_init
      use modFACE_precision
       use modFACE_header
      use modFACE_functions
      use modFACE_step
       use modFACE_output
       use modFACE_input
       use modFACE_IO
       use modFACE_error
      implicit none
      integer ::ifile_Tramp=200
    !  real(DP)::K0abs_l(nspc)
      contains
      subroutine initialize()
!     ******************************************************************
!     * initialization of arrays and parameters                        *
!     *                                                                *
!     * Author: Roman D. Smirnov                                       *
!     * E-mail: rsmirnov@ucsd.edu; rosmirnov@yahoo.com                 *
!     *                                                                *
!     ******************************************************************


!
!      integer i, j, k, l, m, n
!     integer seed, val(1:8)
!     integer lnblnk, unt
!     integer ngrd2
!     real*8 rand
!     integer k
!     character*15 name




!     ------------------------------------------------------------------
!      initialization of species parameters
!     ------------------------------------------------------------------
       call alloc_variables
      call init_misc
      call init_time
 !     call init_seed
      call init_grid
      call init_temp
      call flx_update()
      call init_volume_species
      call init_source
      call init_boundary
      call init_reactions

      call print_milestone('initialization done')
      end subroutine initialize

      subroutine init_misc()
      ! default restart filename
      restart_filename=trim(path_folder)//"dsave.rst"


      ! some material constants
      lambda1c=lambda*clng
      lambda2c=lambda*lambda*csrf
      lambda3c=lambda*lambda*lambda*cvlm

      rhocp=rho*cp

!      tmp=rand(seed)

      ! counter for filename numbers
      sfln_voldata=0
      sfln_srfdata=0
      sfln_heatdata=0

      cnt =0

      ! number fo equations to solve
      if (solve_heat_eq .eq."no") then
       neq=nspc*(ngrd+3)
      else
       neq=nspc*(ngrd+3)+ngrd+1
      endif

      call open_timedata_files

      end subroutine init_misc


!      subroutine init_seed()
!      integer i
!      character*10 cdate, ctime, czone
!      call date_and_time (cdate, ctime, czone, val)
!      seed=0
!
!      do i=1,8
!       seed=seed+abs(val(i))
!      enddo
!      end subroutine init_seed

      subroutine init_reactions()
      integer i,j,k,l,m

!     ------------------------------------------------------------------
!      initialization of reaction constatnts
!     ------------------------------------------------------------------
      do k=1,nspc
       do l=1,nspc
        do m=1,nspc
         kbin0(k,l,m)=kbinar(k,l,m)
         ebin (k,l,m)=ebinar(k,l,m)
         kbin (k,l,m)=kbin0 (k,l,m)*exp(-ee*ebin(k,l,m)/(kb*temp(ndt,0)))
        enddo
        nuth0(k,l)=ktherm(k,l)
        eth  (k,l)=etherm(k,l)
        nuth (k,l)=nuth0 (k,l)*exp(-ee*eth(k,l)/(kb*temp(ndt,0)))
       enddo
      enddo
      do i=1,ndt
       do j=0,ngrd
        do k=1,nspc
         rct(i,j,k)=0.d0
         do l=1,nspc
          rct (i,j,k)=rct   (i,j,k)+nuth  (k,l)*dens(i,j,l)*ctherm(i,j,k,l)
          do m=1,l
           rct (i,j,k)=rct   (i,j,k)+kbin  (k,l,m)*dens(i,j,l)*dens(i,j,m)*cbinar(i,j,k,l,m)
          enddo
         enddo
        enddo
       enddo
      enddo
      if (verbose_init) write(iout,*) "Initialization reaction terms: DONE"
      end subroutine init_reactions

      subroutine init_time()
!     ------------------------------------------------------------------
!      time step initialization
!     ------------------------------------------------------------------
      dt_face=dtmin
!      dt=cdt*ttm
      time=start_time
      if (verbose_init) write(iout,*) "Initialization time parameters: DONE"
      end subroutine init_time


      subroutine init_grid()
      integer:: j,ngrd2
      real(DP)::dx0
!      ------------------------------------------------------------------
!      initialization of grid arrays
!     ------------------------------------------------------------------

      x(0)=0.d0
      if (alpha .eq. 1.d0) then
       dx0=length/ngrd
       do j=1,ngrd
        x  (j  )=j*dx0
        dx (j-1)=x(j)-x(j-1)
       enddo
       dx (ngrd)=dx(ngrd-1)
      else
       ngrd2=ngrd/2
       if (mod(ngrd,2) .ne. 0) then
        dx0=length*(alpha-1.d0)/((1.d0+alpha)*alpha**ngrd2-2.d0)
        x (0)=0.d0
        dx(0)=dx0
        do j=1,ngrd2
         x (j)=x(j-1)+dx(j-1)
         dx(j)=dx(j-1)*alpha
        enddo
        do j=ngrd2+1,ngrd
         x (j)=x(j-1)+dx(j-1)
         dx(j)=dx(j-1)/alpha
        enddo
       else
        dx0=0.5d0*length*(alpha-1.d0)/(alpha**ngrd2-1.d0)
        x (0)=0.d0
        dx(0)=dx0
        do j=1,ngrd2-1
         x (j)=x(j-1)+dx(j-1)
         dx(j)=dx(j-1)*alpha
        enddo
        x (ngrd2)=x (ngrd2-1)+dx(ngrd2-1)
        dx(ngrd2)=dx(ngrd2-1)
        do j=ngrd2+1,ngrd
         x (j)=x(j-1)+dx(j-1)
         dx(j)=dx(j-1)/alpha
        enddo
       endif
       dx (ngrd)=dx(ngrd-1)
      endif
      write(iout,*) "Initialization x grid: DONE"
      end subroutine init_grid

      subroutine init_source()

           integer::i,j,k,l
!     ------------------------------------------------------------------
!      Initialization of sources
!     ------------------------------------------------------------------

      if (solve_heat_eq .eq. "yes") then
       do i=1,ndt
        do j=0,ngrd
         do k=1,nspc
          srs (i,j,k)=source(  j,k)
          src (i,j,k)=srs   (i,j,k)*csours(i,j,k)
          do l=1,nspc
           srb (i,j,k,l)=srcbin(  j,k,l)
           src (i,j,k  )=src   (i,j,k  )+srb(i,j,k,l)*dens(i,j,l)*csrbin(i,j,k,l)
          enddo
         enddo
        enddo
       enddo
      elseif (solve_heat_eq .eq. "no") then
       do i=1,ndt
        do j=0,ngrd
         do k=1,nspc
          srs (i,j,k)=source(  j,k)
          src (i,j,k)=srs   (i,j,k)*csours(i,j,k)
          jout(i,  k)=jout  (i,  k)+srs(i,j,k)*(1.d0-csours(i,j,k))*dx(j)
          do l=1,nspc
           srb (i,j,k,l)=srcbin(  j,k,l)
           src (i,j,k  )=src   (i,j,k  )+srb(i,j,k,l)*dens(i,j,l)*csrbin(i,j,k,l)
           jout(i,  k  )=jout  (i,  k  )+srb(i,j,k,l)*dens(i,j,l)*(1.d0-csrbin(i,j,k,l))*dx(j)
          enddo
         enddo
        enddo
       enddo
       else
        write(iout,*) "ERROR: Unknown option for solve_heat_eq:", solve_heat_eq
        stop
      endif
      write(iout,*) "Initialization source terms: DONE"
      end subroutine

      subroutine init_boundary()
      integer ::i,k
      real(DP)::tmp

!     ------------------------------------------------------------------
!      boundary parameters
!     ------------------------------------------------------------------
      if (verbose_init) write(iout,*) "Initialization boundary variables"
      do k=1,nspc
       K0abs_l(k)=1.d0
       K0des_l(k)=nu(k)*lambda2c
       K0b_l(k)=nu(k)
       K0ads_l(k)=nu(k)*lambda1c
       K0abs_r(k)=1.d0
       K0des_r(k)=nu(k)*lambda2c
       K0b_r(k)=nu(k)
       K0ads_r(k)=nu(k)*lambda1c
       if (mass(k)*gas_temp(k) .ne. 0.d0) then
        j0(k)=j0(k)+gas_pressure(k)/sqrt(twopi*mass(k)*ee*gas_temp(k))
       endif
       tmp=1.d2*(dsrfl0(k)/dsrfm(k)-1.d0)
       qchtl(k)=qchl(k)*0.5d0*(1.d0-erf(tmp))

       Kabs_l(k)=j0(k)*K0abs_l(k)*exp(-     ee* echl(k)          /(kb*temp(ndt,   0)))
       Kdes_l(k)=2.d0 *K0des_l(k)*exp(-2.d0*ee*(echl(k)+qchtl(k))/(kb*temp(ndt,   0)))
       Kb_l(k)=      K0b_l(k) *exp(-     ee*(ebl (k)+qchtl(k))/(kb*temp(ndt,   0)))
       Kads_l(k)=      K0ads_l(k) *exp(-     ee*(ebl (k)-esl  (k))/(kb*temp(ndt,   0)))
       tmp=1.d2*(dsrfr0(k)/dsrfm(k)-1.d0)
       qchtr(k)=qchr(k)*0.5d0*(1.d0-erf(tmp))
       Kabs_r(k)=j0(k)*K0abs_r(k) *exp(-     ee* echr(k)          /(kb*temp(ndt,ngrd)))
       Kdes_r(k)=2.d0 *K0des_r(k) *exp(-2.d0*ee*(echr(k)+qchtr(k))/(kb*temp(ndt,ngrd)))
       Kb_r(k)=      K0b_r(k)*exp(-     ee*(ebr (k)+qchtr(k))/(kb*temp(ndt,ngrd)))
       Kads_r(k)=      K0ads_r(k) *exp(-     ee*(ebr (k)-esr  (k))/(kb*temp(ndt,ngrd)))
       do i=1,ndt
        Gsrf_l(i,k)=0.d0
        Gsrf_r(i,k)=0.d0
        if (dsrfl0(k) .eq. 0.d0) then
         dsrfl(i,k)=1.d0
        else
         dsrfl(i,k)=dsrfl0(k)
        endif
        if (dsrfr0(k) .eq. 0.d0) then
         dsrfr(i,k)=1.d0
        else
         dsrfr(i,k)=dsrfr0(k)
        endif
!


        Gabs_l (i,k)=Kabs_l(k)
        Gdes_l (i,k)=Kdes_l(k)
        Gb_l (i,k)=Kb_l(k)
        Gads_l (i,k)=Kads_l(k)
        Gabs_r (i,k)=Kabs_r(k)
        Gdes_r (i,k)=Kdes_r(k)
        Gb_r (i,k)=Kb_r(k)
        Gads_r (i,k)=Kads_r(k)

        call set_cap_factor_surface(k,i)
!
        if (solve_heat_eq .eq. "yes") then
         jout(i,k)=jout(i,k)+Gdes_l(i,k)
        endif
!
       enddo
      enddo
      write(iout,*) "Initialize boundary: DONE"
            end subroutine init_boundary





      subroutine read_Tramp_file
      integer ios,i

        open(unit=ifile_Tramp, file=trim(framp), iostat=ios,action='read',status='old')
        if ( ios /= 0 ) then
        write(iout,*) 'Opening of temperature ramp file "', framp ,'" : FAIL '
        stop
        endif
        write(iout,*) 'Opening of temperature ramp file "', framp ,'" : DONE '

        read (ifile_Tramp, '(i4)') nramp
        allocate(rtime(nramp))
        allocate(rtemp(nramp))
        do i=1,nramp
         read (20, *) rtime(i), rtemp(i)
        enddo
        close(ifile_Tramp)
        write(iout,*) 'Reading of temperature ramp file "', framp ,'" : DONE '
        end subroutine read_Tramp_file


        subroutine init_temp()
            integer::i,j,n

            if (tramp1 .ne. tramp0) then
                dtemp=(temp1-temp0)/(tramp1-tramp0)
            else
                dtemp=0.d0
            endif

            if (framp .ne. 'none') then
                call read_Tramp_file
            endif

            if (solve_heat_eq .eq. "no") then
                if (framp .ne. 'none') then
                    do j=0,ngrd
                        do i=1,ndt
                            if (rtime(1) .gt. time) then
                                temp(i,j)=rtemp(1)
                            elseif (rtime(nramp) .gt. time) then
                                temp(i,j)=rtemp(nramp)
                            else
                                do n=2,nramp-1
                                    if (rtime(n) .gt. time) then
                                        temp(i,j)=rtemp(n-1)+(rtemp(n)-rtemp(n-1))*(time-rtime(n-1))/(rtime(n)-rtime(n-1))
                                        exit
                                    endif
                                enddo
                            endif
                        enddo ! i
                    enddo ! j

                else
                    do j=0,ngrd
                        do i=1,ndt
                            temp(i,j)=temp0
                        enddo
                    enddo
                endif
             !  solving heat equation -> initial linear profile of T in bulk
            elseif(solve_heat_eq .eq. "yes") then
                do j=0,ngrd
                    do i=1,ndt
                        temp(i,j)=temp0+(temp1-temp0)*x(j)/length
                    enddo
                enddo
            else
                write(iout,*) "ERROR: Unknown option for solve_heat_eq:", solve_heat_eq
                stop
            endif
            write(iout,*) 'Initialization of temperature : DONE '
        end subroutine init_temp

      subroutine init_volume_species()

      integer:: i,j,k
      do k=1,nspc
      write(iout,*) 'initial profile of density:',gprof(k)
       do j=0,ngrd
        do i=1,ndt
         flx (i,j,k)=0.d0
         ero_flx (i,j,k)=0.d0
         cdif(i,j,k)=cdif0(k)*exp(-ee*edif(k)/(kb*temp(i,j)))
         rtd (i,j,k)=0.d0
         if (gprof(k) .eq. 'S'.or.gprof(k) .eq. 'F') then
          dens(i,j,k)=dens0(k)
         elseif(gprof(k) .eq. 'G') then
          dens(i,j,k)=dens0(k) *exp(-0.5d0*abs((x(j)-gxmax(k))/gsigm(k))**2.d0)
         else
         write(iout,*) 'ERROR: unknow option for n0_profile (k=',k,'):', gprof(k)
         stop
         endif
         if (dens(i,j,k) .lt. 1.d0) then
          dens(i,j,k)=1.d0
         endif
        enddo
       enddo
      enddo
      do k=1,nspc
       do i=1,ndt
        jout(i,k)=0.d0
       enddo
      enddo

      end subroutine init_volume_species

subroutine init_casename(case_name)
   character(*)::case_name
   casename=trim(case_name)
   if (verbose_init) write(iout,*) 'Casename : ', casename
   casename=trim(casename)
    end subroutine init_casename


   subroutine init_path(path)

   character(*)::path
   character(30)::name
   integer ios,ifile_testpath

   name='test.path'
   path_folder=trim(path)
   if (verbose_init) write(iout,*) "Files will be saved in the folder :", path_folder
   call system('mkdir -p '//trim(path_folder))
   path_folder=trim(path_folder)//'/'
   dat_folder=trim(path_folder)//trim(casename)//'_dat'
   call system('mkdir -p '//dat_folder)

   call set_ifile(ifile_testpath)
    open (unit=ifile_testpath, file=trim(path_folder)//name,status='replace', form='unformatted', iostat=ios)
    if (ios.ne.0) then
    call face_error('cannot write in the folder : ', trim(path_folder)//name)
    endif
    close(ifile_testpath)

    end subroutine init_path




      end module modFACE_init
