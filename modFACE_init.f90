      module modFACE_init
      use modFACE_precision
       use modFACE_header
      use modFACE_functions
      use modFACE_step
       use modFACE_output
       use modFACE_input
       use modFACE_IO
      implicit none
      integer ::ifile_Tramp=200
    !  real(DP)::r1l(nspc)
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
      lambda1c=lambda*clng
      lambda2c=lambda*lambda*csrf
      lambda3c=lambda*lambda*lambda*cvlm

      rhocp=rho*cp

!      tmp=rand(seed)

      sfln_voldata=0
      sfln_srfdata=0
      sfln_heatdata=0
      cnt =0
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
      real(DP)::c1l, c2l, c3l, c4l
      real(DP)::c1r, c2r, c3r, c4r
      real(DP)::tmp

!     ------------------------------------------------------------------
!      boundary parameters
!     ------------------------------------------------------------------
      if (verbose_init) write(iout,*) "Initialization boundary variables"
      do k=1,nspc
       r1l(k)=1.d0
       r2l(k)=nu(k)*lambda2c
       r3l(k)=nu(k)
       r4l(k)=nu(k)*lambda1c
       r1r(k)=1.d0
       r2r(k)=nu(k)*lambda2c
       r3r(k)=nu(k)
       r4r(k)=nu(k)*lambda1c
       if (mass(k)*tg(k) .ne. 0.d0) then
        j0(k)=prg(k)/sqrt(twopi*mass(k)*ee*tg(k))
       else
        j0(k)=0.d0
       endif
       tmp=1.d2*(dsrfl0(k)/dsrfm(k)-1.d0)
       !write(iout,*) "247: tmp",tmp
       !write(iout,*) "249: qcht",qchl
       !write(iout,*) "erf(tmp)",erf(tmp)
       qchtl(k)=qchl(k)*0.5d0*(1.d0-erf(tmp))
       !write(iout,*) "249: qchtl",qchtl
       k1l(k)=j0(k)*r1l(k)*exp(-     ee* echl(k)          /(kb*temp(ndt,   0)))
       k2l(k)=2.d0 *r2l(k)*exp(-2.d0*ee*(echl(k)+qchtl(k))/(kb*temp(ndt,   0)))
       k3l(k)=      r3l(k) *exp(-     ee*(ebl (k)+qchtl(k))/(kb*temp(ndt,   0)))
       k4l(k)=      r4l(k) *exp(-     ee*(ebl (k)-esl  (k))/(kb*temp(ndt,   0)))
       tmp=1.d2*(dsrfr0(k)/dsrfm(k)-1.d0)
       qchtr(k)=qchr(k)*0.5d0*(1.d0-erf(tmp))
       k1r(k)=j0(k)*r1r(k) *exp(-     ee* echr(k)          /(kb*temp(ndt,ngrd)))
       k2r(k)=2.d0 *r2r(k) *exp(-2.d0*ee*(echr(k)+qchtr(k))/(kb*temp(ndt,ngrd)))
       k3r(k)=      r3r(k)*exp(-     ee*(ebr (k)+qchtr(k))/(kb*temp(ndt,ngrd)))
       k4r(k)=      r4r(k) *exp(-     ee*(ebr (k)-esr  (k))/(kb*temp(ndt,ngrd)))
       do i=1,ndt
        rtsl(i,k)=0.d0
        rtsr(i,k)=0.d0
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
         c2l=dsrfl(i,k)*dsrfl(i,k)
        else
         c2l=0.d0
        endif
        if (dsrfr(i,k) .gt. 0.d0) then
         c2r=dsrfr(i,k)*dsrfr(i,k)
        else
         c2r=0.d0
        endif

        if ((dsrfl(i,k) .gt. 0.d0) .and.(dens(i,0,k) .lt. densm(k))) then
         c3l=dsrfl(i,k)*(1.d0-dens(i,0,k)/densm(k))
        else
         c3l=0.d0
        endif
        if ((dsrfr(i,k) .gt. 0.d0) .and.(dens(i,ngrd,k) .lt. densm(k))) then
         c3r=dsrfr(i,k)*(1.d0-dens(i,ngrd,k)/densm(k))

        else
         c3r=0.d0
        endif

        c4l=dens(i,0   ,k)
        c4r=dens(i,ngrd,k)
!
        j1l (i,k)=k1l(k)*c1l
        j2l (i,k)=k2l(k)*c2l
        j3l (i,k)=k3l(k)*c3l

        !write(iout,*) "313: j3l (i,k) is nan",j3l(ndt,k)

        j4l (i,k)=k4l(k)*c4l
        j1r (i,k)=k1r(k)*c1r
        j2r (i,k)=k2r(k)*c2r
        j3r (i,k)=k3r(k)*c3r
        j4r (i,k)=k4r(k)*c4r
!
        if (solve_heat_eq .eq. "yes") then
         jout(i,k)=jout(i,k)+j2l(i,k)
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
         ero (i,j,k)=0.d0
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
   casename=trim(casename)//'_'
    end subroutine init_casename


   subroutine init_path(path)
   character(*)::path
   character(30)::name
   integer ios
   name='test.path'
   path_folder=trim(path)
   if (verbose_init) write(iout,*) "Files will be saved in the folder :", path_folder
   call system('mkdir -p '//trim(path_folder))
   path_folder=trim(path_folder)//'/'
   call system('mkdir -p '//trim(path_folder)//trim(casename)//'dat')
    open (unit=ifile_testpath, file=trim(path_folder)//name,status='replace', form='unformatted', iostat=ios)
    if (ios.ne.0) then
    stop
    write (iout, '(a)') 'ERROR: cannot write in the save/write folder : ', trim(path_folder)//name
    stop 'Exiting FACE...'
    endif
    close(ifile_testpath)
    end subroutine init_path




      end module modFACE_init
